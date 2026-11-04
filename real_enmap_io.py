"""Deadline-aware local EnMAP acquisition with credentials kept out of output."""
from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urljoin, urlparse, urlunparse

import requests


class StopWork(RuntimeError):
    pass


class Deadline:
    def __init__(self, timestamp, root, report_path):
        self.end = datetime.fromisoformat(timestamp).astimezone(timezone.utc)
        self.root = Path(root)
        self.report = self.root / report_path
        self.timer = None

    def remaining(self):
        return (self.end - datetime.now(timezone.utc)).total_seconds()

    def check(self, minimum_seconds=0):
        if self.remaining() <= minimum_seconds:
            raise StopWork('Hard deadline reached or insufficient time remains for the next operation.')

    def arm(self):
        self.check()
        # The inference code uses RF threads, not detached worker processes.
        # os._exit also stops an in-progress blocking read/model at the hard limit.
        def expire():
            try:
                with self.report.open('a', encoding='utf-8') as stream:
                    stream.write('\n\nStopped automatically at the hard deadline, 15:30 Europe/Berlin. '
                                 'See the last completed step above; unfinished work is not claimed.\n')
            finally:
                os._exit(124)
        self.timer = threading.Timer(self.remaining(), expire)
        self.timer.daemon = True
        self.timer.start()

    def close(self):
        if self.timer:
            self.timer.cancel()


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def local_environment(name):
    """Support user variables saved after the desktop app was started."""
    value = None
    if os.name == 'nt':
        import winreg
        try:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, 'Environment') as key:
                value, _ = winreg.QueryValueEx(key, name)
        except FileNotFoundError:
            pass
    return value or os.environ.get(name)


def credentials(settings):
    username = local_environment(settings['username_environment_variable'])
    password = local_environment(settings['password_environment_variable'])
    if not username or not password:
        raise StopWork('Download credentials are not configured. Set the two named local environment variables.')
    suffix = '-cat1distributor'
    if settings.get('use_ips_main_account_for_geoservice') and username.lower().endswith(suffix):
        username = username[:-len(suffix)]
    return username, password


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.part')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


class LoginForms(HTMLParser):
    def __init__(self):
        super().__init__()
        self.forms, self.current = [], None
        self.error_text, self.capture_depth = [], 0

    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if self.capture_depth:
            self.capture_depth += 1
        elif any(word in attributes.get('class', '').lower().split() for word in ('errors', 'error', 'alert-danger')):
            self.capture_depth = 1
        if tag == 'form':
            self.current = {'action': attributes.get('action', ''), 'fields': {}, 'password': False}
            self.forms.append(self.current)
        elif tag == 'input' and self.current is not None and attributes.get('name'):
            self.current['fields'][attributes['name']] = attributes.get('value', '')
            self.current['password'] |= attributes.get('type') == 'password'

    def handle_endtag(self, tag):
        if self.capture_depth:
            self.capture_depth -= 1
        if tag == 'form':
            self.current = None

    def handle_data(self, data):
        if self.capture_depth and data.strip():
            self.error_text.append(data.strip())


def login_submission_url(page_url, form_action):
    """Preserve CAS service when a relative login-form action omits the query.

    An explicit form service is preserved. The target is copied only within
    the same endpoint. URLs and form values are never printed.
    """
    page = urlparse(page_url)
    target = urlparse(urljoin(page_url, form_action))
    if (target.scheme, target.netloc, target.path) != (page.scheme, page.netloc, page.path):
        return urlunparse(target)
    original = dict(parse_qsl(page.query, keep_blank_values=True))
    fields = parse_qsl(target.query, keep_blank_values=True)
    if 'service' in original and 'service' not in dict(fields):
        fields.append(('service', original['service']))
    return urlunparse(target._replace(query=urlencode(fields)))


class Downloader:
    """Use DLR's CAS login form and an in-memory cookie session.

    Credentials are posted only to the observed DLR SSO password form. Cookies,
    form tokens, headers, credential values and login bodies are never logged.
    Every redirect is checked before following it; credentials are never forwarded.
    """
    def __init__(self, cfg, root, deadline):
        self.cfg, self.root, self.deadline = cfg, Path(root), deadline
        self.raw = self.root / cfg['raw_directory']
        self.raw.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.auth = lambda request: request  # No implicit netrc login.
        self.session_identity = None
        history_path = self.raw / 'download_status.json'
        self.history = json.loads(history_path.read_text(encoding='utf-8')) if history_path.exists() else []

    def checked_url(self, url):
        parsed = urlparse(url)
        if (parsed.scheme != 'https' or parsed.hostname not in
                {'download.geoservice.dlr.de', 'sso.eoc.dlr.de'} or parsed.username or parsed.password):
            raise StopWork('SSO redirected outside the approved DLR hosts; stopping before sending credentials.')
        return url

    def follow(self, response):
        for _ in range(10):
            if response.status_code not in (301, 302, 303, 307, 308):
                return response
            url = self.checked_url(urljoin(response.url, response.headers.get('Location', '')))
            response.close()
            self.deadline.check(5)
            response = self.session.get(url, stream=True, allow_redirects=False,
                                        timeout=(15, min(45, self.deadline.remaining())))
        raise StopWork('Too many SSO redirects; no credentials were logged.')

    def authenticated_response(self, url):
        identity = credentials(self.cfg['credentials'])  # Fresh persisted user variables.
        if self.session_identity != identity:
            self.session.cookies.clear()
            self.session_identity = identity
        response = self.session.get(url, stream=True, allow_redirects=False,
                                    timeout=(15, min(45, self.deadline.remaining())))
        response = self.follow(response)
        parsed = urlparse(response.url)
        if parsed.hostname != 'sso.eoc.dlr.de' or response.status_code != 200:
            return response
        parser = LoginForms()
        parser.feed(response.text)
        forms = [form for form in parser.forms if form['password']]
        if len(forms) != 1:
            response.close()
            raise StopWork('SSO password form is missing or ambiguous; no credentials were posted.')
        form = forms[0]
        action = self.checked_url(login_submission_url(response.url, form['action']))
        if urlparse(action).hostname != 'sso.eoc.dlr.de':
            raise StopWork('Password form action is not the DLR SSO host.')
        fields = form['fields']
        fields.update(username=identity[0], password=identity[1])
        response.close()
        self.deadline.check(5)
        response = self.session.post(action, data=fields, stream=True, allow_redirects=False,
                                     headers={'Origin': 'https://sso.eoc.dlr.de', 'Referer': response.url},
                                     timeout=(15, min(45, self.deadline.remaining())))
        fields.clear()
        if response.status_code in (307, 308):
            response.close()
            raise StopWork('SSO requested a credential-preserving redirect; stopping to prevent forwarding secrets.')
        response = self.follow(response)
        if urlparse(response.url).hostname == 'sso.eoc.dlr.de':
            errors = LoginForms()
            errors.feed(response.text)
            explanation = ' '.join(errors.error_text)
            for secret in sorted(identity, key=len, reverse=True):
                explanation = explanation.replace(secret, '[redacted]')
            response.close()
            raise StopWork('SSO did not establish a download session. '
                           + (explanation[:800] if explanation else 'Check the local login or additional account requirements.'))
        return response

    def item(self):
        self.deadline.check(5)
        url = self.cfg['stac_endpoint'].rstrip('/') + '/collections/' + self.cfg['collection'] + '/items/' + self.cfg['scene_id']
        try:
            response = self.session.get(url, params={'f': 'json'}, headers={'Accept': 'application/json'},
                                        timeout=min(60, self.deadline.remaining()))
        except requests.RequestException:
            raise StopWork('Public STAC item request failed; no credential details are included.') from None
        if response.status_code != 200:
            raise StopWork(f'Public STAC item HTTP {response.status_code}.')
        item = response.json()
        if item['id'] != self.cfg['scene_id']:
            raise StopWork('STAC returned another scene; refusing to continue.')
        write_json(self.raw / 'stac_item.json', item)
        return item

    def asset(self, item, key, filename):
        self.deadline.check(5)
        destination = self.raw / filename
        if destination.exists():
            return destination
        url = item['assets'][key]['href']
        parsed = urlparse(url)
        if parsed.scheme != 'https' or parsed.hostname != 'download.geoservice.dlr.de' or parsed.username or parsed.password:
            raise StopWork('Unexpected asset host or credentials embedded in URL; refusing to send credentials.')
        temporary = destination.with_suffix(destination.suffix + '.part')
        for attempt in range(3):
            self.deadline.check(5)
            try:
                response = self.authenticated_response(url)
            except requests.RequestException:
                self.history.append({'asset': key, 'time_utc': utc_now(), 'status': 'network_error'})
                if attempt == 2:
                    raise StopWork('DLR download network request failed after bounded retries.') from None
                time.sleep(3)
                continue
            status = response.status_code
            self.history.append({'asset': key, 'time_utc': utc_now(), 'status': status})
            write_json(self.raw / 'download_status.json', self.history)
            if status in (429, 500, 502, 503, 504):
                response.close()
                if attempt == 2:
                    raise StopWork(f'DLR asset service HTTP {status} after three bounded attempts; no data were substituted.')
                time.sleep(5 * (attempt + 1))
                continue
            if status != 200:
                response.close()
                raise StopWork(f'DLR asset HTTP {status}. Authentication/role or service access must be checked locally; no credentials are printed.')
            content_type = response.headers.get('Content-Type', '').lower()
            if 'html' in content_type:
                response.close()
                raise StopWork('Received an HTML login page instead of the requested asset; stopping.')
            expected = response.headers.get('Content-Length')
            try:
                with temporary.open('wb') as stream:
                    for block in response.iter_content(1024 * 1024):
                        self.deadline.check()
                        stream.write(block)
            except requests.RequestException:
                raise StopWork('Asset transfer interrupted; incomplete local .part file is not treated as data.') from None
            finally:
                response.close()
            if expected and temporary.stat().st_size != int(expected):
                raise StopWork('Asset byte count differs from Content-Length; incomplete file was not accepted.')
            temporary.replace(destination)
            self.history[-1].update(bytes=destination.stat().st_size, sha256=file_hash(destination))
            write_json(self.raw / 'download_status.json', self.history)
            print(f'Local asset complete: {key}, {destination.stat().st_size} bytes', flush=True)
            return destination
        raise StopWork('Asset download did not complete.')
