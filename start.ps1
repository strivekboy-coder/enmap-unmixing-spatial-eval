$projectRoot = $PSScriptRoot
Set-Location -LiteralPath $projectRoot
$env:VIRTUAL_ENV = Join-Path $projectRoot '.venv'
$env:PATH = (Join-Path $env:VIRTUAL_ENV 'Scripts') + ';' + $env:PATH
$env:TEMP = Join-Path $projectRoot '.tmp'
$env:TMP = $env:TEMP
$env:MPLCONFIGDIR = Join-Path $env:TEMP 'mpl'
$env:MPLBACKEND = 'Agg'
New-Item -ItemType Directory -Path $env:TEMP -Force | Out-Null
Write-Host 'Project Python environment active. Run: python run.py --config config.yaml'
