"""Analytical checks for date deduplication and radiometric harmonization."""
import unittest

import numpy as np
from shapely.geometry import box, mapping
from affine import Affine

from temporal import calibrate_reflectance,choose_daily_items,point_cache_signature,BANDS


def item(identifier,date,offset,flag=False):
    return {'id':identifier,'properties':{'datetime':date+'T10:00:00Z','earthsearch:boa_offset_applied':flag},
            'assets':{'red':{'href':'https://example.invalid/red.tif','raster:bands':[{'scale':.0001,'offset':offset}]}}}


class TemporalScientificChecks(unittest.TestCase):
    def test_point_cache_cannot_reuse_another_endmember_selection(self):
        source={'id':'same-scene','assets':{band:{'href':'https://example.invalid/'+band+'.tif'} for band in ['scl']+BANDS}}
        tx=Affine(30,0,0,0,-30,120)
        original=point_cache_signature([source],[0,1],[0,1],tx,(4,4))
        copied=point_cache_signature([source],np.array([0,1]),np.array([0,1]),tx,(4,4))
        changed=point_cache_signature([source],[0,2],[0,1],tx,(4,4))
        self.assertEqual(original,copied)
        self.assertNotEqual(original,changed)

    def test_date_specific_offset_prevents_false_temporal_change(self):
        items=[item('new','2018-05-02',-.1),item('old','2018-05-01',0)]
        # Same reflectance represented by different raw DN conventions. Reverse
        # item order also checks date alignment rather than positional matching.
        raw=np.array([[[2000]],[[3000]]],dtype='uint16')
        values,audit=calibrate_reflectance(raw,['2018-05-01','2018-05-02'],items,['red'])
        np.testing.assert_allclose(values[:,0,0],[.2,.2],atol=1e-6)
        self.assertEqual(audit[0]['item'],'old')

    def test_audited_override_prevents_double_applied_offset(self):
        items=[item('audited','2018-05-01',-.1,True)]
        values,audit=calibrate_reflectance(np.array([[[2000]]]),['2018-05-01'],items,['red'],['audited'])
        np.testing.assert_allclose(values,[ [[.2]] ],atol=1e-6)
        self.assertEqual(audit[0]['declared_offset'],-.1)
        self.assertEqual(audit[0]['offset'],0)
        self.assertTrue(audit[0]['audited_override'])

    def test_override_rejects_changed_provider_flag(self):
        with self.assertRaises(ValueError):
            calibrate_reflectance(np.array([[[2000]]]),['2018-05-01'],
                                  [item('changed','2018-05-01',-.1,False)],['red'],['changed'])

    def test_daily_selection_prefers_coverage_then_clearer_full_item(self):
        def scene(name,date,footprint,cloud):
            return {'id':name,'geometry':mapping(footprint),
                    'properties':{'datetime':date+'T10:00:00Z','eo:cloud_cover':cloud}}
        items=[scene('partial-clear','2018-05-01',box(0,0,1,1),0),
               scene('full-cloudy','2018-05-01',box(0,0,2,2),80),
               scene('full-clearer','2018-05-01',box(0,0,2,2),20),
               scene('second-date','2018-05-02',box(0,0,2,2),50)]
        chosen=choose_daily_items(items,(0,0,2,2))
        self.assertEqual([x['id'] for x in chosen],['full-clearer','second-date'])


if __name__=='__main__':
    unittest.main()
