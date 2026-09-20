"""Independent checks of source-linked zero counts and variant identity."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[2]/'scripts'))
import mdk_source
import reproduce_genomic_plots as gp

class ConfirmedRecords(unittest.TestCase):
 def test_zero_thresholds_and_baselines(self):
  d=mdk_source.read_records('Figure 3').set_index('source_excel_row')
  for row,baseline,hour in [(62,2.1e8,7),(146,2.25e8,3)]:
   r=d.loc[row]
   self.assertTrue(pd.isna(r.Count));self.assertTrue(pd.isna(r.Dilution))
   self.assertEqual(r.count_effective,0);self.assertEqual(r.factor_effective,10)
   self.assertAlmostEqual(r.TimeZero,baseline,delta=1e-6);self.assertEqual(r.Time,hour)
   self.assertAlmostEqual(r.detection_limit,1/(10*.02*baseline),places=20)
   self.assertTrue(r.below_detection);self.assertTrue(pd.isna(r.fraction_observed))
  self.assertEqual(d.below_detection.sum(),2);self.assertEqual(d.count_effective.isna().sum(),0)
  self.assertEqual(d[d.possible_handling_failure].index.tolist(),[146])
 def test_censored_summary_and_omission(self):
  d=mdk_source.read_records('Figure 3');g=d[(d.Name=='P1')&(d.Time==3)]
  positives=g.fraction_observed.dropna().to_numpy()
  for censored in np.linspace(0,1/(10*.02*2.25e8),101):
   values=np.r_[positives,censored]
   self.assertAlmostEqual(np.median(values),4.807692307692308e-7,places=18)
   self.assertAlmostEqual(median_abs_deviation(values,scale='normal'),1.1840226050565571e-7,places=18)
  self.assertAlmostEqual(np.median(positives),5e-7,places=18)
  self.assertAlmostEqual(median_abs_deviation(positives,scale='normal'),8.236678991697786e-8,places=18)
 def test_two_sites_with_same_label_do_not_average(self):
  d=pd.DataFrame({'Pop':[1,1],'label':['cpxA indel']*2,'coordinate':['4104926','4104944'],
                  'mutation_id':['4104926|A','4104944|T'],'last_freq':[.4,0.]})
  with tempfile.TemporaryDirectory() as tmp,patch.object(gp,'OUT',Path(tmp)),patch.object(gp,'CROSSWALK',[]),patch.object(gp,'CHANGES',[]):
   m,labels=gp.matrix(d,'fixture',['cpxA indel'],{})
   self.assertEqual(m.shape,(2,1));self.assertEqual(m[1].tolist(),[.4,0.]);self.assertEqual(len(set(labels)),2)
 def test_s9_source_positions(self):
  d=gp.table('combined_lineages')
  for pop,position,value in [('PLAC_01',4104926,.4),('PLAC_06',3326144,.47),('PLAC_06',3326602,0.),('PLA_07',396591,.66),('PLA_07',396624,.13)]:
   r=d[(d.population==pop)&(d.coordinate==str(position))]
   self.assertEqual(len(r),1);self.assertEqual(r.iloc[0].last_freq,value)

if __name__=='__main__':unittest.main()
