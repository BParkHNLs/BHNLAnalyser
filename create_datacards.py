import os
from os import path

from computeYields import ComputeYields
from common import PlottingTools
from samples import data_samples, signal_samples, signal_samples_limits

import sys
sys.path.append('../../../../BHNLGen/CMSSW_10_2_15/src/HNLsGen/python/.')
from my_common import getVV

class DatacardsMaker(object):
  def __init__(self, data_file='', signal_file='', outdirlabel=''): # signal_files instead?
    self.data_file = data_file
    self.signal_file = signal_file 
    self.outputdir = './datacards/{}'.format(outdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))


  def getCouplingLabel(self, v2): # move to common?
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2) #.replace('.', 'p').replace('-', 'm')


  def process(self):
    signal_mass = self.signal_file.mass
    signal_ctau = self.signal_file.ctau
    signal_v2 = getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

    #outputdir = PlottingTools().getOutDir('./datacards/{}'.format(self.outdirlabel)

    #label = 'bhnl_m_{}_ctau_{}_v2_{}'.format(signal_mass, signal_ctau, self.getCouplingLabel(signal_v2))
    label = 'bhnl_m_{}_v2_{}'.format(signal_mass, self.getCouplingLabel(signal_v2))
    datacard_name = 'datacard_{}.txt'.format(label)
    #print datacard_name

    signal_yields = ComputeYields(signal_file=self.signal_file, selection='ismatched==1').computeSignalYields() 
    #print signal_yields

    baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<5.5'
    background_yields = ComputeYields(data_file=self.data_file, signal_file=self.signal_file, selection=baseline_selection).computeBkgYieldsFromABCDData()[0] 
    #print background_yields

    datacard = open('{}/{}'.format(self.outputdir, datacard_name), 'w+')
    datacard.write(
'''\
# This is the datacard for signal with mass {m} ctau {ct} and v2 {v2}

imax 1 number of bins
jmax 1 number of backgrounds
kmax * number of nuisance parameters
--------------------------------------------------------------------------------------------------------------------------------------------
bin               {lbl}
observation       {obs}
--------------------------------------------------------------------------------------------------------------------------------------------
bin                                                      {lbl}                          {lbl}             
process                                                  signal                         qcd              
process                                                  -1                             1                   
rate                                                     {sig_yields}                   {bkg_yields}    
--------------------------------------------------------------------------------------------------------------------------------------------
lumi                                       lnN           1.025                          -   
syst_sig_{lbl}                             lnN           1.3                            -    
syst_bkg_{lbl}                             lnN           -                              1.3   
--------------------------------------------------------------------------------------------------------------------------------------------
{lbl} autoMCStats 0 0 1
'''.format(
              m = signal_mass,
              ct = signal_ctau, 
              v2 = self.getCouplingLabel(signal_v2),
              lbl = label,
              obs =  -1, # for the moment, we only look at blinded data
              sig_yields = signal_yields,
              bkg_yields = background_yields,
          )
        )

    datacard.close()

    print '--> {}/{} created'.format(self.outputdir, datacard_name)






if __name__ == '__main__':

  data_file = data_samples[0]
  signal_V20emu = signal_samples[1] 
  signal_V21 = signal_samples[0] 
  signal_V26 = signal_samples[2] 

  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V20emu, outdirlabel='test').process() 
  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V21, outdirlabel='test').process() 
  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V26, outdirlabel='test').process() 

  from samples import signal_samples_limits_m1
  from samples import signal_samples_limits_m3
  from samples import signal_samples_limits_m4p5
  signal_files = signal_samples_limits_m1
  for signal_file in signal_files:
    datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, outdirlabel='test_bkg').process() 

  signal_files = signal_samples_limits_m3
  for signal_file in signal_files:
    datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, outdirlabel='test_bkg').process() 

  signal_files = signal_samples_limits_m4p5
  for signal_file in signal_files:
    datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, outdirlabel='test_bkg').process() 
