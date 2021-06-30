import os
from os import path

from computeYields import ComputeYields
from common import PlottingTools
from samples import data_samples, signal_samples, signal_samples_limits

import sys
sys.path.append('../../../../BHNLGen/CMSSW_10_2_15/src/HNLsGen/python/.')
from my_common import getVV

class DatacardsMaker(object):
  def __init__(self, data_file='', signal_file='', do_categories=True, categories=None, outdirlabel=''): # signal_files instead?
    self.data_file = data_file
    self.signal_file = signal_file 
    self.do_categories = do_categories
    self.categories = categories
    if do_categories and categories == None:
      raise RuntimeError('Please fill the dictionnary with the different categories')
    self.outputdir = './datacards/{}'.format(outdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2) #.replace('.', 'p').replace('-', 'm')


  def writeCard(self, signal_mass='', signal_coupling='', category=''):
    if not self.do_categories:
      label = 'bhnl_m_{}_v2_{}'.format(signal_mass, signal_coupling)
    else:
      label = 'bhnl_cat_{}_m_{}_v2_{}'.format(category, signal_mass, signal_coupling)
    datacard_name = 'datacard_{}.txt'.format(label)

    baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<5.5'

    signal_selection = 'ismatched==1 && {}'.format(baseline_selection)
    if self.do_categories:
      category_selection = categories[category]
      signal_selection += ' && {}'.format(category_selection)

    signal_yields = ComputeYields(signal_file=self.signal_file, selection=signal_selection).computeSignalYields() 
    #print signal_yields

    background_selection = baseline_selection
    if self.do_categories:
      category_selection = categories[category]
      background_selection += ' && {}'.format(category_selection)

    background_yields = ComputeYields(data_file=self.data_file, signal_file=self.signal_file, selection=background_selection).computeBkgYieldsFromABCDData()[0] 
    if background_yields == 0.: background_yields = 1e-9
    #print background_yields

    datacard = open('{}/{}'.format(self.outputdir, datacard_name), 'w+')
    datacard.write(
'''\
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
            lbl = label,
            obs =  -1, # for the moment, we only look at blinded data
            sig_yields = signal_yields,
            bkg_yields = background_yields,
        )
      )

    datacard.close()

    print '--> {}/{} created'.format(self.outputdir, datacard_name)


  def process(self):
    signal_mass = self.signal_file.mass
    signal_ctau = self.signal_file.ctau
    signal_v2 = getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
    signal_coupling = self.getCouplingLabel(signal_v2)

    if not self.do_categories:
      self.writeCard(signal_mass=signal_mass, signal_coupling=signal_coupling)
    else:
      for category in categories:
        self.writeCard(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)




if __name__ == '__main__':

  data_file = data_samples[0]
  signal_V20emu = signal_samples[1] 
  signal_V21 = signal_samples[0] 
  signal_V26 = signal_samples[2] 

  categories = {}
  categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.7'
  categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_pi_mass<4.5 '
  categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_pi_mass<4.5'
  categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.7'
  categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'
  categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'
  #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge'
  #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge '
  #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge'
  #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge'
  #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge'
  #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge'
  '''
  disp3_SS: 
  b_mass > 1.5
  deltaR(mu, pi) < 1.5
  deltaR(trgmu, pi) < 1
  pi pt > 2
  dimu invmass > 0.5
  trgmupi invmass < 4.5

  disp3_OS: 
  b_mass > 1.5
  deltaR(mu, pi) < 1.5
  deltaR(trgmu, pi) < 1
  pi pt > 2
  trgmupi invmass < 4.5


  disp2_SS: 
  b_mass > 1.7
  pi pt > 1
  dimu invmass > 0.5
  trgmupi invmass < 4.5

  disp2_OS: 
  b_mass > 1.7
  pi pt > 1
  trgmupi invmass < 4.5


  disp1_SS: 
  b_mass > 2.7

  disp1_OS:
  '''

  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V20emu, do_categories=True, categories=categories, outdirlabel='test_categories').process() 
  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V21, outdirlabel='test').process() 
  #datacards = DatacardsMaker(data_file=data_file, signal_file=signal_V26, outdirlabel='test').process() 

  from samples import signal_samples_limits_m1
  from samples import signal_samples_limits_m3
  from samples import signal_samples_limits_m4p5
  signal_files = signal_samples_limits_m1
  #for signal_file in signal_files:
  #  datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, do_categories=True, categories=categories, outdirlabel='test_categories_selection').process() 

  signal_files = signal_samples_limits_m3
  for signal_file in signal_files:
    datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, do_categories=False, categories=categories, outdirlabel='test_categories_code').process() 

  signal_files = signal_samples_limits_m4p5
  #for signal_file in signal_files:
  #  datacards = DatacardsMaker(data_file=data_file, signal_file=signal_file, do_categories=True, categories=categories, outdirlabel='test_categories_selection').process() 
