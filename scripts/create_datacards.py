import os
from os import path

import ROOT
import math

from compute_yields import ComputeYields
from tools import Tools

import sys
sys.path.append('../objects')
from samples import signal_samples, data_samples, qcd_samples
from categories import categories
from baseline_selection import selection
from ABCD_regions import ABCD_regions
from qcd_white_list import white_list


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to produce the datacards', add_help=True)
  parser.add_argument('--outdirlabel'     , type=str, dest='outdirlabel'     , help='name of the outdir'                                            , default=None)
  parser.add_argument('--subdirlabel'     , type=str, dest='subdirlabel'     , help='name of the subdir'                                            , default=None)
  #parser.add_argument('--cardlabel'       , type=str, dest='cardlabel'       , help='label of the datacard'                                         , default=None)
  parser.add_argument('--data_label'      , type=str, dest='data_label'      , help='which data samples to consider?'                               , default='V07_18Aug21')
  parser.add_argument('--qcd_label'       , type=str, dest='qcd_label'       , help='which qcd samples to consider?'                                , default='V07_18Aug21')
  parser.add_argument('--signal_label'    , type=str, dest='signal_label'    , help='which signal samples to consider?'                             , default='private')
  parser.add_argument('--selection_label' , type=str, dest='selection_label' , help='apply a baseline selection_label?'                             , default='standard')
  parser.add_argument('--categories_label', type=str, dest='categories_label', help='label of the list of categories'                               , default='standard')
  parser.add_argument('--category_label'  , type=str, dest='category_label'  , help='label of a given category within this list'                    , default=None)
  parser.add_argument('--ABCD_label'      , type=str, dest='ABCD_label'      , help='which ABCD regions?'                                           , default='cos2d_svprob')
  parser.add_argument('--lumi_target'     , type=str, dest='lumi_target'     , help='which luminosity should the yields be normalised to?'          , default='41.599')
  parser.add_argument('--sigma_B'         , type=str, dest='sigma_B'         , help='which value of the B cross section?'                           , default='472.8e9')
  parser.add_argument('--sigma_mult'      , type=str, dest='sigma_mult'      , help='size n*sigma of the window around a given mass'                , default='20')
  parser.add_argument('--weight_hlt'      , type=str, dest='weight_hlt'      , help='name of the branch of hlt weight'                              , default='weight_hlt_A1')
  parser.add_argument('--qcd_white_list ' , type=str, dest='qcd_white_list'  , help='pthat range to consider for qcd samples'                       , default='20to300')
  parser.add_argument('--add_weight_hlt'  ,           dest='add_weight_hlt'  , help='add hlt weight'                           , action='store_true', default=False)
  parser.add_argument('--do_ABCD'         ,           dest='do_ABCD'         , help='compute yields with the ABCD method'      , action='store_true', default=False)
  parser.add_argument('--do_ABCDHybrid'   ,           dest='do_ABCDHybrid'   , help='compute yields with the ABCDHybrid method', action='store_true', default=False)
  parser.add_argument('--do_TF'           ,           dest='do_TF'           , help='compute yields with the TF method'        , action='store_true', default=False)
  parser.add_argument('--do_categories'   ,           dest='do_categories'   , help='compute yields in categories'             , action='store_true', default=False)
  parser.add_argument('--add_Bc'          ,           dest='add_Bc'          , help='add the Bc samples'                       , action='store_true', default=False)
  #parser.add_argument('--submit_batch', dest='submit_batch', help='submit on the batch?', action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):
  '''
    This function determines whether the parser was used by checking the outdirlabel against its default
  '''
  used_parser = (opt.outdirlabel != None)
  if(opt.do_categories == True and opt.categories_label == 'inclusive'):
    raise RuntimeError('The option "--do_categories" is on. Please indicate which categorisation to use (other than "inclusive"')

  return used_parser


def printInfo(opt):
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '                 Running the datacards producer                          '
  print '\n'
  print ' data samples:    {}'.format(opt.data_label)
  print ' signal samples:  {}'.format(opt.signal_label)
  print '\n'
  print ' categorisation:  {}'.format(opt.categories_label)
  print ' category:        {}'.format(opt.category_label)
  print ' selection:       {}'.format(opt.selection_label)
  print '\n'
  print ' outdir label:    {}'.format(opt.outdirlabel)
  print ' subdir label:    {}'.format(opt.subdirlabel)
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '\n'

class DatacardsMaker(Tools):
  def __init__(self, data_files='', signal_files='', qcd_files='', white_list='', baseline_selection='', ABCD_regions='', do_ABCD=True, do_ABCDHybrid=False, do_TF=False, do_categories=True, categories=None, category_label=None, lumi_target=None, sigma_B=None, sigma_mult=None, weight_hlt=None, add_weight_hlt=True, add_Bc=False, outdirlabel='', subdirlabel=''):
    self.tools = Tools()
    self.data_files = data_files
    self.signal_files = signal_files 
    self.qcd_files = qcd_files
    self.white_list = white_list
    self.baseline_selection = baseline_selection
    self.ABCD_regions = ABCD_regions
    self.do_ABCD = do_ABCD 
    self.do_ABCDHybrid = do_ABCDHybrid 
    self.do_TF = do_TF 
    self.do_categories = do_categories
    self.categories = categories
    if do_categories and categories == None:
      raise RuntimeError('Please indicate which categories dictionnary to use')
    self.category_label = category_label
    self.lumi_target = float(lumi_target)
    self.sigma_B = float(sigma_B)
    self.sigma_mult = float(sigma_mult)
    self.weight_hlt = weight_hlt
    self.add_weight_hlt = add_weight_hlt
    self.add_Bc = add_Bc
    self.outputdir = './outputs/{}/datacards/{}'.format(outdirlabel, subdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)


  def getWindowList(self):
    masses = []
    resolutions = []
    windows = []
    for signal_file in signal_files:
      window = {}
      if signal_file.mass not in masses and signal_file.resolution not in resolutions: 
        window['mass'] = signal_file.mass
        window['resolution'] = signal_file.resolution
        masses.append(signal_file.mass)
        resolutions.append(signal_file.resolution)
        windows.append(window)
    return windows


  def getBackgroundYields(self, mass, resolution, category=''):
    background_selection = self.baseline_selection
    if self.do_categories:
      category_selection = category.definition_flat + ' && ' + category.cutbased_selection
      background_selection += ' && {}'.format(category_selection)

    if self.do_ABCD:
      background_yields = ComputeYields(data_files=self.data_files, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions, sigma_mult_window=self.sigma_mult)[0] 
    elif self.do_ABCDHybrid:
      background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions, sigma_mult_window=self.sigma_mult)[0]
    elif self.do_TF:
      background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromMC(mass=mass, resolution=resolution, sigma_mult_window=self.sigma_mult)[0]

    if background_yields == 0.: background_yields = 1e-9

    lumi_true = 0. #0.774 #data_file.lumi
    for data_file in self.data_files:
      lumi_true += data_file.lumi
    if background_yields != 1e-9: background_yields = background_yields * self.lumi_target/lumi_true

    return background_yields


  def getSignalYields(self, signal_file, category=''):
    signal_selection = 'ismatched==1 && hnl_charge==0 && {}'.format(self.baseline_selection) # condition on the charge added in the context of the dsa study
    if self.do_categories:
      category_selection = category.definition_flat + ' && ' + category.cutbased_selection
      signal_selection += ' && {}'.format(category_selection)

    signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=self.lumi_target, sigma_B=self.sigma_B, add_weight_hlt=self.add_weight_hlt, weight_hlt=self.weight_hlt, isBc=False)[0]
    if self.add_Bc and signal_file.filename_Bc != '':
      signal_yields += ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=self.lumi_target, sigma_B=self.sigma_B, add_weight_hlt=self.add_weight_hlt, weight_hlt=self.weight_hlt, isBc=True)[0]

      #print 'yields {} + {} = {}'.format(ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=False)[0], ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=True), signal_yields)[0]

    return signal_yields


  def getSignalMassCoupling(self, signal_point):
    signal_mass = signal_point.mass
    signal_ctau = signal_point.ctau
    signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
    signal_coupling = self.getCouplingLabel(signal_v2)

    return signal_mass, signal_coupling


  def getLabel(self, signal_mass='', signal_coupling='', category=''):
    if not self.do_categories:
      label = 'bhnl_incl_m_{}_v2_{}'.format(signal_mass, signal_coupling)
    else:
      label = 'bhnl_cat_{}_m_{}_v2_{}'.format(category.label, signal_mass, signal_coupling)
    return label


  def writeCard(self, label, signal_yields, background_yields):
    datacard_name = 'datacard_{}.txt'.format(label)
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
process                                                  sig                            qcd              
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
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<5.5'
    #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'
    #baseline_selection = 'trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))'
    # for bkg TF method
    #baseline_selection = 'mu_isdsa!=1 && b_mass<6.4 && hnl_cos2d>0.993 && sv_prob>0.05 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'

    #baseline_selection = 'hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))' # used for sensitivity
    #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'
    #baseline_selection = 'hnl_charge==0 && b_mass<6.4'

    for category in categories:
      if self.do_categories and 'incl' in category.label: continue
      if not self.do_categories and 'incl' not in category.label: continue

      if self.category_label != None and category.label != self.category_label: continue # needed for category parallelisation on the batch

      # loop on the different mass windows
      for window in self.getWindowList():

        # get the background yields
        background_yields = self.getBackgroundYields(mass=window['mass'], resolution=window['resolution'], category=category)

        # loop on the signal points
        for signal_file in signal_files:
          if signal_file.mass != window['mass']: continue

          # get the signal mass/coupling
          signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)

          # get the process label
          label = self.getLabel(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)

          # get the signal yields
          signal_yields = self.getSignalYields(signal_file=signal_file, category=category)

          # create the datacard
          self.writeCard(label=label, signal_yields=signal_yields, background_yields=background_yields)
          


if __name__ == '__main__':

  opt = getOptions()
  used_parser = checkParser(opt)
  if used_parser:
    '''
      Automatic handling 
    '''
    printInfo(opt)

    data_files = data_samples[opt.data_label]
    signal_files = signal_samples[opt.signal_label]
    qcd_files = qcd_samples[opt.qcd_label]
    
    white_list = white_list[opt.qcd_white_list]

    outdirlabel = opt.outdirlabel
    subdirlabel = opt.subdirlabel

    categories = categories[opt.categories_label]
    baseline_selection = selection[opt.selection_label].flat

    do_ABCD = opt.do_ABCD
    do_ABCDHybrid = opt.do_ABCDHybrid
    do_TF = opt.do_TF
    ABCD_regions = ABCD_regions[opt.ABCD_label]

    lumi_target = opt.lumi_target
    sigma_B = opt.sigma_B

    sigma_mult = opt.sigma_mult

    add_weight_hlt = opt.add_weight_hlt
    weight_hlt = opt.weight_hlt

    do_categories = opt.do_categories
    add_Bc = opt.add_Bc

    DatacardsMaker(
        data_files = data_files, 
        signal_files = signal_files, 
        qcd_files = qcd_files,
        white_list = white_list,
        baseline_selection = baseline_selection, 
        do_categories = do_categories, 
        categories = categories, 
        category_label = opt.category_label, 
        ABCD_regions = ABCD_regions,
        do_ABCD = do_ABCD,
        do_ABCDHybrid = do_ABCDHybrid,
        do_TF = do_TF,
        lumi_target = lumi_target,
        sigma_B = sigma_B,
        sigma_mult = sigma_mult,
        weight_hlt = weight_hlt,
        add_weight_hlt = add_weight_hlt,
        add_Bc = add_Bc, 
        outdirlabel = outdirlabel,
        subdirlabel = subdirlabel,
    ).process() 

  else:
    '''
       ## Interactive board ##

    For an interactive usage of the script, use the code below
    '''
    #data_file = data_samples['V07_18Aug21'][0]
    data_file = data_samples['V08_29Sep21'][0]

    categories = categories['13Oct21_slimmed']
    baseline_selection = selection['13Oct21'].flat
    #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.7'
    #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_pi_mass<4.5 '
    #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_pi_mass<4.5'
    #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.7'
    #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'
    #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'

    # used for sensitivity
    #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2 '
    #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
    #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
    #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
    #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
    #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'

    # used for first dsa study
    #categories['lxy0to1_OS_dsa'] = 'mu_isdsa==1 && sv_lxy<1 && trgmu_charge!=mu_charge '
    #categories['lxy1to5_OS_dsa'] = 'mu_isdsa==1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge'
    #categories['lxygt5_OS_dsa'] = 'mu_isdsa==1 && sv_lxy>=5 && trgmu_charge!=mu_charge'
    #categories['lxy0to1_SS_dsa'] = 'mu_isdsa==1 && sv_lxy<1 && trgmu_charge==mu_charge'
    #categories['lxy1to5_SS_dsa'] = 'mu_isdsa==1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge'
    #categories['lxygt5_SS_dsa'] = 'mu_isdsa==1 && sv_lxy>=5 && trgmu_charge==mu_charge'
    #categories['lxy0to1_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2 '
    #categories['lxy1to5_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
    #categories['lxygt5_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
    #categories['lxy0to1_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
    #categories['lxy1to5_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
    #categories['lxygt5_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'

    #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.7'
    #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && pi_pt>1'
    #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.5 && pi_pt>2'
    #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.7'
    #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && pi_pt>1'
    #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>1.5 && pi_pt>2'
    #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge'
    #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge '
    #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge'
    #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge'
    #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge'
    #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge'


    #dirlabel = 'test_categories_selection_29Jun21_fulllumi_muonid'
    #dirlabel = 'categories_selection_19Aug21_fulllumi_combineddsa'
    dirlabel = 'test'

    ABCD_regions = ABCD_regions['cos2d_svprob'] 
    signal_files = signal_samples['limits_m1_29Sep21']
    datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, baseline_selection=baseline_selection, ABCD_regions=ABCD_regions, do_categories=True, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 

    #signal_files = signal_samples['limits_m3']
    #datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, baseline_selection=baseline_selection, do_categories=True, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 

    #signal_files = signal_samples['limits_m4p5']
    #datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, baseline_selection=baseline_selection, do_categories=True, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 
