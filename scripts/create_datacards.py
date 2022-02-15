import os
from os import path

import ROOT
import math

from compute_yields import ComputeYields
from tools import Tools
from fitter import Fitter
from quantity import Quantity

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
  parser.add_argument('--outdirlabel'           , type=str, dest='outdirlabel'           , help='name of the outdir'                                            , default=None)
  parser.add_argument('--subdirlabel'           , type=str, dest='subdirlabel'           , help='name of the subdir'                                            , default=None)
  #parser.add_argument('--cardlabel'            , type=str, dest='cardlabel'             , help='label of the datacard'                                         , default=None)
  parser.add_argument('--data_label'            , type=str, dest='data_label'            , help='which data samples to consider?'                               , default='V07_18Aug21')
  parser.add_argument('--qcd_label'             , type=str, dest='qcd_label'             , help='which qcd samples to consider?'                                , default='V07_18Aug21')
  parser.add_argument('--signal_label'          , type=str, dest='signal_label'          , help='which signal samples to consider?'                             , default='private')
  parser.add_argument('--selection_label'       , type=str, dest='selection_label'       , help='apply a baseline selection_label?'                             , default='standard')
  parser.add_argument('--categories_label '     , type=str, dest='categories_label'      , help='label of the list of categories'                               , default='standard')
  parser.add_argument('--category_label'        , type=str, dest='category_label'        , help='label of a given category within this list'                    , default=None)
  parser.add_argument('--ABCD_label'            , type=str, dest='ABCD_label'            , help='which ABCD regions?'                                           , default='cos2d_svprob')
  parser.add_argument('--signal_model_label'    , type=str, dest='signal_model_label'    , help='name of the signal pdf'                                        , default='voigtiant')
  parser.add_argument('--background_model_label', type=str, dest='background_model_label', help='name of the background pdf'                                    , default='chebychev')
  parser.add_argument('--mass_window_size'      , type=str, dest='mass_window_size'      , help='sigma multiplier for the mass window'                          , default='2')
  parser.add_argument('--fit_window_size'       , type=str, dest='fit_window_size'       , help='sigma multiplier for the fit window'                           , default='6')
  parser.add_argument('--nbins'                 , type=str, dest='nbins'                 , help='number of bins when using shapes'                              , default='40')
  parser.add_argument('--lumi_target'           , type=str, dest='lumi_target'           , help='which luminosity should the yields be normalised to?'          , default='41.599')
  parser.add_argument('--sigma_B'               , type=str, dest='sigma_B'               , help='which value of the B cross section?'                           , default='472.8e9')
  parser.add_argument('--sigma_mult'            , type=str, dest='sigma_mult'            , help='size n*sigma of the window around a given mass'                , default='20')
  parser.add_argument('--weight_hlt'            , type=str, dest='weight_hlt'            , help='name of the branch of hlt weight'                              , default='weight_hlt_A1')
  parser.add_argument('--qcd_white_list '       , type=str, dest='qcd_white_list'        , help='pthat range to consider for qcd samples'                       , default='20to300')
  parser.add_argument('--CMStag '               , type=str, dest='CMStag'                , help='CMS tag to be added if --add_CMSlabel'                         , default='Preliminary')
  parser.add_argument('--add_weight_hlt'        ,           dest='add_weight_hlt'        , help='add hlt weight'                           , action='store_true', default=False)
  parser.add_argument('--do_ABCD'               ,           dest='do_ABCD'               , help='compute yields with the ABCD method'      , action='store_true', default=False)
  parser.add_argument('--do_ABCDHybrid'         ,           dest='do_ABCDHybrid'         , help='compute yields with the ABCDHybrid method', action='store_true', default=False)
  parser.add_argument('--do_TF'                 ,           dest='do_TF'                 , help='compute yields with the TF method'        , action='store_true', default=False)
  parser.add_argument('--do_realData'           ,           dest='do_realData'           , help='get the number of data yields'            , action='store_true', default=False)
  parser.add_argument('--do_counting'           ,           dest='do_counting'           , help='perform counting experiment'              , action='store_true', default=False)
  parser.add_argument('--do_shape_analysis'     ,           dest='do_shape_analysis'     , help='perform shape-based analysis'             , action='store_true', default=False)
  parser.add_argument('--do_shape_TH1'          ,           dest='do_shape_TH1'          , help='perform shape-based analysis with histo'  , action='store_true', default=False)
  parser.add_argument('--use_discrete_profiling',           dest='use_discrete_profiling', help='use discrete profiling method'            , action='store_true', default=False)
  parser.add_argument('--do_binned_fit'         ,           dest='do_binned_fit'         , help='perform binned fit when shape analysis'   , action='store_true', default=False)
  parser.add_argument('--do_blind'              ,           dest='do_blind'              , help='perform blind fit when shape analysis'    , action='store_true', default=False)
  parser.add_argument('--plot_pulls'            ,           dest='plot_pulls'            , help='plot pull distribution'                   , action='store_true', default=False)
  parser.add_argument('--do_categories'         ,           dest='do_categories'         , help='compute yields in categories'             , action='store_true', default=False)
  parser.add_argument('--add_Bc'                ,           dest='add_Bc'                , help='add the Bc samples'                       , action='store_true', default=False)
  parser.add_argument('--plot_prefit'           ,           dest='plot_prefit'           , help='produce prefit plots'                     , action='store_true', default=False)
  parser.add_argument('--add_CMSlabel'          ,           dest='add_CMSlabel'          , help='add CMS label'                            , action='store_true', default=False)
  parser.add_argument('--add_lumilabel'         ,           dest='add_lumilabel'         , help='add CMS label'                            , action='store_true', default=False)
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
  def __init__(self, data_files='', signal_files='', qcd_files='', white_list='', baseline_selection='', ABCD_regions='', do_ABCD=True, do_ABCDHybrid=False, do_TF=False, do_realData=False, do_counting=False, do_shape_analysis=False, do_shape_TH1=False, use_discrete_profiling=False, signal_model_label='', background_model_label='', do_binned_fit=True, do_blind=False, mass_window_size='', fit_window_size='', nbins='', plot_pulls=False, do_categories=True, categories=None, category_label=None, lumi_target=None, sigma_B=None, sigma_mult=None, weight_hlt=None, add_weight_hlt=True, add_Bc=False, plot_prefit=False, outdirlabel='', subdirlabel='', add_CMSlabel=True, add_lumilabel=True, CMStag=''):
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
    self.do_realData = do_realData 
    self.do_counting = do_counting
    self.do_shape_analysis = do_shape_analysis
    self.do_shape_TH1 = do_shape_TH1
    self.use_discrete_profiling = use_discrete_profiling
    self.signal_model_label = signal_model_label
    self.background_model_label = background_model_label
    self.do_binned_fit = do_binned_fit
    self.do_blind = do_blind
    self.mass_window_size = int(mass_window_size)
    self.fit_window_size = int(fit_window_size)
    self.nbins = int(nbins)
    self.plot_pulls = plot_pulls
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
    self.plot_prefit = plot_prefit
    self.outputdir = './outputs/{}/datacards/{}'.format(outdirlabel, subdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))
    self.add_CMSlabel = add_CMSlabel
    self.add_lumilabel = add_lumilabel
    self.CMStag = CMStag
  
    ROOT.gROOT.SetBatch(True)


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
      background_yields = 1e4 #ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions, sigma_mult_window=self.sigma_mult)[0]
    elif self.do_TF:
      background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromMC(mass=mass, resolution=resolution, sigma_mult_window=self.sigma_mult)[0]
    elif self.do_realData:
      # get the number of (unblinded) data entries in the 2 sigma window
      quantity = Quantity(name_flat='hnl_mass', nbins=self.nbins, bin_min=mass-2*resolution, bin_max=mass+2*resolution)
      treename = 'signal_tree'
      tree_data = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree_data.Add(data_file.filename) 
      hist_name = 'data_hist'
      data_hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)
      branch_name = 'hnl_mass'
      selection_data = self.baseline_selection 
      if self.do_categories: selection_data += ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
      tree_data.Project(hist_name, branch_name, selection_data)
      background_yields = data_hist.Integral()

    if background_yields == 0.: background_yields = 1e-9

    lumi_true = self.tools.getDataLumi(self.data_files)
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
    signal_coupling = self.tools.getCouplingLabel(signal_v2)

    return signal_mass, signal_coupling


  def getCardLabel(self, signal_mass='', signal_coupling='', category=''):
    if not self.do_categories:
      label = 'bhnl_m_{}_v2_{}_incl'.format(signal_mass, signal_coupling)
    else:
      label = 'bhnl_m_{}_v2_{}_cat_{}'.format(signal_mass, signal_coupling, category.label)
    label = label.replace('.', 'p').replace('-', 'm')
    return label


  def getCategoryLabel(self, signal_mass='', category=''):
    if not self.do_categories:
      label = 'bhnl_m_{}_incl'.format(signal_mass)
    else:
      label = 'bhnl_m_{}_cat_{}'.format(signal_mass, category.label)
    label = label.replace('.', 'p').replace('-', 'm')
    return label


  def runFitter(self, process='', signal_file=None, mass=None, resolution=None, category='', label=''):
    '''
      Run the parametric shapes and extract the yields
    '''
    if process not in ['signal', 'background', 'data_obs']:
      raise RuntimeError("[create_datacards] Unkown process '{}'. Please choose among ['signal', 'background', 'data_obs']")

    # initialise the fitter
    selection = self.baseline_selection
    if self.do_categories:
      selection += ' && ' + category.definition_flat + ' && ' + category.cutbased_selection

    if process == 'signal':
      fitter = Fitter(signal_file=signal_file, data_files=self.data_files, selection=selection, signal_model_label=self.signal_model_label, background_model_label=self.background_model_label, do_blind=self.do_blind, do_binned_fit=self.do_binned_fit, lumi_target=self.lumi_target, sigma_B=self.sigma_B, mass_window_size=self.mass_window_size, fit_window_size=self.fit_window_size, nbins=self.nbins, outputdir=self.outputdir, category_label=category.label, plot_pulls=self.plot_pulls, add_CMSlabel=self.add_CMSlabel, add_lumilabel=self.add_lumilabel, CMStag=self.CMStag)

      # perform the fits and write the workspaces
      fitter.process_signal(label=label)

      # get the signal yields directly from the histogram
      yields = fitter.getSignalYieldsFromHist()

      # produce prefit plots
      if self.plot_prefit:
        fitter.producePrefitPlot(label=label)

    elif process == 'background':
      fitter = Fitter(data_files=self.data_files, mass=mass, resolution=resolution, selection=selection, background_model_label=self.background_model_label, do_blind=self.do_blind, do_binned_fit=self.do_binned_fit, lumi_target=self.lumi_target, mass_window_size=self.mass_window_size, fit_window_size=self.fit_window_size, nbins=self.nbins, outputdir=self.outputdir, category_label=category.label, plot_pulls=self.plot_pulls, add_CMSlabel=self.add_CMSlabel, add_lumilabel=self.add_lumilabel, CMStag=self.CMStag)

      # perform the fits and write the workspaces
      fitter.process_background(label=label)

      # extract the background yields from the fit and normalise to lumi
      background_yields = fitter.getBackgroundYieldsFromFit()
      lumi_true = self.tools.getDataLumi(self.data_files)
      yields = background_yields * self.lumi_target/lumi_true

    elif process == 'data_obs':
      fitter = Fitter(data_files=self.data_files, mass=mass, resolution=resolution, selection=selection, background_model_label=self.background_model_label, do_blind=self.do_blind, do_binned_fit=self.do_binned_fit, lumi_target=self.lumi_target, mass_window_size=self.mass_window_size, fit_window_size=self.fit_window_size, nbins=self.nbins, outputdir=self.outputdir, category_label=category.label, plot_pulls=self.plot_pulls, add_CMSlabel=self.add_CMSlabel, add_lumilabel=self.add_lumilabel, CMStag=self.CMStag)

      # perform the fits and write the workspaces
      fitter.process_data_obs(label=label)

      yields = -1 if self.do_blind else -1 #FIXME

    return yields


  def runFTestRoutine(self, mass, resolution, window_size, category):
    # define selection
    selection = self.baseline_selection
    if self.do_categories:
      selection += ' && ' + category.definition_flat + ' && ' + category.cutbased_selection

    # produce input workspace 
    fitter = Fitter(data_files=self.data_files, mass=mass, resolution=resolution, selection=selection, fit_window_size=self.fit_window_size, nbins=self.nbins, outputdir=self.outputdir, category_label=category.label, lumi_target=self.lumi_target)
    fitter.createFTestInputWorkspace()

    # run the F-test and save the output multipdf in a workspace
    command_ftest = './flashgg_plugin/bin/fTest -i {inws} --saveMultiPdf {outws} -D {outdir} --category_label {cat} --mN {m} --mN_label {ml} --resolution {rsl} --fit_window_size {fws} --mass_window_size {mws} --nbins {nbins}'.format(
        inws = '{}/input_workspace_fTest_m_{}_cat_{}.root'.format(self.outputdir, mass, category.label),
        outws = '{}/workspace_background_multipdf_bhnl_m_{}_cat_{}.root'.format(self.outputdir, mass, category.label),
        outdir = self.outputdir + '/fTest',
        cat = category.label,
        m = mass,
        ml = str(mass).replace('.', 'p'),
        rsl = resolution, 
        fws = self.fit_window_size, 
        mws = self.mass_window_size,
        nbins = self.nbins,
        )
    if self.do_blind: command_ftest += ' --blind'

    os.system(command_ftest)

    background_yields = 1. # the rate in the datacard is set to 1 as the background normalisation is contained in the workspace

    return background_yields


  def createSigHisto(self, category, signal_file, signal_yields, label):
    signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)

    rootfile_name = 'shape_{}.root'.format(label)
    root_file = ROOT.TFile.Open('{}/{}'.format(self.outputdir, rootfile_name), 'RECREATE')  
    root_file.cd()

    sigma = signal_file.resolution
    quantity = Quantity(name_flat='hnl_mass', nbins=self.nbins, bin_min=signal_mass-2*sigma, bin_max=signal_mass+2*sigma)

    treename = 'signal_tree'
    f_signal = self.tools.getRootFile(signal_file.filename)
    tree_signal = self.getTree(f_signal, treename)

    weight_ctau = self.tools.getCtauWeight(signal_file)
    weight_signal = self.tools.getSignalWeight(signal_file=signal_file, sigma_B=self.sigma_B, lumi=self.lumi_target)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)

    signal_selection = 'ismatched==1 && {}'.format(self.baseline_selection)
    if self.do_categories:
      category_selection = category.definition_flat + ' && ' + category.cutbased_selection
      signal_selection += ' && {}'.format(category_selection)
    hist_sig = self.tools.createHisto(tree_signal, quantity, hist_name='sig', branchname='flat', selection=signal_selection, weight=weight_sig)

    root_file.cd() 
    hist_sig.Write()
    root_file.Close()


  def createBkgHisto(self, category, mass, resolution, background_yields, label):
    rootfile_name = 'shape_{}.root'.format(label)
    root_file = ROOT.TFile.Open('{}/{}'.format(self.outputdir, rootfile_name), 'UPDATE')  
    root_file.cd()

    quantity = Quantity(name_flat='hnl_mass', nbins=self.nbins, bin_min=mass-2*resolution, bin_max=mass+2*resolution)

    use_data=True
    if use_data:
      # background shape, taken from data for now
      treename = 'signal_tree'
      tree_bkg = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree_bkg.Add(data_file.filename) 

      background_selection = self.baseline_selection
      if self.do_categories:
        category_selection = category.definition_flat + ' && ' + category.cutbased_selection
        background_selection += ' && {}'.format(category_selection)
      hist_bkg = self.tools.createHisto(tree_bkg, quantity, hist_name='qcd', branchname='flat', selection=background_selection)

      hist_bkg.Scale(background_yields/hist_bkg.Integral())

    else:
      # take shape from qcd mc
      int_mc_tot = 0.
      hist_bkg = ROOT.TH1D('qcd', 'qcd', quantity.nbins, quantity.bin_min, quantity.bin_max)
      hist_bkg.Sumw2()

      background_selection = baseline_selection
      if self.do_categories:
        category_selection = category.definition_flat + ' && ' + category.cutbased_selection
        background_selection += ' && {}'.format(category_selection)

      for ifile, qcd_file in enumerate(self.qcd_files):
        qcd_file_pthatrange = self.tools.getPthatRange(qcd_file.label)
        if qcd_file_pthatrange not in self.white_list: continue

        f_qcd = ROOT.TFile.Open(qcd_file.filename, 'READ')
        tree_qcd = self.getTree(f_qcd, 'signal_tree')
        tree_run = self.getTree(f_qcd, 'run_tree')

        weight_qcd = self.tools.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
        weight_qcd = '({})'.format(weight_qcd)
        hist_qcd = self.tools.createHisto(tree_qcd, quantity, hist_name='hist_qcd', branchname='flat', selection=background_selection, weight=weight_qcd) 

        int_mc_tot += hist_qcd.Integral()

        hist_bkg.Add(hist_qcd)

      hist_bkg.Scale(background_yields/int_mc_tot)

    root_file.cd() 
    hist_bkg.Write()
    root_file.Close()


  def createDataObsHisto(self, category, mass, resolution, label):
    rootfile_name = 'shape_{}.root'.format(label)
    root_file = ROOT.TFile.Open('{}/{}'.format(self.outputdir, rootfile_name), 'UPDATE')  
    root_file.cd()

    quantity = Quantity(name_flat='hnl_mass', nbins=self.nbins, bin_min=mass-2*resolution, bin_max=mass+2*resolution)

    treename = 'signal_tree'
    tree_data = ROOT.TChain(treename)
    for data_file in self.data_files:
      tree_data.Add(data_file.filename) 
    hist_name = 'data_obs'
    data_hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)
    branch_name = 'hnl_mass'
    selection_data = self.baseline_selection 
    if self.do_categories: selection_data += ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
    tree_data.Project(hist_name, branch_name, selection_data)
    data_hist.Scale(self.lumi_target/self.tools.getDataLumi(self.data_files))

    root_file.cd()
    data_hist.Write()
    root_file.Close()

    print '--> {}/{} created'.format(self.outputdir, rootfile_name)


  def writeCard(self, card_label, cat_label, signal_yields, background_yields):
    datacard_name = 'datacard_{}.txt'.format(card_label)
    if self.do_shape_analysis and not self.use_discrete_profiling:
      shape_line = '\n'.join([
          'shapes     sig        {lbl}   workspace_signal_{cardlbl}.root      workspace:sig'.format(lbl=cat_label, cardlbl=card_label),
          'shapes     qcd        {lbl}   workspace_background_{lbl}.root      workspace:qcd'.format(lbl = cat_label),
          'shapes     data_obs   {lbl}   workspace_data_obs_{lbl}.root        workspace:data_obs'.format(lbl=cat_label),
          ])
      norm_line = 'qcd_norm_{lbl}      rateParam   {lbl}   qcd   1.'.format(
          lbl = cat_label,
          )
      index_line = ''
      autostat_line = ''
      #param_line = 'a0{lbl}  param   {val} {err}'.format(
      #    lbl = label, 
      #    val = 0.1,
      #    err = 0.1,
      #    )
      bkg_syst_line = ''
    elif self.do_shape_analysis and self.use_discrete_profiling:
      shape_line = '\n'.join([
          'shapes     sig        {lbl}   workspace_signal_{cardlbl}.root           workspace:sig'.format(lbl=cat_label, cardlbl=card_label),
          'shapes     qcd        {lbl}   workspace_background_multipdf_{lbl}.root  workspace:qcd_multipdf'.format(lbl = cat_label),
          'shapes     data_obs   {lbl}   workspace_data_obs_{lbl}.root             workspace:data_obs'.format(lbl=cat_label),
          ])
      norm_line = ''.format(
          lbl = cat_label,
          )
      index_line = 'pdfindex_{lbl}    discrete'.format(lbl=cat_label)
      autostat_line = ''
      bkg_syst_line = ''
    elif self.do_shape_TH1:
      shape_line = 'shapes *          {lbl}   shape_{lbl}.root   $PROCESS $PROCESS_$SYSTEMATIC'.format(
          lbl = cat_label,
          )
      norm_line = ''
      index_line = ''
      autostat_line = '{lbl} autoMCStats 0 0 1'.format(
          lbl = cat_label,
          )
      bkg_syst_line = 'syst_bkg_{lbl}                             lnN           -                              1.3 '.format(
          lbl = cat_label,
          )
    else:
      shape_line = ''
      norm_line = ''
      index_line = ''
      autostat_line = '{lbl} autoMCStats 0 0 1'.format(
          lbl = cat_label,
          )
      bkg_syst_line = 'syst_bkg_{lbl}                             lnN           -                              1.3 '.format(
          lbl = cat_label,
          )

    datacard = open('{}/{}'.format(self.outputdir, datacard_name), 'w+')
    datacard.write(
'''\
imax 1 number of bins
jmax 1 number of backgrounds
kmax * number of nuisance parameters
--------------------------------------------------------------------------------------------------------------------------------------------
{shape_line} 
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
{bkg_syst_line}   
--------------------------------------------------------------------------------------------------------------------------------------------
{norm_line}
{index_line}
{autostat_line}
'''.format(
            shape_line = shape_line,
            lbl = cat_label,
            obs =  -1, # for the moment, we only look at blinded data #TODO implement not blind option
            sig_yields = signal_yields,
            bkg_yields = background_yields,
            bkg_syst_line = bkg_syst_line,
            norm_line = norm_line,
            index_line = index_line,
            autostat_line = autostat_line,
        )
      )

    datacard.close()
    print '--> {}/{} created'.format(self.outputdir, datacard_name)


  def writeYieldsForPlots(self, label, signal_yields, background_yields):
    summary_name = 'summary_yields_{}.txt'.format(label)
    summary = open('{}/{}'.format(self.outputdir, summary_name), 'w+')
    summary.write(
'''\
sig {sig_yields}
bkg {bkg_yields}
'''.format(
            sig_yields = signal_yields,
            bkg_yields = background_yields,
        )
      )

    summary.close()
    print '--> {}/{} created'.format(self.outputdir, summary_name)


  def process(self):
    for category in categories:
      if self.do_categories and 'incl' in category.label: continue
      if not self.do_categories and 'incl' not in category.label: continue

      if self.category_label != None and category.label != self.category_label: continue # needed for category parallelisation on the batch

      #if category.label != 'lxy1to5_SS' and category.label != 'lxy0to1_SS': continue

      # loop on the different mass windows
      for window in self.getWindowList():
        # get the category label
        cat_label = self.getCategoryLabel(signal_mass=window['mass'], category=category)

        # get the background yields
        if self.do_counting or self.do_shape_TH1:
          background_yields = self.getBackgroundYields(mass=window['mass'], resolution=window['resolution'], category=category)

        elif self.do_shape_analysis and not self.use_discrete_profiling:
          background_yields = self.runFitter(process='background', mass=window['mass'], resolution=window['resolution'], category=category, label=cat_label) 
          data_obs = self.runFitter(process='data_obs', mass=window['mass'], resolution=window['resolution'], category=category, label=cat_label)

        elif self.do_shape_analysis and self.use_discrete_profiling:
          background_yields = self.runFTestRoutine(mass=window['mass'], resolution=window['resolution'], window_size=self.fit_window_size, category=category)
          data_obs = self.runFitter(process='data_obs', mass=window['mass'], resolution=window['resolution'], category=category, label=cat_label)

        # loop on the signal points
        for signal_file in signal_files:
          if signal_file.mass != window['mass']: continue

          #if signal_file.ctau != 0.1: continue # and signal_file.ctau != 10.: continue
          #if signal_file.ctau != 100.: continue # and signal_file.ctau != 10.: continue

          # get the signal mass/coupling
          signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)

          # get the process label
          card_label = self.getCardLabel(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)

          # get the signal yields (if not shape analysis)
          if self.do_counting or self.do_shape_TH1:
            signal_yields = self.getSignalYields(signal_file=signal_file, category=category)
      
          # get the model shape and yields for shape analysis
          if self.do_shape_analysis:
            signal_yields = self.runFitter(process='signal', signal_file=signal_file, category=category, label=card_label)

          # create histograme for non-parametric shape strategy
          if self.do_shape_TH1:
            self.createSigHisto(category=category, signal_file=signal_file, signal_yields=signal_yields, label=card_label)
            self.createBkgHisto(category=category, mass=window['mass'], resolution=window['resolution'], background_yields=background_yields, label=cat_label)
            self.createDataObsHisto(category=category, mass=window['mass'], resolution=window['resolution'], label=cat_label)

          # create the datacard
          self.writeCard(card_label=card_label, cat_label=cat_label, signal_yields=signal_yields, background_yields=background_yields)

          # save yields summary
          self.writeYieldsForPlots(label=card_label, signal_yields=signal_yields, background_yields=background_yields)



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
    do_realData = opt.do_realData
    ABCD_regions = ABCD_regions[opt.ABCD_label]

    signal_model_label = opt.signal_model_label
    background_model_label = opt.background_model_label
    do_binned_fit = opt.do_binned_fit
    mass_window_size = opt.mass_window_size
    fit_window_size = opt.fit_window_size
    nbins = opt.nbins
    do_blind = opt.do_blind
    plot_pulls = opt.plot_pulls

    lumi_target = opt.lumi_target
    sigma_B = opt.sigma_B

    sigma_mult = opt.sigma_mult

    add_weight_hlt = opt.add_weight_hlt
    weight_hlt = opt.weight_hlt

    do_counting = opt.do_counting
    do_shape_analysis = opt.do_shape_analysis
    do_shape_TH1 = opt.do_shape_TH1
    use_discrete_profiling = opt.use_discrete_profiling
    do_categories = opt.do_categories
    add_Bc = opt.add_Bc

    plot_prefit = opt.plot_prefit
    add_CMSlabel = opt.add_CMSlabel
    add_lumilabel = opt.add_lumilabel
    CMStag = opt.CMStag

    DatacardsMaker(
        data_files = data_files, 
        signal_files = signal_files, 
        qcd_files = qcd_files,
        white_list = white_list,
        baseline_selection = baseline_selection, 
        do_counting = do_counting,
        do_shape_analysis = do_shape_analysis,
        do_shape_TH1 = do_shape_TH1,
        use_discrete_profiling = use_discrete_profiling,
        do_categories = do_categories, 
        categories = categories, 
        category_label = opt.category_label, 
        ABCD_regions = ABCD_regions,
        do_ABCD = do_ABCD,
        do_ABCDHybrid = do_ABCDHybrid,
        do_TF = do_TF,
        do_realData = do_realData,
        signal_model_label = signal_model_label,
        background_model_label = background_model_label,
        do_binned_fit = do_binned_fit,
        mass_window_size = mass_window_size,
        fit_window_size = fit_window_size, 
        nbins = nbins,
        do_blind = do_blind,
        plot_pulls = plot_pulls,
        lumi_target = lumi_target,
        sigma_B = sigma_B,
        sigma_mult = sigma_mult,
        weight_hlt = weight_hlt,
        add_weight_hlt = add_weight_hlt,
        add_Bc = add_Bc, 
        plot_prefit = plot_prefit,
        outdirlabel = outdirlabel,
        subdirlabel = subdirlabel,
        add_CMSlabel = add_CMSlabel,
        add_lumilabel = add_lumilabel,
        CMStag = CMStag,
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
