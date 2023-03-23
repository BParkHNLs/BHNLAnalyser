import os
from os import path
import ROOT
from ROOT import RooFit
import matplotlib.pyplot as plt # this library is used to load the correct version of libpng

from tools import Tools
from mva_tools import MVATools
from samples import signal_samples, data_samples


class Fitter(Tools, MVATools):
  def __init__(self, signal_label=None, data_files='', selection='', mass=None, ctau=None, resolution_p0=None, resolution_p1=None, do_cutbased=False, do_mva=False, training_label=None, do_parametric=False, reweighting_strategy=None, signal_model_label=None, background_model_label=None, do_binned_fit=False, do_blind=False, lumi_target=41.6, lhe_efficiency=0.08244, sigma_B=472.8e9, is_bc=False, file_type='flat', mass_window_size='', fit_window_size='',do_veto_SM=None, veto_SM=None, nbins=250, title=' ', outputdir='', outdirlabel='', category_label='', category_title='', plot_pulls=True, add_weight_hlt=False, add_weight_pu=False, add_weight_muid=False, weight_hlt=None, weight_pusig=None, weight_mu0id=None, weight_muid=None, add_CMSlabel=True, add_lumilabel=True, CMStag='', do_tdrstyle=False):
    self.tools = Tools()
    self.mva_tools = MVATools()
    self.signal_label = signal_label
    self.data_files = data_files
    self.selection = selection
    self.signal_mass = mass
    self.signal_ctau = ctau
    self.do_cutbased = do_cutbased
    self.do_mva = do_mva
    self.training_label = training_label
    self.do_parametric = do_parametric
    self.reweighting_strategy = reweighting_strategy
    self.resolution = resolution_p0 + resolution_p1 * self.signal_mass
    self.signal_model_label = signal_model_label
    self.background_model_label = background_model_label
    self.do_binned_fit = do_binned_fit
    self.do_blind = do_blind
    self.lumi_target = lumi_target
    self.lumi_true = self.tools.getDataLumi(self.data_files)
    if self.lumi_true == 42.001: self.lumi_true = 41.6 
    self.sigma_B = sigma_B
    self.lhe_efficiency = lhe_efficiency
    self.is_bc = is_bc
    self.file_type = file_type
    self.do_veto_SM = do_veto_SM
    self.veto_SM = veto_SM
    self.nbins = int(nbins)
    self.title = title
    self.plot_pulls = plot_pulls
    self.add_weight_hlt = add_weight_hlt
    self.add_weight_pu = add_weight_pu
    self.add_weight_muid = add_weight_muid
    self.weight_hlt = weight_hlt
    self.weight_pusig = weight_pusig
    self.weight_mu0id = weight_mu0id
    self.weight_muid = weight_muid
    self.workspacedir = outputdir #TODO adapt
    if self.workspacedir == '':
      self.workspacedir = './'
    self.outputdir = outputdir
    if self.outputdir != '': self.outputdir = self.outputdir + '/fits'
    else: self.outputdir = './myPlots/fits/' + outdirlabel
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))
    self.category_label = category_label
    self.category_title = category_title
    self.add_CMSlabel = add_CMSlabel
    self.add_lumilabel = add_lumilabel
    self.CMStag = CMStag
    self.do_tdrstyle = do_tdrstyle
    # define window sizes (multiples of sigma)
    self.mass_window_size = mass_window_size
    self.fit_window_size = fit_window_size

    self.do_fixed_shape = True #FIXME configure option


    signal_model_list = ['doubleCB', 'doubleCBPlusGaussian', 'voigtian']
    if self.signal_model_label != None and self.signal_model_label not in signal_model_list:
      raise RuntimeError('Unrecognised signal model "{}". Please choose among {}'.format(self.signal_model_label, signal_model_list))

    background_model_list = ['chebychev']
    if self.background_model_label != None and self.background_model_label not in background_model_list:
      raise RuntimeError('Unrecognised background model "{}". Please choose among {}'.format(self.background_model_label, background_model_list))


    #TODO fit the signal in fit_window, extract the sigma from there, and build the fit_window in the final workspace from that? 
    # (maybe not so optimal since would mean to rerun background for each signal. Otherwise save resolutions in a dict in objects?)


  def getRegion(self, nsigma=2):
    bin_min = self.signal_mass - nsigma * self.resolution
    bin_max = self.signal_mass + nsigma * self.resolution
    return bin_min, bin_max


  def getSignalLabel(self):
    label = 'm_{}_ctau_{}_{}'.format(self.signal_mass, self.signal_ctau, self.signal_model_label).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    return label


  def getBackgroundLabel(self):
    label = 'm_{}_{}'.format(self.signal_mass, self.background_model_label).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    if self.do_blind: label += '_blind'
    return label


  def getMCWeight(self, signal_files, lumi, is_bc=False):
    weight_sig_list = ['gen_hnl_ct']
    weight_ctau = self.tools.getCtauWeight(signal_files=signal_files, ctau=self.signal_ctau, is_bc=is_bc)
    weight_signal = self.tools.getSignalWeight(signal_files=signal_files, mass=self.signal_mass, ctau=self.signal_ctau, sigma_B=self.sigma_B, lumi=lumi, lhe_efficiency=self.lhe_efficiency, is_bc=is_bc)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
    if self.add_weight_hlt: 
      weight_sig += ' * ({})'.format(self.weight_hlt)
      weight_sig_list.append(self.weight_hlt)
    if self.add_weight_pu: 
      weight_sig += ' * ({})'.format(self.weight_pusig)
      weight_sig_list.append(self.weight_pusig)
    if self.add_weight_muid: 
      weight_sig += ' * ({}) * ({})'.format(self.weight_mu0id, self.weight_muid)
      weight_sig_list.append(self.weight_mu0id)
      weight_sig_list.append(self.weight_muid)

    return weight_sig, weight_sig_list


  def getQuantitySet(self):
    '''
      If performing unbinned fit, make sure that all variables used in the selection cuts
      are defined here
    '''
    quantities = [
      ROOT.RooRealVar('ismatched', 'ismatched', -2, 2),
      ROOT.RooRealVar('hnl_pt', 'hnl_pt', 0., 13000.),
      ROOT.RooRealVar('sv_lxy', 'sv_lxy', 0., 13000.),
      ROOT.RooRealVar('trgmu_charge', 'trgmu_charge', -2, 2),
      ROOT.RooRealVar('mu_charge', 'mu_charge', -2, 2),
      ROOT.RooRealVar('trgmu_softid', 'trgmu_softid', -2, 2),
      ROOT.RooRealVar('mu_looseid', 'mu_looseid', -2, 2),
      ROOT.RooRealVar('pi_packedcandhashighpurity', 'pi_packedcandhashighpurity', -2, 2),
      ROOT.RooRealVar('trgmu_mu_mass', 'trgmu_mu_mass', 0., 13000.),
      ROOT.RooRealVar('hnl_charge', 'hnl_charge', -2, 2),
      ROOT.RooRealVar('pi_pt', 'pi_pt', 0., 13000.),
      ROOT.RooRealVar('sv_lxysig', 'sv_lxysig', 0., 13000.),
      ROOT.RooRealVar('mu_dxysig', 'mu_dxysig', -13000., 13000.),
      ROOT.RooRealVar('pi_dxysig', 'pi_dxysig', -13000., 13000.),
    ]

    quantity_set = ROOT.RooArgSet()
    for quantity in quantities:
      ROOT.SetOwnership(quantity, False)
      quantity_set.add(quantity)

    return quantity_set
      
    
  def createFitModels(self, process='', label='', do_recreate=True):
    print ' --- Creating Fit Models --- '

    # get the signal region
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size)

    self.hnl_mass = ROOT.RooRealVar("hnl_mass","hnl_mass", bin_min, bin_max)
    #self.hnl_mass.setBins(self.nbins)

    ### Signal Model ###

    if process == 'signal' or process == 'both':
      if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
        # energy scale correction applied
        self.mean_CB  = ROOT.RooRealVar("mean_CB_"+self.category_label,"mean_CB_"+self.category_label, self.signal_mass, self.signal_mass-0.001*self.signal_mass, self.signal_mass+0.001*self.signal_mass)

        if not self.do_fixed_shape:
          # parameters to be kept floating
          self.sigma_CB = ROOT.RooRealVar("sigma_CB_"+self.category_label, "sigma_CB_"+self.category_label, self.resolution, 0.005, 0.15)
          self.alpha_1 = ROOT.RooRealVar("alpha_1_"+self.category_label, "alpha_1_"+self.category_label, -2, -5, 5)
          self.n_1 = ROOT.RooRealVar("n_1_"+self.category_label, "n_1_"+self.category_label, 0, 5)
          self.alpha_2 = ROOT.RooRealVar("alpha_2_"+self.category_label, "alpha_2_"+self.category_label, 2, -5, 5)
          self.n_2 = ROOT.RooRealVar("n_2_"+self.category_label, "n_2_"+self.category_label, 0, 5)
          self.sigfrac_CB = ROOT.RooRealVar("sigfrac_CB_"+self.category_label,"sigfrac_CB_"+self.category_label, 0.5, 0.0 ,1.0)
        else:
          # parameters fixed
          self.sigma_CB = ROOT.RooRealVar("sigma_CB_"+self.category_label, "sigma_CB_"+self.category_label, self.resolution)
          self.alpha_1 = ROOT.RooRealVar("alpha_1_"+self.category_label, "alpha_1_"+self.category_label, -1.)
          self.n_1 = ROOT.RooRealVar("n_1_"+self.category_label, "n_1_"+self.category_label, 4.)
          self.alpha_2 = ROOT.RooRealVar("alpha_2_"+self.category_label, "alpha_2_"+self.category_label, 1.)
          self.n_2 = ROOT.RooRealVar("n_2_"+self.category_label, "n_2_"+self.category_label, 4.)
          # defines the relative importance of the two CBs
          self.sigfrac_CB = ROOT.RooRealVar("sigfrac_CB_"+self.category_label,"sigfrac_CB_"+self.category_label, 0.5)

        self.CBpdf_1 = ROOT.RooCBShape("CBpdf_1_"+self.category_label, "CBpdf_1_"+self.category_label, self.hnl_mass, self.mean_CB, self.sigma_CB, self.alpha_1, self.n_1)
        self.CBpdf_2 = ROOT.RooCBShape("CBpdf_2_"+self.category_label, "CBpdf_2_"+self.category_label, self.hnl_mass, self.mean_CB, self.sigma_CB, self.alpha_2, self.n_2)

        if self.signal_model_label == 'doubleCB':
          self.signal_model = ROOT.RooAddPdf("sig", "sig", self.CBpdf_1, self.CBpdf_2, self.sigfrac_CB)

        if self.signal_model_label == 'doubleCBPlusGaussian':
          self.sigma_gauss = ROOT.RooRealVar("sigma_gauss", "sigma_gauss", 0.01, 0.005, 0.15)
          self.gaussian = ROOT.RooGaussian('gaussian', 'gaussian', self.hnl_mass, self.mean_CB, self.sigma_gauss)

          # defines the relative importance of gaussian wrt doubleCB
          self.sigfrac_gauss = ROOT.RooRealVar("sigfrac_gauss","sigfrac_gauss", 0.5, 0.0 ,1.0)

          self.doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", self.CBpdf_1, self.CBpdf_2, self.sigfrac_CB)
          self.signal_model = ROOT.RooAddPdf('sig', 'sig', self.doubleCBpdf, self.gaussian, self.sigfrac_gauss) # make sure that the model has the same name as in datacard 

      elif self.signal_model_label == 'voigtian':
        self.mean_voigtian  = ROOT.RooRealVar("mean_voigtian","mean_voigtian", self.signal_mass, self.signal_mass-0.001*self.signal_mass, self.signal_mass+0.001*self.signal_mass)
        if not self.do_fixed_shape:
          self.gamma_voigtian = ROOT.RooRealVar("gamma_voigtian_"+self.category_label, "gamma_voigtian_"+self.category_label, 0.01, 0., 5.)
          self.sigma_voigtian = ROOT.RooRealVar("sigma_voigtian_"+self.category_label, "sigma_voigtian_"+self.category_label, 0.01, 0.005, 0.15)
        else:
          self.gamma_voigtian = ROOT.RooRealVar("gamma_voigtian_"+self.category_label, "gamma_voigtian_"+self.category_label, 0.03)
          self.sigma_voigtian = ROOT.RooRealVar("sigma_voigtian_"+self.category_label, "sigma_voigtian_"+self.category_label, self.resolution)

        self.signal_model = ROOT.RooVoigtian('sig', 'sig', self.hnl_mass, self.mean_voigtian, self.gamma_voigtian, self.sigma_voigtian)

    ### Background Model ###

    if process == 'background' or process == 'both':
      if self.background_model_label == 'chebychev':
        self.n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', 100, 0, 100000)
        self.a0 = ROOT.RooRealVar('a0', 'a0', 0.01, -10, 10)
        self.chebychev = ROOT.RooChebychev('chebychev', 'chebychev', self.hnl_mass, ROOT.RooArgList(self.a0))
        self.background_model = ROOT.RooAddPdf('qcd', 'qcd', ROOT.RooArgList(self.chebychev), ROOT.RooArgList(self.n_bkg))
        #a0 = ROOT.RooRealVar('a0', 'a0', -1.36, -1.52, -1.20)
        #a1 = ROOT.RooRealVar('a1', 'a1', 0.53, 0.29, 0.76)
        #a2 = ROOT.RooRealVar('a2', 'a2', -0.14, -0.32, 0.04)
        #background_model = ROOT.RooChebychev('qcd', 'qcd', hnl_mass, ROOT.RooArgList(a0, a1, a2))

    print '--> Models created'


  def performFit(self, process, label=''):
    '''
      process corresponds to either 'signal' or 'background'
    '''
    print ' --- Running the fits --- '

    if process not in ['signal', 'background', 'both']:
      raise RuntimeError("[fitter] Unrecognised process. Please choose among ['signal', 'background', 'both']")

    # enforce blinding for prefit plots
    if process == 'both': self.do_blind = True

    # get the label
    if label == '' and process == 'signal': label = self.getSignalLabel() 
    if label == '' and process == 'background': label = self.getBackgroundLabel() 
    if label == '' and process == 'both': label = self.getSignalLabel() + '_' + self.getBackgroundLabel()

    # define ranges and binning
    mass_window_min, mass_window_max = self.getRegion(nsigma=self.mass_window_size)
    fit_window_min, fit_window_max = self.getRegion(nsigma=self.fit_window_size)

    # define selection
    cond_sig = 'ismatched==1' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1'
    selection_sig = cond_sig + ' && ' + self.selection
    selection_bkg = self.selection
    #TODO do we want to apply cut on the hnl fit window (10sigma)?
    if self.do_blind:
      selection_bkg += ' && (hnl_mass < {} || hnl_mass > {})'.format(mass_window_min, mass_window_max)
    if self.do_veto_SM:
      selection_bkg += ' && (hnl_mass < {} || hnl_mass > {})'.format(self.veto_SM.range_min, self.veto_SM.range_max)
    print 'sel bkg before: {}'.format(selection_bkg)
    print 'sel sig before: {}'.format(selection_sig)

    # open the file and get the tree
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    if process == 'signal' or process == 'both':
      signal_files = self.tools.getSignalFileList(signal_label=self.signal_label, mass=self.signal_mass, ctau=self.signal_ctau, strategy=self.reweighting_strategy, is_bc=self.is_bc)
      #print '\n\n'
      #print 'list of signal files'
      #for signal_file in signal_files:
      #  print signal_file.filename
      #print '\n\n'

      # define signal weights
      weight_sig, weight_sig_list = self.getMCWeight(signal_files, lumi=self.lumi_true, is_bc=self.is_bc) 
      #TODO do we want to have the prefit plots scaled to the target lumi? In which case, also scale the background

      if self.do_cutbased:
        tree_sig = ROOT.TChain(treename)
        for signal_file in signal_files:
          signal_filename = signal_file.filename if not self.is_bc else signal_file.filename_Bc
          tree_sig.Add(signal_filename)
      elif self.do_mva:
        score_label = self.getSignalLabel()
        filename_sig = self.mva_tools.getFileWithScore(files=signal_files, training_label=self.training_label, do_parametric=self.do_parametric, mass=self.signal_mass, category_label=self.category_label, selection=selection_sig, weights=weight_sig_list, label=score_label, treename=treename, is_bc=self.is_bc) 
        file_sig = self.tools.getRootFile(filename_sig)
        tree_sig = self.tools.getTree(file_sig, treename)

    if process == 'background' or process == 'both':
      if self.do_cutbased:
        tree_data = ROOT.TChain(treename)
        for data_file in self.data_files:
          tree_data.Add(data_file.filename) 
      elif self.do_mva:
        score_label = self.getBackgroundLabel()
        filename_data = self.mva_tools.getFileWithScore(files=self.data_files, training_label=self.training_label, do_parametric=self.do_parametric, mass=self.signal_mass, category_label=self.category_label, selection=selection_bkg, label=score_label, treename=treename) 
        file_data = self.tools.getRootFile(filename_data)
        tree_data = self.tools.getTree(file_data, treename)

    if self.do_binned_fit:
      # build the binned data
      hist_name = 'hist'
      hist = ROOT.TH1D(hist_name, hist_name, self.nbins, fit_window_min, fit_window_max)
      branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'

      if process == 'signal' or process == 'both':
        tree_sig.Project(hist_name, branch_name , '({sel}) * ({wght})'.format(sel=selection_sig if self.do_cutbased else self.mva_tools.getScoreSelection(selection_sig), wght=weight_sig))
        print 'sel sig: ({sel}) * ({wght})'.format(sel=selection_sig if self.do_cutbased else self.mva_tools.getScoreSelection(selection_sig), wght=weight_sig)
        rdh_sig = ROOT.RooDataHist("rdh_sig", "rdh_sig", ROOT.RooArgList(self.hnl_mass), hist)

      if process == 'background' or process == 'both':
        tree_data.Project(hist_name, branch_name , selection_bkg if self.do_cutbased else self.mva_tools.getScoreSelection(selection_bkg))
        print 'sel bkg: {}'.format(selection_bkg if self.do_cutbased else self.mva_tools.getScoreSelection(selection_bkg))
        #hist.Scale(self.lumi_target/self.lumi_true) #FIXME added in the context of prefit plots for the sps, remove?
        rdh_bkg = ROOT.RooDataHist("rdh_bkg", "rdh_bkg", ROOT.RooArgList(self.hnl_mass), hist)
    else:
      # get RooArgSet
      quantity_set = self.getQuantitySet()
      # add hnl_mass to the RooArgSet
      ROOT.SetOwnership(self.hnl_mass, False)
      quantity_set.add(self.hnl_mass)
     
      #TODO correct to have only selection and not selection_sig/bkg below?
      if process == 'signal' or process == 'both':
        print '-> creating unbinned dataset'
        #rds_sig = ROOT.RooDataSet('rds_sig', 'rds_sig', tree_sig, quantity_set, selection)
        rds_sig = ROOT.RooDataSet('rds_sig', 'rds_sig', tree_sig, quantity_set, selection_sig)
        print '-> unbinned dataset created'
      if process == 'background' or process == 'both':
        print '-> creating unbinned dataset'
        #rds_bkg = ROOT.RooDataSet('rds_bkg', 'rds_bkg', tree_bkg, quantity_set, selection)
        rds_bkg = ROOT.RooDataSet('rds_bkg', 'rds_bkg', tree_bkg, quantity_set, selection_bkg)
        print '-> unbinned dataset created'

    # create canvas
    #canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)
    canv = self.tools.createTCanvas(name="canv", dimx=700, dimy=600)

    # and the two pads
    if process == 'both': # remove residuals for prefit plots
      pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.01, 0.99, 0.99)
    else:
      pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    if process != 'both': pad2.Draw()

    frame = self.hnl_mass.frame(ROOT.RooFit.Title(self.title))

    # define fit ranges
    self.hnl_mass.setRange("peak", mass_window_min, mass_window_max)
    self.hnl_mass.setRange('full', fit_window_min, fit_window_max)
    self.hnl_mass.setRange('sideband_left', fit_window_min, mass_window_min)
    self.hnl_mass.setRange('sideband_right', mass_window_max, fit_window_max)

    if process != 'both':
      # plot the data
      if self.do_binned_fit:
        if process == 'signal':
          rdh_sig.plotOn(frame, ROOT.RooFit.Name("data_sig"))
        if process == 'background':
          rdh_bkg.plotOn(frame, ROOT.RooFit.Name("data_bkg"))
      else:
        if process == 'signal':
          rds_sig.plotOn(frame, ROOT.RooFit.Name("data_sig"), ROOT.RooFit.Binning(self.nbins))
        if process == 'background':
          rds_bkg.plotOn(frame, ROOT.RooFit.Name("data_bkg"), ROOT.RooFit.Binning(self.nbins))

      if process == 'signal': # or (process == 'background' and not self.do_blind):
        fit_range = 'peak'
      elif process == 'background' and not self.do_blind:
        fit_range = 'full'
      else:
        fit_range = 'sideband_left,sideband_right'

      # fit the PDF to the data
      if process == 'signal':
        if self.do_binned_fit:
          result_sig = self.signal_model.fitTo(rdh_sig, ROOT.RooFit.Range(fit_range), ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save())
          self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kOrange+7), ROOT.RooFit.Name('sig'), ROOT.RooFit.Components('sig'))
        else:
          result_sig = self.signal_model.fitTo(rds_sig, ROOT.RooFit.Range(fit_range))
      else:
        if self.do_binned_fit:
          result_bkg = self.background_model.fitTo(rdh_bkg, ROOT.RooFit.Range(fit_range))
        else:
          result_bkg = self.background_model.fitTo(rds_bkg, ROOT.RooFit.Range(fit_range))

      # plot the fit 		
      pdf_name = 'sig' if process == 'signal' else 'qcd'
      if process == 'signal':
        if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
          self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1_"+self.category_label),ROOT.RooFit.Components("CBpdf_1_"+self.category_label), ROOT.RooFit.LineStyle(ROOT.kDashed))
          self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2_"+self.category_label),ROOT.RooFit.Components("CBpdf_2_"+self.category_label), ROOT.RooFit.LineStyle(ROOT.kDashed))
          if self.signal_model_label == 'doubleCBPlusGaussian':
            self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))
      else:
        self.background_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))


    else: # fit both the signal and background
      # start with the signal
      # plot the data
      if self.do_binned_fit:
        rdh_sig.plotOn(frame, ROOT.RooFit.Name("data_sig"), ROOT.RooFit.DrawOption('B'), ROOT.RooFit.FillColor(ROOT.kOrange+7), ROOT.RooFit.FillStyle(3003), ROOT.RooFit.DataError(ROOT.RooAbsData.None), ROOT.RooFit.XErrorSize(0))
        #rdh_sig.plotOn(frame, ROOT.RooFit.Name("data_sig"))
      else:
        rds_sig.plotOn(frame, ROOT.RooFit.Name("data_sig"), ROOT.RooFit.Binning(self.nbins), ROOT.RooFit.DrawOption('B'), ROOT.RooFit.FillColor(ROOT.kOrange+7), ROOT.RooFit.FillStyle(3003), ROOT.RooFit.DataError(ROOT.RooAbsData.None), ROOT.RooFit.XErrorSize(0))

      # fit the signal  
      if self.do_binned_fit:
        result_sig = self.signal_model.fitTo(rdh_sig, ROOT.RooFit.Range('peak'), ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save())
      else:
        result_sig = self.signal_model.fitTo(rds_sig, ROOT.RooFit.Range('peak'))

      # plot the fit
      if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
        self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1_"+self.category_label),ROOT.RooFit.Components("CBpdf_1_"+self.category_label), ROOT.RooFit.LineStyle(ROOT.kDashed))
        self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2_"+self.category_label),ROOT.RooFit.Components("CBpdf_2_"+self.category_label), ROOT.RooFit.LineStyle(ROOT.kDashed))
        if self.signal_model_label == 'doubleCBPlusGaussian':
          self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
      self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kOrange+7), ROOT.RooFit.Name('sig'), ROOT.RooFit.Components('sig'))

      # then the background
      # plot the data
      if self.do_binned_fit:
        rdh_bkg.plotOn(frame, ROOT.RooFit.Name("data_bkg"))
      else:
        rds_bkg.plotOn(frame, ROOT.RooFit.Name("data_bkg"), ROOT.RooFit.Binning(self.nbins))

      # fit the background
      if not self.do_blind: fit_range = 'full'
      else: fit_range = 'sideband_left,sideband_right'
      if self.do_binned_fit:
        result_bkg = self.background_model.fitTo(rdh_bkg, ROOT.RooFit.Range(fit_range))
      else:
        result_bkg = self.background_model.fitTo(rds_bkg, ROOT.RooFit.Range(fit_range))

      # plot the fit
      self.background_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name('qcd'), ROOT.RooFit.Components('qcd'))

    # and write the fit parameters
    if process != 'both':
      if process == 'signal':
        self.signal_model.paramOn(frame,   
             ROOT.RooFit.Layout(0.165, 0.4, 0.8),
             ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
             )
      if process == 'background':
        self.background_model.paramOn(frame,   
             ROOT.RooFit.Layout(0.2, 0.4, 0.8),
             ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
             )
      frame.getAttText().SetTextSize(0.03)
      frame.getAttLine().SetLineColorAlpha(0, 0)
      frame.getAttFill().SetFillColorAlpha(0, 0)
    elif process == 'both' and not self.do_tdrstyle:
      self.signal_model.paramOn(frame,   
           ROOT.RooFit.Layout(0.2, 0.3, 0.8),
           ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
           )
      frame.getAttText().SetTextSize(0.03)
      frame.getAttLine().SetLineColorAlpha(0, 0)
      frame.getAttFill().SetFillColorAlpha(0, 0)

      self.background_model.paramOn(frame,   
           ROOT.RooFit.Layout(0.2, 0.3, 0.65),
           ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
           )
      frame.getAttText().SetTextSize(0.03)
      frame.getAttLine().SetLineColorAlpha(0, 0)
      frame.getAttFill().SetFillColorAlpha(0, 0)

    # compute the chisquare
    if process != 'both':
      if process == 'signal':
        chisquare = frame.chiSquare("sig","data_sig")
      else:
        chisquare = frame.chiSquare("qcd","data_bkg")

    # and print it
    if not self.do_tdrstyle:
      label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    else:
      label1 = ROOT.TPaveText(0.2,0.65,0.3,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    if not self.do_tdrstyle:
      label1.SetTextSize(0.03)
    else:
      label1.SetTextSize(0.045)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    if process == 'signal' or process == 'both':
      if not self.do_tdrstyle:
        label_text = 'mass {}GeV, ctau {}mm'.format(self.signal_mass, self.signal_ctau)
      else:
        label_text_mass = 'mass {} GeV'.format(self.signal_mass)
        label_text_ctau = 'ctau {} mm'.format(self.signal_ctau)
    else: 
      label_text = 'Mass window around {}GeV'.format(self.signal_mass)
    if not self.do_tdrstyle:
      label1.AddText(label_text)
    else:
      label1.AddText(label_text_mass)
      label1.AddText(label_text_ctau)
    if process != 'both':
      qte = '#chi^{2}/ndof'
      label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
      print "chisquare = {}".format(chisquare)

    if not self.do_tdrstyle:
      label2 = ROOT.TPaveText(0.62,0.8,0.72,0.86,"brNDC")
    else:
      label2 = ROOT.TPaveText(0.58,0.65,0.68,0.8,"brNDC")
    label2.SetBorderSize(0)
    label2.SetFillColorAlpha(0, 0)
    if not self.do_tdrstyle:
      label2.SetTextSize(0.04)
    else:
      label2.SetTextSize(0.05)
    label2.SetTextFont(22)
    label2.SetTextAlign(11)
    if not self.do_tdrstyle:
      label2.AddText(self.category_label)
    else:
      label2.AddText(self.category_title)

    # We define and plot the pull 		
    if process != 'both':
      if process == 'background' and self.do_blind:
        curve1 = frame.getObject(1)
        curve2 = frame.getObject(2)
        datahist = frame.getHist('data_bkg')
        hpull1 = datahist.makePullHist(curve1, True)
        hpull2 = datahist.makePullHist(curve2, True)
        frame2 = self.hnl_mass.frame(ROOT.RooFit.Title(" "))
        frame2.addPlotable(hpull1,"P")
        frame2.addPlotable(hpull2,"P")
        hpull = frame.pullHist() # needed to get pull distribution
      else:
        hpull = frame.pullHist()
        for i in range(0, frame.GetNbinsX()):
           hpull.SetPointError(i,0,0,0,0)
        frame2 = self.hnl_mass.frame(ROOT.RooFit.Title(" "))
        frame2.addPlotable(hpull,"P")
    else:
      curve1 = frame.getObject(3)
      curve2 = frame.getObject(4)
      datahist = frame.getHist('data_bkg')
      hpull1 = datahist.makePullHist(curve1, True)
      hpull2 = datahist.makePullHist(curve2, True)
      frame2 = self.hnl_mass.frame(ROOT.RooFit.Title(" "))
      frame2.addPlotable(hpull1,"P")
      frame2.addPlotable(hpull2,"P")
      hpull = frame.pullHist() # needed to get pull distribution

    # plot of the curve and the fit
    canv.cd()
    pad1.cd()

    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitle("#mu#pi invariant mass [GeV]")
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()
    label1.Draw()
    label2.Draw()

    # add the legend
    if process != 'both':
      leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.5, xmax=0.8, ymax=0.63, size=0.03, do_alpha=True)
      if process == 'signal':
        if self.signal_model_label == 'doubleCB':
          model_label = 'Double Crystal Ball'
        elif self.signal_model_label == 'doubleCBPlusGaussian':
          model_label = 'Double Crystal Ball + Gaussian'
        elif self.signal_model_label == 'voigtian':
          model_label = 'Voigtian'
        leg.AddEntry(frame.findObject(pdf_name), model_label)
        if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
          leg.AddEntry(frame.findObject('CBpdf_1_'+self.category_label), 'CB_1')
          leg.AddEntry(frame.findObject('CBpdf_2_'+self.category_label), 'CB_2')
          if self.signal_model_label == 'doubleCBPlusGaussian':
            leg.AddEntry(frame.findObject('gaussian'), 'Gaussian')
      else:
        if self.background_model_label == 'chebychev':
          model_label = 'Chebychev'
        leg.AddEntry(frame.findObject(pdf_name), model_label)
      leg.Draw()

    # add the labels
    if self.add_CMSlabel: self.tools.printCMSTag(pad1, self.CMStag, size=0.5, offset=0.1 if process == 'both' else 0.08)
    if self.add_lumilabel and not self.do_tdrstyle: self.tools.printLumiTag(pad1, self.lumi_true, size=0.5, offset=0.50 if process == 'both' else 0.55)
    if self.add_lumilabel and self.do_tdrstyle: self.tools.printLumiTag(pad1, self.lumi_target, size=0.5, offset=0.48 if process == 'both' else 0.55)

    # plot of the residuals
    if process != 'both':
      pad2.cd()
      ROOT.gPad.SetLeftMargin(0.15) 
      ROOT.gPad.SetPad(0.01,0.01,0.99,0.2)

      frame2.GetYaxis().SetNdivisions(3)
      frame2.GetYaxis().SetLabelSize(0.17)
      frame2.GetYaxis().SetTitleSize(0.17)
      frame2.GetYaxis().SetTitleOffset(0.24)
      frame2.GetYaxis().SetRangeUser(-5,5)	
      frame2.GetYaxis().SetTitle("Pulls")	
      frame2.GetXaxis().SetTitle("")	
      frame2.GetXaxis().SetLabelOffset(5)	
      frame2.Draw()

      bin_min = fit_window_min
      bin_max = fit_window_max
      line = ROOT.TLine()
      line.DrawLine(bin_min,0,bin_max,0)
      line.SetLineColor(2)
      line.DrawLine(bin_min,-3,bin_max,-3)
      line.DrawLine(bin_min,3,bin_max,3)

    # save output
    canv.cd()
    if process == 'signal':
      plot_label = self.getSignalLabel() 
    elif process == 'background': 
      plot_label = self.getBackgroundLabel()
    else:
      plot_label = self.getBackgroundLabel() + '_' + self.getSignalLabel()
    if process == 'both': process_label = 'prefit'
    else: process_label = process
    canv.SaveAs("{}/{}_fit_{}.png".format(self.outputdir, process_label, plot_label))
    canv.SaveAs("{}/{}_fit_{}.pdf".format(self.outputdir, process_label, plot_label))

    # additionally, get the pull histogram
    if self.plot_pulls and process != 'both':
      self.getPullDistribution(process, hpull, plot_label)

    print ' --- End Run fit --- '


  def getPullDistribution(self, process, hpull, label):
    '''
      Create pull distribution and fit it with a Gaussian
    '''
    canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)
    hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)

    if process == 'signal': # or (process == 'background' and not self.do_blind):
      bin_min, bin_max = self.getRegion(nsigma=self.mass_window_size)
    else:
      bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size)
   
    for i in range(0, hpull.GetN()):
      x = ROOT.Double()
      point = ROOT.Double()
      hpull.GetPoint(i,x,point) 
      if x<bin_min or x>bin_max: continue
      hist_pull.Fill(point)

    hist_pull.SetTitle("")
    hist_pull.SetLineColor(4)
    hist_pull.SetLineWidth(2)
    hist_pull.Draw()

    Xaxis = hist_pull.GetXaxis()
    Yaxis = hist_pull.GetYaxis()
    Xaxis.SetTitle("pulls")
    Xaxis.SetTitleSize(0.045)
    Xaxis.SetLabelSize(0.045)
    Xaxis.SetTitleOffset(1.1)
    Yaxis.SetTitleSize(0.045)
    Yaxis.SetLabelSize(0.045)
    Yaxis.SetTitleOffset(1.26)
    ROOT.gStyle.SetOptStat(0)
    hist_pull.Draw()

    # and fit it
    fgauss= ROOT.TF1("fgauss", "gaus", -5, 5)
    fgauss.SetLineColor(2)
    hist_pull.Fit("fgauss")
    fgauss.Draw("same")
    ROOT.gStyle.SetOptFit(0011)
   
    if process == 'both': process_label = 'prefit'
    else: process_label = process
    canv_pull.SaveAs("{}/{}_pulls_{}.png".format(self.outputdir, process_label, label))
    canv_pull.SaveAs("{}/{}_pulls_{}.pdf".format(self.outputdir, process_label, label))


  def getBackgroundYieldsFromFit(self):
    # one has to make sure that the fit was run beforehand
    n_bkg = self.n_bkg.getVal()
    if n_bkg ==  100.:
      raise RuntimeError('[fitter] It seems like the fit was not performed. Please check. \n-->Aborting')
    # the following not needed if we keep the sidebands
    #if self.do_blind:
    #  n_bkg = n_bkg * self.mass_window_size / self.fit_window_size
    return n_bkg


  #def getBackgroundYieldsFromFit(self):
  #  self.createFitModels()
  #  self.performFit(process='background', label=label)
  #  n_bkg = self.n_bkg.getVal()
  #  # the following not needed if we keep the sidebands
  #  #if self.do_blind:
  #  #  n_bkg = n_bkg * self.mass_window_size / self.fit_window_size
  #  return n_bkg


  def getSignalYieldsFromHist(self):
    '''
      Returns the normalised number of signal yields in the fit window
    '''
    # define ranges and binning
    fit_window_min, fit_window_max = self.getRegion(nsigma=self.fit_window_size)

    # get the signal files
    signal_files = self.tools.getSignalFileList(signal_label=self.signal_label, mass=self.signal_mass, ctau=self.signal_ctau, strategy=self.reweighting_strategy, is_bc=self.is_bc)
    #TODO do we want to cancel the reweighting for generated points?

    # define selection
    cond_sig = 'ismatched==1' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1'
    selection_sig = cond_sig + ' && ' + self.selection #TODO in case we apply gen-matching condition, do it here /self.efficiency_genmatching
    
    # define signal weights
    weight_sig, weight_sig_list = self.getMCWeight(signal_files, lumi=self.lumi_target, is_bc=self.is_bc)
    print 'sel sig before: {}'.format(selection_sig)

    # get the tree
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    if self.do_cutbased:
      tree_sig = ROOT.TChain(treename)
      for signal_file in signal_files:
        filename = signal_file.filename if not self.is_bc else signal_file.filename_Bc
        tree_sig.Add(filename)
    elif self.do_mva:
      score_label = self.getSignalLabel()
      filename_sig = self.mva_tools.getFileWithScore(files=signal_files, training_label=self.training_label, do_parametric=self.do_parametric, mass=self.signal_mass, category_label=self.category_label, selection=selection_sig, weights=weight_sig_list, label=score_label, treename=treename, is_bc=self.is_bc) 
      file_sig = self.tools.getRootFile(filename_sig)
      tree_sig = self.tools.getTree(file_sig, treename)

    # create histogram
    hist_name = 'hist_signal'
    hist = ROOT.TH1D(hist_name, hist_name, self.nbins, fit_window_min, fit_window_max)
    branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    tree_sig.Project(hist_name, branch_name , '({sel}) * ({wght})'.format(sel=selection_sig if self.do_cutbased else self.mva_tools.getScoreSelection(selection_sig), wght=weight_sig))
    print 'sel sig: ({sel}) * ({wght})'.format(sel=selection_sig if self.do_cutbased else self.mva_tools.getScoreSelection(selection_sig), wght=weight_sig)

    # get the number of yields
    n_sig = hist.Integral()

    print 'n_sig = {}'.format(n_sig)

    return n_sig


  def getSignalYields(self):
    n_sig = self.getSignalYieldsFromHist()

    return n_sig


  def createSignalWorkspace(self, label=''):
    print ' --- Creating Signal Workspace --- '

    if label == '': label = self.getSignalLabel()

    # create workspace
    workspace_filename = '{}/workspace_signal_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')

    # import model
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size) # sidebands are included
    workspace.factory('hnl_mass[{}, {}]'.format(bin_min, bin_max))
    getattr(workspace, 'import')(self.signal_model)

    # make sure that the shape parameters are fixed
    it = workspace.allVars().createIterator() 
    all_vars = [it.Next() for _ in range( workspace.allVars().getSize())] 
    for var in all_vars: 
      # keep the mean floating (energy scale correction)
      if var.GetName() != 'hnl_mass' and 'mean' not in var.GetName():
        var.setConstant()

    workspace.Write()
    workspace.Print()
    output_file.Close()

    print '--> {} created'.format(workspace_filename)


  def createBackgroundWorkspace(self, label=''):
    #NOTE not used when using the discrete profiling; veto not applied 

    print ' --- Creating Background Workspace --- '

    if label == '': label = self.getBackgroundLabel()

    # create workspace
    workspace_filename = '{}/workspace_background_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')

    # create factory
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size) # sidebands are included
    workspace.factory('hnl_mass[{}, {}]'.format(bin_min, bin_max)) #TODO instead use the sigma from the fit

    if self.background_model_label == 'chebychev':
      workspace.factory('RooChebychev::qcd(hnl_mass, a0{lbl}[{ini}, {down}, {up}])'.format(
      #workspace.factory('RooChebychev::qcd(hnl_mass, a0[{ini}, {down}, {up}])'.format(
            lbl = label,
            ini = self.a0.getVal(),
            down = self.a0.getVal() - 2*self.a0.getError(),
            up = self.a0.getVal() + 2*self.a0.getError(),
            )
          )

    it = workspace.allVars().createIterator() 
    all_vars = [it.Next() for _ in range( workspace.allVars().getSize())] 
    for var in all_vars: 
      var.setBins(self.nbins)

    workspace.Write()
    workspace.Print()
    output_file.Close()

    print '--> {} created'.format(workspace_filename)


  def createDataObsWorkspace(self, label=''):
    print ' --- Creating Data Obs Workspace --- '
    #NOTE if using mva as a selection, the tree is going to be the same as for the background. That is, the signal region
    # is going to be blinded, if running blinded. Reason for that is that the tree is going to be taken the same as for
    # the background. To change that, modify the score_label below

    if label == '': label = self.getBackgroundLabel()

    selection_dataobs = self.selection
    if self.do_veto_SM:
      selection_dataobs += ' && (hnl_mass < {} || hnl_mass > {})'.format(self.veto_SM.range_min, self.veto_SM.range_max)

    # getting the observed data
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    if self.do_cutbased:
      tree = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree.Add(data_file.filename) 
    elif self.do_mva:
      score_label = label
      filename = self.mva_tools.getFileWithScore(files=self.data_files, training_label=self.training_label, do_parametric=self.do_parametric, mass=self.signal_mass, category_label=self.category_label, selection=selection_dataobs, label=score_label, treename=treename) 
      file_data = self.tools.getRootFile(filename)
      tree = self.tools.getTree(file_data, treename)

    # create binned dataset
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size)
    hnl_mass = ROOT.RooRealVar("hnl_mass","hnl_mass", bin_min, bin_max)
    hist_name = 'hist'
    hist = ROOT.TH1D(hist_name, hist_name, self.nbins, bin_min, bin_max)
    branch_name = 'hnl_mass'
    tree.Project(hist_name, branch_name, selection_dataobs if self.do_cutbased else self.mva_tools.getScoreSelection(selection_dataobs))

    # normalise to the target luminosity
    hist.Scale(self.lumi_target / self.tools.getDataLumi(self.data_files))
    
    # get the number of yields
    data_obs_yields = hist.Integral()

    # create the workspace
    data_obs_name = 'data_obs_bhnl_m_{}_cat_{}'.format(str(self.signal_mass).replace('.', 'p'), self.category_label)
    data_obs = ROOT.RooDataHist(data_obs_name, data_obs_name, ROOT.RooArgList(hnl_mass), hist)

    # create unbinned dataset
    #quantity_set = self.getQuantitySet()
    # add hnl_mass to the RooArgSet
    #quantity_set.add(hnl_mass)
    #print '-> creating unbinned dataset'
    #data_obs = ROOT.RooDataSet('data_obs', 'data_obs', tree, quantity_set, selection_dataobs)
    #print '-> unbinned dataset created'
    
    # create workspace
    workspace_filename = '{}/workspace_data_obs_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')

    getattr(workspace, 'import')(data_obs)
    workspace.Write()
    workspace.Print()
    output_file.Close()

    print '--> {} created'.format(workspace_filename)

    return data_obs_yields


  def createFTestInputWorkspace(self, label=''):
    '''
      Create workspace containing the RooDataHist used as an input to the F-test
    '''
    print ' --- Creating Input Workspace to F-Test --- '

    treename = 'signal_tree'
    branch_name = 'hnl_mass'
    hist_name = 'hist'

    # define selection
    selection_bkg = self.selection
    if self.do_veto_SM:
      selection_bkg += ' && (hnl_mass < {} || hnl_mass > {})'.format(self.veto_SM.range_min, self.veto_SM.range_max)
    
    # get the tree
    if self.do_cutbased:
      tree = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree.Add(data_file.filename) 
    elif self.do_mva:
      if label == '': label = self.getBackgroundLabel()
      score_label = label
      filename = self.mva_tools.getFileWithScore(files=self.data_files, training_label=self.training_label, do_parametric=self.do_parametric, mass=self.signal_mass, category_label=self.category_label, selection=selection_bkg, label=score_label, treename=treename) 
      file_data = self.tools.getRootFile(filename)
      tree = self.tools.getTree(file_data, treename)

    # create the histogram
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size)
    hist = ROOT.TH1D(hist_name, hist_name, self.nbins, bin_min, bin_max)
    tree.Project(hist_name, branch_name, selection_bkg if self.do_cutbased else self.mva_tools.getScoreSelection(selection_bkg))
    print 'sel bkg: {}'.format(selection_bkg if self.do_cutbased else self.mva_tools.getScoreSelection(selection_bkg))

    # normalise to the target luminosity
    hist.Scale(self.lumi_target / self.tools.getDataLumi(self.data_files))

    # create the binned dataset
    hnl_mass = ROOT.RooRealVar("hnl_mass","hnl_mass", bin_min, bin_max)
    hnl_mass_rdh_name = 'hnl_mass_rdh_bhnl_m_{}_cat_{}'.format(str(self.signal_mass).replace('.', 'p'), self.category_label)
    hnl_mass_rdh = ROOT.RooDataHist(hnl_mass_rdh_name, hnl_mass_rdh_name, ROOT.RooArgList(hnl_mass), hist)

    # import the model in a workspace
    label = 'm_{}_cat_{}'.format(self.signal_mass, self.category_label)
    workspace_filename = '{}/input_workspace_fTest_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE')
    workspace = ROOT.RooWorkspace('fTest_workspace', 'fTest_workspace')

    getattr(workspace, 'import')(hnl_mass_rdh)
    workspace.Write()
    workspace.Print()
    output_file.Close()

    print '--> {} created'.format(workspace_filename)


  def producePrefitPlot(self, label=''):
    self.createFitModels(process='both')
    self.performFit(process='both', label=label)


  def process(self, label=''):
    self.createFitModels()
    self.performFit(process='signal', label=label)
    self.performFit(process='background', label=label)
    self.createWorkspace(label=label)


  def process_signal(self, label=''):
    self.createFitModels(process='signal')
    self.performFit(process='signal', label=label)
    self.createSignalWorkspace(label=label)


  def process_background(self, label=''):
    self.createFitModels(process='background')
    self.performFit(process='background', label=label)
    self.createBackgroundWorkspace(label=label)


  def process_data_obs(self, label=''):
    data_obs_yields = self.createDataObsWorkspace(label=label)
    return data_obs_yields


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  outdirlabel = 'test_code'
  plot_pulls = False
  #selection = 'sv_lxy>5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20 '
  #selection = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge==mu_charge && mu0_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((mu0_charge!=mu_charge && (mu0_mu_mass < 2.9 || mu0_mu_mass > 3.3)) || (mu0_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25'
  selection = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge==mu_charge && mu0_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((mu0_charge!=mu_charge && (mu0_mu_mass < 2.9 || mu0_mu_mass > 3.3)) || (mu0_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25 && score>0.'
  signal_model_label = 'doubleCB'
  signal_files = signal_samples['V10_30Dec21_m1p5'] 
  signal_label = 'V12_08Aug22_m3'
  background_model_label = 'chebychev'
  data_files = data_samples['V10_30Dec21_small']
  do_blind = True
  do_binned_fit = True
  nbins = 100
  lumi_target = 41.6 #5.302
  sigma_B = 472.8e9
  add_CMSlabel = True
  add_lumilabel = True
  CMStag = 'Preliminary'
  resolution_p0 = 0.0004758
  resolution_p1 = 0.007316

  mass_window_size = 3
  fit_window_size = 10

  reweighting_strategy = 'exclusive_fromlargerctau'

  do_cutbased = False
  do_mva = True
  training_label = 'mva/outputs/test_20Aug2022_13h40m03s'

  label = 'test'
  #for signal_file in signal_files:
  #  if signal_file.ctau != 10: continue
  #  category_label = 'lxy1to5_SS'
  #  category_title = '(1<l_{xy}<=5)cm, SS'
  #  fitter = Fitter(signal_file=signal_file, data_files=data_files, selection=selection, signal_model_label=signal_model_label, background_model_label=background_model_label, do_binned_fit=do_binned_fit, do_blind=do_blind, lumi_target=lumi_target, sigma_B=sigma_B, mass_window_size=mass_window_size, fit_window_size=fit_window_size, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls, add_CMSlabel=add_CMSlabel, add_lumilabel=add_lumilabel, CMStag=CMStag, mass=signal_file.mass, resolution=signal_file.resolution, category_label=category_label, category_title=category_title)
  #  #fitter.writeSignalModel(label='test')
  #  #fitter.writeBackgroundModel(label='test')
  #  #fitter.createFitWorkspace()
  #  #fitter.writeFitModels(label='test')
  #  #fitter.performFit(process='signal', label='test')
  #  #fitter.performFit(process='background', label='test')
  #  #fitter.process_signal()
  #  #fitter.process_background()
  #  #fitter.producePrefitPlot()
  #  #print fitter.getSignalYieldsFromHist()

  mass = 3.
  ctau = 100.
  category_label = 'lxy1to5_SS'
  category_title = '(1<l_{xy}<=5)cm, SS'
  fitter = Fitter(signal_label=signal_label, data_files=data_files, selection=selection, do_cutbased=do_cutbased, do_mva=do_mva, training_label=training_label, mass=mass, ctau=ctau, signal_model_label=signal_model_label, background_model_label=background_model_label, do_binned_fit=do_binned_fit, do_blind=do_blind, lumi_target=lumi_target, sigma_B=sigma_B, mass_window_size=mass_window_size, fit_window_size=fit_window_size, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls, add_CMSlabel=add_CMSlabel, add_lumilabel=add_lumilabel, CMStag=CMStag, category_label=category_label, category_title=category_title, reweighting_strategy=reweighting_strategy, resolution_p0=resolution_p0, resolution_p1=resolution_p1)
  fitter.process_signal()
  #fitter.getSignalYields()




