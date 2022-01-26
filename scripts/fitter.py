import os
from os import path
import ROOT
from ROOT import RooFit

from tools import Tools
from samples import signal_samples, data_samples

class Fitter(Tools):
  def __init__(self, signal_file='', data_files='', selection='', signal_model_label=None, background_model_label=None, do_binned_fit=False, do_blind=False, file_type='flat', nbins=250, title=' ', outputdir='', outdirlabel='', category_label='', plot_pulls=True):
    self.tools = Tools()
    self.signal_file = signal_file
    self.data_files = data_files
    self.selection = selection
    self.signal_model_label = signal_model_label
    self.background_model_label = background_model_label
    self.do_binned_fit = do_binned_fit
    self.do_blind = do_blind
    self.file_type = file_type
    self.nbins = nbins
    self.title = title
    self.plot_pulls = plot_pulls
    self.workspacedir = outputdir #TODO adapt
    if self.workspacedir == '':
      self.workspacedir = './'
    self.outputdir = outputdir
    if self.outputdir != '': self.outputdir = self.outputdir + '/fits'
    else: self.outputdir = './myPlots/fits/' + outdirlabel
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))
    self.category_label = category_label
    # define window sizes (multiples of sigma)
    self.mass_window_size = 2
    self.fit_window_size = 4


    signal_model_list = ['doubleCB', 'doubleCBPlusGaussian', 'voigtian']
    if self.signal_model_label != None and self.signal_model_label not in signal_model_list:
      raise RuntimeError('Unrecognised signal model "{}". Please choose among {}'.format(self.signal_model_label, signal_model_list))

    background_model_list = ['chebychev']
    if self.background_model_label != None and self.background_model_label not in background_model_list:
      raise RuntimeError('Unrecognised background model "{}". Please choose among {}'.format(self.background_model_label, background_model_list))

    #TODO there are no weights applied so far (incl. ctau, hlt, pu weights)
    #TODO make sure that the data is normalised to the yields inserted in the datacard
    #TODO make sure that the pdfs are named according to the process name in the datacard

  def getRegion(self, nsigma=2):
    signal_mass = self.signal_file.mass
    signal_resolution = self.signal_file.resolution
    bin_min = signal_mass - nsigma * signal_resolution
    bin_max = signal_mass + nsigma * signal_resolution
    return bin_min, bin_max


  def getSignalLabel(self):
    label = 'm{}_ctau_{}_{}'.format(self.signal_file.mass, self.signal_file.ctau, self.signal_model_label).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    return label


  def getBackgroundLabel(self):
    label = 'm{}_{}'.format(self.signal_file.mass, self.background_model_label).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    if self.do_blind: label += '_blind'
    return label


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
      
    
  def createFitModels(self, label='', do_recreate=True):
    print ' --- Creating Fit Models --- '

    # get the signal region
    signal_mass = self.signal_file.mass
    bin_min, bin_max = self.getRegion(nsigma=self.fit_window_size)
    nbins = self.nbins

    self.mupi_invmass = ROOT.RooRealVar("hnl_mass","hnl_mass", bin_min, bin_max)
    #mupi_invmass.setBins(nbins)

    ### Signal Model ###

    if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
      self.mean_CB  = ROOT.RooRealVar("mean_CB","mean_CB", signal_mass)
      self.sigma_CB = ROOT.RooRealVar("sigma_CB", "sigma_CB", 0.01, 0.005, 0.15)

      self.alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -2, -5, 5)
      self.n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 5)
      self.alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 2, -5, 5)
      self.n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 5)

      self.CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", self.mupi_invmass, self.mean_CB, self.sigma_CB, self.alpha_1, self.n_1)
      self.CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", self.mupi_invmass, self.mean_CB, self.sigma_CB, self.alpha_2, self.n_2)

      # defines the relative importance of the two CBs
      self.sigfrac_CB = ROOT.RooRealVar("sigfrac_CB","sigfrac_CB", 0.5, 0.0 ,1.0)

      if self.signal_model_label == 'doubleCB':
        self.signal_model = ROOT.RooAddPdf("sig", "sig", self.CBpdf_1, self.CBpdf_2, self.sigfrac_CB)

      if self.signal_model_label == 'doubleCBPlusGaussian':
        self.sigma_gauss = ROOT.RooRealVar("sigma_gauss", "sigma_gauss", 0.01, 0.005, 0.15)
        self.gaussian = ROOT.RooGaussian('gaussian', 'gaussian', self.mupi_invmass, self.mean_CB, self.sigma_gauss)

        # defines the relative importance of gaussian wrt doubleCB
        self.sigfrac_gauss = ROOT.RooRealVar("sigfrac_gauss","sigfrac_gauss", 0.5, 0.0 ,1.0)

        self.doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", self.CBpdf_1, self.CBpdf_2, self.sigfrac_CB)
        self.signal_model = ROOT.RooAddPdf('sig', 'sig', self.doubleCBpdf, self.gaussian, self.sigfrac_gauss) # make sure that the model has the same name as in datacard 

    elif self.signal_model_label == 'voigtian':
      self.mean_voigtian  = ROOT.RooRealVar("mean_voigtian","mean_voigtian", signal_mass)
      self.gamma_voigtian = ROOT.RooRealVar("gamma_voigtian", "gamma_voigtian", 0.01, 0., 5.)
      self.sigma_voigtian = ROOT.RooRealVar("sigma_voigtian", "sigma_voigtian", 0.01, 0.005, 0.15)
      self.signal_model = ROOT.RooVoigtian('sig', 'sig', self.mupi_invmass, self.mean_voigtian, self.gamma_voigtian, self.sigma_voigtian)

    ### Background Model ###

    # Define the background model 
    if self.background_model_label == 'chebychev':
      self.n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', 100, 0, 100000)
      self.a0 = ROOT.RooRealVar('a0', 'a0', 0.01, -10, 10)
      self.chebychev = ROOT.RooChebychev('chebychev', 'chebychev', self.mupi_invmass, ROOT.RooArgList(self.a0))
      self.background_model = ROOT.RooAddPdf('qcd', 'qcd', ROOT.RooArgList(self.chebychev), ROOT.RooArgList(self.n_bkg))
      #a0 = ROOT.RooRealVar('a0', 'a0', -1.36, -1.52, -1.20)
      #a1 = ROOT.RooRealVar('a1', 'a1', 0.53, 0.29, 0.76)
      #a2 = ROOT.RooRealVar('a2', 'a2', -0.14, -0.32, 0.04)
      #background_model = ROOT.RooChebychev('qcd', 'qcd', mupi_invmass, ROOT.RooArgList(a0, a1, a2))

    print '--> Models created'


  #def createFitWorkspace(self, label=''):
  #  if label == '': label = self.getSignalLabel()

  #  # get the models
  #  signal_model = self.getSignalModel()
  #  print signal_model
  #  print 'got signal model'
  #  #background_model = self.getBackgroundModel()

  #  # import the model in a workspace
  #  workspace = ROOT.RooWorkspace('workspace', 'workspace')
  #  getattr(workspace, 'import')(signal_model)
  #  #getattr(workspace, 'import')(background_model)
  #  workspace.Print()
  #  workspace.writeToFile('workspace_{}.root'.format(label))

  #  # workspace will remain in memory after macro finishes
  #  #ROOT.gDirectory.Add(workspace)


  #def performFitFromWorkspace(self, process, label=''):
  #  '''
  #    Performs the fit by loading the fit workspace
  #    process corresponds to either 'signal' or 'background'
  #  '''
  #  print ' --- Running the fits --- '

  #  if process not in ['signal', 'background']:
  #    raise RuntimeError("[fitter] Unrecognised process. Please choose among ['signal', 'background']")

  #  # get the label
  #  if label == '' and process == 'signal': label = self.getSignalLabel() 
  #  if label == '' and process == 'background': label = self.getBackgroundLabel() 

  #  # load the workspace
  #  workspace_filename = '{}/workspace_{}.root'.format(self.workspacedir, label)
  #  ws_file = self.tools.getRootFile(workspace_filename, with_ext=False)
  #  ws = ws_file.Get('workspace')
  #  #ws.Print()
  #  pdf_name = 'sig' if process == 'signal' else 'qcd'

  #  # open the file and get the tree
  #  treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
  #  if process == 'signal':
  #    inputfile = self.tools.getRootFile(self.signal_file.filename, with_ext=False)
  #    tree = self.tools.getTree(inputfile, treename)
  #  else:
  #    tree = ROOT.TChain(treename)
  #    for data_file in self.data_files:
  #      tree.Add(data_file.filename) 

  #  # get signal mass
  #  signal_mass = self.signal_file.mass

  #  # declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
  #  bin_min, bin_max = self.getRegion()
  #  nbins = self.nbins

  #  # build the binned data
  #  hist_name = 'hist'
  #  hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
  #  branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
  #  if process == 'signal':
  #    cond = 'ismatched==1 && mu_isdsa==0' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx]==0'
  #  else:
  #    cond = 'hnl_pt > 0' # dummy condition
  #  selection = cond + ' && ' + self.selection
  #  tree.Project(hist_name, branch_name , selection)
  #  rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(ws.var('mupi_invmass')), hist)

  #  # create canvas
  #  canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

  #  # and the two pads
  #  pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
  #  pad1.SetLeftMargin(0.15)
  #  pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
  #  pad2.SetLeftMargin(0.15)
  #  pad1.Draw()
  #  pad2.Draw()

  #  frame = ws.var('mupi_invmass').frame(ROOT.RooFit.Title(self.title))

  #  # plot the data
  #  rdh.plotOn(frame, ROOT.RooFit.Name("data"))

  #  # fit the PDF to the data
  #  fit_range_min = bin_min #signal_mass - 0.1*signal_mass
  #  fit_range_max = bin_max #signal_mass + 0.1*signal_mass
  #  ws.var('mupi_invmass').setRange("peak", fit_range_min, fit_range_max)
  #  result = ws.pdf(pdf_name).fitTo(rdh, ROOT.RooFit.Range("peak"))

  #  # plot the fit 		
  #  if process == 'signal':
  #    if self.signal_model == 'doubleCB' or self.signal_model == 'doubleCBPlusGaussian':
  #      ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1"),ROOT.RooFit.Components("CBpdf_1"), ROOT.RooFit.LineStyle(ROOT.kDashed))
  #      ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2"),ROOT.RooFit.Components("CBpdf_2"), ROOT.RooFit.LineStyle(ROOT.kDashed))
  #      if self.signal_model == 'doubleCBPlusGaussian':
  #        ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
  #  ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))

  #  # and write the fit parameters
  #  ws.pdf(pdf_name).paramOn(frame,   
  #       ROOT.RooFit.Layout(0.2, 0.4, 0.8),
  #       ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
  #       )

  #  frame.getAttText().SetTextSize(0.03)
  #  frame.getAttLine().SetLineColorAlpha(0, 0)
  #  frame.getAttFill().SetFillColorAlpha(0, 0)

  #  # compute the chisquare
  #  chisquare = frame.chiSquare(pdf_name,"data")

  #  # and print it
  #  label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
  #  label1.SetBorderSize(0)
  #  label1.SetFillColor(ROOT.kWhite)
  #  label1.SetTextSize(0.03)
  #  label1.SetTextFont(42)
  #  label1.SetTextAlign(11)
  #  if process == 'signal':
  #    label_text = 'mass {}GeV, ctau {}mm'.format(signal_mass, self.signal_file.ctau)
  #  else: 
  #    label_text = 'Mass window around {}GeV'.format(signal_mass)
  #  label1.AddText(label_text)
  #  qte = '#chi^{2}/ndof'
  #  label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
  #  print "chisquare = {}".format(chisquare)

  #  # We define and plot the pull 		
  #  hpull = frame.pullHist()
  #  for i in range(0, frame.GetNbinsX()):
  #     hpull.SetPointError(i,0,0,0,0)

  #  # create a new frame to draw the pull distribution and add the distribution to the frame
  #  frame2 = ws.var('mupi_invmass').frame(ROOT.RooFit.Title(" "))
  #  frame2.addPlotable(hpull,"P")#,"E3")

  #  # plot of the curve and the fit
  #  canv.cd()
  #  pad1.cd()

  #  frame.GetXaxis().SetTitleSize(0.04)
  #  frame.GetXaxis().SetTitle("#mu#pi invariant mass [GeV]")
  #  frame.GetYaxis().SetTitleSize(0.04)
  #  frame.GetYaxis().SetTitleOffset(1.1)
  #  frame.Draw()
  #  label1.Draw()

  #  # add the legend
  #  leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.4, xmax=0.8, ymax=0.6, size=0.03, do_alpha=True)
  #  if process == 'signal':
  #    if self.signal_model == 'doubleCB':
  #      model_label = 'Double Crystal Ball'
  #    elif self.signal_model == 'doubleCBPlusGaussian':
  #      model_label = 'Double Crystal Ball + Gaussian'
  #    elif self.signal_model == 'voigtian':
  #      model_label = 'Voigtian'
  #    leg.AddEntry(frame.findObject(pdf_name), model_label)
  #    if self.signal_model == 'doubleCB' or self.signal_model == 'doubleCBPlusGaussian':
  #      leg.AddEntry(frame.findObject('CBpdf_1'), 'CB_1')
  #      leg.AddEntry(frame.findObject('CBpdf_2'), 'CB_2')
  #      if self.signal_model == 'doubleCBPlusGaussian':
  #        leg.AddEntry(frame.findObject('gaussian'), 'Gaussian')
  #  else:
  #    if self.background_model == 'chebychev':
  #      model_label = 'Chebychev'
  #    leg.AddEntry(frame.findObject(pdf_name), model_label)
  #  leg.Draw()

  #  # plot of the residuals
  #  pad2.cd()
  #  ROOT.gPad.SetLeftMargin(0.15) 
  #  ROOT.gPad.SetPad(0.01,0.01,0.99,0.2)

  #  frame2.GetYaxis().SetNdivisions(3)
  #  frame2.GetYaxis().SetLabelSize(0.17)
  #  frame2.GetYaxis().SetTitleSize(0.17)
  #  frame2.GetYaxis().SetTitleOffset(0.24)
  #  frame2.GetYaxis().SetRangeUser(-5,5)	
  #  frame2.GetYaxis().SetTitle("Pulls")	
  #  frame2.GetXaxis().SetTitle("")	
  #  frame2.GetXaxis().SetLabelOffset(5)	
  #  frame2.Draw()

  #  line = ROOT.TLine()
  #  line.DrawLine(bin_min,0,bin_max,0)
  #  line.SetLineColor(2)
  #  line.DrawLine(bin_min,-3,bin_max,-3)
  #  line.DrawLine(bin_min,3,bin_max,3)

  #  # save output
  #  canv.cd()
  #  plot_label = self.getSignalLabel() if process == 'signal' else self.getBackgroundLabel()
  #  canv.SaveAs("{}/{}_fit_{}.png".format(self.outputdir, process, plot_label))
  #  canv.SaveAs("{}/{}_fit_{}.pdf".format(self.outputdir, process, plot_label))

  #  # additionally, get the pull histogram
  #  if self.plot_pulls:
  #    self.getPullDistribution(process, hpull, plot_label)

  #  print ' --- End Run fit --- '


  def performFit(self, process, label=''):
    '''
      process corresponds to either 'signal' or 'background'
    '''
    print ' --- Running the fits --- '

    if process not in ['signal', 'background']:
      raise RuntimeError("[fitter] Unrecognised process. Please choose among ['signal', 'background']")

    # get the label
    if label == '' and process == 'signal': label = self.getSignalLabel() 
    if label == '' and process == 'background': label = self.getBackgroundLabel() 

    # get signal mass
    signal_mass = self.signal_file.mass

    # open the file and get the tree
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    if process == 'signal':
      inputfile = self.tools.getRootFile(self.signal_file.filename, with_ext=False)
      tree = self.tools.getTree(inputfile, treename)
    else:
      tree = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree.Add(data_file.filename) 

    # define ranges and binning
    mass_window_min, mass_window_max = self.getRegion(nsigma=self.mass_window_size)
    if self.do_blind:
      fit_window_min, fit_window_max = self.getRegion(nsigma=self.fit_window_size)
    else:
      fit_window_min, fit_window_max = self.getRegion(nsigma=self.mass_window_size)

    # define selection
    if process == 'signal':
      cond = 'ismatched==1' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1'
    else:
      cond = 'hnl_pt > 0' # dummy condition
    selection = cond + ' && ' + self.selection
    if self.do_blind:
      selection += ' && (hnl_mass < {} || hnl_mass > {})'.format(mass_window_min, mass_window_max)

    if self.do_binned_fit:
      # build the binned data
      hist_name = 'hist'
      #hist = ROOT.TH1D(hist_name, hist_name, self.nbins, bin_min, bin_max)
      hist = ROOT.TH1D(hist_name, hist_name, self.nbins, fit_window_min, fit_window_max)
      branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
      tree.Project(hist_name, branch_name , selection)
      rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(self.mupi_invmass), hist)

    else:
      # get RooArgSet
      quantity_set = self.getQuantitySet()
      # add mupi_invmass to the RooArgSet
      ROOT.SetOwnership(self.mupi_invmass, False)
      quantity_set.add(self.mupi_invmass)
     
      print '-> creating unbinned dataset'
      rds = ROOT.RooDataSet('rds', 'rds', tree, quantity_set, selection)
      print '-> unbinned dataset created'

    # create canvas
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    pad2.Draw()

    frame = self.mupi_invmass.frame(ROOT.RooFit.Title(self.title))

    # plot the data
    if self.do_binned_fit:
      rdh.plotOn(frame, ROOT.RooFit.Name("data"))
    else:
      rds.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(self.nbins))

    # define fit ranges
    self.mupi_invmass.setRange("peak", mass_window_min, mass_window_max)
    self.mupi_invmass.setRange('sideband_left', fit_window_min, mass_window_min)
    self.mupi_invmass.setRange('sideband_right', mass_window_max, fit_window_max)
    if process == 'signal' or (process == 'background' and not self.do_blind):
      fit_range = 'peak'
    else:
      fit_range = 'sideband_left,sideband_right'

    # fit the PDF to the data
    if process == 'signal':
      if self.do_binned_fit:
        result = self.signal_model.fitTo(rdh, ROOT.RooFit.Range(fit_range))
      else:
        result = self.signal_model.fitTo(rds, ROOT.RooFit.Range(fit_range))
    else:
      if self.do_binned_fit:
        result = self.background_model.fitTo(rdh, ROOT.RooFit.Range(fit_range))
      else:
        result = self.background_model.fitTo(rds, ROOT.RooFit.Range(fit_range))

    # plot the fit 		
    pdf_name = 'sig' if process == 'signal' else 'qcd'
    if process == 'signal':
      if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
        self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1"),ROOT.RooFit.Components("CBpdf_1"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2"),ROOT.RooFit.Components("CBpdf_2"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        if self.signal_model_label == 'doubleCBPlusGaussian':
          self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
      self.signal_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))
    else:
      self.background_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))

    # and write the fit parameters
    if process == 'signal':
      self.signal_model.paramOn(frame,   
           ROOT.RooFit.Layout(0.2, 0.4, 0.8),
           ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
           )
    else:
      self.background_model.paramOn(frame,   
           ROOT.RooFit.Layout(0.2, 0.4, 0.8),
           ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
           )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # compute the chisquare
    if process == 'signal':
      chisquare = frame.chiSquare("sig","data")
    else:
      chisquare = frame.chiSquare("qcd","data")

    # and print it
    label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    label1.SetTextSize(0.03)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    if process == 'signal':
      label_text = 'mass {}GeV, ctau {}mm'.format(signal_mass, self.signal_file.ctau)
    else: 
      label_text = 'Mass window around {}GeV'.format(signal_mass)
    label1.AddText(label_text)
    qte = '#chi^{2}/ndof'
    label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    print "chisquare = {}".format(chisquare)

    # We define and plot the pull 		
    if self.do_blind:
      curve1 = frame.getObject(1)
      curve2 = frame.getObject(2)
      datahist = frame.getHist('data')
      hpull1 = datahist.makePullHist(curve1, True)
      hpull2 = datahist.makePullHist(curve2, True)
      frame2 = self.mupi_invmass.frame(ROOT.RooFit.Title(" "))
      frame2.addPlotable(hpull1,"P")
      frame2.addPlotable(hpull2,"P")
      hpull = frame.pullHist() # needed to get pull distribution
    else:
      hpull = frame.pullHist()
      for i in range(0, frame.GetNbinsX()):
         hpull.SetPointError(i,0,0,0,0)
      frame2 = self.mupi_invmass.frame(ROOT.RooFit.Title(" "))
      frame2.addPlotable(hpull,"P")

    # plot of the curve and the fit
    canv.cd()
    pad1.cd()

    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitle("#mu#pi invariant mass [GeV]")
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()
    label1.Draw()

    # add the legend
    leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.4, xmax=0.8, ymax=0.6, size=0.03, do_alpha=True)
    if process == 'signal':
      if self.signal_model_label == 'doubleCB':
        model_label = 'Double Crystal Ball'
      elif self.signal_model_label == 'doubleCBPlusGaussian':
        model_label = 'Double Crystal Ball + Gaussian'
      elif self.signal_model_label == 'voigtian':
        model_label = 'Voigtian'
      leg.AddEntry(frame.findObject(pdf_name), model_label)
      if self.signal_model_label == 'doubleCB' or self.signal_model_label == 'doubleCBPlusGaussian':
        leg.AddEntry(frame.findObject('CBpdf_1'), 'CB_1')
        leg.AddEntry(frame.findObject('CBpdf_2'), 'CB_2')
        if self.signal_model_label == 'doubleCBPlusGaussian':
          leg.AddEntry(frame.findObject('gaussian'), 'Gaussian')
    else:
      if self.background_model_label == 'chebychev':
        model_label = 'Chebychev'
      leg.AddEntry(frame.findObject(pdf_name), model_label)
    leg.Draw()

    # plot of the residuals
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

    if process == 'signal' or (process == 'background' and not self.do_blind):
      bin_min = mass_window_min
      bin_max = mass_window_max
    else:
      bin_min = fit_window_min
      bin_max = fit_window_max
    line = ROOT.TLine()
    line.DrawLine(bin_min,0,bin_max,0)
    line.SetLineColor(2)
    line.DrawLine(bin_min,-3,bin_max,-3)
    line.DrawLine(bin_min,3,bin_max,3)

    # save output
    canv.cd()
    plot_label = self.getSignalLabel() if process == 'signal' else self.getBackgroundLabel()
    canv.SaveAs("{}/{}_fit_{}.png".format(self.outputdir, process, plot_label))
    canv.SaveAs("{}/{}_fit_{}.pdf".format(self.outputdir, process, plot_label))

    # additionally, get the pull histogram
    if self.plot_pulls:
      self.getPullDistribution(process, hpull, plot_label)

    print ' --- End Run fit --- '


  def getPullDistribution(self, process, hpull, label):
    '''
      Create pull distribution and fit it with a Gaussian
    '''
    canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)
    hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)

    if process == 'signal' or (process == 'background' and not self.do_blind):
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
   
    canv_pull.SaveAs("{}/{}_pulls_{}.png".format(self.outputdir, process, label))
    canv_pull.SaveAs("{}/{}_pulls_{}.pdf".format(self.outputdir, process, label))


  def getBackgroundYieldsFromFit(self):
    # one has to make sure that the fit was run beforehand
    n_bkg = self.n_bkg.getVal()
    if n_bkg ==  100.:
      raise RuntimeError('[fitter] It seems like the fit was not performed. Please check. \n-->Aborting')
    if self.do_blind:
      n_bkg = n_bkg * self.mass_window_size / self.fit_window_size
    return n_bkg


  def createWorkspace(self, label=''):
    if label == '': label = self.getSignalLabel()

    # getting the observed data #TODO create function
    #treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    #tree = ROOT.TChain(treename)
    #for data_file in self.data_files:
    #  tree.Add(data_file.filename) 

    #bin_min, bin_max = self.getRegion()
    #nbins = self.nbins
    #hist_name = 'hist'
    #hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
    #branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    #selection = self.selection
    #tree.Project(hist_name, branch_name , selection)
    #data_obs = ROOT.RooDataHist("data_obs", "data_obs", ROOT.RooArgList(mupi_invmass), hist)
    data_obs = ROOT.RooDataSet("data_obs", "data_obs", ROOT.RooArgSet(self.mupi_invmass)) #TODO for the moment does not contain sensible data
    
    # import the model in a workspace
    workspace_filename = '{}/workspace_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')

    # create factory
    bin_min, bin_max = self.getRegion()
    workspace.factory('mupi_invmass[{}, {}]'.format(bin_min, bin_max)) #TODO instead use the sigma from the fit
    if self.signal_model_label == 'voigtian':
      workspace.factory('RooVoigtian::sig(mupi_invmass, mean_voigtian[{m}], gamma_voigtian[{g}], sigma_voigtian[{s}])'.format(
            m = self.mean_voigtian.getVal(),
            g = self.gamma_voigtian.getVal(),
            s = self.sigma_voigtian.getVal(),
            )
          )
    if self.background_model_label == 'chebychev':
      workspace.factory('RooChebychev::qcd(mupi_invmass, a0[{ini}, {down}, {up}])'.format(
            ini = self.a0.getVal(),
            down = self.a0.getVal() - 2*self.a0.getError(),
            up = self.a0.getVal() + 2*self.a0.getError(),
            )
          )

    getattr(workspace, 'import')(data_obs)
    workspace.Write()
    workspace.Print()
    output_file.Close()

    print '--> {} created'.format(workspace_filename)


  def process(self, label=''):
    #self.defineVariables()
    self.createFitModels()
    #self.performFit(process='signal', label=label)
    self.performFit(process='background', label=label)
    #self.createWorkspace(label=label)



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  outdirlabel = 'V10_30Dec21_samples_sel'
  plot_pulls = True
  selection = 'sv_lxy>5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20 '
  signal_model_label = 'voigtian'
  signal_files = signal_samples['V10_30Dec21_m3'] 
  background_model_label = 'chebychev'
  data_files = data_samples['V10_30Dec21']
  do_blind = True
  do_binned_fit = True
  nbins = 50

  label = 'test'
  for signal_file in signal_files:
    fitter = Fitter(signal_file=signal_file, data_files=data_files, selection=selection, signal_model_label=signal_model_label, background_model_label=background_model_label, do_binned_fit=do_binned_fit, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
    #fitter.writeSignalModel(label='test')
    #fitter.writeBackgroundModel(label='test')
    #fitter.createFitWorkspace()
    #fitter.writeFitModels(label='test')
    #fitter.performFit(process='signal', label='test')
    #fitter.performFit(process='background', label='test')
  fitter.process()

  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)

  #selection = 'sv_lxy<=1 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10'
  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)

  #selection = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25'
  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)



