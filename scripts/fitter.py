import os
from os import path
import ROOT
from ROOT import RooFit

from tools import Tools
from samples import signal_samples, data_samples

class Fitter(Tools):
  def __init__(self, signal_file='', data_files='', selection='', signal_model=None, background_model=None, do_blind=False, file_type='flat', nbins=250, title=' ', outputdir='', outdirlabel='', category_label='', plot_pulls=True):
    self.tools = Tools()
    self.signal_file = signal_file
    self.data_files = data_files
    self.selection = selection
    self.signal_model = signal_model
    self.background_model = background_model
    self.do_blind = do_blind
    self.file_type = file_type
    self.nbins = nbins
    self.title = title
    self.plot_pulls = plot_pulls
    self.workspacedir = outputdir #TODO adapt
    self.outputdir = outputdir
    if self.outputdir != '': self.outputdir = self.outputdir + '/fits'
    else: self.outputdir = './myPlots/fits/' + outdirlabel
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))
    self.category_label = category_label


    signal_model_list = ['doubleCB', 'doubleCBPlusGaussian']
    if self.signal_model != None and self.signal_model not in signal_model_list:
      raise RuntimeError('Unrecognised signal model "{}". Please choose among {}'.format(self.signal_model, signal_model_list))

    background_model_list = ['chebychev']
    if self.background_model != None and self.background_model not in background_model_list:
      raise RuntimeError('Unrecognised background model "{}". Please choose among {}'.format(self.background_model, background_model_list))

    #TODO there are no weights applied so far (incl. ctau, hlt, pu weights)
    #TODO make sure that the data is normalised to the yields inserted in the datacard
    #TODO make sure that the pdfs are named according to the process name in the datacard

  def getFitRegion(self, nsigma=2):
    signal_mass = self.signal_file.mass
    #bin_min = signal_mass - 0.15 * signal_mass
    #bin_max = signal_mass + 0.15 * signal_mass
    signal_resolution = self.signal_file.resolution
    bin_min = signal_mass - nsigma * signal_resolution
    bin_max = signal_mass + nsigma * signal_resolution
    return bin_min, bin_max


  def getSignalLabel(self):
    label = 'm{}_ctau_{}_{}'.format(self.signal_file.mass, self.signal_file.ctau, self.signal_model).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    return label


  def getBackgroundLabel(self):
    label = 'm{}_{}'.format(self.signal_file.mass, self.background_model).replace('.', 'p') 
    if self.category_label != '': label += '_{}'.format(self.category_label)
    if self.do_blind: label += '_blind'
    return label


  #def writeSignalModel(self, label='', do_recreate=True):
  #  '''
  #    Write the signal model in a workspace
  #  '''
  #  # get the label
  #  if label == '': label = self.getSignalLabel() 

  #  # get the signal region
  #  signal_mass = self.signal_file.mass
  #  bin_min, bin_max = self.getFitRegion()
  #  nbins = self.nbins

  #  mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", bin_min, bin_max)
  #  mupi_invmass.setBins(nbins)

  #  # Define the signal model 
  #  mean_init = signal_mass - 0.01*signal_mass
  #  mean_min = signal_mass - 0.05*signal_mass
  #  mean_max = signal_mass + 0.05*signal_mass
  #  mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

  #  sigma = ROOT.RooRealVar("sigma", "sigma", 0.01, 0.005, 0.15)

  #  alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -2, -5, 5)
  #  n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 5)
  #  alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 2, -5, 5)
  #  n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 5)

  #  CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mupi_invmass, mean, sigma, alpha_1, n_1)
  #  CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mupi_invmass, mean, sigma, alpha_2, n_2)

  #  # defines the relative importance of the two CBs
  #  sigfrac_CB = ROOT.RooRealVar("sigfrac_CB","sigfrac_CB", 0.5, 0.0 ,1.0)

  #  if self.signal_model == 'doubleCBPlusGaussian':
  #    sigma_gauss = ROOT.RooRealVar("sigma_gauss", "sigma_gauss", 0.01, 0.005, 0.15)
  #    gaussian = ROOT.RooGaussian('gaussian', 'gaussian', mupi_invmass, mean, sigma_gauss)

  #    # defines the relative importance of gaussian wrt doubleCB
  #    sigfrac_gauss = ROOT.RooRealVar("sigfrac_gauss","sigfrac_gauss", 0.5, 0.0 ,1.0)

  #  # build the signal model
  #  if self.signal_model == 'doubleCB':
  #    signal_model = ROOT.RooAddPdf("sig", "sig", CBpdf_1, CBpdf_2, sigfrac_CB)

  #  elif self.signal_model == 'doubleCBPlusGaussian':
  #    doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", CBpdf_1, CBpdf_2, sigfrac_CB)
  #    signal_model = ROOT.RooAddPdf('sig', 'sig', doubleCBpdf, gaussian, sigfrac_gauss) # make sure that the model has the same name as in datacard 

  #  # import the model in a workspace
  #  output_file = ROOT.TFile('workspace_{}.root'.format(label), 'RECREATE' if do_recreate else 'UPDATE')
  #  workspace = ROOT.RooWorkspace('workspace', 'workspace')
  #  getattr(workspace, 'import')(signal_model)
  #  workspace.Write()
  #  workspace.Print()
  #  #workspace.writeToFile('workspace_{}.root'.format(label))
  #  output_file.Close()

  #  ## workspace will remain in memory after macro finishes
  #  ##ROOT.gDirectory.Add(workspace)
  #  #ROOT.gDirectory.Add(signal_model)

  #  #return signal_model


  #def writeBackgroundModel(self, label='', do_recreate=False):
  #  '''
  #    Write the background model in a workspace
  #  '''
  #  # get the label
  #  if label == '': label = self.getBackgroundLabel(mass=mass) 

  #  mass = self.signal_file.mass

  #  # declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
  #  bin_min, bin_max = self.getFitRegion()
  #  nbins = self.nbins

  #  mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", bin_min, bin_max)
  #  mupi_invmass.setBins(nbins)

  #  # Define the background model 
  #  #if self.background_model == 'chebychev':
  #  a0 = ROOT.RooRealVar('a0', 'a0', 0.01, -10, 10)
  #  background_model = ROOT.RooChebychev('qcd', 'qcd', mupi_invmass, ROOT.RooArgList(a0))

  #  # import the model in a workspace
  #  output_file = ROOT.TFile('workspace_{}.root'.format(label), 'RECREATE' if do_recreate else 'UPDATE')
  #  workspace = ROOT.RooWorkspace('workspace', 'workspace')
  #  getattr(workspace, 'import')(background_model)
  #  workspace.Write()
  #  workspace.Print()
  #  #workspace.writeToFile('workspace_{}.root'.format(label))
  #  output_file.Close()

  #  ## workspace will remain in memory after macro finishes
  #  ##ROOT.gDirectory.Add(workspace)

  #  #return background_model


  def writeFitModels(self, label='', do_recreate=True):
    '''
      Write the signal model in a workspace
    '''
    # get the label
    if label == '': label = self.getSignalLabel() 

    # get the signal region
    signal_mass = self.signal_file.mass
    bin_min, bin_max = self.getFitRegion()
    nbins = self.nbins

    mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", bin_min, bin_max)
    mupi_invmass.setBins(nbins)

    # Define the signal model 
    mean_init = signal_mass - 0.01*signal_mass
    mean_min = signal_mass - 0.05*signal_mass
    mean_max = signal_mass + 0.05*signal_mass
    mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

    sigma = ROOT.RooRealVar("sigma", "sigma", 0.01, 0.005, 0.15)

    alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -2, -5, 5)
    n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 5)
    alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 2, -5, 5)
    n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 5)

    CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mupi_invmass, mean, sigma, alpha_1, n_1)
    CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mupi_invmass, mean, sigma, alpha_2, n_2)

    # defines the relative importance of the two CBs
    sigfrac_CB = ROOT.RooRealVar("sigfrac_CB","sigfrac_CB", 0.5, 0.0 ,1.0)

    if self.signal_model == 'doubleCBPlusGaussian':
      sigma_gauss = ROOT.RooRealVar("sigma_gauss", "sigma_gauss", 0.01, 0.005, 0.15)
      gaussian = ROOT.RooGaussian('gaussian', 'gaussian', mupi_invmass, mean, sigma_gauss)

      # defines the relative importance of gaussian wrt doubleCB
      sigfrac_gauss = ROOT.RooRealVar("sigfrac_gauss","sigfrac_gauss", 0.5, 0.0 ,1.0)

    # build the signal model
    if self.signal_model == 'doubleCB':
      signal_model = ROOT.RooAddPdf("sig", "sig", CBpdf_1, CBpdf_2, sigfrac_CB)

    elif self.signal_model == 'doubleCBPlusGaussian':
      doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", CBpdf_1, CBpdf_2, sigfrac_CB)
      signal_model = ROOT.RooAddPdf('sig', 'sig', doubleCBpdf, gaussian, sigfrac_gauss) # make sure that the model has the same name as in datacard 

    # Define the background model 
    if self.background_model == 'chebychev':
      a0 = ROOT.RooRealVar('a0', 'a0', 0.01, -10, 10)
      background_model = ROOT.RooChebychev('qcd', 'qcd', mupi_invmass, ROOT.RooArgList(a0))

    # getting the observed data #TODO create function
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    tree = ROOT.TChain(treename)
    for data_file in self.data_files:
      tree.Add(data_file.filename) 

    bin_min, bin_max = self.getFitRegion()
    nbins = self.nbins
    hist_name = 'hist'
    hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
    branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    selection = self.selection
    tree.Project(hist_name, branch_name , selection)
    data_obs = ROOT.RooDataHist("data_obs", "data_obs", ROOT.RooArgList(mupi_invmass), hist)

    # import the model in a workspace
    #workspace_filename = '{}/workspace_{}.root'.format(self.outputdir, label)
    workspace_filename = '{}/workspace_{}.root'.format(self.workspacedir, label)
    output_file = ROOT.TFile(workspace_filename, 'RECREATE') # if do_recreate else 'UPDATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')
    getattr(workspace, 'import')(signal_model)
    getattr(workspace, 'import')(background_model)
    getattr(workspace, 'import')(data_obs)
    workspace.Write()
    #workspace.Print()
    #workspace.writeToFile('workspace_{}.root'.format(label))
    output_file.Close()

    print '--> {} created'.format(workspace_filename)

    ## workspace will remain in memory after macro finishes
    ##ROOT.gDirectory.Add(workspace)


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


  def performFit(self, process, label=''):
    '''
      Performs the fit by loading the fit workspace
      process corresponds to either 'signal' or 'background'
    '''

    if process not in ['signal', 'background']:
      raise RuntimeError("[fitter] Unrecognised process. Please choose among ['signal', 'background']")

    # get the label
    if label == '' and process == 'signal': label = self.getSignalLabel() 
    if label == '' and process == 'background': label = self.getBackgroundLabel() 

    # load the workspace
    workspace_filename = '{}/workspace_{}.root'.format(self.workspacedir, label)
    ws_file = self.tools.getRootFile(workspace_filename, with_ext=False)
    ws = ws_file.Get('workspace')
    ws.Print()

    # open the file and get the tree
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    if process == 'signal':
      inputfile = self.tools.getRootFile(self.signal_file.filename, with_ext=False)
      tree = self.tools.getTree(inputfile, treename)
    else:
      tree = ROOT.TChain(treename)
      for data_file in self.data_files:
        tree.Add(data_file.filename) 

    # get signal mass
    signal_mass = self.signal_file.mass

    # declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
    bin_min, bin_max = self.getFitRegion()
    nbins = self.nbins

    # build the binned data
    hist_name = 'hist'
    hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
    branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    if process == 'signal':
      cond = 'ismatched==1 && mu_isdsa==0' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx]==0'
    else:
      cond = 'hnl_pt > 0' # dummy condition
    selection = cond + ' && ' + self.selection
    tree.Project(hist_name, branch_name , selection)
    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(ws.var('mupi_invmass')), hist)

    # create canvas
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    pad2.Draw()

    frame = ws.var('mupi_invmass').frame(ROOT.RooFit.Title(self.title))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    fit_range_min = bin_min #signal_mass - 0.1*signal_mass
    fit_range_max = bin_max #signal_mass + 0.1*signal_mass
    ws.var('mupi_invmass').setRange("peak", fit_range_min, fit_range_max)
    result = ws.pdf('sig').fitTo(rdh, ROOT.RooFit.Range("peak"))

    # plot the fit 		
    pdf_name = 'sig' if process == 'signal' else 'qcd'
    if process == 'signal':
      ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1"),ROOT.RooFit.Components("CBpdf_1"), ROOT.RooFit.LineStyle(ROOT.kDashed))
      ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2"),ROOT.RooFit.Components("CBpdf_2"), ROOT.RooFit.LineStyle(ROOT.kDashed))
      if self.signal_model == 'doubleCBPlusGaussian':
        ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    ws.pdf(pdf_name).plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name(pdf_name), ROOT.RooFit.Components(pdf_name))

    # and write the fit parameters
    ws.pdf(pdf_name).paramOn(frame,   
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # compute the chisquare
    chisquare = frame.chiSquare(pdf_name,"data")

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
    hpull = frame.pullHist()
    for i in range(0, frame.GetNbinsX()):
       hpull.SetPointError(i,0,0,0,0)

    # create a new frame to draw the pull distribution and add the distribution to the frame
    frame2 = ws.var('mupi_invmass').frame(ROOT.RooFit.Title(" "))
    frame2.addPlotable(hpull,"P")#,"E3")

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
      if self.signal_model == 'doubleCB':
        model_label = 'Double Crystal Ball'
      elif self.signal_model == 'doubleCBPlusGaussian':
        model_label = 'Double Crystal Ball + Gaussian'
      leg.AddEntry(frame.findObject(pdf_name), model_label)
      leg.AddEntry(frame.findObject('CBpdf_1'), 'CB_1')
      leg.AddEntry(frame.findObject('CBpdf_2'), 'CB_2')
      if self.signal_model == 'doubleCBPlusGaussian':
        leg.AddEntry(frame.findObject('gaussian'), 'Gaussian')
    else:
      if self.background_model == 'chebychev':
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


  def getPullDistribution(self, process, hpull, label):
    '''
      Create pull distribution and fit it with a Gaussian
    '''
    canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)
    hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)

    bin_min, bin_max = self.getFitRegion()
   
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



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  outdirlabel = 'V10_30Dec21_samples_sel'
  plot_pulls = True
  selection = 'sv_lxy>5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20 '
  signal_model = 'doubleCB'
  signal_files = signal_samples['V10_30Dec21_m3'] 
  background_model = 'chebychev'
  data_files = data_samples['V10_30Dec21']
  do_blind = False
  nbins = 80

  label = 'test'
  for signal_file in signal_files:
    fitter = Fitter(signal_file=signal_file, data_files=data_files, selection=selection, signal_model=signal_model, background_model=background_model, do_blind=do_blind, nbins=150, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
    #fitter.writeSignalModel(label='test')
    #fitter.writeBackgroundModel(label='test')
    #fitter.createFitWorkspace()
    fitter.writeFitModels(label='test')
    #fitter.performFit(process='signal', label='test')
    #fitter.performFit(process='background', label='test')

  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)

  #selection = 'sv_lxy<=1 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10'
  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)

  #selection = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25'
  #fitter = Fitter(data_files=data_files, selection=selection, background_model=background_model, do_blind=do_blind, nbins=nbins, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
  #fitter.performBackgroundFit(mass=3, resolution=0.023)



