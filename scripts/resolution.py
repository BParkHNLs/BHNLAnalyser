import sys
import ROOT
from ROOT import RooFit
import math
from collections import OrderedDict

from tools import Tools
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples
from baseline_selection import selection
from categories import categories


class Fitter(Tools):
  def __init__(self, signal_labels='', baseline_selection='', nbins=250, outdirlabel='testing'):
    self.tools = Tools()
    self.signal_labels = signal_labels
    self.signal_files = []
    for signal_label in self.signal_labels:
      self.signal_files += signal_samples[signal_label]
    self.baseline_selection = baseline_selection
    self.nbins = nbins
    self.outdirlabel = outdirlabel


  def getMassList(self):
    '''
      Get the list of signal masses used in the training
    '''
    masses = []
    for signal_file in self.signal_files:
      mass = signal_file.mass
      if mass not in masses:
        masses.append(mass)

    masses.sort()

    return masses


  def performFit(self, signal_file, category=None, model_name='voigtian', sigma_ini=0.01, alpha_1_ini=-1., alpha_2_ini=1., n_1_ini=3., n_2_ini=3., gamma_ini=0.01, do_fixed=False):
    if model_name not in ['doubleCB', 'voigtian']:
      raise RuntimeError('Model undefined. Choose among (doubleCB, voigtian)')

    # open the file and get the tree
    inputfile = ROOT.TFile.Open(signal_file.filename)
    treename = 'signal_tree'
    tree = self.tools.getTree(inputfile, treename)

    # get signal mass
    signal_mass = signal_file.mass

    # we declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
    binMin = signal_mass - 0.15*signal_mass
    binMax = signal_mass + 0.15*signal_mass
    nbins = self.nbins

    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    branch_name = 'hnl_mass'
    selection = self.baseline_selection
    if category != None:
      selection += ' && {}'.format(category.definition_flat)
    cond = 'ismatched==1 && {}'.format(selection)
    tree.Draw("{}>>hist".format(branch_name), cond)

    mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", binMin, binMax)
    mupi_invmass.setBins(nbins)

    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(mupi_invmass), hist)

    # Define the PDF to fit: 
    # Double sided crystal ball
    # we declare all the parameters needed for the fits 

    # let the mean float by 1 permille (energy scale correction)
    mean_ini = signal_mass
    if not do_fixed:
      mean_min = signal_mass - 0.001 * signal_mass
      mean_max = signal_mass + 0.001 * signal_mass
      mean  = ROOT.RooRealVar("mean","mean", mean_ini, mean_min, mean_max)
    else:
      mean  = ROOT.RooRealVar("mean","mean", mean_ini)

    #if not do_fixed:
    sigma_min = sigma_ini - 0.8 * sigma_ini 
    sigma_max = sigma_ini + 0.8 * sigma_ini 
    sigma = ROOT.RooRealVar("sigma","sigma", sigma_ini, sigma_min, sigma_max)
    #else:
    #  sigma = ROOT.RooRealVar("sigma","sigma", sigma_ini)

    #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -3, -10, 0)
    #n_1 = ROOT.RooRealVar("n_1", "n_1", 4, 0, 10)
    #alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 3, 0, 30)
    #n_2 = n_1 #ROOT.RooRealVar("n_2", "n_2", 4, 0, 10)

    if not do_fixed:
      alpha_1_min = 0
      alpha_1_max = 30
      alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", alpha_1_ini, alpha_1_min, alpha_1_max)
    else:
      alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", alpha_1_ini)

    if not do_fixed:
      alpha_2_min = 0
      alpha_2_max = 30
      alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", alpha_2_ini, alpha_2_min, alpha_2_max)
    else:
      alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", alpha_2_ini)

    if not do_fixed:
      n_1_min = 0
      n_1_max = 30
      n_1 = ROOT.RooRealVar("n_1", "n_1", n_1_ini, n_1_min, n_1_max)
    else:
      n_1 = ROOT.RooRealVar("n_1", "n_1", n_1_ini)

    if not do_fixed:
      n_2_min = 0
      n_2_max = 30
      n_2 = ROOT.RooRealVar("n_2", "n_2", n_2_ini, n_2_min, n_2_max)
    else:
      n_2 = ROOT.RooRealVar("n_2", "n_2", n_2_ini)
      #n_2 = n_1

    if not do_fixed:
      gamma_min = gamma_ini - gamma_ini
      gamma_max = gamma_ini + gamma_ini
      gamma_voigtian = ROOT.RooRealVar("gamma_voigtian", "gamma_voigtian", gamma_ini, gamma_min, gamma_max)
    else:
      gamma_voigtian = ROOT.RooRealVar("gamma_voigtian", "gamma_voigtian", gamma_ini)

    CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mupi_invmass, mean, sigma, alpha_1, n_1)
    CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mupi_invmass, mean, sigma, alpha_2, n_2)

    # defines the relative importance of the two CBs
    if not do_fixed:
      sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.01, 0., 1.)
    else:
      sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5)

    # we add the two CB pdfs together
    if model_name == 'doubleCB':
      #model = ROOT.RooAddPdf("model", "model", CBpdf_1, CBpdf_2, sigfrac)
      model = ROOT.RooDoubleCBFast("model", "model", mupi_invmass, mean, sigma, alpha_1, n_1, alpha_2, n_2)
    elif model_name == 'voigtian':
      model = ROOT.RooVoigtian('model', 'model', mupi_invmass, mean, gamma_voigtian, sigma)

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    pad2.Draw()

    frame = mupi_invmass.frame(ROOT.RooFit.Title(' '))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    fit_range_min = signal_mass - 0.1*signal_mass
    fit_range_max = signal_mass + 0.1*signal_mass
    mupi_invmass.setRange("peak", fit_range_min, fit_range_max)
    result = model.fitTo(rdh, ROOT.RooFit.Range("peak"))

    # plot the fit    
    #if model_name == 'doubleCB':
    #  model.plotOn(frame, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("CBpdf_1"), ROOT.RooFit.Components("CBpdf_1"))
    #  model.plotOn(frame, ROOT.RooFit.LineColor(8), ROOT.RooFit.Name("CBpdf_2"), ROOT.RooFit.Components("CBpdf_2"))
    model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("model"), ROOT.RooFit.Components("model"))

    # and write the fit parameters
    model.paramOn(frame,   
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # we compute the chisquare
    chisquare = frame.chiSquare("model","data")

    # and print it
    label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    label1.SetTextSize(0.03)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    #label1->AddText("Double-Sided CrystalBall PDF")
    #chi2 = to_string(chisquare)
    #label1.AddText('#chi^{2}/ndof = {}'.format(chisquare))
    qte = '#chi^{2}/ndof'
    if model_name == 'doubleCB':
      name = 'Double sided Crystal Ball'
    elif model_name == 'voigtian':
      name = 'Voigtian'
    label1.AddText(name)
    label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    label1.AddText('mass {} GeV, ctau {} mm'.format(signal_file.mass, signal_file.ctau))
    label1.AddText(category.title)
    print "chisquare = {}".format(chisquare)
    print 'sigma = {} +- {}'.format(sigma.getVal(), (sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2.)

    # We define and plot the residuals    
    # construct a histogram with the pulls of the data w.r.t the curve
    hpull = frame.pullHist()
    for i in range(0, frame.GetNbinsX()):
       hpull.SetPointError(i,0,0,0,0)

    # create a new frame to draw the pull distribution and add the distribution to the frame
    frame2 = mupi_invmass.frame(ROOT.RooFit.Title(" "))
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
    line.DrawLine(binMin,0,binMax,0)
    line.SetLineColor(2)
    line.DrawLine(binMin,-3,binMax,-3)
    line.DrawLine(binMin,3,binMax,3)

    # save output
    canv.cd()
    outputdir = self.tools.getOutDir('./myPlots/fits', self.outdirlabel)
    name = 'fit_m{}_ctau{}'.format(str(signal_file.mass).replace('.', 'p'), str(signal_file.ctau).replace('.', 'p'))
    if model_name == 'doubleCB':
      name += '_doubleCB'
    elif model_name == 'voigtian':
      name += '_voigtian'
    if category != None:
      name += '_{}'.format(category.label)
    canv.SaveAs("{}/{}.png".format(outputdir, name))
    canv.SaveAs("{}/{}.pdf".format(outputdir, name))

    parameters = []
    parameters_error = []

    parameters.append(sigma.getVal())
    parameters_error.append((sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2.)
    if model_name == 'doubleCB':
      parameters.append(alpha_1.getVal())
      parameters_error.append((alpha_1.getAsymErrorHi()+alpha_1.getAsymErrorLo())/2.)
      parameters.append(alpha_2.getVal())
      parameters_error.append((alpha_2.getAsymErrorHi()+alpha_2.getAsymErrorLo())/2.)
      parameters.append(n_1.getVal())
      parameters_error.append((n_1.getAsymErrorHi()+n_1.getAsymErrorLo())/2.)
      parameters.append(n_2.getVal())
      parameters_error.append((n_2.getAsymErrorHi()+n_2.getAsymErrorLo())/2.)
    elif model_name == 'voigtian':
      parameters.append(gamma_voigtian.getVal())
      parameters_error.append((gamma_voigtian.getAsymErrorHi()+gamma_voigtian.getAsymErrorLo())/2.)

    return parameters, parameters_error
    #return sigma.getVal(), (sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2.
    #return sigma.getVal(), (sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2., alpha_1.getVal(), (alpha_1.getAsymErrorHi()+alpha_1.getAsymErrorLo())/2.,  alpha_2.getVal(), (alpha_2.getAsymErrorHi()+alpha_2.getAsymErrorLo())/2.,  n_1.getVal(), (n_1.getAsymErrorHi()+n_1.getAsymErrorLo())/2., n_2.getVal(), (n_2.getAsymErrorHi()+n_2.getAsymErrorLo())/2., sigfrac.getVal(), (sigfrac.getAsymErrorHi()+sigfrac.getAsymErrorLo())/2.
    #return sigma.getVal(), (sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2., gamma_voigtian.getVal(), (gamma_voigtian.getAsymErrorHi()+gamma_voigtian.getAsymErrorLo())/2.
    

  def getResolutionGraph(self):
    ROOT.gStyle.SetPadLeftMargin(0.15) 
    # create graph resolution vs mass
    masses = self.getMassList()
    
    ctaus = [0.1, 1, 10, 100, 1000]
    colours = [ROOT.kBlack, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kGreen-2]
    graphs = []

    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    axis_min = -99
    axis_max = -99

    for ictau, ctau in enumerate(ctaus):
      resolution = dict.fromkeys(masses, -1)
      resolution_err = dict.fromkeys(masses, -1)
      graph = ROOT.TGraphAsymmErrors()
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue

        sigma = self.performFit(signal_file=signal_file)[0]
        sigma_err = self.performFit(signal_file=signal_file)[1]
        resolution[signal_file.mass] = sigma
        resolution_err[signal_file.mass] = sigma_err

        if axis_min == -99: axis_min = sigma
        if axis_min != -99 and sigma < axis_min: axis_min = sigma
        if axis_max == -99: axis_max = sigma
        if axis_max != -99 and sigma > axis_max: axis_max = sigma

      for mass in masses:
        point = graph.GetN()
        graph.SetPoint(point, mass, resolution[mass])
        graph.SetPointError(point, 0, 0, resolution_err[mass], resolution_err[mass]) 
      graph.SetLineColor(colours[ictau])
      graph.SetMarkerColor(colours[ictau])
      graph.SetMarkerStyle(20)
      graph.GetXaxis().SetTitle('Signal mass [GeV]')
      graph.GetXaxis().SetLabelSize(0.037)
      graph.GetXaxis().SetTitleSize(0.042)
      graph.GetXaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetTitle('Resolution [GeV]')
      graph.GetYaxis().SetLabelSize(0.037)
      graph.GetYaxis().SetTitleSize(0.042)
      graph.GetYaxis().SetTitleOffset(1.5)
      graphs.append(graph)
      leg.AddEntry(graph, 'c#tau = {} mm'.format(ctau))

    graphs[0].GetYaxis().SetRangeUser(axis_min-(0.2*axis_min), axis_max+(0.2*axis_max))

    canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    canv.cd()

    for igraph, graph in enumerate(graphs):
      if igraph == 0:
        graph.Draw('AP')
      else:
        graph.Draw('P same')
    leg.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots', 'resolution')
    canv.SaveAs("{}/graph.png".format(outputdir))
      

  def getLifetimeDependencyGraph(self, category=None):
    ROOT.gStyle.SetPadLeftMargin(0.15) 
    # create graph resolution vs mass
    masses = self.getMassList()
    
    ctaus = [0.1, 1, 10, 100, 1000]
    colours = [ROOT.kBlack, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kGreen-2]
    graphs = []

    #canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    axis_min = -99
    axis_max = -99

    # first, compute the average resolution per mass

    resolution_average = dict.fromkeys(masses, 0)
    resolution_average_err = dict.fromkeys(masses, 0)
    n_ctaus = dict.fromkeys(masses, 0) 

    for ictau, ctau in enumerate(ctaus):
      graph = ROOT.TGraphAsymmErrors()
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue

        n_ctaus[signal_file.mass] += 1

        sigma, sigma_err = self.performFit(signal_file=signal_file, category=category)
        resolution_average[signal_file.mass] += sigma
        resolution_average_err[signal_file.mass] += sigma_err

        if axis_min == -99: axis_min = sigma
        if axis_min != -99 and sigma < axis_min: axis_min = sigma
        if axis_max == -99: axis_max = sigma
        if axis_max != -99 and sigma > axis_max: axis_max = sigma

    for mass in masses:
      # for a given mass, compute the average sigma over lifetime
      resolution_average[mass] = resolution_average[mass] / float(n_ctaus[mass]) 
      resolution_average_err[mass] = resolution_average_err[mass] / float(n_ctaus[mass]) 

      point = graph.GetN()
      graph.SetPoint(point, mass, resolution_average[mass])
      graph.SetPointError(point, 0, 0, resolution_average_err[mass], resolution_average_err[mass]) 

    graph.SetLineColor(ROOT.kBlack)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetMarkerStyle(20)
    graph.GetXaxis().SetTitle('Signal mass [GeV]')
    graph.GetXaxis().SetLabelSize(0.037)
    graph.GetXaxis().SetTitleSize(0.042)
    graph.GetXaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitle('Resolution [GeV]')
    graph.GetYaxis().SetLabelSize(0.037)
    graph.GetYaxis().SetTitleSize(0.042)
    graph.GetYaxis().SetTitleOffset(1.5)

    graph.GetYaxis().SetRangeUser(axis_min-(0.2*axis_min), axis_max+(0.2*axis_max))

    graph.Fit('pol1')
    ROOT.gStyle.SetOptFit() 

    canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    canv.cd()

    graph.Draw('AP')

    outputdir = self.tools.getOutDir('./myPlots', 'resolution')
    canv.SaveAs("{}/graph_average.png".format(outputdir))

    # then compute the relative difference between a point and the average
    resolution_diff = dict.fromkeys(masses, 0)
    for ictau, ctau in enumerate(ctaus):
      #resolution_diff = dict.fromkeys(masses, 0)
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue

        sigma = self.performFit(signal_file=signal_file)[0]
        sigma_err = self.performFit(signal_file=signal_file)[1]
        #resolution_diff[signal_file.mass] = abs(sigma - resolution_average[signal_file.mass]) / sigma 
        resolution_diff[signal_file.mass] += ((sigma - resolution_average[signal_file.mass]) * (sigma - resolution_average[signal_file.mass]))

    for mass in masses:
      resolution_diff[mass] = math.sqrt(1./float(n_ctaus[mass]) * resolution_diff[mass])

    filename = './myPlots/signal_parametrisation/resolution_diff_{}'.format(category.label)
    f = open(filename, 'w+')
    f.write(resolution_diff)
    f.close()

    print resolution_diff
      

  def getStandardDeviation(self, category=None):
    # create graph resolution vs mass
    masses = self.getMassList()
    
    ctaus = [0.1, 1, 10, 100, 1000]
    n_ctaus = dict.fromkeys(masses, 0) 

    resolution_diff = dict.fromkeys(masses, 0)
    for ictau, ctau in enumerate(ctaus):
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue
        n_ctaus[signal_file.mass] += 1

        sigma, sigma_err = self.performFit(signal_file=signal_file, category=category)
        sigma_predict = 0.0004758 + 0.007316 * signal_file.mass
        resolution_diff[signal_file.mass] += ((sigma - sigma_predict) * (sigma - sigma_predict))

    for mass in masses:
      resolution_diff[mass] = math.sqrt(1./float(n_ctaus[mass]) * resolution_diff[mass])

    filename = './myPlots/signal_parametrisation/resolution_diff_{}.txt'.format(category.label)
    f = open(filename, 'w+')
    f.write('{}'.format(resolution_diff))
    f.close()

    print resolution_diff
      

  def getResolutionFit(self, category=None):
    ROOT.gStyle.SetPadLeftMargin(0.15) 
    # create graph resolution vs mass
    masses = self.getMassList()
    
    #ctaus = [0.1, 1, 10, 100, 1000]
    ctaus = [100]

    axis_min = -99
    axis_max = -99
      
    graph = ROOT.TGraphAsymmErrors()

    for ictau, ctau in enumerate(ctaus):
      resolution = dict.fromkeys(masses, -99)
      resolution_err = dict.fromkeys(masses, -99)
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue

        sigma, sigma_err = self.performFit(signal_file=signal_file)
        resolution[signal_file.mass] = sigma
        resolution_err[signal_file.mass] = sigma_err

        if axis_min == -99: axis_min = sigma
        if axis_min != -99 and sigma < axis_min: axis_min = sigma
        if axis_max == -99: axis_max = sigma
        if axis_max != -99 and sigma > axis_max: axis_max = sigma

      for mass in masses:
        point = graph.GetN()
        graph.SetPoint(point, mass, resolution[mass])
        graph.SetPointError(point, 0, 0, resolution_err[mass], resolution_err[mass]) 

    graph.SetLineColor(ROOT.kBlack)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetMarkerStyle(20)
    graph.GetXaxis().SetTitle('Signal mass [GeV]')
    graph.GetXaxis().SetLabelSize(0.037)
    graph.GetXaxis().SetTitleSize(0.042)
    graph.GetXaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitle('Resolution')
    graph.GetYaxis().SetLabelSize(0.037)
    graph.GetYaxis().SetTitleSize(0.042)
    graph.GetYaxis().SetTitleOffset(1.5)
    graph.GetYaxis().SetRangeUser(axis_min-(0.2*axis_min), axis_max+(0.2*axis_max))

    graph.Fit('pol1')
    ROOT.gStyle.SetOptFit() 
    #f = ROOT.TF1('f', '0.0002747 + 0.008302*x', 0, 5)
    #f.SetLineColor(ROOT.kRed)
    #f.SetLineWidth(2)

    canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    canv.cd()
    graph.Draw('AP')
    #f.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots', 'resolution')
    name = 'graph_fit'
    if category != None:
      name += '_{}'.format(category.label)
    canv.SaveAs("{}/{}.png".format(outputdir, name))


  def createGraph(self, val, val_err, name, title, bin_min, bin_max, category=None):
    ROOT.gStyle.SetPadLeftMargin(0.15) 

    masses = self.getMassList()
    ctaus = [0.1, 1, 10, 100, 1000]
    colours = [ROOT.kBlack, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kGreen-2]

    graphs = []
    graph_forfit = ROOT.TGraphAsymmErrors()
    axis_min = -99
    axis_max = -99

    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)

    for ictau, ctau in enumerate(ctaus):
      graph = ROOT.TGraphAsymmErrors()
      for mass in masses:
        try:
          point = graph.GetN()
          graph.SetPoint(point, mass, val[(mass, ctau)])
          graph.SetPointError(point, 0, 0, val_err[(mass, ctau)], val_err[(mass, ctau)]) 
          point_fit = graph_forfit.GetN()
          graph_forfit.SetPoint(point, mass, val[(mass, ctau)])
          graph_forfit.SetPointError(point, 0, 0, val_err[(mass, ctau)], val_err[(mass, ctau)]) 
        except:
          continue

        if axis_min == -99: axis_min = val[(mass, ctau)]
        if axis_min != -99 and val[(mass, ctau)] < axis_min: axis_min = val[(mass, ctau)]
        if axis_max == -99: axis_max = val[(mass, ctau)]
        if axis_max != -99 and val[(mass, ctau)] > axis_max: axis_max = val[(mass, ctau)]

      graph.SetLineColor(colours[ictau])
      graph.SetMarkerColor(colours[ictau])
      graph.SetMarkerStyle(20)
      graphs.append(graph)
      leg.AddEntry(graph, 'c#tau = {} mm'.format(ctau))

    graph_frame = ROOT.TGraphAsymmErrors()
    graph_frame.SetPoint(0, 0.9, bin_min)
    graph_frame.SetPoint(1, 4.6, bin_max)
    graph_frame.GetXaxis().SetTitle('Signal mass [GeV]')
    graph_frame.GetXaxis().SetLabelSize(0.037)
    graph_frame.GetXaxis().SetTitleSize(0.042)
    graph_frame.GetXaxis().SetTitleOffset(1.1)
    graph_frame.GetYaxis().SetTitle('{}'.format(title))
    graph_frame.GetYaxis().SetLabelSize(0.037)
    graph_frame.GetYaxis().SetTitleSize(0.042)
    graph_frame.GetYaxis().SetTitleOffset(1.5)

    #graph_forfit.Fit('pol1')
    #ROOT.gStyle.SetOptFit() 
    f = ROOT.TF1('f', '0.0004758 + 0.007316*x', 0, 5)
    f.SetLineColor(ROOT.kRed)
    f.SetLineWidth(2)

    canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    canv.cd()
    
    graph_frame.Draw('AP')
    graph_forfit.Draw('P same')
    for igraph, graph in enumerate(graphs):
      graph.Draw('P same')
    leg.Draw('same')
    f.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')
    name = 'graph_{}'.format(name)
    if category != None:
      name += '_{}'.format(category.label)
    canv.SaveAs("{}/{}.png".format(outputdir, name))


  def getParameterGraph(self, category=None):
    ROOT.gStyle.SetPadLeftMargin(0.15) 
    # create graph resolution vs mass
    masses = self.getMassList()
    
    #ctaus = [0.1, 1, 10, 100, 1000]
    ctaus = [100]
    colours = [ROOT.kBlack, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kGreen-2]
    graphs_resolution = []

    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    axis_min = -99
    axis_max = -99

    resolution = OrderedDict()
    resolution_err = OrderedDict()
    alpha1 = OrderedDict()
    alpha1_err = OrderedDict()
    alpha2 = OrderedDict()
    alpha2_err = OrderedDict()
    n1 = OrderedDict()
    n1_err = OrderedDict()
    n2 = OrderedDict()
    n2_err = OrderedDict()
    frac = OrderedDict()
    frac_err = OrderedDict()
    gamma = OrderedDict()
    gamma_err = OrderedDict()

    for ictau, ctau in enumerate(ctaus):
      for signal_file in self.signal_files:
        if signal_file.ctau != ctau: continue

        #sigma, sigma_err, alpha_1, alpha_1_err, alpha_2, alpha_2_err, n_1, n_1_err, n_2, n_2_err, sigfrac, sigfrac_err = self.performFit(signal_file=signal_file, category=category)
        #resolution[(signal_file.mass, signal_file.ctau)] = sigma
        #resolution_err[(signal_file.mass, signal_file.ctau)] = sigma_err
        #alpha1[(signal_file.mass, signal_file.ctau)] = alpha_1
        #alpha1_err[(signal_file.mass, signal_file.ctau)] = alpha_1_err
        #alpha2[(signal_file.mass, signal_file.ctau)] = alpha_2
        #alpha2_err[(signal_file.mass, signal_file.ctau)] = alpha_2_err
        #n1[(signal_file.mass, signal_file.ctau)] = n_1
        #n1_err[(signal_file.mass, signal_file.ctau)] = n_1_err
        #n2[(signal_file.mass, signal_file.ctau)] = n_2
        #n2_err[(signal_file.mass, signal_file.ctau)] = n_2_err
        #frac[(signal_file.mass, signal_file.ctau)] = sigfrac
        #frac_err[(signal_file.mass, signal_file.ctau)] = sigfrac_err

        #sigma, sigma_err, gamma_voigtian, gamma_voigtian_err = self.performFit(signal_file=signal_file, category=category)
        #resolution[(signal_file.mass, signal_file.ctau)] = sigma
        #resolution_err[(signal_file.mass, signal_file.ctau)] = sigma_err
        #gamma[(signal_file.mass, signal_file.ctau)] = gamma_voigtian
        #gamma_err[(signal_file.mass, signal_file.ctau)] = gamma_voigtian_err

        sigma, sigma_err = self.performFit(signal_file=signal_file, category=category)
        resolution[(signal_file.mass, signal_file.ctau)] = sigma
        resolution_err[(signal_file.mass, signal_file.ctau)] = sigma_err

    self.createGraph(val=resolution, val_err=resolution_err, name='resolution', title='Resolution [GeV]', bin_min=0.005, bin_max=0.045, category=category)
    #self.createGraph(val=gamma, val_err=gamma_err, name='gamma', title='#gamma', bin_min=0., bin_max=0.1, category=category)
    #self.createGraph(val=alpha1, val_err=alpha1_err, name='alpha1', title='#alpha_{1}', bin_min=-3, bin_max=3, category=category)
    #self.createGraph(val=alpha2, val_err=alpha2_err, name='alpha2', title='#alpha_{2}', bin_min=-3, bin_max=3, category=category)
    #self.createGraph(val=n1, val_err=n1_err, name='n1', title='n_{1}', bin_min=0, bin_max=10, category=category)
    #self.createGraph(val=n2, val_err=n2_err, name='n2', title='n_{2}', bin_min=0, bin_max=10, category=category)
    #self.createGraph(val=frac, val_err=frac_err, name='frac', title='CB fraction', bin_min=0, bin_max=1, category=category)


  def getCategoryGraph(self, categories):

    hists_resolution = []
    hists_alpha1 = []
    hists_alpha2 = []
    hists_n1 = []
    hists_n2 = []

    resolution_min = 99.
    resolution_max = -99.

    leg_resolution = self.tools.getRootTLegend(xmin=0.2, ymin=0.6, xmax=0.45, ymax=0.9, size=0.045)
    leg_alpha1 = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    leg_alpha2 = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    leg_n1 = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    leg_n2 = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.45, ymax=0.85, size=0.045)
    colours = [ROOT.kBlue, ROOT.kRed, ROOT.kMagenta, ROOT.kOrange]

    for ictau, signal_file in enumerate(self.signal_files):
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau
      #if signal_ctau != 10.: continue

      #if signal_mass < 3:
      hist_resolution = ROOT.TH1D('hist_resolution', 'hist_resolution', len(categories)-6, 0, len(categories)-6)
      hist_alpha1 = ROOT.TH1D('hist_alpha1', 'hist_alpha1', len(categories)-6, 0, len(categories)-6)
      hist_alpha2 = ROOT.TH1D('hist_alpha2', 'hist_alpha2', len(categories)-6, 0, len(categories)-6)
      hist_n1 = ROOT.TH1D('hist_n1', 'hist_n1', len(categories)-6, 0, len(categories)-6)
      hist_n2 = ROOT.TH1D('hist_n2', 'hist_n2', len(categories)-6, 0, len(categories)-6)
      #else:
      #  hist_resolution = ROOT.TH1D('hist_resolution', 'hist_resolution', len(categories), 0, len(categories))
      #  hist_alpha1 = ROOT.TH1D('hist_alpha1', 'hist_alpha1', len(categories), 0, len(categories))
      #  hist_alpha2 = ROOT.TH1D('hist_alpha2', 'hist_alpha2', len(categories), 0, len(categories))
      #  hist_n1 = ROOT.TH1D('hist_n1', 'hist_n1', len(categories), 0, len(categories))
      #  hist_n2 = ROOT.TH1D('hist_n2', 'hist_n2', len(categories), 0, len(categories))

      for icat, category in enumerate(categories):
        #if 'Bc' in category.label and signal_mass < 3: continue
        if 'Bc' in category.label: continue
        #parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='voigtian', sigma_ini=0.001, gamma_ini=0.007, do_fixed=True)
        sigma_ini = 0.007316 * signal_mass + 0.0004758 
        sigma_ini = sigma_ini - 0.35 * sigma_ini
        parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='doubleCB', sigma_ini=sigma_ini, alpha_1_ini=1.5, alpha_2_ini=1.5, n_1_ini=3., n_2_ini=3., do_fixed=True)
        resolution = parameters[0]
        resolution_err = parameters_err[0]
        alpha1 = parameters[1]
        alpha1_err = parameters_err[1]
        alpha2 = parameters[2]
        alpha2_err = parameters_err[2]
        n1 = parameters[3]
        n1_err = parameters_err[3]
        n2 = parameters[4]
        n2_err = parameters_err[4]

        if resolution < resolution_min: resolution_min = resolution
        if resolution > resolution_max: resolution_max = resolution

        hist_resolution.Fill(category.label, resolution)
        hist_resolution.SetBinError(icat+1, resolution_err)
        hist_resolution.SetMarkerStyle(20)
        hist_resolution.SetMarkerColor(colours[ictau])
        hist_resolution.SetLineColor(colours[ictau])
        hist_resolution.SetTitle('#sigma')
        hist_resolution.GetYaxis().SetRangeUser(resolution_min - 0.5*resolution_min, resolution_max + 0.5*resolution_max)
        hists_resolution.append(hist_resolution)
        if icat == 0:
          leg_resolution.AddEntry(hist_resolution, 'c#tau = {} mm'.format(signal_ctau))

        hist_alpha1.Fill(category.label, alpha1)
        hist_alpha1.SetBinError(icat+1, alpha1_err)
        hist_alpha1.SetMarkerStyle(20)
        hist_alpha1.SetMarkerColor(colours[ictau])
        hist_alpha1.SetLineColor(colours[ictau])
        hist_alpha1.SetTitle('#alpha1')
        hists_alpha1.append(hist_alpha1)
        if icat == 0:
          leg_alpha1.AddEntry(hist_alpha1, 'c#tau = {} mm'.format(signal_ctau))

        hist_alpha2.Fill(category.label, alpha2)
        hist_alpha2.SetBinError(icat+1, alpha2_err)
        hist_alpha2.SetMarkerStyle(20)
        hist_alpha2.SetMarkerColor(colours[ictau])
        hist_alpha2.SetLineColor(colours[ictau])
        hist_alpha2.SetTitle('#alpha2')
        hists_alpha2.append(hist_alpha2)
        if icat == 0:
          leg_alpha2.AddEntry(hist_alpha2, 'c#tau = {} mm'.format(signal_ctau))

        hist_n1.Fill(category.label, n1)
        hist_n1.SetBinError(icat+1, n1_err)
        hist_n1.SetMarkerStyle(20)
        hist_n1.SetMarkerColor(colours[ictau])
        hist_n1.SetLineColor(colours[ictau])
        hist_n1.SetTitle('n1')
        hists_n1.append(hist_n1)
        if icat == 0:
          leg_n1.AddEntry(hist_n1, 'c#tau = {} mm'.format(signal_ctau))

        hist_n2.Fill(category.label, n2)
        hist_n2.SetBinError(icat+1, n2_err)
        hist_n2.SetMarkerStyle(20)
        hist_n2.SetMarkerColor(colours[ictau])
        hist_n2.SetLineColor(colours[ictau])
        hist_n2.SetTitle('n2')
        hists_n2.append(hist_n2)
        if icat == 0:
          leg_n2.AddEntry(hist_n2, 'c#tau = {} mm'.format(signal_ctau))

    ROOT.gStyle.SetOptStat(0)
    outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')

    canv = self.tools.createTCanvas(name="canv_resolution", dimx=900, dimy=800)
    for ihist, hist_resolution in enumerate(hists_resolution):
      if ihist == 0:
        hist_resolution.Draw('PE')
      else:
        hist_resolution.Draw('PE same')
    leg_resolution.Draw()
    name = 'graph_resolution_category_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))

    canv = self.tools.createTCanvas(name="canv_alpha1", dimx=900, dimy=800)
    for ihist, hist_alpha1 in enumerate(hists_alpha1):
      if ihist == 0:
        hist_alpha1.Draw('PE')
      else:
        hist_alpha1.Draw('PE same')
    leg_alpha1.Draw()
    name = 'graph_alpha1_category_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))

    canv = self.tools.createTCanvas(name="canv_alpha2", dimx=900, dimy=800)
    for ihist, hist_alpha2 in enumerate(hists_alpha2):
      if ihist == 0:
        hist_alpha2.Draw('PE')
      else:
        hist_alpha2.Draw('PE same')
    leg_alpha2.Draw()
    name = 'graph_alpha2_category_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))

    canv = self.tools.createTCanvas(name="canv_n1", dimx=900, dimy=800)
    for ihist, hist_n1 in enumerate(hists_n1):
      if ihist == 0:
        hist_n1.Draw('PE')
      else:
        hist_n1.Draw('PE same')
    leg_n1.Draw()
    name = 'graph_n1_category_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))
      
    canv = self.tools.createTCanvas(name="canv_n2", dimx=900, dimy=800)
    for ihist, hist_n2 in enumerate(hists_n2):
      if ihist == 0:
        hist_n2.Draw('PE')
      else:
        hist_n2.Draw('PE same')
    leg_n2.Draw()
    name = 'graph_n2_category_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))

        
  def getCategoryAverageCtau(self, categories):
    ROOT.gStyle.SetPadLeftMargin(0.15) 

    hist_resolution = ROOT.TH1D('hist_resolution', 'hist_resolution', len(categories)-7, 0, len(categories)-7)
    graph_resolution = ROOT.TGraphAsymmErrors()

    for icat, category in enumerate(categories):
      if 'incl' in category.label: continue
      if 'Bc' in category.label: continue
      
      resolution = 0
      resolution_err = 0

      for ictau, signal_file in enumerate(self.signal_files):
        signal_mass = signal_file.mass
        signal_ctau = signal_file.ctau
        #if signal_ctau != 10.: continue

        #parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='voigtian', sigma_ini=0.001, gamma_ini=0.007, do_fixed=True)
        sigma_ini = 0.007316 * signal_mass + 0.0004758 
        sigma_ini = sigma_ini - 0.35 * sigma_ini
        parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='doubleCB', sigma_ini=sigma_ini, alpha_1_ini=1.5, alpha_2_ini=1.5, n_1_ini=3., n_2_ini=3., do_fixed=True)
        resolution += parameters[0]
        resolution_err += parameters_err[0]
        print resolution

      print len(self.signal_files)
      resolution = resolution / float(len(self.signal_files))
      resolution_err = resolution_err / float(len(self.signal_files))
      point = graph_resolution.GetN()
      graph_resolution.SetPoint(point, icat, resolution)
      graph_resolution.SetPointError(point, 0, 0, resolution_err, resolution_err)
      graph_resolution.SetMarkerStyle(20)
      graph_resolution.SetMarkerColor(ROOT.kBlack)
      graph_resolution.GetXaxis().SetTitle('category number')
      graph_resolution.GetYaxis().SetTitle('Resolution (averaged on ctau) [GeV]')

    f = ROOT.TF1("f","[0] + 0*x",0,6);
    graph_resolution.Fit('f')

    ROOT.gStyle.SetOptStat(0)
    outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')

    canv = self.tools.createTCanvas(name="canv_resolution", dimx=900, dimy=800)
    #hist_resolution.Draw('PE')
    graph_resolution.Draw('AP')
    name = 'graph_resolution_category_average_ctau_m{}'.format(str(signal_mass).replace('.', 'p'))
    canv.SaveAs("{}/{}.png".format(outputdir, name))


  def getFittedResolutionGraph(self):
    ROOT.gStyle.SetPadLeftMargin(0.15) 

    masses = [1.0, 1.5, 2.0, 3.0, 4.5]
    #resolutions = [8.29939e-03, 1.21634e-02, 1.55410e-02, 2.31232e-02, 3.56260e-02]
    #errors = [2.07059e-05, 2.80112e-05, 4.49008e-05, 9.31097e-05, 2.19555e-04]
    resolutions = [8.75174e-03, 1.25745e-02, 1.63083e-02, 2.42540e-02, 3.72958e-02]
    errors = [1.65701e-05, 2.76289e-05, 4.16503e-05, 8.36789e-05, 1.97366e-04]

    graph = ROOT.TGraphAsymmErrors()
    for i, mass in enumerate(masses):
      point = graph.GetN()
      graph.SetPoint(point, masses[i], resolutions[i])
      graph.SetPointError(point, 0, 0, errors[i], errors[i])
      graph.SetMarkerStyle(20)
      graph.SetMarkerColor(ROOT.kBlack)
      graph.GetXaxis().SetTitle('Signal mass (GeV)')
      graph.GetYaxis().SetTitle('Resolution (GeV)')

    graph.Fit('pol1')
    #ROOT.gStyle.SetOptFit() 

    outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')
    canv = self.tools.createTCanvas(name="canv_resolution", dimx=900, dimy=800)
    graph.Draw('AP')
    name = 'graph_fitted_resolution_mass'
    canv.SaveAs("{}/{}.png".format(outputdir, name))

    
  def getCategoryResolutionGraph(self, categories):
    ROOT.gStyle.SetPadLeftMargin(0.15) 

    masses = self.getMassList()
    colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kOrange+8, ROOT.kGreen-2]

    graphs = []
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.55, xmax=0.45, ymax=0.85, size=0.045)

    resolution_min = 99.
    resolution_max = -99.

    for icat, category in enumerate(categories):
      if 'incl' not in category.label: continue
      if 'Bc' in category.label: continue
      #if '0to50' not in category.label: continue

      #resolution_dict = OrderedDict()
      #resolution_err_dict = OrderedDict()
      graph_resolution = ROOT.TGraphAsymmErrors()

      for mass in masses:
      
        resolution = 0
        resolution_err = 0

        n_files = 0

        for ictau, signal_file in enumerate(self.signal_files):
          signal_mass = signal_file.mass
          if signal_mass != mass: continue
          signal_ctau = signal_file.ctau
          n_files = n_files + 1
          #if signal_ctau != 10.: continue

          #parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='voigtian', sigma_ini=0.001, gamma_ini=0.007, do_fixed=True)
          sigma_ini = 0.007316 * signal_mass + 0.0004758 
          sigma_ini = sigma_ini - 0.35 * sigma_ini
          parameters, parameters_err = self.performFit(signal_file=signal_file, category=category, model_name='doubleCB', sigma_ini=sigma_ini, alpha_1_ini=1.5, alpha_2_ini=1.5, n_1_ini=3., n_2_ini=3., do_fixed=True)
          resolution += parameters[0]
          resolution_err += parameters_err[0]

        #resolution_dict[mass] = resolution / float(len(self.signal_files))
        #resolution_err_dict[mass] = resolution_err / float(len(self.signal_files))
        resolution = resolution / float(n_files)
        resolution_err = resolution_err / float(n_files)
        if resolution < resolution_min: resolution_min = resolution
        if resolution > resolution_max: resolution_max = resolution

        print '\n\n\n resolution'
        print resolution

        point = graph_resolution.GetN()
        graph_resolution.SetPoint(point, mass, resolution)
        graph_resolution.SetPointError(point, 0, 0, resolution_err, resolution_err)
        graph_resolution.SetMarkerStyle(20)
        graph_resolution.SetMarkerColor(colours[icat])
        graph_resolution.GetXaxis().SetTitle('Signal mass (GeV)')
        graph_resolution.GetYaxis().SetTitle('Resolution (averaged on ctau) [GeV]')
        graph_resolution.GetYaxis().SetRangeUser(resolution_min-0.3*resolution_min, resolution_max+0.3*resolution_max)

      graphs.append(graph_resolution)
      leg.AddEntry(graph_resolution, category.title)
  
      #f = ROOT.TF1("f","[0] + 0*x",0,6);
      #graph_resolution.Fit('f')

    ROOT.gStyle.SetOptStat(0)
    outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')

    canv = self.tools.createTCanvas(name="canv_resolution", dimx=900, dimy=800)
    for igraph, graph in enumerate(graphs):
      graph.Fit('pol1')
      if igraph == 0:
        graph.Draw('AP')
      else:
        graph.Draw('P same')
    leg.Draw()
    name = 'graph_resolution_category'
    canv.SaveAs("{}/{}.png".format(outputdir, name))


  def studyYieldsParametrisation(self, categories):
    ROOT.gStyle.SetPadLeftMargin(0.15) 

    masses = self.getMassList()
    #colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-1, ROOT.kOrange+8, ROOT.kGreen-2]

    ctaus = [10.]

    strategy = 'inclusive'
    #weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable'
    weight_hlt = 'weight_hlt_D1'
    weight_pusig = 'weight_pu_sig_D'
    weight_mu0id = 'weight_mu0_softid'
    weight_muid = 'weight_mu_looseid'
    add_weight_hlt = True
    add_weight_pu = True
    add_weight_muid = True

    #graphs = []
    #leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.55, xmax=0.45, ymax=0.85, size=0.045)

    for icat, category in enumerate(categories):
      if 'incl' in category.label: continue
      #if 'Bc' in category.label: continue

      for ctau in ctaus:
        graph = ROOT.TGraph()
        for mass in masses:
          for signal_label in self.signal_labels:
            for signal_file in signal_samples[signal_label]:
              signal_mass = signal_file.mass
              signal_ctau = signal_file.ctau
              if signal_mass != mass: continue
              if signal_ctau != ctau: continue
              if signal_mass < 3 and '_Bc' in category.label: continue

              signal_selection = self.baseline_selection + ' && ismatched==1 &&' + category.definition_flat
              #signal_yields, err = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=41.6, sigma_B=472.8e9, is_bc=category.is_bc, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_pusig=weight_pusig, weight_mu0id=weight_mu0id, weight_muid=weight_muid, strategy=strategy) 
              signal_efficiency = self.tools.getSignalEfficiency(signal_files=[signal_file], selection=signal_selection, is_bc=category.is_bc) 
              print '{} {} {} {}'.format(category.label, mass, ctau, signal_efficiency)
              point = graph.GetN()
              #graph.SetPoint(point, mass, signal_yields)
              graph.SetPoint(point, mass, signal_efficiency)
              graph.SetMarkerStyle(20)


        ROOT.gStyle.SetOptStat(0)
        outputdir = self.tools.getOutDir('./myPlots', 'signal_parametrisation')

        canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)
        graph.Draw('AP')
        #leg.Draw()
        name = 'signal_yields_{}'.format(category.label)
        canv.SaveAs("{}/{}.png".format(outputdir, name))
        canv.SaveAs("{}/{}.C".format(outputdir, name))




if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  baseline_selection = selection['baseline_08Aug22'].flat

  signal_labels = [
    #'V12_08Aug22_training_large',
    'V13_06Feb23_m1',
    'V13_06Feb23_m1p5',
    'V13_06Feb23_m2',
    'V13_06Feb23_m3',
    'V13_06Feb23_m4p5',
    ##'V42_06Feb23_m1p02',
    #'V42_06Feb23_m1p04',
    ##'V42_06Feb23_m1p06',
    ##'V42_06Feb23_m1p08',
    #'V42_06Feb23_m1p1',
    ##'V42_06Feb23_m1p12',
    ###'V42_06Feb23_m1p14',
    ##'V42_06Feb23_m1p16',
    ##'V42_06Feb23_m1p18',
    #'V42_06Feb23_m1p2',
    ##'V42_06Feb23_m1p22',
    ##'V42_06Feb23_m1p24',
    ##'V42_06Feb23_m1p26',
    ###'V42_06Feb23_m1p28',
    #'V42_06Feb23_m1p3',
    ##'V42_06Feb23_m1p32',
    ##'V42_06Feb23_m1p34',
    ##'V42_06Feb23_m1p36',
    ##'V42_06Feb23_m1p38',
    #'V42_06Feb23_m1p4',
    ##'V42_06Feb23_m1p42',
    ##'V42_06Feb23_m1p44',
    ##'V42_06Feb23_m1p46',
    ##'V42_06Feb23_m1p48',
    ##'V42_06Feb23_m1p53',
    ##'V42_06Feb23_m1p56',
    #'V42_06Feb23_m1p59',
    ###'V42_06Feb23_m1p62',
    ##'V42_06Feb23_m1p65',
    ###'V42_06Feb23_m1p68',
    #'V42_06Feb23_m1p71',
    ##'V42_06Feb23_m1p74',
    ##'V42_06Feb23_m1p77',
    #'V42_06Feb23_m1p8',
    ##'V42_06Feb23_m1p83',
    ###'V42_06Feb23_m1p86',
    ###'V42_06Feb23_m1p89',
    #'V42_06Feb23_m1p92',
    ##'V42_06Feb23_m1p95',
    ###'V42_06Feb23_m1p98',
    ###'V42_06Feb23_m2p05',
    #'V42_06Feb23_m2p1',
    ##'V42_06Feb23_m2p15',
    ###'V42_06Feb23_m2p2',
    ##'V42_06Feb23_m2p25',
    #'V42_06Feb23_m2p3',
    ##'V42_06Feb23_m2p35',
    ###'V42_06Feb23_m2p4',
    ##'V42_06Feb23_m2p45',
    #'V42_06Feb23_m2p5',
    ##'V42_06Feb23_m2p55',
    ##'V42_06Feb23_m2p6',
    ##'V42_06Feb23_m2p65',
    #'V42_06Feb23_m2p7',
    ###'V42_06Feb23_m2p75',
    ##'V42_06Feb23_m2p8',
    ##'V42_06Feb23_m2p85',
    ##'V42_06Feb23_m2p9',
    ##'V42_06Feb23_m2p95',
    ##'V42_06Feb23_m3p05',
    #'V42_06Feb23_m3p1',
    ##'V42_06Feb23_m3p15',
    #'V42_06Feb23_m3p2',
    ##'V42_06Feb23_m3p25',
    #'V42_06Feb23_m3p3',
    ##'V42_06Feb23_m3p35',
    #'V42_06Feb23_m3p4',
    ##'V42_06Feb23_m3p45',
    #'V42_06Feb23_m3p5',
    ##'V42_06Feb23_m3p55',
    #'V42_06Feb23_m3p6',
    ##'V42_06Feb23_m3p65',
    #'V42_06Feb23_m3p7',
    ##'V42_06Feb23_m3p75',
    #'V42_06Feb23_m3p8',
    ##'V42_06Feb23_m3p85',
    #'V42_06Feb23_m3p9',
    ##'V42_06Feb23_m3p95',
    #'V42_06Feb23_m4p0',
    #'V42_06Feb23_m4p1',
    #'V42_06Feb23_m4p2',
    ###'V42_06Feb23_m4p3',
    ###'V42_06Feb23_m4p4',
    ##'V42_06Feb23_m4p5',
    #'V42_06Feb23_m4p6',
    ###'V42_06Feb23_m4p7',
    ###'V42_06Feb23_m4p8',
  ]
  outdirlabel = 'study_signal_parametrisation_V13_06Feb23_category_doubleCBFast'

  fitter = Fitter(signal_labels=signal_labels, baseline_selection=baseline_selection, nbins=150, outdirlabel=outdirlabel)
  #fitter.getResolutionGraph()
  #fitter.getResolutionFit()
  #fitter.getLifetimeDependencyGraph()
  #fitter.getParameterGraph()

  categories = categories['categories_0_50_150']
  #fitter.getCategoryGraph(categories=categories)
  #fitter.getCategoryAverageCtau(categories=categories)
  #fitter.getFittedResolutionGraph()
  fitter.getCategoryResolutionGraph(categories=categories)
  #fitter.studyYieldsParametrisation(categories=categories)
  

  #for category in categories:
  #  #fitter.getResolutionFit(category=category)
    #fitter.getParameterGraph(category=category)
    #fitter.getStandardDeviation(category=category)
    

  signal_labels = [
    #'V42_06Feb23_m0p5',
    #'V42_06Feb23_m0p6',
    #'V42_06Feb23_m0p7',
    #'V42_06Feb23_m0p8',
    #'V42_06Feb23_m0p9',
    'V42_06Feb23_m1p02',
    'V42_06Feb23_m1p04',
    'V42_06Feb23_m1p06',
    'V42_06Feb23_m1p08',
    'V42_06Feb23_m1p1',
    'V42_06Feb23_m1p12',
    #'V42_06Feb23_m1p14',
    'V42_06Feb23_m1p16',
    'V42_06Feb23_m1p18',
    'V42_06Feb23_m1p2',
    'V42_06Feb23_m1p22',
    'V42_06Feb23_m1p24',
    'V42_06Feb23_m1p26',
    #'V42_06Feb23_m1p28',
    'V42_06Feb23_m1p3',
    'V42_06Feb23_m1p32',
    'V42_06Feb23_m1p34',
    'V42_06Feb23_m1p36',
    'V42_06Feb23_m1p38',
    'V42_06Feb23_m1p4',
    'V42_06Feb23_m1p42',
    'V42_06Feb23_m1p44',
    'V42_06Feb23_m1p46',
    'V42_06Feb23_m1p48',
    'V42_06Feb23_m1p53',
    'V42_06Feb23_m1p56',
    'V42_06Feb23_m1p59',
    #'V42_06Feb23_m1p62',
    'V42_06Feb23_m1p65',
    #'V42_06Feb23_m1p68',
    'V42_06Feb23_m1p71',
    'V42_06Feb23_m1p74',
    'V42_06Feb23_m1p77',
    'V42_06Feb23_m1p8',
    'V42_06Feb23_m1p83',
    #'V42_06Feb23_m1p86',
    #'V42_06Feb23_m1p89',
    'V42_06Feb23_m1p92',
    'V42_06Feb23_m1p95',
    #'V42_06Feb23_m1p98',
    #'V42_06Feb23_m2p05',
    'V42_06Feb23_m2p1',
    'V42_06Feb23_m2p15',
    #'V42_06Feb23_m2p2',
    'V42_06Feb23_m2p25',
    'V42_06Feb23_m2p3',
    'V42_06Feb23_m2p35',
    #'V42_06Feb23_m2p4',
    'V42_06Feb23_m2p45',
    'V42_06Feb23_m2p5',
    'V42_06Feb23_m2p55',
    'V42_06Feb23_m2p6',
    'V42_06Feb23_m2p65',
    'V42_06Feb23_m2p7',
    #'V42_06Feb23_m2p75',
    'V42_06Feb23_m2p8',
    'V42_06Feb23_m2p85',
    'V42_06Feb23_m2p9',
    'V42_06Feb23_m2p95',
    'V42_06Feb23_m3p05',
    'V42_06Feb23_m3p1',
    'V42_06Feb23_m3p15',
    'V42_06Feb23_m3p2',
    'V42_06Feb23_m3p25',
    'V42_06Feb23_m3p3',
    'V42_06Feb23_m3p35',
    'V42_06Feb23_m3p4',
    'V42_06Feb23_m3p45',
    'V42_06Feb23_m3p5',
    'V42_06Feb23_m3p55',
    'V42_06Feb23_m3p6',
    'V42_06Feb23_m3p65',
    'V42_06Feb23_m3p7',
    'V42_06Feb23_m3p75',
    'V42_06Feb23_m3p8',
    'V42_06Feb23_m3p85',
    'V42_06Feb23_m3p9',
    'V42_06Feb23_m3p95',
    'V42_06Feb23_m4p0',
    'V42_06Feb23_m4p1',
    'V42_06Feb23_m4p2',
    #'V42_06Feb23_m4p3',
    #'V42_06Feb23_m4p4',
    'V42_06Feb23_m4p5',
    'V42_06Feb23_m4p6',
    #'V42_06Feb23_m4p7',
    #'V42_06Feb23_m4p8',
  ]

  #fitter = Fitter(signal_labels=signal_labels, baseline_selection=baseline_selection, nbins=150, outdirlabel=outdirlabel)
  ##fitter.getResolutionGraph()
  ##fitter.getResolutionFit()



