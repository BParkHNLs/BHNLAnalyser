import sys
import ROOT
from ROOT import RooFit
import math

from tools import Tools
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


  def performFit(self, signal_file):
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
    cond = 'ismatched==1 && {}'.format(baseline_selection)
    tree.Draw("{}>>hist".format(branch_name), cond)

    mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", binMin, binMax)
    mupi_invmass.setBins(nbins)

    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(mupi_invmass), hist)

    # Define the PDF to fit: 
    # Double sided crystal ball
    # we declare all the parameters needed for the fits 
    mean_init = signal_mass - 0.01*signal_mass
    mean_min = signal_mass - 0.05*signal_mass
    mean_max = signal_mass + 0.05*signal_mass
    mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

    sigma = ROOT.RooRealVar("sigma","sigma", 0.01, 0.005, 0.15)

    alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -1, -3, 3)
    n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 3)
    alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 1, -3, 3)
    n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 3)

    CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mupi_invmass, mean, sigma, alpha_1, n_1)
    CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mupi_invmass, mean, sigma, alpha_2, n_2)

    # defines the relative importance of the two CBs
    sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5, 0.0 ,1.0)

    # we add the two CB pdfs together
    doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", CBpdf_1, CBpdf_2, sigfrac)

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
    result = doubleCBpdf.fitTo(rdh, ROOT.RooFit.Range("peak"))

    # plot the fit    
    doubleCBpdf.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("doubleCBpdf"), ROOT.RooFit.Components("doubleCBpdf"))

    # and write the fit parameters
    doubleCBpdf.paramOn(frame,   
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # we compute the chisquare
    chisquare = frame.chiSquare("doubleCBpdf","data")

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
    label1.AddText('Double sided Crystal Ball')
    label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    label1.AddText('mass {} GeV, ctau {} mm'.format(signal_file.mass, signal_file.ctau))
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
    canv.SaveAs("{}/fit_m{}_ctau{}.png".format(outputdir, str(signal_file.mass).replace('.', 'p'), str(signal_file.ctau).replace('.', 'p')))
    canv.SaveAs("{}/fit_m{}_ctau{}.pdf".format(outputdir, str(signal_file.mass).replace('.', 'p'), str(signal_file.ctau).replace('.', 'p')))

    return sigma.getVal(), (sigma.getAsymErrorHi()+sigma.getAsymErrorLo())/2.


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
      graph.SetLineColor(colours[ictau])
      graph.SetMarkerColor(colours[ictau])
      graph.SetMarkerStyle(20)
      graph.GetXaxis().SetTitle('Signal mass [GeV]')
      graph.GetXaxis().SetLabelSize(0.037)
      graph.GetXaxis().SetTitleSize(0.042)
      graph.GetXaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetTitle('Resolution')
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
      

  def getLifetimeDependencyGraph(self):
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

        sigma, sigma_err = self.performFit(signal_file=signal_file)
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
    graph.GetYaxis().SetTitle('Resolution')
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

        sigma, sigma_err = self.performFit(signal_file=signal_file)
        #resolution_diff[signal_file.mass] = abs(sigma - resolution_average[signal_file.mass]) / sigma 
        resolution_diff[signal_file.mass] += ((sigma - resolution_average[signal_file.mass]) * (sigma - resolution_average[signal_file.mass]))

    for mass in masses:
      resolution_diff[mass] = math.sqrt(1./float(n_ctaus[mass]) * resolution_diff[mass])

    print resolution_diff
      

  def getResolutionFit(self):
    ROOT.gStyle.SetPadLeftMargin(0.15) 
    # create graph resolution vs mass
    masses = self.getMassList()
    
    ctaus = [0.1, 1, 10, 100, 1000]

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

    #graph.Fit('pol1')
    #ROOT.gStyle.SetOptFit() 
    f = ROOT.TF1('f', '0.0002747 + 0.008302*x', 0, 5)
    f.SetLineColor(ROOT.kRed)
    f.SetLineWidth(2)

    canv = self.tools.createTCanvas(name='canv', dimx=800, dimy=700)
    canv.cd()
    graph.Draw('AP')
    f.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots', 'resolution')
    canv.SaveAs("{}/graph_fit.png".format(outputdir))


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  baseline_selection = selection['baseline_08Aug22'].flat

  signal_labels = [
    'V12_08Aug22_training_large',
  ]
  outdirlabel = 'study_resolution_central_V12_08Aug22'

  fitter = Fitter(signal_labels=signal_labels, baseline_selection=baseline_selection, nbins=150, outdirlabel=outdirlabel)
  #fitter.getResolutionGraph()
  #fitter.getResolutionFit()
  #fitter.getLifetimeDependencyGraph()


  signal_labels = [
    'V42_08Aug22_m0p5',
    'V42_08Aug22_m0p6',
    'V42_08Aug22_m0p7',
    'V42_08Aug22_m0p8',
    'V42_08Aug22_m0p9',
    'V42_08Aug22_m1p02',
    'V42_08Aug22_m1p04',
    'V42_08Aug22_m1p06',
    'V42_08Aug22_m1p08',
    'V42_08Aug22_m1p1',
    'V42_08Aug22_m1p12',
    #'V42_08Aug22_m1p14',
    'V42_08Aug22_m1p16',
    'V42_08Aug22_m1p18',
    'V42_08Aug22_m1p2',
    'V42_08Aug22_m1p22',
    'V42_08Aug22_m1p24',
    'V42_08Aug22_m1p26',
    #'V42_08Aug22_m1p28',
    'V42_08Aug22_m1p3',
    'V42_08Aug22_m1p32',
    'V42_08Aug22_m1p34',
    'V42_08Aug22_m1p36',
    'V42_08Aug22_m1p38',
    'V42_08Aug22_m1p4',
    'V42_08Aug22_m1p42',
    'V42_08Aug22_m1p44',
    'V42_08Aug22_m1p46',
    'V42_08Aug22_m1p48',
    'V42_08Aug22_m1p53',
    'V42_08Aug22_m1p56',
    'V42_08Aug22_m1p59',
    #'V42_08Aug22_m1p62',
    'V42_08Aug22_m1p65',
    #'V42_08Aug22_m1p68',
    'V42_08Aug22_m1p71',
    'V42_08Aug22_m1p74',
    'V42_08Aug22_m1p77',
    'V42_08Aug22_m1p8',
    'V42_08Aug22_m1p83',
    #'V42_08Aug22_m1p86',
    #'V42_08Aug22_m1p89',
    'V42_08Aug22_m1p92',
    'V42_08Aug22_m1p95',
    #'V42_08Aug22_m1p98',
    #'V42_08Aug22_m2p05',
    'V42_08Aug22_m2p1',
    'V42_08Aug22_m2p15',
    #'V42_08Aug22_m2p2',
    'V42_08Aug22_m2p25',
    'V42_08Aug22_m2p3',
    'V42_08Aug22_m2p35',
    #'V42_08Aug22_m2p4',
    'V42_08Aug22_m2p45',
    'V42_08Aug22_m2p5',
    'V42_08Aug22_m2p55',
    'V42_08Aug22_m2p6',
    'V42_08Aug22_m2p65',
    'V42_08Aug22_m2p7',
    #'V42_08Aug22_m2p75',
    'V42_08Aug22_m2p8',
    'V42_08Aug22_m2p85',
    'V42_08Aug22_m2p9',
    'V42_08Aug22_m2p95',
    'V42_08Aug22_m3p05',
    'V42_08Aug22_m3p1',
    'V42_08Aug22_m3p15',
    'V42_08Aug22_m3p2',
    'V42_08Aug22_m3p25',
    'V42_08Aug22_m3p3',
    'V42_08Aug22_m3p35',
    'V42_08Aug22_m3p4',
    'V42_08Aug22_m3p45',
    'V42_08Aug22_m3p5',
    'V42_08Aug22_m3p55',
    'V42_08Aug22_m3p6',
    'V42_08Aug22_m3p65',
    'V42_08Aug22_m3p7',
    'V42_08Aug22_m3p75',
    'V42_08Aug22_m3p8',
    'V42_08Aug22_m3p85',
    'V42_08Aug22_m3p9',
    'V42_08Aug22_m3p95',
    'V42_08Aug22_m4p0',
    'V42_08Aug22_m4p1',
    'V42_08Aug22_m4p2',
    #'V42_08Aug22_m4p3',
    #'V42_08Aug22_m4p4',
    'V42_08Aug22_m4p5',
    'V42_08Aug22_m4p6',
    #'V42_08Aug22_m4p7',
    #'V42_08Aug22_m4p8',
  ]
  outdirlabel = 'study_resolution_private_V42_08Aug22'

  fitter = Fitter(signal_labels=signal_labels, baseline_selection=baseline_selection, nbins=150, outdirlabel=outdirlabel)
  #fitter.getResolutionGraph()
  fitter.getResolutionFit()



