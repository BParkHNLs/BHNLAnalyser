import ROOT
from ROOT import RooFit

from tools import Tools
from samples import signal_samples

class Fitter(Tools):
  def __init__(self, signal_file='', selection='', signal_model='', file_type='flat', nbins=250, title=' ', outdirlabel='testing', plot_pulls=True):
    self.tools = Tools()
    self.signal_file = signal_file
    self.selection = selection
    self.signal_model = signal_model
    self.file_type = file_type
    self.nbins = nbins
    self.title = title
    self.outdirlabel = outdirlabel
    self.plot_pulls = plot_pulls

    self.label = 'm{}_ctau_{}'.format(self.signal_file.mass, self.signal_file.ctau).replace('.', 'p') 

    signal_model_list = ['doubleCB', 'doubleCBPlusGaussian']
    if self.signal_model not in signal_model_list:
      raise RuntimeError('Unrecognised signal model "{}". Please choose among {}'.format(self.signal_model, signal_model_list))

    #TODO there are no weights applied so far (incl. ctau, hlt, pu weights)

  def performFit(self):
    # open the file and get the tree
    inputfile = self.tools.getRootFile(self.signal_file.filename, with_ext=False)
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    tree = self.tools.getTree(inputfile, treename)

    # get signal mass
    for entry in tree:
      signal_mass = entry.gen_hnl_mass
      break 
    if signal_mass != self.signal_file.mass:
      raise RuntimeError('The signal mass indicated in Samples does not match the actual signal mass. \n-->Aborting')

    # declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
    binMin = signal_mass - 0.15*signal_mass
    binMax = signal_mass + 0.15*signal_mass
    nbins = self.nbins

    hist_name = 'hist'
    hist = ROOT.TH1D(hist_name, hist_name, nbins, binMin, binMax)
    branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    cond = 'ismatched==1 && mu_isdsa==0' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx]==0'
    signal_selection = cond + ' && ' + self.selection
    tree.Project(hist_name, branch_name , signal_selection)

    mupi_invmass = ROOT.RooRealVar("mupi_invmass","mupi_invmass", binMin, binMax)
    mupi_invmass.setBins(nbins)

    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(mupi_invmass), hist)

    # Define the PDF to fit: 
    # Double sided crystal ball
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
      signal_model = ROOT.RooAddPdf("signal_model", "signal_model", CBpdf_1, CBpdf_2, sigfrac_CB)

    elif self.signal_model == 'doubleCBPlusGaussian':
      doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", CBpdf_1, CBpdf_2, sigfrac_CB)
      signal_model = ROOT.RooAddPdf('signal_model', 'signal_model', doubleCBpdf, gaussian, sigfrac_gauss) 

    # define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    pad2.Draw()

    frame = mupi_invmass.frame(ROOT.RooFit.Title(self.title))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    fit_range_min = signal_mass - 0.1*signal_mass
    fit_range_max = signal_mass + 0.1*signal_mass
    mupi_invmass.setRange("peak", fit_range_min, fit_range_max)
    result = signal_model.fitTo(rdh, ROOT.RooFit.Range("peak"))

    # plot the fit 		
    signal_model.plotOn(frame, ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("CBpdf_1"),ROOT.RooFit.Components("CBpdf_1"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    signal_model.plotOn(frame, ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("CBpdf_2"),ROOT.RooFit.Components("CBpdf_2"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    if self.signal_model == 'doubleCBPlusGaussian':
      signal_model.plotOn(frame, ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("gaussian"),ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    signal_model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("signal_model"), ROOT.RooFit.Components("signal_model"))

    # and write the fit parameters
    signal_model.paramOn(frame,   
         #ROOT.RooFit.Layout(0.67, 0.87, 0.85),
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # compute the chisquare
    chisquare = frame.chiSquare("signal_model","data")

    # and print it
    label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    label1.SetTextSize(0.03)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    label1.AddText('mass {}GeV, ctau {}mm'.format(signal_mass, self.signal_file.ctau))
    qte = '#chi^{2}/ndof'
    label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    print "chisquare = {}".format(chisquare)

    # We define and plot the residuals 		
    # construct a histogram with the pulls of the data w.r.t the curve
    #RooHist* hpull = frame->residHist()
    # pull is same as residuals, but is furthermore divided by the low y error in each bin
    #RooHist* hpull = frame->pullHist()
    hpull = frame.pullHist()
    for i in range(0, frame.GetNbinsX()):
       hpull.SetPointError(i,0,0,0,0)
    #for(int i(0) i<frame->GetNbinsX() ++i)
    #{
    #    hpull->SetPointError(i,0,0,0,0)
    #}

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

    # add the legend
    leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.4, xmax=0.8, ymax=0.6, size=0.03, do_alpha=True)
    if self.signal_model == 'doubleCB':
      model_label = 'Double Crystal Ball'
    elif self.signal_model == 'doubleCBPlusGaussian':
      model_label = 'Double Crystal Ball + Gaussian'
    leg.AddEntry(frame.findObject('signal_model'), model_label)
    leg.AddEntry(frame.findObject('CBpdf_1'), 'CB_1')
    leg.AddEntry(frame.findObject('CBpdf_2'), 'CB_2')
    if self.signal_model == 'doubleCBPlusGaussian':
      leg.AddEntry(frame.findObject('gaussian'), 'Gaussian')
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
    line.DrawLine(binMin,0,binMax,0)
    line.SetLineColor(2)
    line.DrawLine(binMin,-3,binMax,-3)
    line.DrawLine(binMin,3,binMax,3)

    # save output
    canv.cd()
    outputdir = self.tools.getOutDir('./myPlots/fits', self.outdirlabel)
    canv.SaveAs("{}/fit_{}.png".format(outputdir, self.label))
    canv.SaveAs("{}/fit_{}.pdf".format(outputdir, self.label))

    #delete hist
    #delete canv

    if self.plot_pulls:
      # additionally, get the pull histogram
      canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)
     
      hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)
     
      for i in range(0, hpull.GetN()):
        x = ROOT.Double()
        point = ROOT.Double()
        hpull.GetPoint(i,x,point) 
        if x<binMin or x>binMax: continue
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
     
      canv_pull.SaveAs("{}/pulls_{}.png".format(outputdir, self.label))
      canv_pull.SaveAs("{}/pulls_{}.pdf".format(outputdir, self.label))


if __name__ == '__main__':
  #ROOT.gROOT.SetBatch(True)

  outdirlabel = 'V10_30Dec21_samples_sel'
  plot_pulls = False
  selection = 'sv_lxy>5 && trgmu_charge!=mu_charge && trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge)) && hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20 '
  signal_model = 'doubleCBPlusGaussian'
  signal_files = signal_samples['V10_30Dec21_m3'] 

  for signal_file in signal_files:
    fitter = Fitter(signal_file=signal_file, selection=selection, signal_model=signal_model, nbins=150, outdirlabel=outdirlabel, plot_pulls=plot_pulls)
    fitter.performFit()

