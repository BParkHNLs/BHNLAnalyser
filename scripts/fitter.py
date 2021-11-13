import ROOT
from ROOT import RooFit

from tools import Tools

class Fitter(Tools):
  def __init__(self, filename='', file_type='flat', nbins=250, title=' ', outdirlabel='testing', label=''):
    self.tools = Tools()
    self.filename = filename
    self.file_type = file_type
    self.nbins = nbins
    self.title = title
    self.outdirlabel = outdirlabel
    self.label = label

  def performFit(self):
    # open the file and get the tree
    inputfile = ROOT.TFile.Open(self.filename)
    treename = 'signal_tree' if self.file_type == 'flat' else 'Events'
    tree = self.tools.getTree(inputfile, treename)

    # get signal mass
    for entry in tree:
      signal_mass = entry.gen_hnl_mass
      #signal_mass = 3.0
      break 

    # we declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
    binMin = signal_mass - 0.15*signal_mass
    binMax = signal_mass + 0.15*signal_mass
    nbins = self.nbins

    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    #tree->Draw("TMath::Sqrt(E1*E1 - px1*px1 - py1*py1 - pz1*pz1 + E2*E2 - px2*px2 - py2*py2 - pz2*pz2 + 2*E1*E2 - 2*px1*px2 -2*py1*py2 -2*pz1*pz2)>>hist", "pt1>13 && pt2>11")
    branch_name = 'hnl_mass' if self.file_type == 'flat' else 'BToMuMuPi_hnl_mass'
    cond = 'ismatched==1 && mu_isdsa==0' if self.file_type == 'flat' else 'BToMuMuPi_isMatched==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx]==0'
    tree.Draw("{}>>hist".format(branch_name), cond, 'goff')

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

    alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -2, -5, 5)
    n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 5)
    alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 2, -5, 5)
    n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 5)

    CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mupi_invmass, mean, sigma, alpha_1, n_1)
    CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mupi_invmass, mean, sigma, alpha_2, n_2)

    # defines the relative importance of the two CBs
    sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5, 0.0 ,1.0)

    # we add the two CB pdfs together
    #RooAddPdf *signalPdf= new RooAddPdf("signalPdf", "signalPdf", RooArgList(*CBpdf_1, *CBpdf_2), *sigfrac) 
    doubleCBpdf = ROOT.RooAddPdf("doubleCBpdf", "doubleCBpdf", CBpdf_1, CBpdf_2, sigfrac)

    # background
    #RooRealVar *tau = new RooRealVar("tau", "tau", -0.039, -0.045, 0.)
    #tau->setConstant()
    #RooExponential *ExpPdf = new RooExponential("ExpPdf", "ExpPdf", *mupi_invmass, *tau)

    #RooRealVar *a0 = new RooRealVar("a1", "a1", 100, 90, 110)
    #RooRealVar *a1 = new RooRealVar("a0", "a0", -0.05, -10, 0.)
    #RooRealVar *a2 = new RooRealVar("a2", "a2", -0.09, -100, 0.)
    #RooPolynomial *PolyPdf = new RooPolynomial("PolyPdf","PolyPdf", *mupi_invmass ,RooArgList(*a0,*a1))
    #RooRealVar *bkgfrac  = new RooRealVar("bkgfrac","bkgfrac", 0.0 ,1.0)
    #RooAddPdf *ExpPdf= new RooAddPdf("ExpPdf", "ExpPdf", RooArgList(*PolyPdf, *ExpPdf1), *bkgfrac) 

    # model (signal + background)
    #RooRealVar* nsig = new RooRealVar("nsig", "nsig", 200, 0, 1000000)
    #RooRealVar* nbkg = new RooRealVar("nbkg", "nbkg", 500, 0, 1000000)

    #RooRealVar* fsig = new RooRealVar("fsig", "fsig", 0, 1)

    #RooAddPdf* model = new RooAddPdf("model","model" ,RooArgList(*signalPdf,*ExpPdf), *fsig)

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    #pad1->SetLogx()
    #pad1->SetLogy()
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    #pad2->SetLogx()
    pad1.Draw()
    pad2.Draw()

    #RooPlot *frame = mupi_invmass->frame(Title(""))
    frame = mupi_invmass.frame(ROOT.RooFit.Title(self.title))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    #RooFitResult *result
    #result = doubleCBpdf.fitTo(rdh)
    fit_range_min = signal_mass - 0.1*signal_mass
    fit_range_max = signal_mass + 0.1*signal_mass
    mupi_invmass.setRange("peak", fit_range_min, fit_range_max)
    result = doubleCBpdf.fitTo(rdh, ROOT.RooFit.Range("peak"))

    # plot the fit 		
    #model->plotOn(frame,LineColor(2),RooFit::Name("CBpdf_1"),Components("CBpdf_1"))
    #model->plotOn(frame,LineColor(3),RooFit::Name("CBpdf_2"),Components("CBpdf_2"))
    #model->plotOn(frame,LineColor(6),RooFit::Name("signalPdf"),Components("signalPdf"))
    #model->plotOn(frame,LineColor(2),RooFit::Name("ExpPdf"),Components("ExpPdf"), LineStyle(kDashed))
    #model->plotOn(frame,LineColor(4),RooFit::Name("model"), Components("model"))
    doubleCBpdf.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("doubleCBpdf"), ROOT.RooFit.Components("doubleCBpdf"))

    # and write the fit parameters
    #model->paramOn(frame,   
    doubleCBpdf.paramOn(frame,   
         #ROOT.RooFit.Layout(0.67, 0.87, 0.85),
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # we compute the chisquare
    chisquare = frame.chiSquare("doubleCBpdf","data")

    # and print it
    #label1 = ROOT.TPaveText(0.2,0.75,0.3,0.85,"brNDC")
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
    print "chisquare = {}".format(chisquare)

    # we add the legend
    '''
    TLegend* leg = new TLegend(0.2, 0.15, 0.3, 0.35)
    #leg->AddEntry(frame->findObject("CBpdf_1"), "CB_1")
    #leg->AddEntry(frame->findObject("CBpdf_2"), "CB_2")
    #leg->AddEntry(frame->findObject("signalPdf"), "signal")
    leg->AddEntry(frame->findObject("ExpPdf"), "bkg")
    leg->AddEntry(frame->findObject("model"), "signal+bkg")
    leg->SetTextSize(0.03)
    leg->SetTextFont(42)
    leg->SetTextAlign(11)
    leg->SetLineColor(0)
    leg->SetFillColorAlpha(0, 0)
    leg->SetBorderSize(0)
    '''

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
    #RooPlot* frame2 = mupi_invmass->frame(Title(" "))
    frame2 = mupi_invmass.frame(ROOT.RooFit.Title(" "))
    #frame2->addPlotable(hpull,"P")#,"E3")
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
    #leg.Draw()


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
    canv.SaveAs("{}/fit_{}.png".format(outputdir, label))
    canv.SaveAs("{}/fit_{}.pdf".format(outputdir, label))

    #delete hist
    #delete canv

    # additionally, get the pull histogram
    canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)
   
    hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)
   
    for i in range(0, hpull.GetN()):
    #for(Int_t i=0 i<hpull->GetN() i++) {
 
      #Double_t x,point
      #hpull->GetPoint(i,x,point) 
      #x = 0.
      x = ROOT.Double()
      #point = 0.
      point = ROOT.Double()
      hpull.GetPoint(i,x,point) 
      if x<binMin or x>binMax: continue

      #if (x<binMin || x>binMax) continue 

      #cout << x << " " << point << endl
      #hist_pull->Fill(point)
      hist_pull.Fill(point)
   #}

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
    #TF1* fgauss= new TF1("fgauss", "gaus", -5, 5)
    fgauss= ROOT.TF1("fgauss", "gaus", -5, 5)
    fgauss.SetLineColor(2)
    hist_pull.Fit("fgauss")
    fgauss.Draw("same")
    ROOT.gStyle.SetOptFit(0011)
   
    #canv_pull.SaveAs("test.png")  
    canv_pull.SaveAs("{}/pulls_{}.png".format(outputdir, label))
    canv_pull.SaveAs("{}/pulls_{}.pdf".format(outputdir, label))


if __name__ == '__main__':
  #ROOT.gROOT.SetBatch(True)

  #inputfile = ROOT.TFile.Open("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_standardgenmatching.root")
  #TFile* inputfile = TFile::Open("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_updatedgenmatching_mu_0p1_0p25_pi_0p15_0p5_massreldiff_0p1.root")

  #outdirlabel = 'dsamuon_study'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_dsamuons_dr0p25.root'
  #title = 'DSA Muons'
  #label = 'dsa_muons'
  #fitter = Fitter(filename=filename, file_type='nano', nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter = Fitter(filename=filename, file_type='nano', nbins=90, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  outdirlabel = 'Bc_29Jun21'
  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root'
  label = 'm3'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root'
  label = 'm4p5'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  outdirlabel = 'dsa_study'
  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_unresolved_fittedmass_looseselection_originalmatching.root'
  label = 'unresolved_fittedmass_looseselection_originalmatching_nodsa'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_resolved_motherpdgid_unfittedmass.root'
  label = 'resolved_nodsa'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  outdirlabel = 'matching_study_resolved_vs_unresolved'
  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root'
  label = 'resolved_nodsa'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root'
  label = 'unresolved_nodsa'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  outdirlabel = 'check_central_samples'
  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V01_selected/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root'
  label = 'm3_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V01_selected/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root'
  label = 'm4p5_ctau1'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()


  outdirlabel = 'central_samples'

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1_ctau1000'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1_ctau10'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1p5_ctau1000'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1p5_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm1p5_ctau10'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm2_ctau1000'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm2_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm2_ctau10'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm3_ctau1000'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm3_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm3_ctau10'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm3_ctau1'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm4p5_ctau100'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm4p5_ctau10'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm4p5_ctau1'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root'
  label = 'm4p5_ctau0p1'
  fitter = Fitter(filename=filename, nbins=150, outdirlabel=outdirlabel, label=label)
  fitter.performFit()

  #outdirlabel = 'genmatching_comparison'

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root'
  #title = 'Updated gen-matching'
  #label = 'V20_emu_updatedgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_initialgenmatching.root'
  #title = 'Initial gen-matching'
  #label = 'V20_emu_initialgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root'
  #title = 'Updated gen-matching'
  #label = 'V21_updatedgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_initialgenmatching.root'
  #title = 'Initial gen-matching'
  #label = 'V21_initialgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root'
  #title = 'Updated gen-matching'
  #label = 'V26_updatedgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_initialgenmatching.root'
  #title = 'Initial gen-matching'
  #label = 'V26_initialgenmatching'
  #fitter = Fitter(filename=filename, nbins=150, title=title, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()


