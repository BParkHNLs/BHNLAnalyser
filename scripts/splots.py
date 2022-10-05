import ROOT
from ROOT import RooFit

from tools import Tools
import sys
sys.path.append('../')
from samples import data_samples
from baseline_selection import selection


class SPlotter(Tools):
  def __init__(self, dirname):
    self.tools = Tools()
    self.dirname = dirname


  def getQuantitySet(self):
    '''
      If performing unbinned fit, make sure that all variables used in the selection cuts
      are defined here
    '''
    quantities = [
      ROOT.RooRealVar('mass', 'mass', 0., 13000.),
      ROOT.RooRealVar('pt', 'pt', 0., 13000.),
      ROOT.RooRealVar('probe_pt', 'probe_pt', 0., 13000.),
      ROOT.RooRealVar('lxy', 'lxy', 0., 13000.),
      ROOT.RooRealVar('lxy_sig', 'lxy_sig', 0., 13000.),
      ROOT.RooRealVar('probe_dxy_bs', 'probe_dxy_bs', -13000., 13000.),
      ROOT.RooRealVar('probe_dxy_sig_bs', 'probe_dxy_sig_bs', -13000., 13000.),
      ROOT.RooRealVar('probe_isloose', 'probe_isloose', 0, 1),
    ]

    quantity_set = ROOT.RooArgSet()
    for quantity in quantities:
      ROOT.SetOwnership(quantity, False)
      quantity_set.add(quantity)

    return quantity_set
      

  def computeSWeights(self):
    # get the data file
    #filename = data_samples['tag_and_probe'][0].filename 
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6.root' 
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_nj1.root' 
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_smearing.root' 
    f = ROOT.TFile.Open(filename)
    tree_data = f.Get('tree')

    # define window
    signal_mass = 3.097 

    binMin = 2.9
    binMax = 3.3
    nbins = 80

    # define object to fit
    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    branch_name = 'mass'
    cond = selection['tag_and_probe'].flat + ' && probe_dxy_bs > 0.015' 
    tree_data.Draw("{}>>hist".format(branch_name), cond)

    mass = ROOT.RooRealVar("mass","mass", binMin, binMax)
    mass.setBins(nbins)

    #rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(mass), hist)
    quantity_set = self.getQuantitySet()
    fulldata = ROOT.RooDataSet('fulldata', 'fulldata', tree_data, quantity_set, cond)

    # Define the PDF to fit: 
    # signal: double sided crystal ball
    mean_init = signal_mass - 0.01*signal_mass
    mean_min = signal_mass - 0.05*signal_mass
    mean_max = signal_mass + 0.05*signal_mass
    mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

    sigma = ROOT.RooRealVar("sigma","sigma", 0.024, 0.005, 0.15)

    alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -0.9, -3, 3)
    n_1 = ROOT.RooRealVar("n_1", "n_1", 100)
    alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 0.7, -3, 3)
    n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 100)

    CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mass, mean, sigma, alpha_1, n_1)
    CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mass, mean, sigma, alpha_2, n_2)

    # defines the relative importance of the two CBs
    sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5, 0.0 ,1.0)

    # we add the two CB pdfs together
    double_CB = ROOT.RooAddPdf("double_CB", "double_CB", CBpdf_1, CBpdf_2, sigfrac)

    # background: chebychev
    a0 = ROOT.RooRealVar('a0', 'a0', -0.45, -10, 10)
    chebychev = ROOT.RooChebychev('chebychev', 'chebychev', mass, ROOT.RooArgList(a0))

    # build the full model
    frac_bkg = ROOT.RooRealVar("frac_bkg","frac_bkg", 0.8, 0.0 ,1.0)
    model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(double_CB, chebychev), ROOT.RooArgList(frac_bkg))

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    frame = mass.frame(ROOT.RooFit.Title(""))

    # plot the data
    #rdh.plotOn(frame, ROOT.RooFit.Name("data"))
    fulldata.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(nbins))

    # fit the PDF to the data
    #result = model.fitTo(rdh)
    result = model.fitTo(fulldata)

    # plot the fit    
    model.plotOn(frame, ROOT.RooFit.LineColor(6), ROOT.RooFit.Name("model"), ROOT.RooFit.Components("model"))
    model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name('chebychev'), ROOT.RooFit.Components('chebychev'))

    # and write the fit parameters
    model.paramOn(frame,   
         ROOT.RooFit.Layout(0.1, 0.4, 0.8),
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
    qte = '#chi^{2}/ndof'
    #label1.AddText('Double sided Crystal Ball + Chebychev')
    #label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    print "chisquare = {}".format(chisquare)

    canv.cd()

    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitle("#mu#mu invariant mass [GeV]")
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()
    label1.Draw()

    # save output
    canv.cd()
    outputdir = self.tools.getOutDir('./myPlots/splots', self.dirname)
    canv.SaveAs("{}/fit.png".format(outputdir))
    canv.SaveAs("{}/fit.pdf".format(outputdir))

    ## now going unbinned...
    ##quantity_set = self.getQuantitySet()
    ##fulldata = ROOT.RooDataSet('fulldata', 'fulldata', tree_data, quantity_set, cond)
   
    nsig = ROOT.RooRealVar('nsig', 'nsig', 7000, 0., 10e6)
    nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 800, 0., 10e6)  

    fit_function_extended = ROOT.RooAddPdf('fit_function_extended', 'fit_function_extended', ROOT.RooArgList(double_CB, chebychev), ROOT.RooArgList(nbkg, nsig)) 
    results_data_extended = fit_function_extended.fitTo(fulldata, ROOT.RooFit.Extended(True)) 

    sData = ROOT.RooStats.SPlot('sData', 'An SPlot', fulldata, fit_function_extended, ROOT.RooArgList(nsig, nbkg))

    output = ROOT.TFile.Open('splots.root', 'recreate')
    tree = fulldata.GetClonedTree()
    tree.SetName('tree')
    tree.Write()
    output.Close()



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)

  dirname = 'study_v0'
  splotter = SPlotter(dirname=dirname)
  splotter.computeSWeights() 



