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
      #ROOT.RooRealVar('mass', 'mass', 0., 13000.),
      #ROOT.RooRealVar('pt', 'pt', 0., 13000.),
      #ROOT.RooRealVar('probe_pt', 'probe_pt', 0., 13000.),
      #ROOT.RooRealVar('lxy', 'lxy', 0., 13000.),
      #ROOT.RooRealVar('lxy_sig', 'lxy_sig', 0., 13000.),
      #ROOT.RooRealVar('probe_dxy_bs', 'probe_dxy_bs', -13000., 13000.),
      #ROOT.RooRealVar('probe_dxy_sig_bs', 'probe_dxy_sig_bs', -13000., 13000.),
      #ROOT.RooRealVar('probe_isloose', 'probe_isloose', 0, 1),
      ROOT.RooRealVar('b_mass', 'b_mass', 5., 6.),
      ROOT.RooRealVar('b_pt', 'b_pt', 0., 13000.),
      ROOT.RooRealVar('b_cos2d', 'b_cos2d', 0., 13000.),
      ROOT.RooRealVar('dimu_mass', 'dimu_mass', 0., 13000.),
      ROOT.RooRealVar('sv_lxy', 'sv_lxy', 0., 13000.),
      ROOT.RooRealVar('sv_prob', 'sv_prob', 0., 13000.),
      ROOT.RooRealVar('k_pt', 'k_pt', 0., 13000.),
      ROOT.RooRealVar('k_eta', 'k_eta', -13000., 13000.),
      ROOT.RooRealVar('l2_pt', 'l2_pt', 0., 13000.),
      ROOT.RooRealVar('l2_eta', 'l2_eta', -13000., 13000.),
      ROOT.RooRealVar('l1_softid', 'l1_softid', 0., 13000.),
      ROOT.RooRealVar('sv_lxysig', 'sv_lxysig', 0., 13000.),
      ROOT.RooRealVar('pi_pt', 'pi_pt', 0., 13000.),
      ROOT.RooRealVar('mu_pt', 'mu_pt', 0., 13000.),
      ROOT.RooRealVar('mu0_pt', 'mu0_pt', 0., 13000.),
      ROOT.RooRealVar('hnl_cos2d', 'hnl_cos2d', 0., 13000.),
      ROOT.RooRealVar('mu_pi_mass', 'mu_pi_mass', 0., 13000.),
      ROOT.RooRealVar('pi_dcasig', 'pi_dcasig', 0., 13000.),
      ROOT.RooRealVar('sv_prob', 'sv_prob', 0., 13000.),
      ROOT.RooRealVar('mu0_mu_mass', 'mu0_mu_mass', 0., 13000.),
      ROOT.RooRealVar('mu0_pi_mass', 'mu0_pi_mass', 0., 13000.),
      ROOT.RooRealVar('deltar_mu0_mu', 'deltar_mu0_mu', 0., 13000.),
      ROOT.RooRealVar('deltar_mu0_pi', 'deltar_mu0_pi', 0., 13000.),
      ROOT.RooRealVar('mu0_pfiso03_rel', 'mu0_pfiso03_rel', 0., 13000.),
      ROOT.RooRealVar('mu_pfiso03_rel', 'mu_pfiso03_rel', 0., 13000.),
      ROOT.RooRealVar('pi_numberofpixellayers', 'pi_numberofpixellayers', 0., 13000.),
      ROOT.RooRealVar('pi_numberoftrackerlayers', 'pi_numberoftrackerlayers', 0., 13000.),
      ROOT.RooRealVar('mu_numberofpixellayers', 'mu_numberofpixellayers', 0., 13000.),
      ROOT.RooRealVar('mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', 0., 13000.),
      ROOT.RooRealVar('mu0_numberofpixellayers', 'mu0_numberofpixellayers', 0., 13000.),
      ROOT.RooRealVar('mu0_numberoftrackerlayers', 'mu0_numberoftrackerlayers', 0., 13000.),
    ]

    ['pi_pt', 'mu_pt', 'mu0_pt', 'hnl_cos2d', 'sv_lxysig', 'pi_dcasig', 'sv_prob', 'mu0_mu_mass', 'mu0_pi_mass', 'b_mass', 'deltar_mu0_mu', 'deltar_mu0_pi', 'mu0_pfiso03_rel', 'mu_pfiso03_rel', 'pi_numberofpixellayers', 'pi_numberoftrackerlayers', 'mu_numberofpixellayers', 'mu_numberoftrackerlayers', 'mu0_numberofpixellayers', 'mu0_numberoftrackerlayers']

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
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_smearing.root' 
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_control_pNN.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/Chunk0_n500/flat/flat_bparknano_control_pNN.root'
    f = ROOT.TFile.Open(filename)
    tree_data = f.Get('control_tree')

    # define window
    signal_mass = 5.27934 

    binMin = 5.
    binMax = 6.
    nbins = 100

    # define object to fit
    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    branch_name = 'b_mass'
    #cond = selection['control'].flat #+ ' && probe_dxy_bs > 0.015' 
    cond = 'b_cos2d > 0.9995 && abs(dimu_mass-3.097) < 0.05  && sv_lxy > 0.035  && sv_prob > 0.08  && k_pt > 1.0  && abs(k_eta) < 1.6  && l2_pt > 2.  && abs(l2_eta) < 1.8  && l1_softid==1'
    #cond += ' && sv_lxysig<50'
    #cond += ' && sv_lxysig>50 && sv_lxysig<150'
    cond += ' && sv_lxysig>150'
    tree_data.Draw("{}>>hist".format(branch_name), cond)

    b_mass = ROOT.RooRealVar("b_mass","b_mass", binMin, binMax)
    b_mass.setBins(nbins)

    #rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(mass), hist)
    quantity_set = self.getQuantitySet()
    fulldata = ROOT.RooDataSet('fulldata', 'fulldata', tree_data, quantity_set, cond)

    ## Define the PDF to fit: 
    ## signal: double sided crystal ball
    #mean_init = signal_mass - 0.01*signal_mass
    #mean_min = signal_mass - 0.05*signal_mass
    #mean_max = signal_mass + 0.05*signal_mass
    #mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

    #sigma = ROOT.RooRealVar("sigma","sigma", 0.024, 0.005, 0.15)

    #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -0.9, -3, 3)
    #n_1 = ROOT.RooRealVar("n_1", "n_1", 100)
    #alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 0.7, -3, 3)
    #n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 100)

    #CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", mass, mean, sigma, alpha_1, n_1)
    #CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", mass, mean, sigma, alpha_2, n_2)

    ## defines the relative importance of the two CBs
    #sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5, 0.0 ,1.0)

    ## we add the two CB pdfs together
    #double_CB = ROOT.RooAddPdf("double_CB", "double_CB", CBpdf_1, CBpdf_2, sigfrac)

    ## background: chebychev
    #a0 = ROOT.RooRealVar('a0', 'a0', -0.45, -10, 10)
    #chebychev = ROOT.RooChebychev('chebychev', 'chebychev', mass, ROOT.RooArgList(a0))

    ## build the full model
    #frac_bkg = ROOT.RooRealVar("frac_bkg","frac_bkg", 0.8, 0.0 ,1.0)
    #model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(double_CB, chebychev), ROOT.RooArgList(frac_bkg))

    ## bkg
    mult = 1.
    nbkg_comb = ROOT.RooRealVar('nbkg_comb','number of bkg events (combinatorial)',            1000*mult,0,1.0E07) 
    nbkg_prec = ROOT.RooRealVar('nbkg_prec','number of bkg events (partially reconstructed)', 40000*mult,0,1.0E07) 
    
    # crystal ball
    # initial (a la MG) values
    #cb_mean   = ROOT.RooRealVar('cb_mean','cb_mean',  5.0, 4.9, 5.1)
    #cb_sigma  = ROOT.RooRealVar('cb_sigma','cb_sigma', 0.1, 0., 0.3) 
    #cb_alpha  = ROOT.RooRealVar('cb_alpha', 'cb_alpha', 0.5, 0., 2.) # alpha: where the Gaussian and exponential connect, at alpha * sigma of the Gaussian (pos)
    #cb_n      = ROOT.RooRealVar('cb_n', 'cb_n', 1., 0., 10.)         # the slope of the exponential 
    #cb        = ROOT.RooCBShape('cb', 'cb', mass, cb_mean, cb_sigma, cb_alpha, cb_n)

    cb_mean   = ROOT.RooRealVar('cb_mean','cb_mean',  5.0, 4.9, 5.1)
    cb_sigma  = ROOT.RooRealVar('cb_sigma','cb_sigma', 0.5, 0.01, 0.1) 
    cb_alpha  = ROOT.RooRealVar('cb_alpha', 'cb_alpha', 0.001, 0., 2.) # alpha: where the Gaussian and exponential connect, at alpha * sigma of the Gaussian (pos)
    cb_n      = ROOT.RooRealVar('cb_n', 'cb_n', 0.1, 0., 10.)         # the slope of the exponential 
    cb        = ROOT.RooCBShape('cb', 'cb', b_mass, cb_mean, cb_sigma, cb_alpha, cb_n)
    
    ### 2) combinatorial
    #### exponential
    #exp_a = ROOT.RooRealVar('exp_a', 'exp_a', -0.8, -0.9, -0.7)
    exp_a = ROOT.RooRealVar('exp_a', 'exp_a', -0.0001, -5.,0.)
    exp = ROOT.RooExponential('exp', 'Exponential', b_mass, exp_a)
    
    ## signal
    nsig  = ROOT.RooRealVar('nsig','number of signal events', 20000*mult,0,1E07)
   
    ### RooVoigtian is an efficient implementation of the convolution of a Breit-Wigner with a Gaussian
    #voigt_mean  = ROOT.RooRealVar('voigt_mean', 'voigt_mean', 5.28, 5.26, 5.3)
    voigt_mean  = ROOT.RooRealVar('voigt_mean', 'voigt_mean', 5.27934, 5.26, 5.29)
    #voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.02, 0., 0.06)
    voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.015, 0.005, 0.03)
    #voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.002, 0., 0.06)
    voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.028, 0.005, 0.05)
    voigt       = ROOT.RooVoigtian('voigt', 'voigt', b_mass, voigt_mean, voigt_sigma, voigt_width)
    
    # choose signal model
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','J/#psi K signal',ROOT.RooArgList(voigt),ROOT.RooArgList(nsig))
    
    # choose bkg model for the full fit
    fitmodel_bkg_comb = ROOT.RooAddPdf('fitmodel_bkg_comb', 'Combinatorial bkg', ROOT.RooArgList(exp), ROOT.RooArgList(nbkg_comb))
    fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(cb), ROOT.RooArgList(nbkg_prec))
        
    # model
    #model = ROOT.RooAddPdf('model', 'signal + bkg(comb) + bkg(prec)', ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb,fitmodel_bkg_prec), ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec)) 
    #model = ROOT.RooAddPdf('model','model',ROOT.RooArgList(voigt),ROOT.RooArgList(nsig))
    #model = ROOT.RooAddPdf('model', 'Combinatorial bkg', ROOT.RooArgList(exp), ROOT.RooArgList(nbkg_comb))
    #model = ROOT.RooAddPdf('model', 'signal + bkg(comb) + bkg(prec)', ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb), ROOT.RooArgList(nsig,nbkg_comb)) 
    #model = ROOT.RooAddPdf('model', 'Partially reco bkg', ROOT.RooArgList(cb), ROOT.RooArgList(nbkg_prec))
    #model = ROOT.RooAddPdf('model', 'signal + bkg(comb) + bkg(prec)', ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb,fitmodel_bkg_prec), ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec)) 
    model = ROOT.RooAddPdf('model',
                              'signal + bkg(comb) + bkg(prec)',
                              ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb,fitmodel_bkg_prec),
                              ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec)) 

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    frame = b_mass.frame(ROOT.RooFit.Title(" "))

    # plot the data
    #rdh.plotOn(frame, ROOT.RooFit.Name("data"))
    fulldata.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(nbins))

    # fit the PDF to the data
    #result = model.fitTo(rdh)
    #result = model.fitTo(fulldata)
    results = model.fitTo(fulldata, ROOT.RooFit.Extended(True), ROOT.RooFit.Save()) 
    #results = model.fitTo(fulldata, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(), ROOT.RooFit.Verbose(True)) 

    # plot the fit    
    #model.plotOn(frame, ROOT.RooFit.LineColor(6), ROOT.RooFit.Name("model"), ROOT.RooFit.Components("model"))
    #model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name('chebychev'), ROOT.RooFit.Components('chebychev'))

    model.plotOn(frame, ROOT.RooFit.Components('fitmodel_signal'),ROOT.RooFit.LineColor(ROOT.kOrange+7), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name('signal'))
    model.plotOn(frame, ROOT.RooFit.Components('fitmodel_bkg_comb'),ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name('bkg_comb'))
    model.plotOn(frame, ROOT.RooFit.Components('fitmodel_bkg_prec'),ROOT.RooFit.LineColor(ROOT.kGreen-2), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name('bkg_prec'))
    #model.plotOn(frame, ROOT.RooFit.Components('model'), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name('model'))

    model.plotOn(frame, ROOT.RooFit.Components('model'), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name('model'))

    # and write the fit parameters
    model.paramOn(frame,   
         ROOT.RooFit.Layout(0.95, 0.4, 0.8),
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
    #frame.GetXaxis().SetTitle("#mu#mu invariant mass [GeV]")
    frame.GetXaxis().SetTitle("K#mu#mu invariant mass (GeV)")
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
   
    #nsig = ROOT.RooRealVar('nsig', 'nsig', 7000, 0., 10e6)
    #nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 800, 0., 10e6)  

    #fit_function_extended = ROOT.RooAddPdf('fit_function_extended', 'fit_function_extended', ROOT.RooArgList(double_CB, chebychev), ROOT.RooArgList(nbkg, nsig)) 
    #results_data_extended = fit_function_extended.fitTo(fulldata, ROOT.RooFit.Extended(True)) 

    #sData = ROOT.RooStats.SPlot('sData', 'An SPlot', fulldata, fit_function_extended, ROOT.RooArgList(nsig, nbkg))
    sData = ROOT.RooStats.SPlot('sData', 'An SPlot', fulldata, model, ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec)) 

    output = ROOT.TFile.Open('splots.root', 'recreate')
    tree = fulldata.GetClonedTree()
    tree.SetName('control_tree')
    tree.Write()
    output.Close()



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)

  dirname = 'study_control'
  splotter = SPlotter(dirname=dirname)
  splotter.computeSWeights() 



