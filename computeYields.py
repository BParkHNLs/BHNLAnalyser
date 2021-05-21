import os
import os.path
from os import path
import ROOT
import math
from common import PlottingTools
from samples import data_samples, data_samples_small, qcd_samples, qcd_samples, signal_samples
from quantity import Quantity
#from decays import Decays, HNLDecays


class ComputeYields(PlottingTools):
  def __init__(self, data_file='', qcd_files='', signal_file='', selection='', white_list=''):
    self.data_file = data_file
    self.qcd_files = qcd_files
    self.signal_file = signal_file
    self.selection = selection
    self.white_list = white_list


  #def getHistoMC(self, quantity, selection):
  #  hist_mc_tot = ROOT.TH1D('hist_mc_tot', 'hist_mc_tot', quantity.nbins, quantity.bin_min, quantity.bin_max)
  #  hist_mc_tot.Sumw2()

  #  for ifile, qcd_file in enumerate(self.qcd_files):
  #    if qcd_file.label not in self.white_list: continue

  #    f_mc = ROOT.TFile.Open(qcd_file.filename, 'READ')
      
  #    #weight_mc = PlottingTools.computeQCDMCWeight(self, PlottingTools.getTree(self, f_mc, 'signal_tree'), qcd_file.cross_section, qcd_file.filter_efficiency)
  #    weight_mc = PlottingTools.computeQCDMCWeight(self, f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
  #    hist_mc = PlottingTools.createHisto(self, f_mc, 'signal_tree', quantity, branchname='flat', selection=selection, weight=weight_mc) 
  #    hist_mc.Sumw2()
    
  #    hist_mc_tot.Add(hist_mc)

  #  return hist_mc_tot

  '''
  def gamma_total(self, mass,vv): # taken from common.py
      gamma_total =  HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['tot']   # GeV
      return gamma_total


  def ctau_from_gamma(self, gamma): # taken from common.py
      tau_natural = 1. / gamma                  # 1/GeV
      tau = tau_natural * const_hbar            # s
      ctau = tau * const_c * 1000               # mm
      return ctau


  def getVV(self, mass=-99.,ctau=-99.,ismaj=True): # taken from common.py  
    mult = 2. if ismaj else 1. 
    ref_m = 1. # GeV 
    ref_vv = 1.  
    ref_ctau = ctau_from_gamma(mult*gamma_total(mass=ref_m,vv=ref_vv)) 
 
    k = ref_ctau * np.power(ref_m,5) * ref_vv 
     
    return k/(np.power(mass, 5) * ctau) 


  #def getHNLBRs(self, mass):
  #    #HNL BR tot = sum_B_species(HNL BR)

  #  dec = Decays(mass=mass, mixing_angle_square=1)
  #'''


  def computeCrossRatioFromQCDMC(self):
    '''
      Study the correlation of two variables computing the ratio
      r = N_A * ND / N_B * N_C
    '''
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=5.4)
    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge==0')
    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && hnl_charge==0')
    #hist_mc_tot_A = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge==0')
    hist_mc_tot_A = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && sv_prob>0.05')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge!=0')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && hnl_charge!=0')
    #hist_mc_tot_B = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge!=0')
    hist_mc_tot_B = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && sv_prob<0.05')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge==0')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && hnl_charge==0')
    #hist_mc_tot_C = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge==0')
    hist_mc_tot_C = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && sv_prob>0.05')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge!=0')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && hnl_charge!=0')
    #hist_mc_tot_D = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge!=0')
    hist_mc_tot_D = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && sv_prob<0.05')

    n_obs_mc_A = hist_mc_tot_A.GetBinContent(1)
    n_obs_mc_B = hist_mc_tot_B.GetBinContent(1)
    n_obs_mc_C = hist_mc_tot_C.GetBinContent(1)
    n_obs_mc_D = hist_mc_tot_D.GetBinContent(1)

    n_err_mc_A = hist_mc_tot_A.GetBinError(1) #math.sqrt(n_obs_mc_A) 
    n_err_mc_B = hist_mc_tot_B.GetBinError(1) #math.sqrt(n_obs_mc_B) 
    n_err_mc_C = hist_mc_tot_C.GetBinError(1) #math.sqrt(n_obs_mc_C) 
    n_err_mc_D = hist_mc_tot_D.GetBinError(1) #math.sqrt(n_obs_mc_D) 

    ratio = n_obs_mc_A*n_obs_mc_D/(n_obs_mc_B*n_obs_mc_C)
    err = ratio * (n_err_mc_A/n_obs_mc_A + n_err_mc_B/n_obs_mc_B + n_err_mc_C/n_obs_mc_C + n_err_mc_D/n_obs_mc_D)

    return ratio, err


  def computeApproxWeightQCDMCtoData(self, quantity, selection):
    '''
      weight = lumi_data / lumi_mc = N_data * sigma_mc / (N_mc * sigma_data) estimated as N_data / N_mc
    '''

    f_data = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+self.data_file.filename, 'READ')
    hist_data = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=selection)
    hist_data.Sumw2()
    n_obs_data = hist_data.GetBinContent(1)
    n_err_data = math.sqrt(n_obs_data) #hist_data.GetBinError(1)

    hist_mc_tot = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=selection)
    n_obs_mc = hist_mc_tot.GetBinContent(1)
    n_err_mc = math.sqrt(n_obs_mc) #hist_mc_tot.GetBinError(1)

    weight = n_obs_data / n_obs_mc

    if n_obs_data != 0 and n_obs_mc != 0:
      err = weight* (n_err_data / n_obs_data + n_err_mc / n_obs_mc)
    else: 
      err = 0

    return weight, err


  def computeWeightQCDMCtoData(self, lumi_data):
    '''
      weight = lumi_data / lumi_mc
    '''

    # get lumi mc
    lumi_mc = 0

    for ifile, qcd_file in enumerate(self.qcd_files):
      if qcd_file.label not in self.white_list: continue

      f_mc = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+qcd_file.filename, 'READ')
      
      weight = PlottingTools.computeQCDMCWeight(self, PlottingTools.getTree(self, f_mc, 'signal_tree'), qcd_file.cross_section, qcd_file.filter_efficiency)
      weight_mc = 1./weight
      lumi_mc += weight_mc

    weight = lumi_data / lumi_mc 
    err = 0.

    return weight, err


  def computeBkgYieldsFromABCDData(self):
    '''
      Estimate background yields from data using the ABCD method
      A = b_mass < 6.27 && hnl_charge == 0 (SR)
      B = b_mass < 6.27 && hnl_charge != 0 
      C = b_mass > 6.27 && hnl_charge == 0 
      D = b_mass > 6.27 && hnl_charge != 0 

      N_A = N_B * N_C/N_D
    '''
    f_data = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+self.data_file.filename, 'READ')
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=1000)
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=2.83268, bin_max=3.13932)

    bin_selection = self.selection

    #hist_data_B = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && b_mass<6.27 && hnl_charge!=0')
    #hist_data_C = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && b_mass>6.27 && hnl_charge==0')
    #hist_data_D = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && b_mass>6.27 && hnl_charge!=0')
    hist_data_B = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && hnl_cos2d>0.993 && sv_prob<0.05')
    hist_data_C = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && hnl_cos2d<0.993 && sv_prob>0.05')
    hist_data_D = PlottingTools.createHisto(self, f_data, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && hnl_cos2d<0.993 && sv_prob<0.05')

    #hist_data.Sumw2()
    n_obs_data_B = hist_data_B.GetBinContent(1)
    n_obs_data_C = hist_data_C.GetBinContent(1)
    n_obs_data_D = hist_data_D.GetBinContent(1)

    n_err_data_B = math.sqrt(hist_data_B.GetBinContent(1))
    n_err_data_C = math.sqrt(hist_data_C.GetBinContent(1))
    n_err_data_D = math.sqrt(hist_data_D.GetBinContent(1))

    n_obs_data_A = n_obs_data_B * (n_obs_data_C / n_obs_data_D) 
    n_err_data_A = n_obs_data_A * (n_err_data_B / n_obs_data_B + n_err_data_C / n_obs_data_C + n_err_data_D / n_obs_data_D)

    return int(n_obs_data_A), int(n_err_data_A)


  def computeBkgYieldsFromABCDQCDMC(self):
    '''
      Estimate background yields from data using the ABCD method
      A = b_mass < 6.27 && hnl_charge == 0 (SR)
      B = b_mass < 6.27 && hnl_charge != 0 
      C = b_mass > 6.27 && hnl_charge == 0 
      D = b_mass > 6.27 && hnl_charge != 0 

      N_A = N_B * N_C/N_D
    '''
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=1000)

    bin_selection = self.selection
    hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge==0')
    hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge!=0')
    hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge==0')
    hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge!=0')

    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection='b_mass<6.27 && hnl_charge==0')
    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection='sv_prob>0.05 && hnl_charge==0')
    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection='b_mass<6.27 && hnl_mass<5.4')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection='b_mass<6.27 && hnl_charge!=0')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection='sv_prob>0.05 && hnl_charge!=0')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection='b_mass<6.27 && hnl_mass>5.4')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection='b_mass>6.27 && hnl_charge==0')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection='sv_prob<0.05 && hnl_charge==0')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection='b_mass>6.27 && hnl_mass<5.4')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection='b_mass>6.27 && hnl_charge!=0')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection='sv_prob<0.05 && hnl_charge!=0')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection='b_mass>6.27 && hnl_mass>5.4')

    n_obs_mc_A = hist_mc_tot_A.GetBinContent(1)
    n_obs_mc_B = hist_mc_tot_B.GetBinContent(1)
    n_obs_mc_C = hist_mc_tot_C.GetBinContent(1)
    n_obs_mc_D = hist_mc_tot_D.GetBinContent(1)

    n_err_mc_A = hist_mc_tot_A.GetBinError(1) #math.sqrt(n_obs_mc_A) 
    n_err_mc_B = hist_mc_tot_B.GetBinError(1) #math.sqrt(n_obs_mc_B) 
    n_err_mc_C = hist_mc_tot_C.GetBinError(1) #math.sqrt(n_obs_mc_C) 
    n_err_mc_D = hist_mc_tot_D.GetBinError(1) #math.sqrt(n_obs_mc_D) 

    n_obs_mc_A = n_obs_mc_B * (n_obs_mc_C / n_obs_mc_D) 
    n_err_mc_A = n_obs_mc_A * (n_err_mc_B / n_obs_mc_B + n_err_mc_C / n_obs_mc_C + n_err_mc_D / n_obs_mc_D)

    n_real_mc_A = hist_mc_tot_A.GetBinContent(1)
    n_realerr_mc_A = hist_mc_tot_A.GetBinError(1)

    print 'ABCD estimation SR: {} +- {}'.format(int(n_obs_mc_A), int(n_err_mc_A))
    print 'True SR: {} +- {}'.format(int(n_real_mc_A), int(n_realerr_mc_A))

    #return int(n_obs_mc_A), int(n_err_mc_A)


  def validateABCDOnQCDMC(self, plot_name='closure', label=''):
    '''
      Get mupi mass ditribution from data using the ABCD method
      A = b_mass < 6.27 && hnl_charge == 0 (SR)
      B = b_mass < 6.27 && hnl_charge != 0 
      C = b_mass > 6.27 && hnl_charge == 0 
      D = b_mass > 6.27 && hnl_charge != 0 

      and compare it to the actual distribution in the SR
    '''
    #quantity = Quantity(name_flat='hnl_mass', nbins=100, bin_min=0, bin_max=5.6)
    quantity = Quantity(name_flat='hnl_mass', nbins=50, bin_min=0, bin_max=5.37)

    bin_selection = self.selection
    #hist_mc_tot_A = self.getHistoMC(quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && sv_prob>0.05')
    #hist_mc_tot_B = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge!=0')
    #hist_mc_tot_C = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge==0')
    #hist_mc_tot_D = self.getHistoMC(quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge!=0')

    #hist_mc_tot_A = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge==0')
    #hist_mc_tot_B = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass<6.27 && hnl_charge!=0')
    #hist_mc_tot_C = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge==0')
    #hist_mc_tot_D = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && b_mass>6.27 && hnl_charge!=0')

    hist_mc_tot_A = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && sv_prob>0.05')
    hist_mc_tot_B = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d>0.993 && sv_prob<0.05')
    hist_mc_tot_C = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && sv_prob>0.05')
    hist_mc_tot_D = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && hnl_cos2d<0.993 && sv_prob<0.05')

    #hist_estimate = ROOT.TH1D('hist_estimate', 'hist_estimate', 100, 0, 5.6)
    hist_estimate = ROOT.TH1D('hist_estimate', 'hist_estimate', 50, 0, 5.37)
    for ibin in range(1, hist_mc_tot_B.GetNbinsX()+1):
      n_obs_mc_B = hist_mc_tot_B.GetBinContent(ibin)
      n_obs_mc_C = hist_mc_tot_C.GetBinContent(ibin)
      n_obs_mc_D = hist_mc_tot_D.GetBinContent(ibin)
      n_err_mc_B = hist_mc_tot_B.GetBinError(ibin) 
      n_err_mc_C = hist_mc_tot_C.GetBinError(ibin) 
      n_err_mc_D = hist_mc_tot_D.GetBinError(ibin) 
      content = n_obs_mc_B * (n_obs_mc_C / n_obs_mc_D) if n_obs_mc_C != 0 and n_obs_mc_D != 0 else 0
      err = content * (n_err_mc_B / n_obs_mc_B + n_err_mc_C / n_obs_mc_C + n_err_mc_D / n_obs_mc_D) if n_obs_mc_B!=0 and n_obs_mc_C!=0 and n_obs_mc_D!=0 else 0 
      hist_estimate.SetBinContent(ibin, content)
      hist_estimate.SetBinError(ibin, err)

    ROOT.gStyle.SetOptStat(0)

    canv = PlottingTools.getTCanvas(self, 'canv', 800, 900)

    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
    pad_up.SetBottomMargin(0.03)
    pad_up.Draw()
    canv.cd()
    pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
    pad_down.SetBottomMargin(0.25)
    pad_down.Draw()

    hist_estimate.SetLineWidth(3)
    hist_estimate.SetLineColor(ROOT.kOrange)
    #hist_estimate.SetTitle('Closure of the ABCD method in the Signal Region')
    hist_estimate.SetTitle('')
    #hist_estimate.GetXaxis().SetTitle('#mu#pi invariant mass [GeV]')
    hist_estimate.GetXaxis().SetLabelSize(0.0)
    hist_estimate.GetXaxis().SetTitleSize(0.0)
    #hist_estimate.GetXaxis().SetTitleOffset(1.1)
    hist_estimate.GetYaxis().SetTitle('Entries')
    hist_estimate.GetYaxis().SetLabelSize(0.037)
    hist_estimate.GetYaxis().SetTitleSize(0.042)
    hist_estimate.GetYaxis().SetTitleOffset(1.1)
    #hist_estimate.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_estimate, hist_mc_stack, do_log))

    hist_mc_tot_A.SetLineWidth(3)
    hist_mc_tot_A.SetLineColor(ROOT.kBlue)

    int_estimate = hist_estimate.Integral()
    hist_estimate.Scale(1/int_estimate)
    int_A = hist_mc_tot_A.Integral()
    hist_mc_tot_A.Scale(1/int_A)

    pad_up.cd()
    hist_estimate.Draw('histo')
    hist_mc_tot_A.Draw('histo same')

    legend = PlottingTools.getRootTLegend(self, xmin=0.55, ymin=0.65, xmax=0.8, ymax=0.85, size=0.043)
    legend.AddEntry(hist_estimate, 'QCD ABCD estimate')
    legend.AddEntry(hist_mc_tot_A, 'QCD true')
    legend.Draw()

    label_box = ROOT.TPaveText(0.65,0.4,0.8,0.6, "brNDC")
    label_box.SetBorderSize(0)
    label_box.SetFillColor(ROOT.kWhite)
    label_box.SetTextSize(0.042)
    label_box.SetTextAlign(11)
    label_box.SetTextFont(42)
    label_box.AddText(label)
    label_box.Draw()

    pad_down.cd()

    hist_ratio = PlottingTools.getRatioHistogram(self, hist_estimate, hist_mc_tot_A)
    hist_ratio.Sumw2()

    for ibin in range(0, hist_ratio.GetNbinsX()+1):
      if hist_estimate.GetBinContent(ibin) != 0 and hist_mc_tot_A.GetBinContent(ibin) != 0:
        #err = hist_ratio.GetBinContent(ibin) * (math.sqrt(hist_data.GetBinContent(ibin))/hist_data.GetBinContent(ibin) + math.sqrt(hist_mc_tot.GetBinContent(ibin))/hist_mc_tot.GetBinContent(ibin))
        #err = math.sqrt((math.sqrt(hist_data.GetBinContent(ibin))/hist_mc_tot.GetBinContent(ibin))**2 + (math.sqrt(hist_mc_tot.GetBinContent(ibin))*hist_data.GetBinContent(ibin)/(hist_mc_tot.GetBinContent(ibin))**2)**2)
        #err = math.sqrt((hist_data.GetBinError(ibin)/hist_mc_tot.GetBinContent(ibin))**2 + (hist_mc_tot.GetBinError(ibin)*hist_data.GetBinContent(ibin)/(hist_mc_tot.GetBinContent(ibin))**2)**2)
        err = hist_ratio.GetBinContent(ibin) * (hist_estimate.GetBinError(ibin) / hist_estimate.GetBinContent(ibin) + hist_mc_tot_A.GetBinError(ibin) / hist_mc_tot_A.GetBinContent(ibin))
      else: 
        err = 0
      if hist_ratio.GetBinContent(ibin) != 0: hist_ratio.SetBinError(ibin, err)

    hist_ratio.SetLineWidth(2)
    hist_ratio.SetLineColor(ROOT.kBlack)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetTitle('')
    hist_ratio.GetXaxis().SetTitle('#mu#pi invariant mass [GeV]')
    hist_ratio.GetXaxis().SetLabelSize(0.1)
    hist_ratio.GetXaxis().SetTitleSize(0.15)
    hist_ratio.GetXaxis().SetTitleOffset(0.73)
    hist_ratio.GetYaxis().SetTitle('Estimate/True')
    hist_ratio.GetYaxis().SetLabelSize(0.1)
    hist_ratio.GetYaxis().SetTitleSize(0.13)
    hist_ratio.GetYaxis().SetTitleOffset(0.345)
    #val_min = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
    #val_max = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
    hist_ratio.GetYaxis().SetRangeUser(-1, 3)

    hist_ratio.Draw('PE')

    line = ROOT.TLine(0, 1, 5.6, 1)
    line.SetLineColor(4)
    line.SetLineWidth(2)
    line.Draw('same')

    
    canv.cd()
    if not path.exists('./myPlots/ABCDClosure'):
      os.system('mkdir -p ./myPlots/ABCDClosure')
    canv.SaveAs('./myPlots/ABCDClosure/{}.png'.format(plot_name))
    canv.SaveAs('./myPlots/ABCDClosure/{}.pdf'.format(plot_name))
    canv.SaveAs('./myPlots/ABCDClosure/{}.C'.format(plot_name))


  def computeBkgYieldsFromMC(self):
    '''
      QCD MC (background) yields are computed as N_exp(SR) = N_obs(SR) * weight(CR)
    '''
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=1000)
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=2.83268, bin_max=3.13932)

    selection_extra = self.selection

    # we compute the weight in the control region
    weight, err_weight = self.computeApproxWeightQCDMCtoData(quantity, selection=('hnl_charge!=0' if selection_extra=='' else 'hnl_charge!=0 &&' + selection_extra))
    #weight, err_weight = self.computeWeightQCDMCtoData(lumi_data=774)

    #hist_mc_tot = self.getHistoMC(quantity=quantity, selection='hnl_charge==0' if selection_extra=='' else 'hnl_charge==0 &&' + selection_extra)
    hist_mc_tot = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection='hnl_charge==0' if selection_extra=='' else 'hnl_charge==0 &&' + selection_extra)
    n_obs_mc = hist_mc_tot.GetBinContent(1)
    n_err_mc = math.sqrt(n_obs_mc) #hist_mc_tot.GetBinError(1)

    n_exp_mc = n_obs_mc * weight
    err = n_exp_mc * (n_err_mc / n_obs_mc + err_weight / weight)

    return int(n_exp_mc), int(err), weight


  def getSignalEfficiency(self):
    '''
      eff(bin) = N_flat(bin) / N_gen
      N_gen = N_reco / filter_efficiency
    '''

    # get number of generated events
    f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+self.signal_file.filename, 'READ')
    n_reco = PlottingTools.getNminiAODEvts(self, f)
    n_gen = n_reco / self.signal_file.filter_efficiency

    # get number of selected reco events
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=1000)
    hist_flat_bin = PlottingTools.createHisto(self, f, 'signal_tree', quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection)
    n_selected_bin = hist_flat_bin.GetBinContent(1)
    #print 'n_selected_bin',n_selected_bin

    efficiency = n_selected_bin / n_gen
    #print 'efficiency',efficiency

    return efficiency


  #def computeHNLBRs(self):
  #  from decays import Decays 

  #  # first, the production BRs
  #  dec = Decays(mass=self.signal_file.mass, mixing_angle_square=1) # we factorise the mixing angle 
  #  # in our model, we only consider the case where the hnl is produced with a muon
  #  BR_prod = 1. * dec.BR_tot_mu  + 0. * dec.BR_tot_el # Bc excluded
  #  print BR_prod

  #  # then, the decay BR = Gamma(N->mupi)/Gamma_tot
  #  BR_NToMuPi = 0.0119363120843 # from our computation in HNLsGen common.py
  #  BR_decay = 0.5 * BR_NToMuPi # in our samples, the hnl decays only half of the time into muons
  #  #gamma_tot = HNLDecays(mass=self.signal_file.mass).decay_rate['tot'] # the total decay rate is computed with 0.5 mu and 0.5 el
  #  print BR_decay


  #def computeCouplingSquare(self):
  #  '''
  #    V^2 = 1/tau * 1/Gamma_tot = c/ctau * 1/Gamma_tot
  #  '''

  #  tau = self.signal_file.ctau * 1e-3 / 299792458.0 
  #  print 'tau ',tau
  #  from decays import HNLDecays
  #  gamma_tot = HNLDecays(mass=self.signal_file.mass).decay_rate['tot'] # the total decay rate is computed with 0.5 mu and 0.5 el
  #  print gamma_tot

  #  coupling_square = 1.0/(tau * gamma_tot)
  #  print coupling_square
  #  return coupling_square


  #def computeSignalYields(self, mass, lumi, sigma_B_tot, filter_efficiency):
  def computeSignalYields(self):
    '''
      signal yields computed as sigma_HNL * lumi * efficiency
    '''

    # we get the sigma_HNL
    #sigma_Bcharged = 15.3 * 1e9  # (taken from  https://www.sciencedirect.com/science/article/pii/S0370269317304379)
    sigma_Bcharged = 1e6 * 1e6  # (taken from  https://www.researchgate.net/figure/Standard-model-cross-sections-at-the-Tevatron-and-LHC-colliders-calculated-at_fig2_254469235)
    frag_Bcharged = 0.4
    
    sigma_Btot_approx = sigma_Bcharged / frag_Bcharged

    ## coupling for mass 3 and ctau 184mm, 0.5 mu 0.5 el
    v_square = 1.17457274924e-05 # obtained with HNLsGen common.py

    ## total production branching ratio
    from decays import Decays 
    dec = Decays(mass=self.signal_file.mass, mixing_angle_square=1) # we factorise the mixing angle 
    BR_prod = 1. * dec.BR_tot_mu  + 0. * dec.BR_tot_el # in our model, we only consider the case where the hnl is produced with a muon # Bc excluded
    print 'BR_prod ',BR_prod

    ## total decay branching ratio
    BR_NToMuPi = 0.0119363120843 # obtained with HNLsGen common.py
    BR_decay = 0.5 * BR_NToMuPi # in our samples, the hnl decays only half of the time into muons
    print 'BR decay, ',BR_decay

    sigma_HNL = sigma_Btot_approx * BR_prod * BR_decay * v_square
    print 'sigma_hnl ',sigma_HNL

    # lumi corresponding to the A1 period
    lumi_A1 = 0.774

    # efficiency in the mu-channel
    efficiency = self.getSignalEfficiency()
    print 'efficiency ',efficiency

    signal_yields = sigma_HNL * lumi_A1 * efficiency

    return signal_yields




if __name__ == '__main__':

  #white_list_20to300 = ['QCD_pt20to30 (V04)', 'QCD_pt30to50 (V04)', 'QCD_pt50to80 (V04)', 'QCD_pt80to120 (V04)', 'QCD_pt80to120_ext (V04)', 'QCD_pt120to170 (V04)', 'QCD_pt120to170_ext (V04)', 'QCD_pt170to300 (V04)']
  white_list_15to300 = ['QCD_pt15to20 (V05)', 'QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
  white_list_20to300 = ['QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
  white_list_20to30 = ['QCD_pt20to30 (V05)', 'QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)']

  yields = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxysig>20', white_list=white_list_15to300)

  #yields.validateABCDOnQCDMC()


  #cross_ratio = yields.computeCrossRatioFromQCDMC()
  #print '{} +- {}'.format(cross_ratio[0], cross_ratio[1])
  #bkg_yields_mc_weight = yields.computeBkgYieldsFromMC()[2]
  #print bkg_yields_mc_weight

  #yields.computeBkgYieldsFromABCDQCDMC()

  do_computeCrossRatio = False
  do_computeBkgYieldsTF = False
  do_computeBkgYieldsABCD = True
  do_testClosureABCD = True
  do_computeSigYields = False


  flag = True
  if flag:
    #baseline_selection = 'b_mass<5.3 && fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20 && deltaeta_pi_fit_pi<0.015 && deltaphi_pi_fit_pi<0.03' # if we want to categorise on cos2D, we dont want this cut to be too high
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_cos2d>0.993 && sv_prob>0.05 && trgmu_mu_mass<5.4'
    baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<6.3'
    #baseline_selection = 'sv_lxysig>20' #fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20' # if we want to categorise on cos2D, we dont want this cut to be too high
    yields_incl = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection, white_list=white_list_15to300)
    yields_disp1_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy<1 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    yields_disp2_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    #yields_disp3_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+'sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    yields_disp3_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy>5 && trgmu_charge==mu_charge', white_list=white_list_15to300)

    yields_disp1_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy<1 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    yields_disp2_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    #yields_disp3_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+'sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    yields_disp3_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && sv_lxy>5 && trgmu_charge!=mu_charge', white_list=white_list_15to300)


    if do_computeCrossRatio:
      cross_ratio_incl = yields_incl.computeCrossRatioFromQCDMC()

      cross_ratio_disp1_SS = yields_disp1_SS.computeCrossRatioFromQCDMC()
      cross_ratio_disp2_SS = yields_disp2_SS.computeCrossRatioFromQCDMC()
      cross_ratio_disp3_SS = yields_disp3_SS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp4_SS = yields_disp4_SS.computeCrossRatioFromQCDMC()

      cross_ratio_disp1_OS = yields_disp1_OS.computeCrossRatioFromQCDMC()
      cross_ratio_disp2_OS = yields_disp2_OS.computeCrossRatioFromQCDMC()
      cross_ratio_disp3_OS = yields_disp3_OS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp4_OS = yields_disp4_OS.computeCrossRatioFromQCDMC()

      print 'cross ratio, inclusive, {} +- {}'.format(round(cross_ratio_incl[0], 3), round(cross_ratio_incl[1], 3))

      print 'cross ratio, disp1, SS, {} +- {}'.format(round(cross_ratio_disp1_SS[0], 3), round(cross_ratio_disp1_SS[1], 3))
      print 'cross ratio, disp2, SS, {} +- {}'.format(round(cross_ratio_disp2_SS[0], 3), round(cross_ratio_disp2_SS[1], 3))
      print 'cross ratio, disp3, SS, {} +- {}'.format(round(cross_ratio_disp3_SS[0], 3), round(cross_ratio_disp3_SS[1], 3))
      #print 'cross ratio, disp4, SS, {} +- {}'.format(cross_ratio_disp4_SS[0], cross_ratio_disp4_SS[1])

      print 'cross ratio, disp1, OS, {} +- {}'.format(round(cross_ratio_disp1_OS[0], 3), round(cross_ratio_disp1_OS[1], 3))
      print 'cross ratio, disp2, OS, {} +- {}'.format(round(cross_ratio_disp2_OS[0], 3), round(cross_ratio_disp2_OS[1], 3))
      print 'cross ratio, disp3, OS, {} +- {}'.format(round(cross_ratio_disp3_OS[0], 3), round(cross_ratio_disp3_OS[1], 3))
      #print 'cross ratio, disp4, OS, {} +- {}'.format(cross_ratio_disp4_OS[0], cross_ratio_disp4_OS[1])


    if do_testClosureABCD:
      yields_incl.validateABCDOnQCDMC(plot_name='closure_incl_v2', label='inclusive region')
      #yields_disp1_SS.validateABCDOnQCDMC(plot_name='closure_disp1_SS_v2', label='sv l_{xy}<1cm, SS')
      #yields_disp2_SS.validateABCDOnQCDMC('closure_disp2_SS_v2', label='(1<sv l_{xy}<5)cm, SS')
      #yields_disp3_SS.validateABCDOnQCDMC('closure_disp3_SS_v2', label='sv l_{xy}>5cm, SS')
      #yields_disp1_OS.validateABCDOnQCDMC('closure_disp1_OS_v2', label='sv l_{xy}<1cm, OS')
      #yields_disp2_OS.validateABCDOnQCDMC('closure_disp2_OS_v2', label='(1<sv l_{xy}<5)cm, OS')
      #yields_disp3_OS.validateABCDOnQCDMC('closure_disp3_OS_v2', label='sv l_{xy}>5cm, OS')


    if do_computeBkgYieldsTF:
      bkg_yields_incl_mc = yields_incl.computeBkgYieldsFromMC()

      bkg_yields_disp1_SS_mc = yields_disp1_SS.computeBkgYieldsFromMC()
      bkg_yields_disp2_SS_mc = yields_disp2_SS.computeBkgYieldsFromMC()
      bkg_yields_disp3_SS_mc = yields_disp3_SS.computeBkgYieldsFromMC()
      #bkg_yields_disp4_SS_mc = yields_disp4_SS.computeBkgYieldsFromMC()

      bkg_yields_disp1_OS_mc = yields_disp1_OS.computeBkgYieldsFromMC()
      bkg_yields_disp2_OS_mc = yields_disp2_OS.computeBkgYieldsFromMC()
      bkg_yields_disp3_OS_mc = yields_disp3_OS.computeBkgYieldsFromMC()

      print 'incl, {} +- {}'.format(bkg_yields_incl_mc[0], bkg_yields_incl_mc[1])

      print 'disp1, SS, {} +- {}'.format(bkg_yields_disp1_SS_mc[0], bkg_yields_disp1_SS_mc[1])
      print 'disp2, SS, {} +- {}'.format(bkg_yields_disp2_SS_mc[0], bkg_yields_disp2_SS_mc[1])
      print 'disp3, SS, {} +- {}'.format(bkg_yields_disp3_SS_mc[0], bkg_yields_disp3_SS_mc[1])

      print 'disp1, OS, {} +- {}'.format(bkg_yields_disp1_OS_mc[0], bkg_yields_disp1_OS_mc[1])
      print 'disp2, OS, {} +- {}'.format(bkg_yields_disp2_OS_mc[0], bkg_yields_disp2_OS_mc[1])
      print 'disp3, OS, {} +- {}'.format(bkg_yields_disp3_OS_mc[0], bkg_yields_disp3_OS_mc[1])


    if do_computeBkgYieldsABCD:
      bkg_yields_incl_data = yields_incl.computeBkgYieldsFromABCDData()

      bkg_yields_disp1_SS_data = yields_disp1_SS.computeBkgYieldsFromABCDData()
      bkg_yields_disp2_SS_data = yields_disp2_SS.computeBkgYieldsFromABCDData()
      bkg_yields_disp3_SS_data = yields_disp3_SS.computeBkgYieldsFromABCDData()

      bkg_yields_disp1_OS_data = yields_disp1_OS.computeBkgYieldsFromABCDData()
      bkg_yields_disp2_OS_data = yields_disp2_OS.computeBkgYieldsFromABCDData()
      bkg_yields_disp3_OS_data = yields_disp3_OS.computeBkgYieldsFromABCDData()

      print 'incl data, {} +- {}'.format(bkg_yields_incl_data[0], bkg_yields_incl_data[1])

      print 'disp1, SS data, {} +- {}'.format(bkg_yields_disp1_SS_data[0], bkg_yields_disp1_SS_data[1])
      print 'disp2, SS data, {} +- {}'.format(bkg_yields_disp2_SS_data[0], bkg_yields_disp2_SS_data[1])
      print 'disp3, SS data, {} +- {}'.format(bkg_yields_disp3_SS_data[0], bkg_yields_disp3_SS_data[1])

      print 'disp1, OS data, {} +- {}'.format(bkg_yields_disp1_OS_data[0], bkg_yields_disp1_OS_data[1])
      print 'disp2, OS data, {} +- {}'.format(bkg_yields_disp2_OS_data[0], bkg_yields_disp2_OS_data[1])
      print 'disp3, OS data, {} +- {}'.format(bkg_yields_disp3_OS_data[0], bkg_yields_disp3_OS_data[1])


  if do_computeSigYields:
    signal_V20emu = signal_samples[1]

    sig_eff = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1').getSignalEfficiency()
    sig_eff_disp1_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge==mu_charge').getSignalEfficiency()
    sig_eff_disp2_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge').getSignalEfficiency()
    sig_eff_disp3_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge==mu_charge').getSignalEfficiency()

    sig_eff_disp1_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge!=mu_charge').getSignalEfficiency()
    sig_eff_disp2_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge').getSignalEfficiency()
    sig_eff_disp3_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge!=mu_charge').getSignalEfficiency()

    #BRs = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1').computeCouplingSquare()

    sig_incl = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1')
    sig_disp1_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge==mu_charge')
    sig_disp2_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge')
    sig_disp3_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge==mu_charge')

    sig_disp1_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge!=mu_charge')
    sig_disp2_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge')
    sig_disp3_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge!=mu_charge')

    sig_yields_incl = sig_incl.computeSignalYields()

    sig_yields_disp1_SS = sig_disp1_SS.computeSignalYields()
    sig_yields_disp2_SS = sig_disp2_SS.computeSignalYields()
    sig_yields_disp3_SS = sig_disp3_SS.computeSignalYields()

    sig_yields_disp1_OS = sig_disp1_OS.computeSignalYields()
    sig_yields_disp2_OS = sig_disp2_OS.computeSignalYields()
    sig_yields_disp3_OS = sig_disp3_OS.computeSignalYields()

    print 'signal yields, incl: {}'.format(sig_yields_incl)

    print 'signal yields, disp1, SS: {}'.format(sig_yields_disp1_SS)
    print 'signal yields, disp2, SS: {}'.format(sig_yields_disp2_SS)
    print 'signal yields, disp2, SS: {}'.format(sig_yields_disp2_SS)

    print 'signal yields, disp1, OS: {}'.format(sig_yields_disp1_OS)
    print 'signal yields, disp2, OS: {}'.format(sig_yields_disp2_OS)
    print 'signal yields, disp2, OS: {}'.format(sig_yields_disp2_OS)


    #print 'MC lxy<1, SS, {} +- {}'.format(mc_yields_disp1_SS[0], mc_yields_disp1_SS[1])
    #print 'data lxy<1, SS, {} +- {}'.format(mc_yields_disp1_SS_data[0], mc_yields_disp1_SS_data[1])

    #print 'lxy<1, SS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy<1 && trgmu_charge==mu_charge', white_list=white_list_20to300))
    #print '1<lxy<5, SS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', white_list=white_list_20to300))
    #print 'lxy>5, SS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>5 && trgmu_charge==mu_charge', white_list=white_list_20to300))

    #print 'lxy<1, OS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy<1 && trgmu_charge!=mu_charge', white_list=white_list_20to300))
    #print '1<lxy<5, OS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', white_list=white_list_20to300))
    #print 'lxy>5, OS, {}'.format(computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>5 && trgmu_charge!=mu_charge', white_list=white_list_20to300))


    #print computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>10', white_list=white_list_20to300)
    #print computeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxy>15', white_list=white_list_20to300)
