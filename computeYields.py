import os
import os.path
from os import path
import ROOT
import math
from common import PlottingTools
from samples import data_samples, data_samples_small, qcd_samples, qcd_samples, signal_samples
from quantity import Quantity
#from decays import Decays, HNLDecays

import sys
sys.path.append('../../../../BHNLGen/CMSSW_10_2_15/src/HNLsGen/python/.')
from my_common import getVV, gamma_partial, gamma_total


class ComputeYields(PlottingTools):
  def __init__(self, data_file='', qcd_files='', signal_file='', selection='', white_list=''):
    self.data_file = data_file
    self.qcd_files = qcd_files
    self.signal_file = signal_file
    self.selection = selection
    self.white_list = white_list


  def computeCrossRatioFromQCDMC(self):
    '''
      Study the correlation of two variables computing the ratio
      r = N_A * ND / N_B * N_C
    '''

    # defining mass window
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*sigma, bin_max=mass+2*sigma)

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

    if n_obs_mc != 0:
      weight = n_obs_data / n_obs_mc
    else:
      weight = 0.

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
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=2.83268, bin_max=3.13932)

    # defining mass window
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*sigma, bin_max=mass+2*sigma)

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

    if n_obs_data_D != 0:
      n_obs_data_A = n_obs_data_B * (n_obs_data_C / n_obs_data_D) 
    else:
      n_obs_data_A = 1e-9 # choose an arbitrarily number for those pathological cases

    if n_obs_data_B != 0 and n_obs_data_C != 0 and n_obs_data_D != 0:
      n_err_data_A = n_obs_data_A * (n_err_data_B / n_obs_data_B + n_err_data_C / n_obs_data_C + n_err_data_D / n_obs_data_D)
    else:
      n_err_data_A = 0

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
    # defining mass window
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    bin_min = mass-2*sigma
    bin_max = mass+2*sigma
    quantity = Quantity(name_flat='hnl_mass', nbins=50, bin_min=bin_min, bin_max=bin_max)

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
    hist_estimate = ROOT.TH1D('hist_estimate', 'hist_estimate', 50, bin_min, bin_max)
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

    canv = PlottingTools.createTCanvas(self, 'canv', 800, 900)

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
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=2.83268, bin_max=3.13932)
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=2.9514, bin_max=3.0486)

    # defining mass window
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*sigma, bin_max=mass+2*sigma)

    selection_extra = self.selection

    # we compute the weight in the control region
    weight, err_weight = self.computeApproxWeightQCDMCtoData(quantity, selection=('hnl_charge!=0' if selection_extra=='' else 'hnl_charge!=0 &&' + selection_extra))
    #weight, err_weight = self.computeWeightQCDMCtoData(lumi_data=774)

    #hist_mc_tot = self.getHistoMC(quantity=quantity, selection='hnl_charge==0' if selection_extra=='' else 'hnl_charge==0 &&' + selection_extra)
    hist_mc_tot = PlottingTools.createWeightedHistoQCDMC(self, self.qcd_files, self.white_list, quantity=quantity, selection='hnl_charge==0' if selection_extra=='' else 'hnl_charge==0 &&' + selection_extra)
    n_obs_mc = hist_mc_tot.GetBinContent(1)
    n_err_mc = math.sqrt(n_obs_mc) #hist_mc_tot.GetBinError(1)

    n_exp_mc = n_obs_mc * weight
    if n_obs_mc != 0 and weight != 0:
      err = n_exp_mc * (n_err_mc / n_obs_mc + err_weight / weight)
    else:
      err = 0.

    return int(n_exp_mc), int(err), weight


  def getCtauWeight(self): # move to common
    filename = self.signal_file.filename # might need to be modified later on for Bc
    original_ctau = filename[filename.find('ctau')+4:filename.find('/', filename.find('ctau')+1)]
    target_ctau = self.signal_file.ctau

    ctau_weight = -99
    if float(original_ctau) != float(target_ctau):
      ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

    return ctau_weight


  def getSignalEfficiency(self, isBc=False):
    '''
      eff(bin) = N_flat(bin) / N_gen
      N_gen = N_reco / filter_efficiency
    '''

    # get number of generated events
    filename = self.signal_file.filename if not isBc else self.signal_file.filename_Bc
    f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+filename, 'READ')
    n_reco = PlottingTools.getNminiAODEvts(self, f)
    filter_efficiency = self.signal_file.filter_efficiency if not isBc else self.signal_file.filter_efficiency_Bc
    n_gen = n_reco / filter_efficiency
    if self.signal_file.mass == 3.0 and not isBc: 
      n_gen = 0.5*n_gen

    # get number of selected reco events
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=13000) # go to 13TeV
    weight = self.getCtauWeight()
    #print 'begin'
    hist_flat_bin = PlottingTools.createHisto(self, f, 'signal_tree', quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection, weight=weight)
    #print 'end'
    #n_selected_bin = hist_flat_bin.Integral() 
    bin_err = ROOT.double(0.)
    n_selected_bin = hist_flat_bin.IntegralAndError(0, 13000, bin_err)
    
    #print 'n_selected_bin',n_selected_bin

    efficiency = n_selected_bin / n_gen

    # uncertainty
    err_filter_efficiency = 0.15
    err_n_selected_bin = bin_err
    err_efficiency = efficiency * (err_filter_efficiency / filter_efficiency + err_n_selected_bin / n_selected_bin) if n_selected_bin!=0 else 0.
    #print 'efficiency',efficiency

    return efficiency, err_efficiency


  def computeSignalYields(self, lumi=0.774, sigma_B=327.0e9, isBc=False):
    '''
      signal yields computed as sigma_HNL * lumi * efficiency
    '''

    # we get the sigma_HNL
    #sigma_Bcharged = 15.3 * 1e9  # (taken from  https://www.sciencedirect.com/science/article/pii/S0370269317304379)
    #sigma_Bcharged = 1e6 * 1e6  # (taken from  https://www.researchgate.net/figure/Standard-model-cross-sections-at-the-Tevatron-and-LHC-colliders-calculated-at_fig2_254469235)
    #sigma_Bcharged = 327.0e9 # fb, from MG measurement
    #sigma_Bcharged = 15.3 * 2 
    sigma_Bcharged = sigma_B
    frag_Bcharged = 0.4
    
    sigma_Btot_approx = sigma_Bcharged / frag_Bcharged

    ## coupling for mass 3 and ctau 184mm, 0.5 mu 0.5 el
    v_square = getVV(mass=self.signal_file.mass, ctau=self.signal_file.ctau, ismaj=True) #1.17457274924e-05 # obtained with HNLsGen common.py
    #print 'v2: ',v_square

    ## total production branching ratio
    from decays import Decays 
    dec = Decays(mass=self.signal_file.mass, mixing_angle_square=1) # we factorise the mixing angle 
    #BR_prod = 1. * dec.BR_tot_mu  + 0. * dec.BR_tot_el # in our model, we only consider the case where the hnl is produced with a muon # Bc excluded

    if not isBc:
      BR_prod = dec.BR_tot_mu 
    else:
      BR_prod = dec.BR_Bc_mu 
    #print 'BR_prod ',BR_prod

    # using same fashion as in gen analysis
    dec_SM = Decays(mass=0., mixing_angle_square=1)
    B_w =        (dec.B_to_uHNL.BR+       \
                  dec.B_to_D0uHNL.BR+     \
                  dec.B_to_D0staruHNL.BR+ \
                  dec.B_to_pi0uHNL.BR+    \
                  dec.B_to_rho0uHNL.BR)  

    B0_w =       (dec.B0_to_DuHNL.BR+     \
                  dec.B0_to_DstaruHNL.BR+ \
                  dec.B0_to_piuHNL.BR+    \
                  dec.B0_to_rhouHNL.BR)   
   
    Bs_w =       (dec.Bs_to_DsuHNL.BR+ \
                  dec.Bs_to_DsstaruHNL.BR+ \
                  dec.Bs_to_KuHNL.BR+ \
                  dec.Bs_to_KstaruHNL.BR)  
    
    Bc_w =     (dec.Bc_to_uHNL.BR) 

    Bspecies = ['B', 'B0', 'Bs', 'Bc']
    Bweights = [B_w, B0_w, Bs_w, Bc_w]
    Bfracs =   [0.4, 0.4, 0.1, 0.001]
    #Bspecies = ['B', 'B0', 'Bs']
    #Bweights = [B_w, B0_w, Bs_w]
    #Bfracs =   [0.4, 0.4, 0.1]
    BR_prod_2 = 0. 
    for ib,b in enumerate(Bspecies):
      BR_prod_2 += Bfracs[ib] * Bweights[ib]
    #print 'BR_prod_1 ',BR_prod 
    #print 'BR_prod_2 ',BR_prod_2 


    ## total decay branching ratio
    BR_NToMuPi = gamma_partial(mass=self.signal_file.mass, vv=v_square) / gamma_total(mass=self.signal_file.mass, vv=v_square) #0.0119363120843 # obtained with HNLsGen common.py
    #print 'BR_NToMuPi ',BR_NToMuPi

    #NToMuPi_fraction = 0.5 if self.signal_file.mass == 3. else 1.
    NToMuPi_fraction = 1.
    BR_decay = NToMuPi_fraction * BR_NToMuPi # in our samples, the hnl decays only half of the time into muons
    #print 'BR decay, ',BR_decay

    sigma_HNL = sigma_Btot_approx * BR_prod * BR_decay * v_square
    #print 'sigma_hnl ',sigma_HNL

    # lumi corresponding to the A1 period
    #lumi_A1 = 0.774 # fb-1
    #lumi_A1 = 41.6 * 1e9
    lumi_A1 = lumi

    # efficiency in the mu-channel
    #efficiency = 7.35236765292e-05 #self.getSignalEfficiency()
    efficiency, err_efficiency = self.getSignalEfficiency(isBc=isBc)
    #print 'efficiency ',efficiency #efficiency/self.signal_file.filter_efficiency

    #norm = 15.3 * 2 * 41.6 * 1e9 * self.signal_file.filter_efficiency
    norm = sigma_B * lumi * self.signal_file.filter_efficiency
    #rest = 1./0.4 * efficiency/self.signal_file.filter_efficiency * BR_prod_2 * BR_decay * v_square 
    #print 'norm ',norm
    #print 'BR_decay * vv ',BR_decay * v_square
    #print 'without acceptance ',sigma_HNL * lumi_A1 * self.signal_file.filter_efficiency
    #print 'acceptance ',efficiency/self.signal_file.filter_efficiency
    #print 'rest ',rest

    signal_yields = sigma_HNL * lumi_A1 * efficiency

    # uncertainty
    err_sigma = 0.15
    err_BR_prod = 0.05
    err_BR_decay = 0.05
    err_v_square = 0.05 
    err_lumi = 0.025
    
    err_signal_yields = 0. #signal_yields * (err_sigma / sigma_Bcharged + err_BR_prod / BR_prod + err_BR_decay / BR_decay + err_v_square / v_square + err_lumi / lumi + err_efficiency / efficiency)
    
    #print 'yields ',signal_yields

    return signal_yields, err_signal_yields







if __name__ == '__main__':

  #white_list_20to300 = ['QCD_pt20to30 (V04)', 'QCD_pt30to50 (V04)', 'QCD_pt50to80 (V04)', 'QCD_pt80to120 (V04)', 'QCD_pt80to120_ext (V04)', 'QCD_pt120to170 (V04)', 'QCD_pt120to170_ext (V04)', 'QCD_pt170to300 (V04)']
  #white_list_15to300 = ['QCD_pt15to20 (V05)', 'QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
  white_list_15to300 = ['QCD_pt15to20 (V06_29Jun21)', 'QCD_pt20to30 (V06_29Jun21)', 'QCD_pt30to50 (V06_29Jun21)', 'QCD_pt50to80 (V06_29Jun21)', 'QCD_pt80to120 (V06_29Jun21)', 'QCD_pt120to170 (V06_29Jun21)', 'QCD_pt170to300 (V06_29Jun21)']
  white_list_20to300 = ['QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
  white_list_20to30 = ['QCD_pt20to30 (V05)', 'QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)']

  yields = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, selection='sv_lxysig>20', white_list=white_list_15to300)
  #vv = yields.getVV(mass=3., ctau=184, ismajorana=True)
  vv = getVV(mass=4.5, ctau=1.2, ismaj=True)

  #yields.validateABCDOnQCDMC()


  #cross_ratio = yields.computeCrossRatioFromQCDMC()
  #print '{} +- {}'.format(cross_ratio[0], cross_ratio[1])
  #bkg_yields_mc_weight = yields.computeBkgYieldsFromMC()[2]
  #print bkg_yields_mc_weight

  #yields.computeBkgYieldsFromABCDQCDMC()

  do_computeCrossRatio = False
  do_computeBkgYieldsTF = False
  do_computeBkgYieldsABCD = False
  do_testClosureABCD = False
  do_computeSigYields = True
  plotEfficiencyScan = False


  flag = True
  if flag:
    #baseline_selection = 'b_mass<5.3 && fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20 && deltaeta_pi_fit_pi<0.015 && deltaphi_pi_fit_pi<0.03' # if we want to categorise on cos2D, we dont want this cut to be too high
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_cos2d>0.993 && sv_prob>0.05 && trgmu_mu_mass<5.4'
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<3.5 && b_mass>1.5 && pi_pt>1.3 && pi_dcasig>10 && pi_dxysig>20 && pi_dzsig>10 && mu_segmentcompatibility>0.25'
    #baseline_selection = 'sv_lxysig>20' #fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20' # if we want to categorise on cos2D, we dont want this cut to be too high
    #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4'
    baseline_selection = 'mu_isdsa!=1 && b_mass<6.4 && hnl_cos2d>0.993 && sv_prob>0.05'

    signal_V20emu = signal_samples[1]
    yields_incl = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection, white_list=white_list_15to300)
    yields_disp1_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy<1 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    yields_disp2_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    #yields_disp3_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+'sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', white_list=white_list_15to300)
    yields_disp3_SS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy>5 && trgmu_charge==mu_charge', white_list=white_list_15to300)

    yields_disp1_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy<1 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    yields_disp2_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    #yields_disp3_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+'sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', white_list=white_list_15to300)
    yields_disp3_OS = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_V20emu, selection=baseline_selection+' && sv_lxy>5 && trgmu_charge!=mu_charge', white_list=white_list_15to300)


    if do_computeCrossRatio:
      cross_ratio_incl = yields_incl.computeCrossRatioFromQCDMC()

      #cross_ratio_disp1_SS = yields_disp1_SS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp2_SS = yields_disp2_SS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp3_SS = yields_disp3_SS.computeCrossRatioFromQCDMC()
      ##cross_ratio_disp4_SS = yields_disp4_SS.computeCrossRatioFromQCDMC()

      #cross_ratio_disp1_OS = yields_disp1_OS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp2_OS = yields_disp2_OS.computeCrossRatioFromQCDMC()
      #cross_ratio_disp3_OS = yields_disp3_OS.computeCrossRatioFromQCDMC()
      ##cross_ratio_disp4_OS = yields_disp4_OS.computeCrossRatioFromQCDMC()

      print 'cross ratio, inclusive, {} +- {}'.format(round(cross_ratio_incl[0], 3), round(cross_ratio_incl[1], 3))

      #print 'cross ratio, disp1, SS, {} +- {}'.format(round(cross_ratio_disp1_SS[0], 3), round(cross_ratio_disp1_SS[1], 3))
      #print 'cross ratio, disp2, SS, {} +- {}'.format(round(cross_ratio_disp2_SS[0], 3), round(cross_ratio_disp2_SS[1], 3))
      #print 'cross ratio, disp3, SS, {} +- {}'.format(round(cross_ratio_disp3_SS[0], 3), round(cross_ratio_disp3_SS[1], 3))
      ##print 'cross ratio, disp4, SS, {} +- {}'.format(cross_ratio_disp4_SS[0], cross_ratio_disp4_SS[1])

      #print 'cross ratio, disp1, OS, {} +- {}'.format(round(cross_ratio_disp1_OS[0], 3), round(cross_ratio_disp1_OS[1], 3))
      #print 'cross ratio, disp2, OS, {} +- {}'.format(round(cross_ratio_disp2_OS[0], 3), round(cross_ratio_disp2_OS[1], 3))
      #print 'cross ratio, disp3, OS, {} +- {}'.format(round(cross_ratio_disp3_OS[0], 3), round(cross_ratio_disp3_OS[1], 3))
      ##print 'cross ratio, disp4, OS, {} +- {}'.format(cross_ratio_disp4_OS[0], cross_ratio_disp4_OS[1])


    if do_testClosureABCD:
      yields_incl.validateABCDOnQCDMC(plot_name='closure_incl_29Jun21', label='inclusive region')
      #yields_disp1_SS.validateABCDOnQCDMC(plot_name='closure_disp1_SS_v2', label='sv l_{xy}<1cm, SS')
      #yields_disp2_SS.validateABCDOnQCDMC('closure_disp2_SS_v2', label='(1<sv l_{xy}<5)cm, SS')
      #yields_disp3_SS.validateABCDOnQCDMC('closure_disp3_SS_v2', label='sv l_{xy}>5cm, SS')
      #yields_disp1_OS.validateABCDOnQCDMC('closure_disp1_OS_v2', label='sv l_{xy}<1cm, OS')
      #yields_disp2_OS.validateABCDOnQCDMC('closure_disp2_OS_v2', label='(1<sv l_{xy}<5)cm, OS')
      #yields_disp3_OS.validateABCDOnQCDMC('closure_disp3_OS_v2', label='sv l_{xy}>5cm, OS')


    if do_computeBkgYieldsTF:
      bkg_yields_incl_mc = yields_incl.computeBkgYieldsFromMC()

      #bkg_yields_disp1_SS_mc = yields_disp1_SS.computeBkgYieldsFromMC()
      #bkg_yields_disp2_SS_mc = yields_disp2_SS.computeBkgYieldsFromMC()
      #bkg_yields_disp3_SS_mc = yields_disp3_SS.computeBkgYieldsFromMC()
      ##bkg_yields_disp4_SS_mc = yields_disp4_SS.computeBkgYieldsFromMC()

      #bkg_yields_disp1_OS_mc = yields_disp1_OS.computeBkgYieldsFromMC()
      #bkg_yields_disp2_OS_mc = yields_disp2_OS.computeBkgYieldsFromMC()
      #bkg_yields_disp3_OS_mc = yields_disp3_OS.computeBkgYieldsFromMC()

      print 'incl, {} +- {}'.format(bkg_yields_incl_mc[0], bkg_yields_incl_mc[1])

      #print 'disp1, SS, {} +- {}'.format(bkg_yields_disp1_SS_mc[0], bkg_yields_disp1_SS_mc[1])
      #print 'disp2, SS, {} +- {}'.format(bkg_yields_disp2_SS_mc[0], bkg_yields_disp2_SS_mc[1])
      #print 'disp3, SS, {} +- {}'.format(bkg_yields_disp3_SS_mc[0], bkg_yields_disp3_SS_mc[1])

      #print 'disp1, OS, {} +- {}'.format(bkg_yields_disp1_OS_mc[0], bkg_yields_disp1_OS_mc[1])
      #print 'disp2, OS, {} +- {}'.format(bkg_yields_disp2_OS_mc[0], bkg_yields_disp2_OS_mc[1])
      #print 'disp3, OS, {} +- {}'.format(bkg_yields_disp3_OS_mc[0], bkg_yields_disp3_OS_mc[1])


    if do_computeBkgYieldsABCD:
      signal_file = signal_samples[0]
      yields_bkg_incl = ComputeYields(data_file=data_samples[0], qcd_files=qcd_samples, signal_file=signal_file, selection=baseline_selection, white_list=white_list_15to300)
      bkg_yields_incl_data = yields_bkg_incl.computeBkgYieldsFromABCDData()

      print 'incl data, {} +- {}'.format(bkg_yields_incl_data[0]*41.6/0.774, bkg_yields_incl_data[1])

      #bkg_yields_disp1_SS_data = yields_disp1_SS.computeBkgYieldsFromABCDData()
      #bkg_yields_disp2_SS_data = yields_disp2_SS.computeBkgYieldsFromABCDData()
      #bkg_yields_disp3_SS_data = yields_disp3_SS.computeBkgYieldsFromABCDData()

      #bkg_yields_disp1_OS_data = yields_disp1_OS.computeBkgYieldsFromABCDData()
      #bkg_yields_disp2_OS_data = yields_disp2_OS.computeBkgYieldsFromABCDData()
      #bkg_yields_disp3_OS_data = yields_disp3_OS.computeBkgYieldsFromABCDData()

      #print 'disp1, SS data, {} +- {}'.format(bkg_yields_disp1_SS_data[0], bkg_yields_disp1_SS_data[1])
      #print 'disp2, SS data, {} +- {}'.format(bkg_yields_disp2_SS_data[0], bkg_yields_disp2_SS_data[1])
      #print 'disp3, SS data, {} +- {}'.format(bkg_yields_disp3_SS_data[0], bkg_yields_disp3_SS_data[1])

      #print 'disp1, OS data, {} +- {}'.format(bkg_yields_disp1_OS_data[0], bkg_yields_disp1_OS_data[1])
      #print 'disp2, OS data, {} +- {}'.format(bkg_yields_disp2_OS_data[0], bkg_yields_disp2_OS_data[1])
      #print 'disp3, OS data, {} +- {}'.format(bkg_yields_disp3_OS_data[0], bkg_yields_disp3_OS_data[1])


  if do_computeSigYields:
    from samples import signal_samples_limits_m3
    #signal_V20emu = signal_samples[1]
    signal_V20emu = signal_samples_limits_m3[2]

    sig_eff = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1').getSignalEfficiency()[0]
    sig_eff_disp1_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge==mu_charge').getSignalEfficiency()[0]
    sig_eff_disp2_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge').getSignalEfficiency()[0]
    sig_eff_disp3_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge==mu_charge').getSignalEfficiency()[0]

    sig_eff_disp1_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge!=mu_charge').getSignalEfficiency()[0]
    sig_eff_disp2_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge').getSignalEfficiency()[0]
    sig_eff_disp3_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge!=mu_charge').getSignalEfficiency()[0]

    #BRs = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1').computeCouplingSquare()

    sig_incl = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1')
    sig_disp1_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge==mu_charge')
    sig_disp2_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge')
    sig_disp3_SS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge==mu_charge')

    sig_disp1_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy<1 && trgmu_charge!=mu_charge')
    sig_disp2_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge')
    sig_disp3_OS = ComputeYields(signal_file=signal_V20emu, selection='ismatched==1 & sv_lxy>5 && trgmu_charge!=mu_charge')

    sig_incl.getCtauWeight()

    sig_yields_incl = sig_incl.computeSignalYields(lumi=41.6)[0]

    #sig_yields_disp1_SS = sig_disp1_SS.computeSignalYields()[0]
    #sig_yields_disp2_SS = sig_disp2_SS.computeSignalYields()[0]
    #sig_yields_disp3_SS = sig_disp3_SS.computeSignalYields()[0]

    #sig_yields_disp1_OS = sig_disp1_OS.computeSignalYields()[0]
    #sig_yields_disp2_OS = sig_disp2_OS.computeSignalYields()[0]
    #sig_yields_disp3_OS = sig_disp3_OS.computeSignalYields()[0]

    print 'signal yields, incl: {}'.format(sig_yields_incl)

    #print 'signal yields, disp1, SS: {}'.format(sig_yields_disp1_SS)
    #print 'signal yields, disp2, SS: {}'.format(sig_yields_disp2_SS)
    #print 'signal yields, disp2, SS: {}'.format(sig_yields_disp2_SS)

    #print 'signal yields, disp1, OS: {}'.format(sig_yields_disp1_OS)
    #print 'signal yields, disp2, OS: {}'.format(sig_yields_disp2_OS)
    #print 'signal yields, disp2, OS: {}'.format(sig_yields_disp2_OS)

 
  if plotEfficiencyScan:
    from samples import signal_samples_effscan_m3 as samples_m3
    from samples import signal_samples_effscan_m1 as samples_m1
    from samples import signal_samples_effscan_m4p5 as samples

    canv = PlottingTools().createTCanvas('canv', 800, 800)
    #canv = ROOT.TCanvas()
    canv.SetLogx()
    canv.SetLogy()
    graph = ROOT.TGraph() 

    for sample in samples:
      sig_incl = ComputeYields(signal_file=sample, selection='ismatched==1')
      eff, err_eff = sig_incl.getSignalEfficiency()/sample.filter_efficiency
      v_square = getVV(mass=sample.mass, ctau=sample.ctau, ismaj=True)
      point = graph.GetN()
      graph.SetPoint(point, v_square, eff)
      print '{} {}'.format(v_square, eff)

    graph.Draw('AP')
    graph.SetMarkerStyle(43)
    graph.SetMarkerSize(4)
    graph.SetMarkerColor(ROOT.kRed+2)
      
    graph.GetYaxis().SetLimits(1e-6, 2e-3)
    graph.GetYaxis().SetLimits(1e-5, 1.5)

    if not path.exists('./myPlots/recoEfficiencyScan'):
      os.system('mkdir -p ./myPlots/recoEfficiencyScan')
    
    canv.SaveAs('./myPlots/recoEfficiencyScan/mass_4p5.png')

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
