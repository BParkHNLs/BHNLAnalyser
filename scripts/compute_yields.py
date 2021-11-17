import os
import os.path
from os import path
import ROOT
import math
from tools import Tools

import sys
sys.path.append('../objects')
from samples import data_samples, qcd_samples, signal_samples
from quantity import Quantity
from categories import categories
from ABCD_regions import ABCD_regions
from baseline_selection import selection
from qcd_white_list import white_list 

#TODO apply weights 

class ComputeYields(Tools):
  def __init__(self, data_files='', qcd_files='', signal_file='', selection='', white_list=''):
    self.tools = Tools()
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_file = signal_file
    self.selection = selection
    self.white_list = white_list


  ###
  # BACKGROUND TF METHOD 
  ###

  def computeApproxWeightQCDMCtoData(self, quantity, selection):
    '''
      weight = lumi_data / lumi_mc = N_data * sigma_mc / (N_mc * sigma_data) estimated as N_data / N_mc
    '''

    hist_data = self.tools.createHisto(self.data_files, 'signal_tree', quantity, branchname='flat', selection=selection)
    hist_data.Sumw2()
    n_obs_data = hist_data.Integral()
    n_err_data = math.sqrt(n_obs_data) #hist_data.GetBinError(1)

    hist_mc_tot = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=selection)
    n_obs_mc = hist_mc_tot.Integral()
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

      f_mc = self.tools.getRootFile(qcd_file.filename)
      
      weight = self.tools.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
      weight_mc = 1./weight
      lumi_mc += weight_mc

    weight = lumi_data / lumi_mc 
    err = 0.

    return weight, err


  def computeBkgYieldsFromMC(self, mass, resolution):
    '''
      QCD MC (background) yields are computed as N_exp(SR) = N_obs(SR) * weight(CR)
    '''

    # defining mass window
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*resolution, bin_max=mass+2*resolution)

    # make sure there is no selection on the hnl charge 
    if 'hnl_charge' in self.selection:
      substr = self.selection[self.selection.find('hnl_charge'):self.selection.find('&&',self.selection.find('hnl_charge')+1)+2]
      selection_extra = self.selection.replace(substr, '')
    else:
      selection_extra = self.selection

    # we compute the weight in the control region
    weight, err_weight = self.computeApproxWeightQCDMCtoData(quantity, selection=('hnl_charge!=0' if selection_extra=='' else 'hnl_charge!=0 &&' + selection_extra))

    hist_mc_tot = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection='hnl_charge==0' if selection_extra=='' else 'hnl_charge==0 &&' + selection_extra)
    n_obs_mc = hist_mc_tot.Integral()
    n_err_mc = math.sqrt(n_obs_mc) #hist_mc_tot.GetBinError(1)

    n_exp_mc = n_obs_mc * weight
    if n_obs_mc != 0 and weight != 0:
      err = n_exp_mc * (n_err_mc / n_obs_mc + err_weight / weight)
    else:
      err = 0.

    return n_exp_mc, err, weight


  ###
  # BACKGROUND ABCD METHOD 
  ###

  def computeBkgYieldsFromABCDData(self, mass, resolution, ABCD_regions):
    '''
      Estimate background yields from data using the ABCD method
      A = b_mass < 6.27 && hnl_charge == 0 (SR)
      B = b_mass < 6.27 && hnl_charge != 0 
      C = b_mass > 6.27 && hnl_charge == 0 
      D = b_mass > 6.27 && hnl_charge != 0 

      N_A = N_B * N_C/N_D
    '''
    # defining mass window
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*resolution, bin_max=mass+2*resolution)

    bin_selection = self.selection

    hist_data_B = self.tools.createHisto(self.data_files, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_B_selection)
    hist_data_C = self.tools.createHisto(self.data_files, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_C_selection)
    hist_data_D = self.tools.createHisto(self.data_files, 'signal_tree', quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_D_selection)

    n_obs_data_B = hist_data_B.Integral()
    n_obs_data_C = hist_data_C.Integral()
    n_obs_data_D = hist_data_D.Integral()

    n_err_data_B = math.sqrt(n_obs_data_B)
    n_err_data_C = math.sqrt(n_obs_data_C)
    n_err_data_D = math.sqrt(n_obs_data_D)

    #print '{} * ({} / {})'.format(n_obs_data_B, n_obs_data_C, n_obs_data_D)

    if n_obs_data_D != 0:
      n_obs_data_A = n_obs_data_B * (n_obs_data_C / n_obs_data_D) 
    else:
      n_obs_data_A = 1e-9 # choose an arbitrarily small number for those pathological cases

    if n_obs_data_B != 0 and n_obs_data_C != 0 and n_obs_data_D != 0:
      n_err_data_A = n_obs_data_A * (n_err_data_B / n_obs_data_B + n_err_data_C / n_obs_data_C + n_err_data_D / n_obs_data_D)
    else:
      n_err_data_A = 0

    return n_obs_data_A, n_err_data_A, n_obs_data_B, n_obs_data_C, n_obs_data_D


  def computeBkgYieldsFromABCDHybrid(self, mass, resolution, ABCD_regions):
    '''
      Estimate background yields using the ABCD method on data
      unless one CR has 0 entries, then compute the background
      yields using the TF method
    '''

    background_yields, err, n_CR_B, n_CR_C, n_CR_D = self.computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=ABCD_regions)
    if n_CR_B == 0 or n_CR_C == 0 or n_CR_D == 0: #TODO correct to use TF if NB or NC = 0?
      background_yields, err, wght = self.computeBkgYieldsFromMC(mass=mass, resolution=resolution)

    return background_yields, err


  ###
  # BACKGROUND ABCD VALIDATION 
  ###

  def computeCrossRatioFromQCDMC(self, ABCD_regions):
    '''
      Study the correlation of two variables computing the ratio
      r = N_A * ND / N_B * N_C
    '''

    # defining mass window
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*sigma, bin_max=mass+2*sigma)
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=6.4)

    hist_mc_tot_A = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.SR_selection)
    hist_mc_tot_B = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_B_selection)
    hist_mc_tot_C = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_C_selection)
    hist_mc_tot_D = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_D_selection)

    n_obs_mc_A = hist_mc_tot_A.Integral()
    n_obs_mc_B = hist_mc_tot_B.Integral()
    n_obs_mc_C = hist_mc_tot_C.Integral()
    n_obs_mc_D = hist_mc_tot_D.Integral()

    n_err_mc_A = hist_mc_tot_A.GetBinError(1) #math.sqrt(n_obs_mc_A) 
    n_err_mc_B = hist_mc_tot_B.GetBinError(1) #math.sqrt(n_obs_mc_B) 
    n_err_mc_C = hist_mc_tot_C.GetBinError(1) #math.sqrt(n_obs_mc_C) 
    n_err_mc_D = hist_mc_tot_D.GetBinError(1) #math.sqrt(n_obs_mc_D) 

    #print '{} * {} / ({} * {})'.format(n_obs_mc_A, n_obs_mc_D, n_obs_mc_B, n_obs_mc_C)

    if n_obs_mc_B != 0 and n_obs_mc_C !=0:
      ratio = n_obs_mc_A*n_obs_mc_D/(n_obs_mc_B*n_obs_mc_C)
    else:
      ratio = 0.

    if n_obs_mc_B != 0 and n_obs_mc_C !=0 and n_obs_mc_A != 0 and n_obs_mc_D != 0:
      err = ratio * (n_err_mc_A/n_obs_mc_A + n_err_mc_B/n_obs_mc_B + n_err_mc_C/n_obs_mc_C + n_err_mc_D/n_obs_mc_D)
    else: err = 0.

    return ratio, err


  def validateABCDOnQCDMC(self, outlabel='', label='', ABCD_regions=''):
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

    hist_mc_tot_A = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.SR_selection)
    hist_mc_tot_B = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_B_selection)
    hist_mc_tot_C = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_C_selection)
    hist_mc_tot_D = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, selection=self.selection+' && '+ABCD_regions.CR_D_selection)

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

    canv = self.tools.createTCanvas('canv', 900, 800)

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
    hist_estimate.GetXaxis().SetLabelSize(0.0)
    hist_estimate.GetXaxis().SetTitleSize(0.0)
    hist_estimate.GetYaxis().SetTitle('Entries')
    hist_estimate.GetYaxis().SetLabelSize(0.037)
    hist_estimate.GetYaxis().SetTitleSize(0.042)
    hist_estimate.GetYaxis().SetTitleOffset(1.1)

    hist_mc_tot_A.SetLineWidth(3)
    hist_mc_tot_A.SetLineColor(ROOT.kBlue)

    int_estimate = hist_estimate.Integral()
    if int_estimate != 0: hist_estimate.Scale(1/int_estimate)
    int_A = hist_mc_tot_A.Integral()
    hist_mc_tot_A.Scale(1/int_A)

    pad_up.cd()
    hist_estimate.Draw('histo')
    hist_mc_tot_A.Draw('histo same')

    legend = self.tools.getRootTLegend(xmin=0.55, ymin=0.65, xmax=0.8, ymax=0.85, size=0.043)
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

    hist_ratio = self.tools.getRatioHistogram(hist_estimate, hist_mc_tot_A)
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

    line = ROOT.TLine(bin_min, 1, bin_max, 1)
    line.SetLineColor(4)
    line.SetLineWidth(2)
    line.Draw('same')
    
    canv.cd()
    outdir = './myPlots/ABCDClosure/{}/{}'.format(outlabel, ABCD_regions.label)
    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))
    canv.SaveAs('{}/closure_{}.png'.format(outdir, label))
    canv.SaveAs('{}/closure_{}.pdf'.format(outdir, label))

  ###
  # SIGNAL 
  ###

  def getSignalEfficiency(self, add_weight_hlt=True, weight_hlt='', isBc=False):
    '''
      eff(bin) = N_flat(bin) / N_gen
      N_gen = N_reco / filter_efficiency
    '''

    # get number of generated events
    filename = self.signal_file.filename if not isBc else self.signal_file.filename_Bc
    f = self.tools.getRootFile(filename)
    n_reco = self.tools.getNminiAODEvts(f)
    filter_efficiency = self.signal_file.filter_efficiency if not isBc else self.signal_file.filter_efficiency_Bc
    n_gen = n_reco / filter_efficiency
    #if self.signal_file.mass == 3.0 and not isBc: 
    #  n_gen = 0.5*n_gen

    # central samples are produced 50% muon and 50% electron
    n_gen = 0.5*n_gen #TODO do not hardcode it?

    # get number of selected reco events
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=13000) # go to 13TeV
    weight = '({})'.format(self.tools.getCtauWeight(self.signal_file))
    if add_weight_hlt : weight += ' * ({})'.format(weight_hlt)
    #TODO add pu weight
    hist_flat_bin = self.tools.createHisto(self.signal_file, 'signal_tree', quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection, weight=weight)
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


  def computeSignalYields(self, lumi=0.774, sigma_B=472.8e9, add_weight_hlt=True, weight_hlt='', isBc=False):
    '''
      signal yields computed as sigma_HNL * lumi * efficiency
    '''

    # we get the sigma_HNL
    ## measured sigma_B_tot from normalisation study
    sigma_Bcharged = sigma_B

    frag_Bcharged = 0.4
    sigma_Btot_approx = sigma_Bcharged / frag_Bcharged

    v_square = self.tools.getVV(mass=self.signal_file.mass, ctau=self.signal_file.ctau, ismaj=True)
    #print 'v2: ',v_square

    ## total production branching ratio
    from decays import Decays 
    dec = Decays(mass=self.signal_file.mass, mixing_angle_square=1) # we factorise the mixing angle 
    if not isBc:
      BR_prod = dec.BR_tot_mu 
    else:
      BR_prod = dec.BR_Bc_mu 
    #print 'BR_prod ',BR_prod

    ## total decay branching ratio
    BR_NToMuPi = self.tools.gamma_partial(mass=self.signal_file.mass, vv=v_square) / self.tools.gamma_total(mass=self.signal_file.mass, vv=v_square)
    #print 'BR_NToMuPi ',BR_NToMuPi

    #NToMuPi_fraction = 0.5 if self.signal_file.mass == 3. else 1. #FIXME understand why this line was commented
    NToMuPi_fraction = 1.
    BR_decay = NToMuPi_fraction * BR_NToMuPi # in our samples, the hnl decays only half of the time into muons
    #print 'BR decay, ',BR_decay

    ## and finally sigma_HNL
    sigma_HNL = sigma_Btot_approx * BR_prod * BR_decay * v_square
    #print 'sigma_hnl ',sigma_HNL

    # lumi
    lumi_A1 = lumi

    # efficiency in the mu-channel
    efficiency, err_efficiency = self.getSignalEfficiency(add_weight_hlt=add_weight_hlt, weight_hlt=weight_hlt, isBc=isBc)
    #print 'efficiency ',efficiency #efficiency/self.signal_file.filter_efficiency

    signal_yields = sigma_HNL * lumi_A1 * efficiency

    # uncertainty
    err_sigma = 0.10
    err_BR_prod = 0.05
    err_BR_decay = 0.05
    err_v_square = 0.05 
    err_lumi = 0.025
    
    err_signal_yields = 0. #signal_yields * (err_sigma / sigma_Bcharged + err_BR_prod / BR_prod + err_BR_decay / BR_decay + err_v_square / v_square + err_lumi / lumi + err_efficiency / efficiency)
    
    #print 'yields ',signal_yields

    return signal_yields, err_signal_yields



if __name__ == '__main__':

  data_samples = data_samples['V09_06Nov21']
  qcd_samples = qcd_samples['V09_06Nov21']
  signal_file = signal_samples['central_V09_06Nov21_benchmark'][1]

  baseline_selection = selection['standard'].flat
  #baseline_selection = 'hnl_pt>0' # && mu_isdsa==1' #selection[''].flat
  categories = categories['inclusive']
  #categories = categories['category_study_combined_dsa']

  ABCD_regions = ABCD_regions['cos2d_svprob']
  #ABCD_regions = ABCD_regions['bmass_hnlcharge']

  #TODO fix white list handling 
  #white_list_20to300 = ['QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
  white_list_20to300 = white_list['20to300']

  do_computeCrossRatio = False
  do_computeBkgYieldsTF = False
  do_computeBkgYieldsABCD = True
  do_computeBkgYieldsABCDHybrid = False
  do_testClosureABCD = False
  do_computeSigYields = False

  #baseline_selection = 'b_mass<5.3 && fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20 && deltaeta_pi_fit_pi<0.015 && deltaphi_pi_fit_pi<0.03' # if we want to categorise on cos2D, we dont want this cut to be too high
  #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_cos2d>0.993 && sv_prob>0.05 && trgmu_mu_mass<5.4'
  #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<3.5 && b_mass>1.5 && pi_pt>1.3 && pi_dcasig>10 && pi_dxysig>20 && pi_dzsig>10 && mu_segmentcompatibility>0.25'
  #baseline_selection = 'sv_lxysig>20' #fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20' # if we want to categorise on cos2D, we dont want this cut to be too high
  #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4'
  #baseline_selection = 'mu_isdsa!=1 && b_mass<6.4 && hnl_cos2d>0.993 && sv_prob>0.05'

  if do_computeCrossRatio:
    for category in categories:
      yields = ComputeYields(data_files=data_samples[0], qcd_files=qcd_samples, signal_file=signal_file, selection=baseline_selection+' && '+category.definition_flat, white_list=white_list_20to300)
      cross_ratio = yields.computeCrossRatioFromQCDMC(ABCD_regions=ABCD_regions)
      print 'cross ratio, {}, {} +- {}'.format(category.label, round(cross_ratio[0], 3), round(cross_ratio[1], 3))


  if do_testClosureABCD:
    ROOT.gROOT.SetBatch(True)
    for category in categories:
      outlabel = 'category_study_combined_dsa' #TODO
      yields = ComputeYields(data_files=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && '+category.definition_flat+' && '+category.cutbased_selection, white_list=white_list_20to300)
      yields.validateABCDOnQCDMC(outlabel=outlabel, label=category.label, ABCD_regions=ABCD_regions)


  if do_computeBkgYieldsTF:
    #baseline_selection = selection['sensitivity_study_July'].flat
    for category in categories:
      yields = ComputeYields(data_files=data_samples, qcd_files=qcd_samples, selection=baseline_selection+' && '+category.definition_flat+' && '+category.cutbased_selection, white_list=white_list_20to300)
      bkg_yields_mc = yields.computeBkgYieldsFromMC(mass=signal_file.mass, resolution=signal_file.resolution)
      print '{} TF, {} +- {}'.format(category.title, bkg_yields_mc[0], bkg_yields_mc[1])


  if do_computeBkgYieldsABCD:
    #baseline_selection = selection['sensitivity_study_July'].flat
    for category in categories:
      yields_bkg_incl = ComputeYields(data_files=data_samples, qcd_files=qcd_samples, signal_file=signal_file, selection=baseline_selection+' && '+category.definition_flat+' && '+category.cutbased_selection, white_list=white_list_20to300)
      bkg_yields_incl_data = yields_bkg_incl.computeBkgYieldsFromABCDData(mass=signal_file.mass, resolution=signal_file.resolution, ABCD_regions=ABCD_regions)
      print '{} ABCD, {} +- {}'.format(category.title, bkg_yields_incl_data[0], bkg_yields_incl_data[1])


  if do_computeBkgYieldsABCDHybrid:
    #baseline_selection = selection['sensitivity_study_July'].flat
    for category in categories:
      yields_bkg_incl = ComputeYields(data_files=data_samples[0], qcd_files=qcd_samples, signal_file=signal_file, selection=baseline_selection+' && '+category.definition_flat+' && '+category.cutbased_selection, white_list=white_list_20to300)
      bkg_yields_incl_data = yields_bkg_incl.computeBkgYieldsFromABCDHybrid(mass=signal_file.mass, resolution=signal_file.resolution, ABCD_regions=ABCD_regions)
      print '{} hybrid ABCD, {} +- {}'.format(category.title, bkg_yields_incl_data[0], bkg_yields_incl_data[1])


  if do_computeSigYields:
    for category in categories:
      sig_eff = ComputeYields(signal_file=signal_file, selection='ismatched==1 && '+category.definition_flat).getSignalEfficiency()[0]
      #BRs = ComputeYields(signal_file=signal_file, selection='ismatched==1').computeCouplingSquare()

      sig = ComputeYields(signal_file=signal_file, selection='ismatched==1 && '+category.definition_flat)
      #sig.getCtauWeight()

      sig_yields = sig.computeSignalYields(lumi=41.6, weight_hlt='weight_hlt_A1')[0]
      print 'signal yields, {}: {}'.format(category.title, sig_yields)
 
