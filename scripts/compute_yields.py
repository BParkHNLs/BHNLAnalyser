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
  def __init__(self, data_files='', qcd_files='', signal_file='', signal_label='', selection='', white_list=''):
    self.tools = Tools()
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_file = signal_file
    self.signal_label = signal_label # this will be needed for the inclusive signal reweighting
    self.selection = selection
    self.white_list = white_list


  ###
  # BACKGROUND TF METHOD 
  ###

  def computeApproxWeightQCDMCtoData(self, quantity, selection):
    '''
      weight = lumi_data / lumi_mc = N_data * sigma_mc / (N_mc * sigma_data) estimated as N_data / N_mc
    '''

    tree = ROOT.TChain('signal_tree')
    for data_file in self.data_files:
      tree.Add(data_file.filename)  

    hist_data = self.tools.createHisto(tree, quantity, branchname='flat', selection=selection)
    hist_data.Sumw2()
    n_obs_data = hist_data.Integral()
    n_err_data = math.sqrt(n_obs_data) #hist_data.GetBinError(1)

    #TODO add weights
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
      qcd_file_pthatrange = self.tools.getPthatRange(qcd_file.label)
      if qcd_file_pthatrange not in self.white_list: continue

      f_mc = self.tools.getRootFile(qcd_file.filename)
      tree_run = self.tools.getTree(f_mc, 'run_tree')
      
      weight = self.tools.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
      weight_mc = 1./weight
      lumi_mc += weight_mc

    weight = lumi_data / lumi_mc 
    err = 0.

    return weight, err


  def computeBkgYieldsFromMC(self, mass, resolution, sigma_mult_window=2):
    '''
      QCD MC (background) yields are computed as N_exp(SR) = N_obs(SR) * weight(CR)
    '''

    # defining mass window
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-sigma_mult_window*resolution, bin_max=mass+sigma_mult_window*resolution)

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

    # rescale yields to the 2 sigma window
    n_exp_mc = n_exp_mc * 2. / sigma_mult_window

    return n_exp_mc, err, weight


  ###
  # BACKGROUND ABCD METHOD 
  ###

  def computeBkgYieldsFromABCDData(self, mass, resolution, ABCD_regions, sigma_mult_window=2):
    '''
      Estimate background yields from data using the ABCD method
      A = b_mass < 6.27 && hnl_charge == 0 (SR)
      B = b_mass < 6.27 && hnl_charge != 0 
      C = b_mass > 6.27 && hnl_charge == 0 
      D = b_mass > 6.27 && hnl_charge != 0 

      N_A = N_B * N_C/N_D
    '''
    # defining mass window of width sigma_mult_window*sigma
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-sigma_mult_window*resolution, bin_max=mass+sigma_mult_window*resolution)
    print '{} {}'.format(mass-sigma_mult_window*resolution, mass+sigma_mult_window*resolution)

    bin_selection = self.selection

    # create tree
    tree = ROOT.TChain('signal_tree')
    for data_file in self.data_files:
      tree.Add(data_file.filename)  

    hist_data_B = self.tools.createHisto(tree, quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_B_selection)
    hist_data_C = self.tools.createHisto(tree, quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_C_selection)
    hist_data_D = self.tools.createHisto(tree, quantity, branchname='flat', selection=bin_selection+' && '+ABCD_regions.CR_D_selection)

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

    # rescale yields to the 2 sigma window
    n_obs_data_A = n_obs_data_A * 2. / sigma_mult_window
    n_obs_data_B = n_obs_data_B * 2. / sigma_mult_window
    n_obs_data_C = n_obs_data_C * 2. / sigma_mult_window
    n_obs_data_D = n_obs_data_D * 2. / sigma_mult_window

    return n_obs_data_A, n_err_data_A, n_obs_data_B, n_obs_data_C, n_obs_data_D


  def computeBkgYieldsFromABCDHybrid(self, mass, resolution, ABCD_regions, sigma_mult_window=2):
    '''
      Estimate background yields using the ABCD method on data
      unless one CR has 0 entries, then compute the background
      yields using the TF method
    '''

    background_yields, err, n_CR_B, n_CR_C, n_CR_D = self.computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=ABCD_regions, sigma_mult_window=sigma_mult_window)
    if n_CR_B == 0 or n_CR_C == 0 or n_CR_D == 0: #TODO correct to use TF if NB or NC = 0?
      background_yields, err, wght = self.computeBkgYieldsFromMC(mass=mass, resolution=resolution, sigma_mult_window=sigma_mult_window)

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

    hist_mc_tot_A = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, hist_name='hist_mc_tot_A', selection=self.selection+' && '+ABCD_regions.SR_selection)
    hist_mc_tot_B = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, hist_name='hist_mc_tot_B', selection=self.selection+' && '+ABCD_regions.CR_B_selection)
    hist_mc_tot_C = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, hist_name='hist_mc_tot_C', selection=self.selection+' && '+ABCD_regions.CR_C_selection)
    hist_mc_tot_D = self.tools.createWeightedHistoQCDMC(self.qcd_files, self.white_list, quantity=quantity, hist_name='hist_mc_tot_D', selection=self.selection+' && '+ABCD_regions.CR_D_selection)

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
    bin_min = 0. #mass-2*sigma
    bin_max = 5.4 #mass+2*sigma
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
    #if int_estimate != 0: hist_estimate.Scale(1/int_estimate)
    int_A = hist_mc_tot_A.Integral()
    #hist_mc_tot_A.Scale(1/int_A)

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

  def getSignalEfficiency(self, add_weight_hlt=False, add_weight_pu=False, weight_hlt='weight_hlt_D', weight_pusig='weight_pu_sig_D', isMixed=False, lhe_efficiency=0.08244, isBc=False):
    '''
      eff(bin) = N_flat(bin) / N_gen
      N_gen = N_reco / filter_efficiency
    '''

    # get the generated that will be used to do the reweighting (all ctau points at a given mass)
    samples_for_reweighting = []
    the_signal_samples = signal_samples[self.signal_label]
    do_exclusive_reweighting_fromlargerctau = True
    do_exclusive_reweighting_fromsmallerctau = False
    do_full_inclusive_reweighting = False
    do_partial_inclusive_reweighting = False
    do_unique_reweighting = False

    the_generated_samples = []
    for the_signal_sample in the_signal_samples:
      if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
      if '*' in the_signal_sample.filename: continue
      the_generated_samples.append(the_signal_sample)

    if do_full_inclusive_reweighting:
      # perform the reweighting with a set of generated samples
      for the_signal_sample in the_generated_samples:
        if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        samples_for_reweighting.append(the_signal_sample)
    elif do_partial_inclusive_reweighting:
      # perform the reweighting with a set of generated samples
      for the_signal_sample in the_generated_samples:
        if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        if the_signal_sample.ctau < self.signal_file.ctau: 
        #  if len(samples_for_reweighting) == 0: # for mass 1, we reweight to even larger ctaus than the ones generated
        #    samples_for_reweighting.append(the_signal_sample)
          samples_for_reweighting.append(the_signal_sample)
          break
        samples_for_reweighting.append(the_signal_sample)
    elif do_exclusive_reweighting_fromlargerctau:
      # only use the generated sample with the closest larger ctau
      for the_signal_sample in the_generated_samples:
        if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        if the_signal_sample.ctau >= self.signal_file.ctau:
          samples_for_reweighting.append(the_signal_sample)
        elif the_signal_sample.ctau < self.signal_file.ctau and len(samples_for_reweighting) == 0: 
          # for mass 1, we reweight to even larger ctaus than the ones generated
          samples_for_reweighting.append(the_signal_sample)

      if len(samples_for_reweighting) > 1:
        # only keep last element of the list
        samples_for_reweighting_tmp = []
        samples_for_reweighting_tmp.append(samples_for_reweighting.pop())
        samples_for_reweighting = samples_for_reweighting_tmp

    elif do_exclusive_reweighting_fromsmallerctau:
      # only use the generated sample with the closest larger ctau
      for the_signal_sample in the_generated_samples:
        if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        if the_signal_sample.ctau <= self.signal_file.ctau:
          samples_for_reweighting.append(the_signal_sample)
        elif the_signal_sample.ctau < self.signal_file.ctau and len(samples_for_reweighting) == 0: 
          # for mass 1, we reweight to even larger ctaus than the ones generated
          samples_for_reweighting.append(the_signal_sample)

      if len(samples_for_reweighting) > 1:
        # only keep last element of the list
        samples_for_reweighting_tmp = []
        samples_for_reweighting_tmp.append(samples_for_reweighting[0])
        samples_for_reweighting = samples_for_reweighting_tmp

      if len(samples_for_reweighting) == 0:
        samples_for_reweighting.append(the_generated_samples.pop())

    elif do_unique_reweighting:
      for the_signal_sample in the_generated_samples:
        if str(self.signal_file.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if 'mN1p0' in the_signal_sample.filename and 'ctau1000p0mm' in the_signal_sample.filename: 
          samples_for_reweighting.append(the_signal_sample)
        if 'mN3p0' in the_signal_sample.filename and 'ctau1p0mm' in the_signal_sample.filename: 
          samples_for_reweighting.append(the_signal_sample)
        if 'mN4p5' in the_signal_sample.filename and 'ctau100p0mm' in the_signal_sample.filename: 
          samples_for_reweighting.append(the_signal_sample)

    #print '\n {}'.format(self.signal_file.ctau)
    #for the_signal_sample in samples_for_reweighting:
    #  print the_signal_sample.filename

    # get the trees
    filename = self.signal_file.filename if not isBc else self.signal_file.filename_Bc
    if '*' not in filename: 
      # no ctau reweighting
      f = self.tools.getRootFile(filename)
      tree_sig = self.getTree(f, 'signal_tree')
      tree_run = self.getTree(f, 'run_tree')
    else: 
      # use inclusive (in ctau) samples
      tree_sig = ROOT.TChain('signal_tree')
      tree_run = ROOT.TChain('run_tree')
      for the_signal_sample in samples_for_reweighting:
        tree_sig.Add(the_signal_sample.filename)
        tree_run.Add(the_signal_sample.filename)

    # get number of gen events (= number of miniAOD events)
    n_miniaod_tot = self.tools.getNminiAODEvts(tree_run)

    # get the filter efficiency
    if '*' not in filename:
      filter_efficiency = self.signal_file.filter_efficiency if not isBc else self.signal_file.filter_efficiency_Bc
    else:
      # take the weighted average of the filter efficiencies of the generated samples
      filter_efficiency = 0.
      for the_signal_sample in samples_for_reweighting:
        the_file = self.tools.getRootFile(the_signal_sample.filename)
        the_tree_run = self.getTree(the_file, 'run_tree')
        n_miniaod = self.tools.getNminiAODEvts(the_tree_run)
        filter_efficiency += n_miniaod * the_signal_sample.filter_efficiency
      filter_efficiency = filter_efficiency / n_miniaod_tot

    # get the number of generated events
    n_generated = n_miniaod_tot / filter_efficiency
    if isMixed: n_generated = n_generated * 0.5

    # get number of selected reco events
    mass = self.signal_file.mass
    sigma = self.signal_file.resolution
    #quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-2*sigma, bin_max=mass+2*sigma)
    quantity = Quantity(name_flat='hnl_mass', nbins=1, bin_min=mass-50*sigma, bin_max=mass+50*sigma)

    # question: is the Ni the miniaod or the generated ones?

    if '*' not in filename:
      weight = '(1)'
    else:
      deno_weight = ''
      for ifile, the_signal_sample in enumerate(samples_for_reweighting):
        the_file = self.tools.getRootFile(the_signal_sample.filename)
        the_tree_run = self.getTree(the_file, 'run_tree')
        n_miniaod = self.tools.getNminiAODEvts(the_tree_run)
        #the_filter_efficiency = the_signal_sample.filter_efficiency
        #the_n_generated = n_miniaod / the_filter_efficiency
        #the_tree_sig = self.getTree(the_file, 'signal_tree')
        #the_hist = self.tools.createHisto(the_tree_sig, quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection, weight='(1.)')
        #n_events = the_hist.Integral()
        #hist_tot = self.tools.createHisto(tree_sig, quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection, weight='(1.)')
        #n_events_tot = hist_tot.Integral()
        if ifile == 0:
          deno_weight += ' {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod,
                #n0 = n_generated,
                #n0 = n_events,
                ctau0 = the_signal_sample.ctau,
                )
        else:
          deno_weight += ' + {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod,
                #n0 = n_generated,
                #n0 = n_events,
                ctau0 = the_signal_sample.ctau,
                )
      weight = '({ntot} / {ctau1} * exp(-gen_hnl_ct / {ctau1}) * (1. / ({deno_weight})))'.format(
          ntot = n_miniaod_tot,
          #ntot = n_generated,
          #ntot = n_events_tot,
          ctau1 = self.signal_file.ctau,
          deno_weight = deno_weight,
          )

    #((365329.0 / 1.5 * exp(-gen_hnl_ct / 1.5)) * (1. / ( 365329.0 / 10.0 * exp(-gen_hnl_ct / 10.0))))
    #(365329.0 / 1.5 * exp(-gen_hnl_ct / 1.5) * (1. / ( 365329.0 / 10.0 * exp(-gen_hnl_ct / 10.0))))
    #(3629644.0 / 1.5 * exp(-gen_hnl_ct / 1.5) * (1. / ( 2128601.0 / 1000.0 * exp(-gen_hnl_ct / 1000.0 + 808459.0 / 100.0 * exp(-gen_hnl_ct / 100.0 + 365329.0 / 10.0 * exp(-gen_hnl_ct / 10.0 + 327255.0 / 1.0 * exp(-gen_hnl_ct / 1.0))))

    if add_weight_hlt : weight += ' * ({})'.format(weight_hlt)
    if add_weight_pu : weight += ' * ({})'.format(weight_pusig)

    # plot the weight distribution
    #canv_weight = ROOT.TCanvas('canv_weight', 'canv_weight', 800, 700)
    #hist = ROOT.TH1D('hist_m{}_ctau{}'.format(self.signal_file.mass, self.signal_file.ctau), 'hist_m{}_ctau{}'.format(self.signal_file.mass, self.signal_file.ctau), 100, 0, 5)
    #tree_sig.Draw('{}>>hist_m{}_ctau{}'.format(weight, self.signal_file.mass, self.signal_file.ctau))
    #hist.Draw('hist')
    #canv_weight.SaveAs('myPlots/lifetime_weights/weight_m{}_ctau{}mm.png'.format(str(self.signal_file.mass).replace('.', 'p'), str(self.signal_file.ctau).replace('.', 'p')))

    # this was to check if a cut on the weight would impact the performance
    #self.selection += ' && {}<30'.format(weight)
    #print self.selection

    hist_flat_bin = self.tools.createHisto(tree_sig, quantity, branchname='flat', selection='ismatched==1' if self.selection=='' else 'ismatched==1 && '+self.selection, weight=weight)
    bin_err = ROOT.double(0.)
    n_selected_bin = hist_flat_bin.IntegralAndError(0, 13000, bin_err)
    #n_selected_bin = hist_flat_bin.Integral() #FIXME error, overflow?

    # compute the efficiency
    efficiency = n_selected_bin / n_generated
    if isBc:
      efficiency = efficiency * lhe_efficiency
    #err_efficiency = 0. #TODO

    err_filter_efficiency = 0. * filter_efficiency
    err_n_selected_bin = math.sqrt(n_selected_bin) #bin_err
    err_efficiency = efficiency * (err_filter_efficiency / filter_efficiency + err_n_selected_bin / n_selected_bin) if n_selected_bin!=0 else 0.
    #print '{} {} {} {} {}%'.format(self.signal_file.mass, self.signal_file.ctau, efficiency, err_efficiency, err_efficiency/efficiency*100)

    return efficiency, err_efficiency


  def computeSignalYields(self, mass='', ctau='', lumi=0.774, sigma_B=472.8e9, add_weight_hlt=False, add_weight_pu=False, weight_hlt='weight_hlt_D', weight_pusig='weight_pu_sig_D', isBc=False):
    '''
      signal yields computed as sigma_HNL * lumi * efficiency
    '''

    # define ranges and binning (full window)
    fit_window_min = mass - 1.5
    fit_window_max = mass + 1.5 

    # get the signal files
    signal_files = self.tools.getSignalFileList(signal_label=self.signal_label, mass=mass, ctau=ctau, strategy='exclusive_fromlargerctau')

    # get the tree
    treename = 'signal_tree'
    tree_sig = ROOT.TChain(treename)
    for signal_file in signal_files:
      filename = signal_file.filename if not isBc else signal_file.filename_Bc
      tree_sig.Add(filename)

    # define selection
    cond_sig = 'ismatched==1'
    selection_sig = cond_sig + ' && ' + self.selection

    # define signal weights
    weight_ctau = self.tools.getCtauWeight(signal_files=signal_files, ctau=ctau)
    lhe_efficiency = 0.08244 #FIXME
    weight_signal = self.tools.getSignalWeight(signal_files=signal_files, mass=mass, ctau=ctau, sigma_B=sigma_B, lumi=lumi, lhe_efficiency=lhe_efficiency)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
    if add_weight_hlt: weight_sig += ' * ({})'.format(weight_hlt)
    if add_weight_pu: weight_sig += ' * ({})'.format(weight_pusig)
    #print 'weight ',weight_sig

    # create histogram
    hist_name = 'hist_signal_{}'.format(isBc)
    hist = ROOT.TH1D(hist_name, hist_name, 100, fit_window_min, fit_window_max)
    branch_name = 'hnl_mass'
    tree_sig.Project(hist_name, branch_name , '({sel}) * ({wght})'.format(sel=selection_sig, wght=weight_sig))

    # get the number of yields
    n_sig = hist.Integral()

    return n_sig



if __name__ == '__main__':

  data_samples = data_samples['V09_06Nov21']
  qcd_samples = qcd_samples['V09_06Nov21']
  signal_file = signal_samples['central_V09_06Nov21_m3'][0]

  baseline_selection = selection['study_Nov21'].flat
  #baseline_selection = 'hnl_pt>0' # && mu_isdsa==1' #selection[''].flat
  categories = categories['inclusive']
  #categories = categories['category_study_combined_dsa']

  ABCD_regions = ABCD_regions['cos2d_svprob_0p996']

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
      print 'cross ratio, {}, {} $\pm$ {} &'.format(category.label, round(cross_ratio[0], 2), round(cross_ratio[1], 2))


  if do_testClosureABCD:
    ROOT.gROOT.SetBatch(True)
    for category in categories:
      outlabel = 'bmass_hnlcharge_0_1_5' #TODO
      #yields = ComputeYields(data_files=data_samples[0], qcd_files=qcd_samples, selection=baseline_selection+' && '+category.definition_flat+' && '+category.cutbased_selection, white_list=white_list_20to300)
      yields = ComputeYields(data_files=data_samples[0], qcd_files=qcd_samples, signal_file=signal_file, selection=baseline_selection+' && '+category.definition_flat, white_list=white_list_20to300)
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
 
