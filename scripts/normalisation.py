import os
import sys
from os import path

import ROOT
import matplotlib.pyplot as plt # this library is used to load the correct version of libpng
import math

from tools import Tools
from mva_tools import MVATools

sys.path.append('../objects')
from samples import signal_samples
from baseline_selection import selection
from categories import categories


class NormalisationStudy(object):
  def __init__(self, mass, ctaus, baseline_selection, categories, resolution_p0, resolution_p1, add_weight_hlt, add_weight_pu, add_weight_muid, weight_hlt, weight_pusig, weight_mu0id, weight_muid):
    self.tools = Tools()
    self.mva_tools = MVATools(path_mva='./mva/outputs')
    self.mass = mass
    self.ctaus = ctaus
    self.baseline_selection = baseline_selection
    self.categories = categories
    self.resolution_p0 = resolution_p0
    self.resolution_p1 = resolution_p1
    #self.resolution = resolution_p0 + self.mass * resolution_p1
    self.add_weight_hlt = add_weight_hlt
    self.add_weight_pu = add_weight_pu
    self.add_weight_muid = add_weight_muid
    self.weight_hlt = weight_hlt
    self.weight_pusig = weight_pusig
    self.weight_mu0id = weight_mu0id
    self.weight_muid = weight_muid
    self.lumi_target = 41.6
    self.sigma_B = 472.8e9
    self.lhe_efficiency = 0.08244 
    self.training_label = 'V13_06Feb23_2023Apr06_14h13m31s'
    self.cut_score = 0.99
    self.outputdir = './myPlots/normalisation'

    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))

    ROOT.gStyle.SetPadLeftMargin(0.15) 


  def getRegion(self, mass, nsigma=2):
    resolution = self.resolution_p0 + mass * self.resolution_p1
    bin_min = mass - nsigma * resolution
    bin_max = mass + nsigma * resolution
    return bin_min, bin_max


  def getSignalWeight(self, signal_files, mass, ctau, sigma_B, lumi, lhe_efficiency=0.08244, is_bu=False, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False):
    '''
      weight = sigma_B * lumi * v_square * BR(B->muNX) * BR(N->mupi) * filter_eff / N_mini 
    '''
    # coupling square
    v_square = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)

    # production branching ratio
    from decays import Decays 
    dec = Decays(mass=mass, mixing_angle_square=1) # we factorise the mixing angle 
    if is_bu:
      BR_prod = dec.BR_B_mu / 0.4
    elif is_bd:
      BR_prod = dec.BR_B0_mu / 0.4
    elif is_bs:
      BR_prod = dec.BR_Bs_mu / 0.1
    elif is_bc:
      BR_prod = dec.BR_Bc_mu
    else:
      BR_prod = dec.BR_tot_mu  #/ 0.4 #FIXME this is to test closure at high mass

    # decay branching ratio
    BR_NToMuPi = self.tools.gamma_partial(mass=mass, vv=1.) / self.tools.gamma_total(mass=mass, vv=1.) # the coupling cancels in the ratio

    # number of generated events (= n_gen / filter_efficiency = n_miniaod / filter_efficiency)
    n_gen_tot = 0
    if not is_exclusive: # sample generated for all Bq at once
      if is_bu:
        for signal_file in signal_files:
          n_gen_tot += signal_file.n_miniaod_Bu
      elif is_bd:
        for signal_file in signal_files:
          n_gen_tot += signal_file.n_miniaod_Bd
      elif is_bs:
        for signal_file in signal_files:
          n_gen_tot += signal_file.n_miniaod_Bs
      elif is_bc:
        tree_run = ROOT.TChain('run_tree') 
        for signal_file in signal_files:
          tree_run.Add(signal_file.filename_Bc)
        n_gen_tot = self.tools.getNminiAODEvts(tree_run)
      else: # previous normalisation strategy
        tree_run = ROOT.TChain('run_tree') 
        for signal_file in signal_files:
          signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
          tree_run.Add(signal_filename)
        n_gen_tot = self.tools.getNminiAODEvts(tree_run)
    else: # sample generated for Bq separately
      tree_run = ROOT.TChain('run_tree') 
      for signal_file in signal_files:
        signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
        tree_run.Add(signal_filename)
      n_gen_tot = self.tools.getNminiAODEvts(tree_run)

    # take the weighted average of the filter efficiencies of the generated samples
    filter_efficiency = 0.
    for signal_file in signal_files:
      if not is_exclusive: # sample generated for all Bq at once
        if is_bu:
          the_filter_efficiency = signal_file.filter_efficiency_Bu
          n_gen = signal_file.n_miniaod_Bu
        elif is_bd:
          the_filter_efficiency = signal_file.filter_efficiency_Bd
          n_gen = signal_file.n_miniaod_Bd
        elif is_bs:
          the_filter_efficiency = signal_file.filter_efficiency_Bs
          n_gen = signal_file.n_miniaod_Bs
        elif is_bc:
          the_filter_efficiency = signal_file.filter_efficiency_Bc
          the_file = self.tools.getRootFile(signal_file.filename_Bc)
          the_tree_run = self.tools.getTree(the_file, 'run_tree')
          n_gen = self.tools.getNminiAODEvts(the_tree_run)
        else: # previous normalisation strategy
          the_filter_efficiency = signal_file.filter_efficiency if not is_bc else signal_file.filter_efficiency_Bc
          signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
          the_file = self.tools.getRootFile(signal_filename)
          the_tree_run = self.tools.getTree(the_file, 'run_tree')
          n_gen = self.tools.getNminiAODEvts(the_tree_run)
      else: # sample generated for Bq separately
        the_filter_efficiency = signal_file.filter_efficiency if not is_bc else signal_file.filter_efficiency_Bc
        signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
        the_file = self.tools.getRootFile(signal_filename)
        the_tree_run = self.tools.getTree(the_file, 'run_tree')
        n_gen = self.tools.getNminiAODEvts(the_tree_run)

      filter_efficiency += n_gen * the_filter_efficiency
    filter_efficiency = filter_efficiency / n_gen_tot
    #print 'filter efficiency {}'.format(filter_efficiency)

    # in the case where the samples were produced with both the electron and muon channels, apply a correction
    if not is_exclusive and not is_bu and not is_bd and not is_bs:
      corr = signal_file.muon_rate # only consider events that were generated in the muon channel 
    else:
      corr = 1. # number of events are fetched for the muon channel only

    efficiency = filter_efficiency if not is_bc else filter_efficiency * lhe_efficiency 
    n_generated = corr * n_gen_tot / efficiency
    print 'n_generated = corr * n_gen_tot / efficiency'
    print 'n_generated = {} * {} / {}'.format(corr, n_gen_tot, efficiency)

    f_u = 0.4 # B fragmentation fraction

    weight = sigma_B / f_u * lumi * v_square * BR_prod * BR_NToMuPi / n_generated 
    print 'sigma_B / f_u * lumi * v_square * BR_prod * BR_NToMuPi / n_generated'
    print 'weight = {} / 0.4 * {} * {} * {} * {} / {}'.format(sigma_B, lumi, v_square, BR_prod, BR_NToMuPi, n_generated)

    return weight, n_gen_tot


  def getCtauWeight(self, signal_files, ctau, strategy='new', is_bu=False, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False):
    if strategy == 'new':
      # get the total number of gen (=miniaod) events
      n_miniaod_tot = 0
      if not is_exclusive:
        if is_bu:
          for signal_file in signal_files:
            n_miniaod_tot += signal_file.n_miniaod_Bu
        elif is_bd:
          for signal_file in signal_files:
            n_miniaod_tot += signal_file.n_miniaod_Bd
        elif is_bs:
          for signal_file in signal_files:
            n_miniaod_tot += signal_file.n_miniaod_Bs
        elif is_bc:
          tree_run = ROOT.TChain('run_tree') 
          for signal_file in signal_files:
            tree_run.Add(signal_file.filename_Bc)
          n_miniaod_tot = self.tools.getNminiAODEvts(tree_run)
        else:
          tree_run = ROOT.TChain('run_tree') 
          for signal_file in signal_files:
            tree_run.Add(signal_file.filename)
          n_miniaod_tot = self.tools.getNminiAODEvts(tree_run)

      else:
        tree_run = ROOT.TChain('run_tree') 
        for signal_file in signal_files:
          signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
          tree_run.Add(signal_filename)
        n_miniaod_tot = self.tools.getNminiAODEvts(tree_run)

      filter_efficiency_avg = 0.
      for signal_file in signal_files:
        if not is_exclusive:
          if is_bu:
            the_filter_efficiency = signal_file.filter_efficiency_Bu
            n_gen = signal_file.n_miniaod_Bu
          elif is_bd:
            the_filter_efficiency = signal_file.filter_efficiency_Bd
            n_gen = signal_file.n_miniaod_Bd
          elif is_bs:
            the_filter_efficiency = signal_file.filter_efficiency_Bs
            n_gen = signal_file.n_miniaod_Bs
          elif is_bc:
            the_filter_efficiency = signal_file.filter_efficiency_Bc
            the_file = self.tools.getRootFile(signal_file.filename_Bc)
            the_tree_run = self.tools.getTree(the_file, 'run_tree')
            n_gen = self.tools.getNminiAODEvts(the_tree_run)
          else:
            the_filter_efficiency = signal_file.filter_efficiency
            the_file = self.tools.getRootFile(signal_file.filename)
            the_tree_run = self.tools.getTree(the_file, 'run_tree')
            n_gen = self.tools.getNminiAODEvts(the_tree_run)

        else:
          the_filter_efficiency = signal_file.filter_efficiency if not is_bc else signal_file.filter_efficiency_Bc
          signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
          the_file = self.tools.getRootFile(signal_filename)
          the_tree_run = self.tools.getTree(the_file, 'run_tree')
          n_gen = self.tools.getNminiAODEvts(the_tree_run)

        filter_efficiency_avg += n_gen * the_filter_efficiency
      filter_efficiency_avg = filter_efficiency_avg / n_miniaod_tot

      deno_weight = ''
      for ifile, signal_file in enumerate(signal_files):
        if not is_exclusive:
          if is_bu:
            the_filter_efficiency = signal_file.filter_efficiency_Bu
            n_miniaod = signal_file.n_miniaod_Bu
          elif is_bd:
            the_filter_efficiency = signal_file.filter_efficiency_Bd
            n_miniaod = signal_file.n_miniaod_Bd
          elif is_bs:
            the_filter_efficiency = signal_file.filter_efficiency_Bs
            n_miniaod = signal_file.n_miniaod_Bs
          elif is_bc:
            the_filter_efficiency = signal_file.filter_efficiency_Bc
            the_file = self.tools.getRootFile(signal_file.filename_Bc)
            the_tree_run = self.tools.getTree(the_file, 'run_tree')
            n_miniaod = self.tools.getNminiAODEvts(the_tree_run)
          else:
            the_filter_efficiency = signal_file.filter_efficiency
            the_file = self.tools.getRootFile(signal_file.filename)
            the_tree_run = self.tools.getTree(the_file, 'run_tree')
            n_miniaod = self.tools.getNminiAODEvts(the_tree_run)
        else:
          the_filter_efficiency = signal_file.filter_efficiency if not is_bc else signal_file.filter_efficiency_Bc
          signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
          the_file = self.tools.getRootFile(signal_filename)
          the_tree_run = self.tools.getTree(the_file, 'run_tree')
          n_miniaod = self.tools.getNminiAODEvts(the_tree_run)

        if ifile == 0:
          deno_weight += ' {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod / the_filter_efficiency,
                ctau0 = signal_file.ctau,
                )
        else:
          deno_weight += ' + {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod / the_filter_efficiency,
                ctau0 = signal_file.ctau,
                )
      weight_ctau = '({ntot} / {ctau1} * exp(-gen_hnl_ct / {ctau1}) * (1. / ({deno_weight})))'.format(
          ntot = n_miniaod_tot / filter_efficiency_avg,
          ctau1 = ctau,
          deno_weight = deno_weight,
          )

    elif strategy == 'old':
      signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
      if signal_file.is_private:
        original_ctau = signal_filename[signal_filename.find('ctau')+4:signal_filename.find('/', signal_filename.find('ctau')+1)]
      else:
        original_ctau = signal_filename[signal_filename.find('ctau')+4:signal_filename.find('mm', signal_filename.find('ctau')+1)]
      original_ctau = original_ctau.replace('p', '.')
      target_ctau = signal_file.ctau

      weight_ctau = 1. #-99
      if float(original_ctau) != float(target_ctau):
        ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

    return weight_ctau


  def getMCWeight(self, signal_files, mass, ctau, lumi, is_bu=False, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False):
    weight_sig_list = ['gen_hnl_ct']
    #weight_ctau = '1.'
    weight_ctau = self.getCtauWeight(signal_files=signal_files, ctau=ctau, is_bu=is_bu, is_bd=is_bd, is_bs=is_bs, is_bc=is_bc, is_exclusive=is_exclusive)
    weight_signal, n_gen = self.getSignalWeight(signal_files=signal_files, mass=mass, ctau=ctau, sigma_B=self.sigma_B, lumi=lumi, lhe_efficiency=self.lhe_efficiency, is_bu=is_bu, is_bd=is_bd, is_bs=is_bs, is_bc=is_bc, is_exclusive=is_exclusive)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
    if self.add_weight_hlt: 
      weight_sig += ' * ({})'.format(self.weight_hlt)
      weight_sig_list.append(self.weight_hlt)
    if self.add_weight_pu: 
      weight_sig += ' * ({})'.format(self.weight_pusig)
      weight_sig_list.append(self.weight_pusig)
    if self.add_weight_muid: 
      weight_sig += ' * ({}) * ({})'.format(self.weight_mu0id, self.weight_muid)
      weight_sig_list.append(self.weight_mu0id)
      weight_sig_list.append(self.weight_muid)

    #return weight_sig, weight_sig_list
    return weight_sig, weight_sig_list, n_gen


  def getSignalYields(self, signal_label, mass, ctau, category, selection, is_bu=False, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False):
    '''
      Returns the normalised number of signal yields in the fit window
    '''
    # define ranges and binning
    fit_window_min, fit_window_max = self.getRegion(mass=mass, nsigma=3)

    # get the signal files
    signal_files = self.tools.getSignalFileList(signal_label=signal_label, mass=mass, ctau=ctau, strategy='inclusive', is_bc=False)

    # define selection
    selection_sig = selection
    
    # define signal weights
    weight_sig, weight_sig_list, n_gen = self.getMCWeight(signal_files=signal_files, mass=mass, ctau=ctau, lumi=self.lumi_target, is_bu=is_bu, is_bd=is_bd, is_bs=is_bs, is_bc=is_bc, is_exclusive=is_exclusive)

    # get the tree
    treename = 'signal_tree'
    score_label = 'normalisation_m_{}_ctau_{}_{}'.format(mass, ctau, category.label)
    filename_sig = self.mva_tools.getFileWithScore(files=signal_files, training_label=self.training_label, do_parametric=True, mass=mass, category_label=category.label, selection=selection_sig, weights=weight_sig_list, label=score_label, treename=treename, is_bc=is_bc, force_overwrite=True) 
    file_sig = self.tools.getRootFile(filename_sig)
    tree_sig = self.tools.getTree(file_sig, treename)

    # create histogram
    hist_name = 'hist_signal'
    hist = ROOT.TH1D(hist_name, hist_name, 100, fit_window_min, fit_window_max)
    branch_name = 'hnl_mass'
    tree_sig.Project(hist_name, branch_name , '({sel}) * ({wght})'.format(sel=self.mva_tools.getScoreSelection(selection_sig), wght=weight_sig))

    # get the number of yields
    print 'efficiency = {}'.format(float(hist.GetEntries()) / float(n_gen)) 
    n_sig = hist.Integral()

    print 'n_sig = {}'.format(n_sig)

    return n_sig

    
  def studyClosure(self):
    '''
      Compare the signal yields between the inclusive and exclusive samples
    '''
    signal_label_inclusive = 'V42_study_normalisation'
    signal_label_exclusive_Bu = 'V44_study_normalisation_Bu' 
    signal_label_exclusive_Bd = 'V44_study_normalisation_Bd' 
    signal_label_exclusive_Bs = 'V44_study_normalisation_Bs' 

    for category in self.categories:
      #if self.mass < 3 and '_Bc' in category.label: continue
      if '_Bc' in category.label: continue
      if 'gt150_OS' not in category.label: continue

      graph_Bu = ROOT.TGraphAsymmErrors()
      graph_Bd = ROOT.TGraphAsymmErrors()
      graph_Bs = ROOT.TGraphAsymmErrors()
      graph_tot = ROOT.TGraphAsymmErrors()

      for ctau in self.ctaus:
        # compute signal yields with exclusive samples
        print '\nexclusive Bu'
        selection_exclusive = 'ismatched==1 && {} && {} && score>{}'.format(self.baseline_selection, category.definition_flat, self.cut_score)
        signal_yields_exclusive_Bu = self.getSignalYields(signal_label_exclusive_Bu, self.mass, ctau, category, selection_exclusive, is_bu=True, is_bd=False, is_bs=False, is_bc=False, is_exclusive=True)

        print '\nexclusive Bd'
        signal_yields_exclusive_Bd = self.getSignalYields(signal_label_exclusive_Bd, self.mass, ctau, category, selection_exclusive, is_bu=False, is_bd=True, is_bs=False, is_bc=False, is_exclusive=True)

        print '\nexclusive Bs'
        signal_yields_exclusive_Bs = self.getSignalYields(signal_label_exclusive_Bs, self.mass, ctau, category, selection_exclusive, is_bu=False, is_bd=False, is_bs=True, is_bc=False, is_exclusive=True)

        signal_yields_exclusive_tot = signal_yields_exclusive_Bu + signal_yields_exclusive_Bd + signal_yields_exclusive_Bs

        # compute signal yields with inclusive sample
        print '\ninclusive Bu'
        selection_inclusive_Bu = 'ismatched==1 && {} && {} && score>{} && isbu==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
        signal_yields_inclusive_Bu = self.getSignalYields(signal_label_inclusive, self.mass, ctau, category, selection_inclusive_Bu, is_bu=True, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False)

        print '\ninclusive Bd'
        selection_inclusive_Bd = 'ismatched==1 && {} && {} && score>{} && isbd==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
        signal_yields_inclusive_Bd = self.getSignalYields(signal_label_inclusive, self.mass, ctau, category, selection_inclusive_Bd, is_bu=False, is_bd=True, is_bs=False, is_bc=False, is_exclusive=False)

        print '\ninclusive Bs'
        selection_inclusive_Bs = 'ismatched==1 && {} && {} && score>{} && isbs==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
        signal_yields_inclusive_Bs = self.getSignalYields(signal_label_inclusive, self.mass, ctau, category, selection_inclusive_Bs, is_bu=False, is_bd=False, is_bs=True, is_bc=False, is_exclusive=False)

        signal_yields_inclusive_tot = signal_yields_inclusive_Bu + signal_yields_inclusive_Bd + signal_yields_inclusive_Bs

        print 'signal_yields_exclusive_Bu = {}'.format(signal_yields_exclusive_Bu)
        print 'signal_yields_inclusive_Bu = {}'.format(signal_yields_inclusive_Bu)
        print 'exclusive/inclusive Bu = {}'.format(signal_yields_exclusive_Bu/signal_yields_inclusive_Bu)

        print '\n'
        print 'signal_yields_exclusive_Bd = {}'.format(signal_yields_exclusive_Bd)
        print 'signal_yields_inclusive_Bd = {}'.format(signal_yields_inclusive_Bd)
        print 'exclusive/inclusive Bd = {}'.format(signal_yields_exclusive_Bd/signal_yields_inclusive_Bd if signal_yields_inclusive_Bd !=0  else 1.)

        print '\n'
        print 'signal_yields_exclusive_Bs = {}'.format(signal_yields_exclusive_Bs)
        print 'signal_yields_inclusive_Bs = {}'.format(signal_yields_inclusive_Bs)
        print 'exclusive/inclusive Bs = {}'.format(signal_yields_exclusive_Bs/signal_yields_inclusive_Bs if signal_yields_inclusive_Bs !=0  else 1.)

        print '\n'
        print 'signal_yields_exclusive_tot = {}'.format(signal_yields_exclusive_tot)
        print 'signal_yields_inclusive_tot = {}'.format(signal_yields_inclusive_tot)
        print 'exclusive/inclusive tot = {}'.format(signal_yields_exclusive_tot/signal_yields_inclusive_tot)

        # impose 25% error on the yields
        error_yields = 0.25

        # fill graphs
        ratio_Bu = signal_yields_exclusive_Bu / signal_yields_inclusive_Bu
        error_Bu = ratio_Bu * math.sqrt(math.pow(error_yields * signal_yields_exclusive_Bu / signal_yields_exclusive_Bu, 2) + math.pow(error_yields * signal_yields_inclusive_Bu / signal_yields_inclusive_Bu, 2)) 
        point_Bu = graph_Bu.GetN()
        graph_Bu.SetPoint(point_Bu, ctau, ratio_Bu)
        graph_Bu.SetPointError(point_Bu, 0, 0, error_Bu, error_Bu)
        graph_Bu.SetLineColor(4)
        graph_Bu.SetMarkerColor(4)
        graph_Bu.SetMarkerStyle(20)

        ratio_Bd = signal_yields_exclusive_Bd / signal_yields_inclusive_Bd if signal_yields_inclusive_Bd != 0 else -99.
        if ratio_Bd != -99.:
          error_Bd = ratio_Bd * math.sqrt(math.pow(error_yields * signal_yields_exclusive_Bd / signal_yields_exclusive_Bd, 2) + math.pow(error_yields * signal_yields_inclusive_Bd / signal_yields_inclusive_Bd, 2)) 
        else: error_Bd = 0.
        point_Bd = graph_Bd.GetN()
        graph_Bd.SetPoint(point_Bd, ctau, ratio_Bd)
        graph_Bd.SetPointError(point_Bd, 0, 0, error_Bd, error_Bd)
        graph_Bd.SetLineColor(4)
        graph_Bd.SetMarkerColor(4)
        graph_Bd.SetMarkerStyle(20)

        ratio_Bs = signal_yields_exclusive_Bs / signal_yields_inclusive_Bs if signal_yields_inclusive_Bs != 0 else -99.
        if ratio_Bs != -99.:
          error_Bs = ratio_Bs * math.sqrt(math.pow(error_yields * signal_yields_exclusive_Bs / signal_yields_exclusive_Bs, 2) + math.pow(error_yields * signal_yields_inclusive_Bs / signal_yields_inclusive_Bs, 2)) 
        else: error_Bs = 0.
        point_Bs = graph_Bs.GetN()
        graph_Bs.SetPoint(point_Bs, ctau, ratio_Bs)
        graph_Bs.SetPointError(point_Bs, 0, 0, error_Bs, error_Bs)
        graph_Bs.SetLineColor(4)
        graph_Bs.SetMarkerColor(4)
        graph_Bs.SetMarkerStyle(20)

        ratio_tot = signal_yields_exclusive_tot / signal_yields_inclusive_tot
        error_yields_tot = error_yields #math.sqrt(3 * math.pow(error_yields, 2))
        error_tot = ratio_tot * math.sqrt(math.pow(error_yields_tot * signal_yields_exclusive_tot / signal_yields_exclusive_tot, 2) + math.pow(error_yields_tot * signal_yields_inclusive_tot / signal_yields_inclusive_tot, 2)) 
        point_tot = graph_tot.GetN()
        graph_tot.SetPoint(point_tot, ctau, ratio_tot)
        graph_tot.SetPointError(point_tot, 0, 0, error_tot, error_tot)
        graph_tot.SetLineColor(4)
        graph_tot.SetMarkerColor(4)
        graph_tot.SetMarkerStyle(20)

      # create canvases
      canv_Bu = self.tools.createTCanvas('canv_Bu_{}'.format(category.label))
      canv_Bu.SetLogx()
      canv_Bd = self.tools.createTCanvas('canv_Bd_{}'.format(category.label))
      canv_Bd.SetLogx()
      canv_Bs = self.tools.createTCanvas('canv_Bs_{}'.format(category.label))
      canv_Bs.SetLogx()
      canv_tot = self.tools.createTCanvas('canv_tot_{}'.format(category.label))
      canv_tot.SetLogx()

      # create frame
      x_min = 1. #6.5
      x_max = 1000. #1050
      frame = ROOT.TGraph()
      frame.SetPoint(0, x_min, 0)
      frame.SetPoint(1, x_max, 1.8)
      frame.SetMarkerStyle(0)
      frame.SetMarkerSize(0)
      frame.SetMarkerColor(0)
      frame.GetXaxis().SetTitle('c#tau (mm)')
      frame.GetXaxis().SetLabelSize(0.037)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetTitle('signal yields exclusive / signal yields inclusive')
      frame.GetYaxis().SetLabelSize(0.037)
      frame.GetYaxis().SetTitleSize(0.042)
      frame.GetYaxis().SetTitleOffset(1.5)

      # legend
      text = 'pNN score > {}'.format(self.cut_score)

      # draw line
      line = ROOT.TLine(x_min, 1, x_max, 1)
      line.SetLineColor(2)
      line.SetLineWidth(2)
      line.SetLineStyle(9)

      # draw graphs
      canv_Bu.cd()
      frame.SetTitle('Bu signal')
      frame.Draw('AP')
      graph_Bu.Draw('P same')
      line.Draw('same')
      self.tools.printLatexBox(0.75, 0.8, category.title, size=0.04, pos='center', font=42)
      self.tools.printLatexBox(0.75, 0.7, text, size=0.04, pos='center', font=42)
      name_Bu = 'graph_Bu_{}_score{}'.format(category.label, str(self.cut_score).replace('.', 'p'))
      canv_Bu.SaveAs('{}/{}.png'.format(self.outputdir, name_Bu))
      canv_Bu.SaveAs('{}/{}.pdf'.format(self.outputdir, name_Bu))

      canv_Bd.cd()
      frame.SetTitle('Bd signal')
      frame.Draw('AP')
      graph_Bd.Draw('P same')
      line.Draw('same')
      self.tools.printLatexBox(0.75, 0.8, category.title, size=0.04, pos='center', font=42)
      self.tools.printLatexBox(0.75, 0.7, text, size=0.04, pos='center', font=42)
      name_Bd = 'graph_Bd_{}_score{}'.format(category.label, str(self.cut_score).replace('.', 'p'))
      canv_Bd.SaveAs('{}/{}.png'.format(self.outputdir, name_Bd))
      canv_Bd.SaveAs('{}/{}.pdf'.format(self.outputdir, name_Bd))

      canv_Bs.cd()
      frame.SetTitle('Bs signal')
      frame.Draw('AP')
      graph_Bs.Draw('P same')
      line.Draw('same')
      self.tools.printLatexBox(0.75, 0.8, category.title, size=0.04, pos='center', font=42)
      self.tools.printLatexBox(0.75, 0.7, text, size=0.04, pos='center', font=42)
      name_Bs = 'graph_Bs_{}_score{}'.format(category.label, str(self.cut_score).replace('.', 'p'))
      canv_Bs.SaveAs('{}/{}.png'.format(self.outputdir, name_Bs))
      canv_Bs.SaveAs('{}/{}.pdf'.format(self.outputdir, name_Bs))

      canv_tot.cd()
      frame.SetTitle('(Bu + Bd + Bs) signal')
      frame.Draw('AP')
      graph_tot.Draw('P same')
      line.Draw('same')
      self.tools.printLatexBox(0.75, 0.8, category.title, size=0.04, pos='center', font=42)
      self.tools.printLatexBox(0.75, 0.7, text, size=0.04, pos='center', font=42)
      name_tot = 'graph_tot_{}_score{}'.format(category.label, str(self.cut_score).replace('.', 'p'))
      canv_tot.SaveAs('{}/{}.png'.format(self.outputdir, name_tot))
      canv_tot.SaveAs('{}/{}.pdf'.format(self.outputdir, name_tot))
    

  def studyYields(self):
    '''
      Compare the signal yields between the old and new normalisation strategies
    '''
    #signal_labels = ['V42_study_normalisation_m4p5'] 
    #masses = [4.5]

    signal_labels = ['V42_study_normalisation_m1', 'V42_study_normalisation_m1p5', 'V42_study_normalisation', 'V42_study_normalisation_m3', 'V42_study_normalisation_m4p5'] 
    masses = [1., 1.5, 2., 3., 4.5]
    colours = [ROOT.kBlue, ROOT.kRed, ROOT.kMagenta, ROOT.kGreen+2, ROOT.kCyan-7]

    # impose cut on the score
    self.cut_score = 0.99

    for category in self.categories:
      if '_Bc' in category.label: continue
      #if 'gt150_SS' not in category.label: continue

      graphs = []
      leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.15, xmax=0.85, ymax=0.37, size=0.04)

      for imass, mass in enumerate(masses):
        graph_ratio = ROOT.TGraph()

        for ctau in self.ctaus:
          # previous normalisation strategy
          print '\nprevious strategy'
          selection_old = 'ismatched==1 && {} && {} && score>{}'.format(self.baseline_selection, category.definition_flat, self.cut_score)
          signal_yields_old = self.getSignalYields(signal_labels[imass], mass, ctau, category, selection_old, is_bu=False, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False)

          # revisited normalisation strategy
          print '\nrevisited strategy - Bu'
          selection_Bu = 'ismatched==1 && {} && {} && score>{} && isbu==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
          signal_yields_Bu = self.getSignalYields(signal_labels[imass], mass, ctau, category, selection_Bu, is_bu=True, is_bd=False, is_bs=False, is_bc=False, is_exclusive=False)

          print '\nrevisited strategy - Bd'
          selection_Bd = 'ismatched==1 && {} && {} && score>{} && isbd==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
          signal_yields_Bd = self.getSignalYields(signal_labels[imass], mass, ctau, category, selection_Bd, is_bu=False, is_bd=True, is_bs=False, is_bc=False, is_exclusive=False)

          print '\nrevisited strategy - Bs'
          selection_Bs = 'ismatched==1 && {} && {} && score>{} && isbs==1'.format(self.baseline_selection, category.definition_flat, self.cut_score)
          signal_yields_Bs = self.getSignalYields(signal_labels[imass], mass, ctau, category, selection_Bs, is_bu=False, is_bd=False, is_bs=True, is_bc=False, is_exclusive=False)

          signal_yields_new  = signal_yields_Bu + signal_yields_Bd + signal_yields_Bs

          ratio_yields = signal_yields_new / signal_yields_old if signal_yields_old != 0 else -99. 

          print '\n'
          print 'previous strategy: {}'.format(signal_yields_old)
          print 'revisited strategy: {}'.format(signal_yields_new)
          print 'Bu + Bd + Bs = {}  + {}  + {}'.format(signal_yields_Bu, signal_yields_Bd, signal_yields_Bs)
          print 'ratio revisited / previous: {}'.format(ratio_yields)

          # fill graph
          point = graph_ratio.GetN()
          graph_ratio.SetPoint(point, ctau, ratio_yields)
          #graph_ratio.SetPointError(point, 0, 0, error_Bu, error_Bu)
          graph_ratio.SetLineColor(colours[imass])
          graph_ratio.SetMarkerColor(colours[imass])
          graph_ratio.SetMarkerStyle(20)
          graphs.append(graph_ratio)
        leg.AddEntry(graph_ratio, '{} GeV'.format(mass))

      canv = self.tools.createTCanvas('canv_{}'.format(category.label))
      canv.SetLogx()
      canv.SetGridy()

      # create frame
      x_min = 1.
      x_max = 1000.
      frame = ROOT.TGraph()
      frame.SetPoint(0, x_min, 0)
      frame.SetPoint(1, x_max, 4.0)
      frame.SetMarkerStyle(0)
      frame.SetMarkerSize(0)
      frame.SetMarkerColor(0)
      frame.GetXaxis().SetTitle('c#tau (mm)')
      frame.GetXaxis().SetLabelSize(0.037)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetTitle('signal yields new / old')
      frame.GetYaxis().SetLabelSize(0.037)
      frame.GetYaxis().SetTitleSize(0.042)
      frame.GetYaxis().SetTitleOffset(1.5)

      # legend
      text = 'pNN score > {}'.format(self.cut_score)

      # draw graph
      canv.cd()
      frame.SetTitle('')
      frame.Draw('AP')
      for graph_ratio in graphs:
        graph_ratio.Draw('PL same')
      self.tools.printLatexBox(0.75, 0.8, category.title, size=0.04, pos='center', font=42)
      self.tools.printLatexBox(0.75, 0.7, text, size=0.04, pos='center', font=42)
      leg.Draw()
      name = 'graph_ratio_yields_{}_score{}'.format(category.label, str(self.cut_score).replace('.', 'p'))
      canv.SaveAs('{}/{}.png'.format(self.outputdir, name))
      canv.SaveAs('{}/{}.pdf'.format(self.outputdir, name))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  do_studyClosure = False
  do_studyYields = True

  resolution_p0 = 6.98338e-04 
  resolution_p1 = 7.78382e-03 

  baseline_selection = selection['baseline_08Aug22'].flat
  categories = categories['categories_0_50_150_Bc']

  add_weight_hlt = True
  add_weight_pu = True
  add_weight_muid = True
  weight_hlt = 'weight_hlt_fullBpark'
  weight_pusig = 'weight_pu_sig_tot'
  weight_mu0id = 'weight_mu0_softid'
  weight_muid = 'weight_mu_looseid'

  mass = 2.
  #ctaus = [10.]
  ctaus = [1., 10., 100., 1000.]
  #ctaus = [7., 10., 15., 20., 30., 40., 50., 70., 100., 150., 200., 300., 400., 500., 700., 1000.]

  analyser = NormalisationStudy(
      mass = mass, 
      ctaus = ctaus,
      baseline_selection = baseline_selection,
      categories = categories,
      resolution_p0 = resolution_p0,
      resolution_p1 = resolution_p1,
      add_weight_hlt = add_weight_hlt,
      add_weight_pu = add_weight_pu,
      add_weight_muid = add_weight_muid,
      weight_hlt = weight_hlt,
      weight_pusig = weight_pusig,
      weight_mu0id = weight_mu0id,
      weight_muid = weight_muid,
      )

  if do_studyClosure:
    analyser.studyClosure()

  if do_studyYields:
    analyser.studyYields()





