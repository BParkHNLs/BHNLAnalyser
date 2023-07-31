import ROOT
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import sys
from glob import glob
import numpy as np
from os import path
import math
from tools import Tools
from mva_tools import MVATools
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples
from categories import categories
from baseline_selection import selection
from quantity import Quantity
from ctau_points import ctau_points

import pandas as pd
from root_pandas import read_root


class Systematics(Tools):
  def __init__(self, signal_labels='', categories='', baseline_selection='', lumi=41.6):
    self.tools = Tools()
    self.mva_tools = MVATools()
    self.signal_labels = signal_labels
    self.categories = categories
    self.baseline_selection = baseline_selection
    self.lumi = lumi

  def studyHLTsyst(self, weight_hlt):
    print '\n',weight_hlt
    for signal_label in self.signal_labels:
      signal_files = signal_samples[signal_label]

      for signal_file in signal_files:
        signal_mass = signal_file.mass
        signal_ctau = signal_file.ctau

        signal_selection = 'ismatched==1' if self.baseline_selection=='' else 'ismatched==1 && {}'.format(self.baseline_selection)
        signal_yields_no_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=False)[0] 
        signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, weight_hlt=weight_hlt)[0]
        signal_yields_weight_plus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, weight_hlt=weight_hlt+'_plus_one_sigma')[0]
        signal_yields_weight_minus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, weight_hlt=weight_hlt+'_minus_one_sigma')[0]

        syst_up = round((signal_yields_weight_plus_one_sigma / signal_yields_weight -1 ) * 100, 2)
        syst_down = round((signal_yields_weight_minus_one_sigma / signal_yields_weight -1) * 100, 2)

        print 'm{}\_ctau{} & {} & {} & {} & {} & {}\% & {}\% \\\ '.format(signal_mass, signal_ctau, int(signal_yields_no_weight), int(signal_yields_weight), int(signal_yields_weight_plus_one_sigma), int(signal_yields_weight_minus_one_sigma), syst_up, syst_down)


  def plotOccupancy(self, muon_type='primary'):
    if muon_type not in ['primary', 'displaced']:
      raise RuntimeError('Unknown muon typei "{}". Choose among ["primary", "displaced"]'.format(muon_type))

    pt_bins = ['6.0', '7.0', '8.0', '8.5', '9.0', '10.0', '10.5', '11.0', '12.0', '20.0', '100.0']
    eta_bins = ['0.0', '0.5', '1.0', '1.5', '2.0']
    dxysig_bins = ['0.0', '3.0', '4.0', '6.0', '8.0', '10.0', '20.0', '500.0']
    
    quantity = Quantity(name_flat='hnl_mass', nbins=30, bin_min=0, bin_max=13000)

    for signal_label in self.signal_labels:
      signal_files = signal_samples[signal_label]

      for signal_file in signal_files:
        #if signal_file.ctau != 1000.: continue
        f = self.tools.getRootFile(signal_file.filename)
        tree = self.tools.getTree(f, 'signal_tree')

        canv_name = 'canv'+signal_file.filename
        canv = self.tools.createTCanvas(canv_name)
        canv.SetRightMargin(0.15)
        
        hist2d_name = 'hist2d' + signal_file.filename
        hist2d = ROOT.TH2D(hist2d_name, hist2d_name, len(pt_bins)-1, 0, len(pt_bins)-1, len(dxysig_bins)-1, 0, len(dxysig_bins)-1) 

        hist_incl_name = 'hist_incl'
        selection_incl = 'ismatched == 1'
        if muon_type == 'primary':
          #selection_incl = 'ismatched == 1 && mu0_istriggering==1'
          selection_incl = 'ismatched == 1'
        elif muon_type == 'displaced':
          #selection_incl = 'ismatched == 1 && mu_istriggering==1'
          selection_incl = 'ismatched == 1'
        weight_incl = '(1)'
        hist_incl = self.tools.createHisto(tree, quantity, hist_incl_name, weight_incl, selection_incl) 
        entries_incl = hist_incl.GetEntries()

        tot = 0

        for ipt, pt_bin in enumerate(pt_bins):
          if ipt == len(pt_bins)-1: continue
          for idxysig, dxysig_bin in enumerate(dxysig_bins):
            if idxysig == len(dxysig_bins)-1: continue
            pt_min = pt_bins[ipt]
            pt_max = pt_bins[ipt+1]
            dxysig_min = dxysig_bins[idxysig]
            dxysig_max = dxysig_bins[idxysig+1]

            hist_name = 'hist' + pt_bin + dxysig_bin
            if muon_type == 'primary':
              #selection = 'ismatched == 1 && mu0_istriggering==1 && mu0_pt > {} && mu0_pt < {} && abs(mu0_dxysig) > {} && abs(mu0_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
              selection = 'ismatched == 1 && mu0_pt > {} && mu0_pt < {} && abs(mu0_dxysig) > {} && abs(mu0_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
            elif muon_type == 'displaced':
              #selection = 'ismatched == 1 && mu_istriggering==1 && mu_pt > {} && mu_pt < {} && abs(mu_dxysig) > {} && abs(mu_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
              selection = 'ismatched == 1 && mu_pt > {} && mu_pt < {} && abs(mu_dxysig) > {} && abs(mu_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
            weight = '(1)'
            hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
            entries = hist.GetEntries()
            tot = tot + entries
            hist2d.Fill(pt_bins[ipt+1], dxysig_bins[idxysig+1], entries/entries_incl)
            #print '{} {} {} {} {}'.format(pt_min, pt_max, dxysig_min, dxysig_max, entries/entries_incl)
  
        print 'check: {}'.format(tot/entries_incl)

        if muon_type == 'primary':
          x_label = 'primary muon p_{T} [GeV]'
          y_label = 'primary muon d_{xy}/#sigma_{xy}'
        elif muon_type == 'displaced':
          x_label = 'displaced muon p_{T} [GeV]'
          y_label = 'displaced muon d_{xy}/#sigma_{xy}'

        hist2d.SetTitle('mass {} GeV, ctau {} mm'.format(signal_file.mass, signal_file.ctau))
        hist2d.GetXaxis().SetTitle(x_label)
        hist2d.GetXaxis().SetLabelSize(0.05)
        hist2d.GetXaxis().SetTitleSize(0.042)
        hist2d.GetXaxis().SetTitleOffset(1.1)
        hist2d.GetYaxis().SetTitle(y_label)
        hist2d.GetYaxis().SetLabelSize(0.05)
        hist2d.GetYaxis().SetTitleSize(0.042)
        hist2d.GetYaxis().SetTitleOffset(1.1)
        hist2d.GetZaxis().SetTitle("Occupancy")
        hist2d.GetZaxis().SetTitleOffset(1.5)
        hist2d.GetZaxis().SetRangeUser(0, 0.05)
        hist2d.SetOption("colztext");
        hist2d.Draw()
        ROOT.gStyle.SetTitleFillColor(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPaintTextFormat(".3f")

        outdir = 'myPlots/systematics/'
        if not path.exists(outdir):
          os.system('mkdir -p {}'.format(outdir))
        name = 'occupancy_m{}_ctau{}_{}_muon'.format(str(signal_file.mass).replace('.', 'p'), str(signal_file.ctau).replace('.', 'p'), muon_type)
        canv.SaveAs(outdir + name + '.png')
        
              
  def plotTriggeringTable(self):

    mu0_bins = ['#mu_{0} triggers', '#mu_{0} !triggers']
    mu_bins = ['#mu triggers', '#mu !triggers']
    
    quantity = Quantity(name_flat='hnl_mass', nbins=30, bin_min=0, bin_max=13000)

    for signal_label in self.signal_labels:
      signal_files = signal_samples[signal_label]

      for signal_file in signal_files:
        #if signal_file.ctau != 1000.: continue
        f = self.tools.getRootFile(signal_file.filename_Bc)
        tree = self.tools.getTree(f, 'signal_tree')

        canv_name = 'canv'+signal_file.filename
        canv = self.tools.createTCanvas(canv_name)
        canv.SetLeftMargin(0.15)
        canv.SetRightMargin(0.15)
        
        hist2d_name = 'hist2d' + signal_file.filename
        hist2d = ROOT.TH2D(hist2d_name, hist2d_name, 2, 0, 2, 2, 0, 2) 

        hist_incl_name = 'hist_incl'
        selection_incl = 'ismatched == 1'
        weight_incl = '(1)'
        hist_incl = self.tools.createHisto(tree, quantity, hist_incl_name, weight_incl, selection_incl) 
        entries_incl = hist_incl.GetEntries()

        tot = 0

        hist_name = 'hist00'
        #selection = 'ismatched == 1 && trgmu_istriggering==1 && mu_istriggering==1'
        selection = 'ismatched == 1 && mu0_istriggering==1 && mu_istriggering==1'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        print entries/entries_incl
        tot = tot + entries
        hist2d.Fill(mu0_bins[0], mu_bins[0], entries/entries_incl)

        hist_name = 'hist01'
        #selection = 'ismatched == 1 && trgmu_istriggering==1 && mu_istriggering==0'
        selection = 'ismatched == 1 && mu0_istriggering==1 && mu_istriggering==0'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        print entries/entries_incl
        tot = tot + entries
        hist2d.Fill(mu0_bins[0], mu_bins[1], entries/entries_incl)

        hist_name = 'hist10'
        #selection = 'ismatched == 1 && trgmu_istriggering==0 && mu_istriggering==1'
        selection = 'ismatched == 1 && mu0_istriggering==0 && mu_istriggering==1'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        print entries/entries_incl
        tot = tot + entries
        hist2d.Fill(mu0_bins[1], mu_bins[0], entries/entries_incl)

        hist_name = 'hist11'
        #selection = 'ismatched == 1 && trgmu_istriggering==0 && mu_istriggering==0'
        selection = 'ismatched == 1 && mu0_istriggering==0 && mu_istriggering==0'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        print entries/entries_incl
        tot = tot + entries
        hist2d.Fill(mu0_bins[1], mu_bins[1], entries/entries_incl)

        print 'check: {}'.format(tot/entries_incl)

        hist2d.SetTitle('mass {} GeV, ctau {} mm'.format(signal_file.mass, signal_file.ctau))
        hist2d.GetXaxis().SetLabelSize(0.05)
        hist2d.GetYaxis().SetLabelSize(0.05)
        hist2d.GetZaxis().SetTitle("Occupancy")
        hist2d.GetZaxis().SetTitleOffset(1.5)
        hist2d.GetZaxis().SetRangeUser(0, 1)
        hist2d.SetOption("colztext");
        hist2d.Draw()
        ROOT.gStyle.SetTitleFillColor(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPaintTextFormat(".3f")

        outdir = 'myPlots/systematics/'
        if not path.exists(outdir):
          os.system('mkdir -p {}'.format(outdir))
        name = 'trigger_table_m{}_ctau{}'.format(str(signal_file.mass).replace('.', 'p'), str(signal_file.ctau).replace('.', 'p'))
        canv.SaveAs(outdir + name + '.png')
        canv.SaveAs(outdir + name + '.pdf')


  def getMassFromResultFilename(self, filename):
    idx1 = filename.rfind('_m_')+3
    idx2 = filename.rfind('_ctau_')
    return filename[idx1:idx2]


  def getCouplingFromResultFilename(self, filename):
    idx1 = filename.rfind('_v2_')+4
    idx2 = filename.rfind('.txt')
    return filename[idx1:idx2]


  def getMedianValue(self, filename):
    try:
      f = open(filename)
      lines = f.readlines()
      median = -99.
      for line in lines:
        if 'Expected 50.0%' not in line: continue
        # Expected 50.0%: r < 345.0000
        idx = line.find('r < ')+4
        median = float(line[idx:])
    except: 
      median = -99.

    return median


  def getResultFile(self, path, filename):
    idx1 = filename.rfind('/')+1
    idx2 = filename.rfind('.txt')+4
    name = filename[idx1:idx2] 

    return path + '/' + name
      

  def studyShapeSyst(self):
    '''
      Study the signal parametrisation systematics based on the median value
    '''

    #path_results = '../outputs/V12_08Aug22/limits/signal_parametrisation_systematics_v1/results'
    #path_results_plusonesigma = '../outputs/V12_08Aug22/limits/signal_parametrisation_systematics_v1_plusonesigmacat/results'
    #path_results_minusonesigma = '../outputs/V12_08Aug22/limits/signal_parametrisation_systematics_v1_minusonesigmacat/results'

    #path_results = '../outputs/V12_08Aug22/limits/study_parametrisation_doubleCBFast_floating/results'
    path_results = '../outputs/V12_08Aug22/limits/study_parametrisation_doubleCBFast_inclusive/results'
    path_results_plusonesigma = '../outputs/V12_08Aug22/limits/study_parametrisation_doubleCBFast_percat/results'
    #path_results_plusonesigma = '../outputs/V12_08Aug22/limits/study_parametrisation_doubleCBFast_percat/results'
    path_results_minusonesigma = path_results_plusonesigma

    files_results = [f for f in glob(path_results + '/result_Majorana_m_*_v2_*.txt')]
    if len(files_results) == 0:
      raise RuntimeError('No limit results found under path "{}"'.format(path_results))

    #files_results_plusonesigma = [f for f in glob(path_results_plusonesigma + '/result_m_*_v2_*.txt')]
    #if len(files_results_plusonesigma) == 0:
    #  raise RuntimeError('No limit results found under path "{}"'.format(path_results_plusonesigma))

    #files_results_minusonesigma = [f for f in glob(path_results_minusonesigma + '/result_m_*_v2_*.txt')]
    #if len(files_results_minusonesigma) == 0:
    #  raise RuntimeError('No limit results found under path "{}"'.format(path_results_minusonesigma))

    for file_results in files_results:
      # get corresponding files for plus/minus one sigma 
      if 'ctau_5000.0' not in file_results and 'ctau_1500.0' not in file_results and 'ctau_300.0' not in file_results and 'ctau_15.0' not in file_results and 'ctau_0.02' not in file_results: continue
      file_results_plusonesigma = self.getResultFile(path=path_results_plusonesigma, filename=file_results)
      file_results_minusonesigma = self.getResultFile(path=path_results_minusonesigma, filename=file_results)

      median = self.getMedianValue(filename=file_results)
      median_plusonesigma = self.getMedianValue(filename=file_results_plusonesigma)
      median_minusonesigma = self.getMedianValue(filename=file_results_minusonesigma)
    
      if median == -99. or median_plusonesigma == -99. or median_minusonesigma == -99.: 
        continue

      mass = self.getMassFromResultFilename(filename=file_results)
      v2 = self.getCouplingFromResultFilename(filename=file_results)
      #if float(mass) not in [1, 1.5, 2, 3., 4.5]: continue
      #if float(mass) <= 3. and float(v2) > 1e-3: continue
      #if float(mass) > 3. and float(v2) < 1e-3: continue

      syst_up = (median_plusonesigma / median -1) * 100. 
      syst_down = (median_minusonesigma / median -1) * 100. 

      #print 'mass {} v2 {}: syst up = {}% syst down = {}%'.format(mass, v2, round(syst_up, 1), round(syst_down, 1))
      print '{} & {} & {}\% \\\ '.format(mass, v2, round(syst_up, 1))


  def studyShapeSystPerCategory(self):
    '''
      Study the signal parametrisation systematics based on the median value
    '''

    from categories import categories
    from ctau_points import ctau_points

    path_datacards_resoincl = '../outputs/V12_08Aug22/datacards/study_parametrisation_v2'
    path_datacards_resopercat = '../outputs/V12_08Aug22/datacards/study_parametrisation_percat_resolution'

    categories = categories['categories_0_50_150']

    masses = [1.0, 1.5, 2.0, 3.0, 4.5]
    ctau_points = ctau_points['close_to_exclusion']

    outputdir = './myPlots/signal_parametrisation/limits'
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))

    for category in categories:
      if 'incl' in category.label: continue

      for mass in masses:
        for ctau_point_list in ctau_points:
          # get the signal mass
          if mass not in ctau_point_list.mass_list: continue
          
          for ctau_point in ctau_point_list.ctau_list:
            signal_v2 = self.tools.getVV(mass=mass, ctau=ctau_point, ismaj=True)
            signal_coupling = self.tools.getCouplingLabel(signal_v2)

            datacard_name = 'datacard_bhnl_m_{}_ctau_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(ctau_point).replace('.', 'p'), str(signal_coupling).replace('.', 'p').replace('-', 'm'), category.label) 
            result_name_resoincl = 'result_resoincl_m_{}_ctau_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(ctau_point).replace('.', 'p'), str(signal_coupling).replace('.', 'p').replace('-', 'm'), category.label) 
            if not path.exists('{}/{}'.format(outputdir, result_name_resoincl)):
              combine_command_resoincl = 'combine -M AsymptoticLimits {}/{}  --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams --run blind &> {}/{}'.format(path_datacards_resoincl, datacard_name, outputdir, result_name_resoincl) 
              #print combine_command_resoincl
              os.system(combine_command_resoincl)

            result_name_resopercat = 'result_resopercat_m_{}_ctau_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(ctau_point).replace('.', 'p'), str(signal_coupling).replace('.', 'p').replace('-', 'm'), category.label) 
            if not path.exists('{}/{}'.format(outputdir, result_name_resopercat)):
              combine_command_resopercat = 'combine -M AsymptoticLimits {}/{}  --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams --run blind &> {}/{}'.format(path_datacards_resopercat, datacard_name, outputdir, result_name_resopercat) 
              #print combine_command_resopercat
              os.system(combine_command_resopercat)

            ## get corresponding files for plus/minus one sigma 
            #file_results_plusonesigma = self.getResultFile(path=path_results_plusonesigma, filename=file_results)
            #file_results_minusonesigma = self.getResultFile(path=path_results_minusonesigma, filename=file_results)

            median_resoincl = self.getMedianValue(filename='{}/{}'.format(outputdir, result_name_resoincl))
            median_resopercat = self.getMedianValue(filename='{}/{}'.format(outputdir, result_name_resopercat))
          
            if median_resoincl == -99. or median_resopercat == -99.: 
              continue

            #mass = self.getMassFromResultFilename(filename=file_results)
            #v2 = self.getCouplingFromResultFilename(filename=file_results)
            #if float(mass) not in [1, 1.5, 2, 3., 4.5]: continue
            #if float(mass) <= 3. and float(v2) > 1e-3: continue
            #if float(mass) > 3. and float(v2) < 1e-3: continue

            syst = (median_resopercat / median_resoincl -1) * 100. 
            #syst_down = (median_minusonesigma / median -1) * 100. 

            print 'mass {} v2 {} cat {}: syst = {}%'.format(mass, signal_coupling, category.label, round(syst, 1))
            #print '{} & {} & {}\% \\\ '.format(mass, v2, round(syst_up, 1))


  def plotMCSyst(self, weight):
    print '\n',weight

    hists_up = []

    syst_min = 1e9
    syst_max = 1e-9

    leg = self.tools.getRootTLegend(xmin=0.45, ymin=0.55, xmax=0.7, ymax=0.85, size=0.035)
    canv = ROOT.TCanvas('canv', 'canv', 800, 700)

    for signal_label in self.signal_labels:
      signal_files = signal_samples[signal_label]

      for signal_file in signal_files:
        signal_mass = signal_file.mass
        signal_ctau = signal_file.ctau
        signal_resolution = 0.0002747 + signal_mass * 0.008302  
        #if signal_mass != 1 and signal_ctau != 1000: continue
      

        #sig_label = 'sig_{}_{}_{}'.format(signal_mass, signal_ctau, category.label)
        #filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=signal_mass, selection=self.baseline_selection, weights=extra_branches+['gen_hnl_ct', self.weight_hlt, self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=False) 
        #file_sig = self.tools.getRootFile(filename_sig)
        #tree_sig = self.tools.getTree(file_sig, 'signal_tree')

        #hist_up = ROOT.TH1D('hist_up', 'hist_up', len(self.categories), 0, len(self.categories)) 
        hist_up = ROOT.TH1D('hist_up', 'hist_up', len(self.categories)-1, 0, len(self.categories)-1) 
        hist_up.SetMarkerSize(2)
        hist_up.SetMarkerStyle(20)
        hist_up.SetMarkerColor(signal_file.colour)
        hist_up.GetXaxis().SetLabelSize(0.033)
        hist_up.GetXaxis().SetTitleSize(0.042)
        hist_up.GetXaxis().SetTitleOffset(1.1)
        hist_up.GetYaxis().SetTitle('|syst up + syst down| / 2')
        hist_up.GetYaxis().SetLabelSize(0.033)
        hist_up.GetYaxis().SetTitleSize(0.042)
        hist_up.GetYaxis().SetTitleOffset(1.1)
        #hist_up.SetDirectory(ROOT.gROOT)
          
        leg.AddEntry(hist_up, 'mass {} GeV c#tau {} mm'.format(signal_mass, signal_ctau))

        for icat, category in enumerate(self.categories):
          if 'incl' in category.label:
            continue
          #print '\n ',category.label
          #if 'incl' in category.label:
          #  continue
          #signal_selection = 'ismatched==1' if self.baseline_selection=='' else 'ismatched==1 && {} && {}'.format(self.baseline_selection, category.definition_flat)

          ### for hlt weight
          ##signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=weight, weight_pusig='weight_pu_sig_D', weight_mu0id='weight_mu0_softid', weight_muid='weight_mu_looseid')[0] 
          ##signal_yields_weight_plus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=weight+'_plus_one_sigma', weight_pusig='weight_pu_sig_D', weight_mu0id='weight_mu0_softid', weight_muid='weight_mu_looseid')[0] 
          ##signal_yields_weight_minus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=weight+'_minus_one_sigma', weight_pusig='weight_pu_sig_D', weight_mu0id='weight_mu0_softid', weight_muid='weight_mu_looseid')[0] 

          ## for mu0 weight
          #signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight, weight_muid='weight_mu_looseid')[0] 
          #signal_yields_weight_plus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight+'_plus_one_sigma', weight_muid='weight_mu_looseid')[0] 
          #signal_yields_weight_minus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight+'_plus_one_sigma', weight_muid='weight_mu_looseid')[0] 

          #sig_label = 'm_{}_ctau_{}_cat_{}'.format(signal_mass, signal_ctau, category.label)
          #filename_sig = '{}.root'.format(sig_label.replace('.', 'p'))
          #file_sig = self.tools.getRootFile(filename_sig)
          #tree_sig = self.tools.getTree(file_sig, 'signal_tree')
          #print tree_sig.GetEntries()

          ## for weight_hlt
          #selection_string_weight = '(score > 0.95) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid)'
          #selection_string_weight_plus_one_sigma = '(score > 0.95) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_plus_one_sigma * weight_mu0_softid * weight_mu_looseid)'
          #selection_string_weight_minus_one_sigma = '(score > 0.95) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_minus_one_sigma * weight_mu0_softid * weight_mu_looseid)'
          #hist_name_weight = 'hist_signal_m_{}_ctau_{}_cat_{}_weight'.format(signal_mass, signal_ctau, category.label)
          #hist_weight = ROOT.TH1D(hist_name_weight, hist_name_weight, 100, signal_mass-10*signal_resolution, signal_mass+10*signal_resolution)
          #tree_sig.Project(hist_name_weight, 'hnl_mass', selection_string_weight)

          #signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight, weight_muid='weight_mu_looseid')[0] 

          ## for mu weight
          ##signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight, weight_muid='weight_mu0_softid')[0] 
          ##signal_yields_weight_plus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight+'_plus_one_sigma', weight_muid='weight_mu0_softid')[0] 
          ##signal_yields_weight_minus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight+'_plus_one_sigma', weight_muid='weight_mu0_softid')[0] 

          #syst_up = (signal_yields_weight_plus_one_sigma / signal_yields_weight -1 ) * 100
          #syst_down = (signal_yields_weight_minus_one_sigma / signal_yields_weight -1) * 100
          ##if syst_up < syst_up_min: syst_up_min = syst_up
          ##if syst_up > syst_up_max: syst_up_max = syst_up
          ##if syst_down < syst_down_min: syst_down_min = syst_down
          ##if syst_down > syst_down_max: syst_down_max = syst_down

          #syst = (abs(syst_up) + abs(syst_down)) / 2.

          #if syst < syst_min: syst_min = syst
          #if syst > syst_max: syst_max = syst

          #hist_up.SetBinContent(icat+1, syst)
          #hist_up.SetDirectory(0)
          #hist_up.GetXaxis().SetBinLabel(icat+1, category.title)
          #hists_up.append(hist_up)
          #hist_up.Draw()

          #print 'm{}\_ctau{} & {}\% & {}\% \\\ '.format(signal_mass, signal_ctau, syst_up, syst_down)
          ##else:

          signal_selection = 'ismatched==1' if self.baseline_selection=='' else 'ismatched==1 && {} && {}'.format(self.baseline_selection, category.definition_flat)

          training_label = './mva/outputs/test_2022Nov29_09h26m28s'

          sig_label = 'm_{}_ctau_{}_cat_{}'.format(signal_mass, signal_ctau, category.label)
          #print 'creating file with score'
          extra_branches = [
            #'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable',
            #'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_plus_one_sigma',
            #'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_minus_one_sigma',
            'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_v2',
            'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_plus_one_sigma_v2',
            'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_minus_one_sigma_v2',
            'weight_mu0_softid',
            'weight_mu0_softid_plus_one_sigma',
            'weight_mu0_softid_minus_one_sigma',
            'weight_mu_looseid',
            'weight_mu_looseid_plus_one_sigma',
            'weight_mu_looseid_minus_one_sigma',
            'weight_pu_sig_D',
            ]
          filename_sig = self.mva_tools.getFileWithScore(files=[signal_file], training_label=training_label, category_label=category.label, do_parametric=True, mass=signal_mass, selection=signal_selection, weights=extra_branches, label=sig_label, treename='signal_tree', force_overwrite=False) 
          filename_sig = '{}.root'.format(sig_label.replace('.', 'p'))
          file_sig = self.tools.getRootFile(filename_sig)
          tree_sig = self.tools.getTree(file_sig, 'signal_tree')
          print tree_sig.GetEntries()

          # for weight_hlt
          selection_string_weight = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_v2 * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D)'
          selection_string_weight_plus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_plus_one_sigma_v2 * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D)'
          selection_string_weight_minus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_minus_one_sigma_v2 * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D)'

          # initial (everything below 6 put in first binning)
          #selection_string_weight = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D)'
          #selection_string_weight_plus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_plus_one_sigma * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D)'
          #selection_string_weight_minus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable_minus_one_sigma * weight_mu0_softid * weight_mu_looseidi * weight_pu_sig_D)'

          # for weight_mu0_softid
          #selection_string_weight = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D )'
          #selection_string_weight_plus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid_plus_one_sigma * weight_mu_looseid * weight_pu_sig_D)'
          #selection_string_weight_minus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid_minus_one_sigma * weight_mu_looseid * weight_pu_sig_D)'

          # for weight_mu_looseid
          #selection_string_weight = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D )'
          #selection_string_weight_plus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid_plus_one_sigma * weight_pu_sig_D)'
          #selection_string_weight_minus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid_minus_one_sigma * weight_pu_sig_D)'

          # for pileup
          #selection_string_weight = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D )'
          #selection_string_weight_plus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D*1.05)'
          #selection_string_weight_minus_one_sigma = '(score > 0.99) * (weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable * weight_mu0_softid * weight_mu_looseid * weight_pu_sig_D*0.95)'

          hist_name_weight = 'hist_signal_m_{}_ctau_{}_cat_{}_weight'.format(signal_mass, signal_ctau, category.label)
          hist_weight = ROOT.TH1D(hist_name_weight, hist_name_weight, 100, signal_mass-10*signal_resolution, signal_mass+10*signal_resolution)
          tree_sig.Project(hist_name_weight, 'hnl_mass', selection_string_weight)

          hist_name_weight_plus_one_sigma = 'hist_signal_m_{}_ctau_{}_cat_{}_weight_plus_one_sigma'.format(signal_mass, signal_ctau, category.label)
          hist_weight_plus_one_sigma = ROOT.TH1D(hist_name_weight_plus_one_sigma, hist_name_weight_plus_one_sigma, 100, signal_mass-10*signal_resolution, signal_mass+10*signal_resolution)
          tree_sig.Project(hist_name_weight_plus_one_sigma, 'hnl_mass', selection_string_weight_plus_one_sigma)

          hist_name_weight_minus_one_sigma = 'hist_signal_m_{}_ctau_{}_cat_{}_weight_minus_one_sigma'.format(signal_mass, signal_ctau, category.label)
          hist_weight_minus_one_sigma = ROOT.TH1D(hist_name_weight_minus_one_sigma, hist_name_weight_minus_one_sigma, 100, signal_mass-10*signal_resolution, signal_mass+10*signal_resolution)
          tree_sig.Project(hist_name_weight_minus_one_sigma, 'hnl_mass', selection_string_weight_minus_one_sigma)

          signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, is_bc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable', weight_pusig='weight_pu_sig_D', weight_mu0id=weight, weight_muid='weight_mu_looseid')[0] 

          # get the number of yields
          n_sig_weight = hist_weight.Integral()
          n_sig_weight_plus_one_sigma = hist_weight_plus_one_sigma.Integral()
          n_sig_weight_minus_one_sigma = hist_weight_minus_one_sigma.Integral()

          syst_up = (n_sig_weight_plus_one_sigma / n_sig_weight -1 ) * 100
          syst_down = (n_sig_weight_minus_one_sigma / n_sig_weight -1 ) * 100

          syst = (abs(syst_up) + abs(syst_down)) / 2.

          if syst < syst_min: syst_min = syst
          if syst > syst_max: syst_max = syst

          canv.cd()

          #hist_up = ROOT.TH1D('hist_up', 'hist_up', len(self.categories), 0, len(self.categories)) 
          #print 'creating hist' 
          #hist_up.SetMarkerSize(2)
          #hist_up.SetMarkerStyle(20)
          #hist_up.SetMarkerColor(signal_file.colour)
          #hist_up.SetDirectory(0)
          ##hist_up = hist_up.Clone('hist_up')

          hist_up.SetBinContent(icat, syst)
          hist_up.GetXaxis().SetBinLabel(icat, category.title)
          hists_up.append(hist_up)
          #hist_up.Draw()

          print 'm{}\_ctau{} & {}\% & {}\% \\\ '.format(signal_mass, signal_ctau, syst_up, syst_down)

    #canv = ROOT.TCanvas('canv', 'canv', 800, 700)
    ROOT.gStyle.SetOptStat(0)

    for ihist, hist_up in enumerate(hists_up):
      if ihist == 0:
        hist_up.GetYaxis().SetRangeUser(syst_min-0.15*syst_min, syst_max+0.35*syst_max)
        hist_up.SetTitle('')
        hist_up.Draw('P')
      else:
        hist_up.Draw('P same')

    leg.Draw()

    canv.SaveAs('test.png')
    canv.SaveAs('test.pdf')


  def getSelection(self, selection):
    '''
      Function that removes the cut on the score if any from the selection string
    '''
    if 'score' not in selection:
      selection_string = selection
    else:
      idx = selection.find('score')
      idx_in = selection.rfind('&&')
      if idx_in < idx:
        selection_string = selection[:idx_in-1]
      else:
        selection_string = selection[:idx-1] + selection[idx_in+2:]
      if 'score' in selection_string:
        idx = selection_string.find('score')
        idx_in = selection_string.rfind('&&')
        if idx_in < idx:
          selection_string = selection_string[:idx_in-1]
        else:
          selection_string = selection_string[:idx-1] + selection_string[idx_in+2:]
    
    return selection_string


  def convertRootToDF(self, sample, training_info, signal_treename, selection, weights):
    '''
      Function that converts root samples to pandas dataframe
      The only kept branches are the hnl mass, score and signal weights
    '''
    if weights == None:
      extra_columns = []
    else:
      extra_columns = weights

    df = read_root(sample, signal_treename, where=self.getSelection(selection), warn_missing_tree=True, columns=training_info.features+extra_columns)

    return df

      
  def createDataframe(self, training_info, do_parametric, mass, samples_filename, selection, weights, signal_treename):
    '''
      Function that returns a dataframe out of a list of samples
    '''
    df = pd.concat([self.convertRootToDF(idt, training_info, signal_treename, selection, weights) for idt in samples_filename], sort=False)

    # remove inf and nan
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

    # re-index
    df = df.reset_index(drop=True)
    
    # add parameters
    if do_parametric:
      df['mass_key'] = mass

    return df


  def studySelectionSyst(self, categories, baseline_selection):
    ROOT.gStyle.SetOptStat(0)
    #ROOT.TH1.SetDefaultSumw2()

    #training_label = 'V13_06Feb23_2023Apr06_14h13m31s'
    training_label = 'V13_06Feb23_2023Jul18_15h36m17s'

    #filename_mc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/Chunk0_n500/flat/flat_bparknano_06Feb23_nj1.root'
    #filename_mc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/Chunk0_n500/flat/flat_bparknano_06Feb23_pNN_nj1.root'
    filename_mc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_06Feb23_pNN.root'
    #filename_data = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_control_pNN.root'
    filename_data = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_control_splot.root'

    mass = 4.5
    resolution_p0 = 6.98338e-04
    resolution_p1 = 7.78382e-03 
    resolution = resolution_p0 + mass * resolution_p1
    window_size = 1000
    treename = 'control_tree'

    mva_tools = MVATools(path_mva='./mva/outputs')

    outputdir = './myPlots/systematics/selection'
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
        
    hist_efficiency_mc = ROOT.TH1D('hist_efficiency_mc', 'hist_efficiency_mc', len(categories)-1, 0, len(categories)-1) 
    hist_efficiency_mc.SetTitle(' ')
    hist_efficiency_mc.SetFillColor(ROOT.kBlue-3)
    hist_efficiency_mc.SetFillStyle(3005)
    hist_efficiency_mc.GetXaxis().SetLabelSize(0.0)
    hist_efficiency_mc.GetXaxis().SetTitleSize(0.0)
    hist_efficiency_mc.GetYaxis().SetTitle('Efficiency')
    hist_efficiency_mc.GetYaxis().SetLabelSize(0.037)
    hist_efficiency_mc.GetYaxis().SetTitleSize(0.043)
    hist_efficiency_mc.GetYaxis().SetTitleOffset(1.1)

    hist_efficiency_data = ROOT.TH1D('hist_efficiency_data', 'hist_efficiency_data', len(categories)-1, 0, len(categories)-1) 
    hist_efficiency_data.SetMarkerStyle(20)
    hist_efficiency_data.SetMarkerSize(2)

    hist_syst = ROOT.TH1D('hist_syst', 'hist_syst', len(categories)-1, 0, len(categories)-1) 
    hist_syst.SetTitle(' ')
    #hist_syst.SetMarkerSize(2)
    #hist_syst.SetMarkerStyle(20)
    hist_syst.SetMarkerColor(2)
    hist_syst.SetLineWidth(2)
    hist_syst.SetMarkerStyle(33)
    hist_syst.SetMarkerSize(2)
    hist_syst.SetTitle('')
    #hist_syst.GetXaxis().SetTitle('pNN score')
    hist_syst.GetXaxis().SetLabelSize(0.16)
    hist_syst.GetXaxis().SetTitleSize(0.17)
    hist_syst.GetXaxis().SetTitleOffset(0.73)
    hist_syst.GetYaxis().SetTitle('Data/MC')
    hist_syst.GetYaxis().SetLabelSize(0.1)
    hist_syst.GetYaxis().SetTitleSize(0.13)
    hist_syst.GetYaxis().SetTitleOffset(0.345)
    hist_syst.GetYaxis().SetRangeUser(0.5, 1.5)
    hist_syst.GetYaxis().SetNdivisions(5)

    leg_syst = self.tools.getRootTLegend(xmin=0.3, ymin=0.55, xmax=0.5, ymax=0.85, size=0.045)

    for icat, category in enumerate(categories): 
      if 'incl' in category.label: continue
      if 'Bc' in category.label: continue
      #if 'SS' in category.label: continue
      #if category.label != 'lxysig0to50_OS': continue
      #if category.label != 'lxysiggt150_SS': continue

      training_info = mva_tools.getTrainingInfo(training_label, category.label)

      selection = baseline_selection + ' && mu_pi_mass > {} && mu_pi_mass < {}'.format(mass-window_size*resolution, mass+window_size*resolution)
      #selection += ' && sv_lxysig<50'

      weights_mc = ['weight_hlt_A1']
      df_mc = self.createDataframe(training_info=training_info, do_parametric=True, mass=mass, samples_filename=[filename_mc], selection=selection, weights=weights_mc, signal_treename=treename)
      weight_hlt = df_mc['weight_hlt_A1']

      weights_data = ['nsig_sw']
      df_data = self.createDataframe(training_info=training_info, do_parametric=True, mass=mass, samples_filename=[filename_data], selection=selection, weights=weights_data, signal_treename=treename)
      sweight = df_data['nsig_sw']

      score_mc = mva_tools.predictScore(training_info=training_info, df=df_mc, do_parametric=True) 
      score_data = mva_tools.predictScore(training_info=training_info, df=df_data, do_parametric=True) 

      canv_name = 'canv_score_{}'.format(category.label)
      canv = self.tools.createTCanvas(canv_name)
      pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
      #pad_up.SetTopMargin(0.13)
      pad_up.SetBottomMargin(0.03)
      pad_up.Draw()
      canv.cd()
      pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
      pad_down.SetBottomMargin(0.25)
      pad_down.Draw()
      leg = self.tools.getRootTLegend(xmin=0.3, ymin=0.55, xmax=0.5, ymax=0.85, size=0.045)

      nbins = 30
      bin_min = 0.05
      bin_max = 1

      score_WP = 0.99

      hist_score_mc_full = ROOT.TH1D('hist_score_mc_full_{}'.format(category.label), 'hist_score_mc_full', nbins, 0, 1)
      hist_score_mc = ROOT.TH1D('hist_score_mc_{}'.format(category.label), 'hist_score_mc', nbins, bin_min, bin_max)
      hist_score_mc_WP = ROOT.TH1D('hist_score_mc_WP_{}'.format(category.label), 'hist_score_mc_WP', nbins, score_WP, 1)
      for idx, i in enumerate(score_mc):
        #hist_score_mc_full.Fill(i[0])
        #hist_score_mc.Fill(i[0])
        #hist_score_mc_WP.Fill(i[0])
        hist_score_mc_full.Fill(i[0], weight_hlt[idx])
        hist_score_mc.Fill(i[0], weight_hlt[idx])
        hist_score_mc_WP.Fill(i[0], weight_hlt[idx])

      hist_score_data_full = ROOT.TH1D('hist_score_data_full_{}'.format(category.label), 'hist_score_data_full', nbins, 0, 1)
      hist_score_data = ROOT.TH1D('hist_score_data_{}'.format(category.label), 'hist_score_data', nbins, bin_min, bin_max)
      hist_score_data_WP = ROOT.TH1D('hist_score_data_WP_{}'.format(category.label), 'hist_score_data_WP', nbins, score_WP, 1)
      for idx, i in enumerate(score_data):
        #hist_score_data_full.Fill(i[0])
        #hist_score_data.Fill(i[0])
        #hist_score_data_WP.Fill(i[0])
        hist_score_data_full.Fill(i[0], sweight[idx])
        hist_score_data.Fill(i[0], sweight[idx])
        hist_score_data_WP.Fill(i[0], sweight[idx])

      print hist_score_mc_full.Integral()
      print hist_score_mc_WP.Integral()
  
      err_mc_full = ROOT.double(0.)
      int_mc_full = hist_score_mc_full.IntegralAndError(0, 10000, err_mc_full)

      err_data_full = ROOT.double(0.)
      int_data_full = hist_score_data_full.IntegralAndError(0, 100000, err_data_full)

      hist_score_mc.Scale(int_data_full/int_mc_full)

      err_mc_WP = ROOT.double(0.)
      int_mc_WP = hist_score_mc_WP.IntegralAndError(hist_score_mc_WP.FindBin(score_WP), 100000, err_mc_WP)

      err_data_WP = ROOT.double(0.)
      int_data_WP = hist_score_data_WP.IntegralAndError(hist_score_data_WP.FindBin(score_WP), 100000, err_data_WP)

      efficiency_mc = int_mc_WP / int_mc_full
      efficiency_data = int_data_WP / int_data_full

      try:
        err_mc = efficiency_mc * math.sqrt(math.pow(err_mc_WP/int_mc_WP, 2) + math.pow(err_mc_full/int_mc_full, 2))
      except:
        err_mc = 0
      try:
        err_data = efficiency_data * math.sqrt(math.pow(err_data_WP/int_data_WP, 2) + math.pow(err_data_full/int_data_full, 2))
      except:
        err_data = 0

      hist_efficiency_mc.SetBinContent(icat, efficiency_mc)
      hist_efficiency_mc.SetBinError(icat, err_mc)
      hist_efficiency_data.SetBinContent(icat, efficiency_data)
      hist_efficiency_data.SetBinError(icat, err_data)

      #syst = int_data_WP / int_mc_WP
      try:
        syst = efficiency_data / efficiency_mc
      except:
        syst = 0
      try: 
        err_syst = syst * math.sqrt(math.pow(err_mc/efficiency_mc, 2) + math.pow(err_data/efficiency_data, 2)) 
      except:
        err_syst = 0
      print '{} syst {} / {} = {}'.format(category.label, efficiency_data, efficiency_mc, syst)
      hist_syst.SetBinContent(icat, syst)
      hist_syst.SetBinError(icat, err_syst)
      hist_syst.GetXaxis().SetBinLabel(icat, category.title)

      #hist_score_mc.Scale(1./hist_score_mc.Integral())
      #hist_score_data.Scale(1./hist_score_data.Integral())

      hist_score_mc.SetTitle(category.title)
      hist_score_mc.SetFillColor(ROOT.kBlue-3)
      hist_score_mc.SetFillStyle(3005)
      hist_score_mc.GetXaxis().SetLabelSize(0.0)
      hist_score_mc.GetXaxis().SetTitleSize(0.0)
      hist_score_mc.GetYaxis().SetTitle('Entries')
      hist_score_mc.GetYaxis().SetLabelSize(0.037)
      hist_score_mc.GetYaxis().SetTitleSize(0.043)
      hist_score_mc.GetYaxis().SetTitleOffset(1.1)
      range_min = min(hist_score_mc.GetMinimum(), hist_score_data.GetMinimum())
      range_max = max(hist_score_mc.GetMaximum(), hist_score_data.GetMaximum())
      hist_score_mc.GetYaxis().SetRangeUser(range_min-0.15*range_min, range_max+0.15*range_max)

      hist_score_data.SetMarkerStyle(20)
      hist_score_data.SetMarkerSize(2)

      hist_score_mc_err = hist_score_mc.Clone('hist_score_mc_err')
      hist_score_mc_err.SetLineWidth(0)
      hist_score_mc_err.SetFillStyle(3244)
      hist_score_mc_err.SetFillColor(ROOT.kGray+2)

      leg.AddEntry(hist_score_data, 'Data')
      leg.AddEntry(hist_score_mc, 'MC')

      pad_up.cd()
      hist_score_mc.Draw('hist')
      hist_score_data.Draw('PE same')
      hist_score_mc_err.Draw('E2 same')

      leg.Draw()

      pad_down.cd()
      hist_ratio = self.tools.getRatioHistogram(hist_score_data, hist_score_mc)
      hist_ratio.SetLineWidth(2)
      hist_ratio.SetMarkerStyle(20)
      hist_ratio.SetTitle('')
      hist_ratio.GetXaxis().SetTitle('pNN score')

      hist_ratio.GetXaxis().SetLabelSize(0.1)
      hist_ratio.GetXaxis().SetTitleSize(0.13)
      hist_ratio.GetXaxis().SetTitleOffset(0.73)
      hist_ratio.GetYaxis().SetTitle('Data/MC')
      hist_ratio.GetYaxis().SetLabelSize(0.1)
      hist_ratio.GetYaxis().SetTitleSize(0.13)
      hist_ratio.GetYaxis().SetTitleOffset(0.345)
      val_min = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
      val_max = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
      hist_ratio.GetYaxis().SetRangeUser(val_min-0.15*val_min, val_max+0.15*val_max)

      hist_ratio.Draw('P')

      ## draw line at ratio = 1
      line = ROOT.TLine(bin_min, 1, bin_max, 1)
      line.SetLineColor(4)
      line.SetLineWidth(2)
      line.Draw('same')

      canv.cd()
      name = 'score_{}'.format(category.label)
      canv.SaveAs('{}/{}.png'.format(outputdir, name))
      canv.SaveAs('{}/{}.pdf'.format(outputdir, name))

    canv_name_syst = 'canv_syst'
    canv_syst = self.tools.createTCanvas(canv_name_syst)
    pad_up_syst = ROOT.TPad("pad_up_syst","pad_up_syst",0,0.25,1,1)
    pad_up_syst.SetBottomMargin(0.03)
    pad_up_syst.Draw()
    canv_syst.cd()
    pad_down_syst = ROOT.TPad("pad_down_syst","pad_down_syst",0,0,1,0.25)
    pad_down_syst.SetBottomMargin(0.35)
    pad_down_syst.Draw()

    pad_up_syst.cd()
    range_max = max(hist_efficiency_mc.GetMaximum(), hist_efficiency_data.GetMaximum())
    #range_min = min(hist_efficiency_mc.GetMinimum(), hist_efficiency_data.GetMinimum())
    print '{} {} {}'.format(hist_efficiency_mc.GetMaximum(), hist_efficiency_data.GetMaximum(), range_max)
    #print range_min
    hist_efficiency_mc.GetYaxis().SetRangeUser(0, range_max + 0.15*range_max)
    hist_efficiency_mc.Draw('hist')
    hist_efficiency_data.Draw('PE same')

    leg_syst.AddEntry(hist_efficiency_data, 'Data')
    leg_syst.AddEntry(hist_efficiency_mc, 'MC')
    leg.Draw()

    hist_efficiency_mc_err = hist_efficiency_mc.Clone('hist_efficiency_mc_err')
    hist_efficiency_mc_err.SetLineWidth(0)
    hist_efficiency_mc_err.SetFillStyle(3244)
    hist_efficiency_mc_err.SetFillColor(ROOT.kGray+2)
    hist_efficiency_mc_err.Draw('E2 same')

    pad_down_syst.cd() 
    hist_syst.Draw('PE')

    ## draw line at ratio = 1
    line = ROOT.TLine(0, 1, len(categories)-1, 1)
    line.SetLineColor(1)
    line.SetLineWidth(2)
    line.Draw('same')
    line_up = ROOT.TLine(0, 1.3, len(categories)-1, 1.3)
    line_up.SetLineColor(1)
    line_up.SetLineWidth(1)
    line_up.Draw('same')
    line_down = ROOT.TLine(0, 0.7, len(categories)-1, 0.7)
    line_down.SetLineColor(1)
    line_down.SetLineWidth(1)
    line_down.Draw('same')

    canv_syst.cd()
    name_syst = 'systematics'
    canv_syst.SaveAs('{}/{}.png'.format(outputdir, name_syst))
    canv_syst.SaveAs('{}/{}.pdf'.format(outputdir, name_syst))
    

      #bins = np.arange(0, 1.1, 0.1)
      #plt.hist(score_mc, bins=bins)
      #plt.show()
      

      #filename_mc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/Chunk0_n500/flat/flat_bparknano_06Feb23_nj1.root'
      #f_mc_original = self.tools.getRootFile(filename=filename_mc)
      #tree_mc_original = self.tools.getTree(f_mc_original, 'control_tree')

      ## to apply pNN on control sample, the branches are to be renamed
      ### create dir
      #ourdir_file = '/scratch/anlyon/systematics'
      #if not path.exists(ourdir_file):
      #  os.system('mkdir -p {}'.format(ourdir_file))

      ### create file
      #outfile_name_mc = '{}/file_mc.root'.format(ourdir_file)
      #outfile_mc = ROOT.TFile(outfile_name_mc, 'RECREATE')
      ##['pi_pt', 'mu_pt', 'mu0_pt', 'b_mass', 'hnl_cos2d', 'sv_lxysig', 'sv_prob', 'sv_chi2', 'b_pt', 'mu0_mu_mass', 'mu0_pi_mass', 'deltar_mu0_mu', 'deltar_mu0_pi', 'mu0_pfiso03_rel', 'mu_pfiso03_rel']

      #tree_mc_renamed = tree_mc_original.CloneTree()
      #tree_mc_renamed.GetBranch('k_pt').SetTitle('pi_pt')
      #tree_mc_renamed.GetBranch('k_pt').SetName('pi_pt')
      #tree_mc_renamed.GetBranch('l1_pt').SetTitle('mu0_pt')
      #tree_mc_renamed.GetBranch('l1_pt').SetName('mu0_pt')
      #tree_mc_renamed.GetBranch('l2_pt').SetTitle('mu_pt')
      #tree_mc_renamed.GetBranch('l2_pt').SetName('mu_pt')
      #tree_mc_renamed.GetBranch('b_cos2d').SetTitle('hnl_cos2d')
      #tree_mc_renamed.GetBranch('b_cos2d').SetName('hnl_cos2d')
      #tree_mc_renamed.GetBranch('dimu_sv_prob').SetTitle('sv_prob') #FIXME sv_chi2 not in tree
      #tree_mc_renamed.GetBranch('dimu_sv_prob').SetName('sv_prob')
      #tree_mc_renamed.GetBranch('k_phi').SetTitle('sv_chi2') #FIXME sv_chi2 not in tree
      #tree_mc_renamed.GetBranch('k_phi').SetName('sv_chi2')
      #tree_mc_renamed.GetBranch('dimu_mass').SetTitle('mu0_mu_mass')
      #tree_mc_renamed.GetBranch('dimu_mass').SetName('mu0_mu_mass')
      #tree_mc_renamed.GetBranch('l1_eta').SetTitle('mu0_pi_mass') #FIXME mu0_pi_mass not in tree
      #tree_mc_renamed.GetBranch('l1_eta').SetName('mu0_pi_mass')
      #tree_mc_renamed.GetBranch('l1_phi').SetTitle('deltar_mu0_mu') #FIXME
      #tree_mc_renamed.GetBranch('l1_phi').SetName('deltar_mu0_mu')
      #tree_mc_renamed.GetBranch('l2_eta').SetTitle('deltar_mu0_pi') #FIXME
      #tree_mc_renamed.GetBranch('l2_eta').SetName('deltar_mu0_pi')
      #tree_mc_renamed.GetBranch('b_eta').SetTitle('mu0_pfiso03_rel') #FIXME
      #tree_mc_renamed.GetBranch('b_eta').SetName('mu0_pfiso03_rel')
      #tree_mc_renamed.GetBranch('k_eta').SetTitle('mu_pfiso03_rel') #FIXME
      #tree_mc_renamed.GetBranch('k_eta').SetName('mu_pfiso03_rel')

      #outfile_mc.cd()
      #tree_mc_renamed.Print()
      #outfile_mc.Write()







if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  do_studyHLTsyst = False
  do_plotOccupancy = False
  do_plotTriggeringTable = False
  do_studyShapeSyst = False
  do_plotMCSyst = False
  do_studySelectionSyst = True

  signal_labels = []
  signal_label_m1 = 'V12_08Aug22_m1'
  signal_label_m1p5 = 'V12_08Aug22_m1p5'
  signal_label_m2 = 'V12_08Aug22_m2'
  signal_label_m3 = 'V12_08Aug22_m3'
  signal_label_m4p5 = 'V12_08Aug22_m4p5'
  signal_labels.append(signal_label_m1)
  signal_labels.append(signal_label_m1p5)
  signal_labels.append(signal_label_m2)
  signal_labels.append(signal_label_m3)
  signal_labels.append(signal_label_m4p5)

  #signal_labels = ['V12_08Aug22_benchmark']

  #categories = categories['inclusive']
  #categories = categories['categories_0_50_150']
  #baseline_selection = selection['baseline_08Aug22'].flat

  if do_studyHLTsyst:
    #weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable' 
    weight_hlt = 'weight_mu_looseid' 
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).studyHLTsyst(weight_hlt=weight_hlt)

  if do_plotOccupancy:
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).plotOccupancy(muon_type='primary')
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).plotOccupancy(muon_type='displaced')

  if do_plotTriggeringTable:
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).plotTriggeringTable()

  if do_studyShapeSyst:
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).studyShapeSyst()
    #Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).studyShapeSystPerCategory()

  if do_plotMCSyst:
    #weight = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable'
    weight = 'weight_mu0_softid'
    #weight = 'weight_mu_looseid'
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).plotMCSyst(weight=weight)

  if do_studySelectionSyst:
    categories = categories['categories_0_50_150']
    baseline_selection = selection['control'].flat
    Systematics(signal_labels=signal_labels, categories=categories, baseline_selection=baseline_selection, lumi=41.6).studySelectionSyst(
        categories=categories, 
        baseline_selection=baseline_selection,
        )



