import os
import sys
import numpy as np
from os import path
import math
from tools import Tools
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples
from categories import categories
from baseline_selection import selection
from quantity import Quantity

import ROOT


class Systematics(Tools):
  def __init__(self, signal_labels='', categories='', selection='', lumi=41.6):
    self.tools = Tools()
    self.signal_labels = signal_labels
    self.categories = categories
    self.selection = selection
    self.lumi = lumi

  def studyHLTsyst(self, weight_hlt):
    print '\n',weight_hlt
    for signal_label in self.signal_labels:
      signal_files = signal_samples[signal_label]

      for signal_file in signal_files:
        signal_mass = signal_file.mass
        signal_ctau = signal_file.ctau

        signal_selection = 'ismatched==1' if self.selection=='' else 'ismatched==1 && {}'.format(self.selection)
        signal_yields_no_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, isBc=False, add_weight_hlt=False) 
        signal_yields_weight = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, weight_hlt=weight_hlt) 
        signal_yields_weight_plus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, weight_hlt=weight_hlt+'_plus_one_sigma') 
        signal_yields_weight_minus_one_sigma = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=self.lumi, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, weight_hlt=weight_hlt+'_minus_one_sigma') 

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
          #selection_incl = 'ismatched == 1 && mu0_istriggering==1)
          selection_incl = 'ismatched == 1 && trgmu_istriggering==1'
        elif muon_type == 'displaced':
          selection_incl = 'ismatched == 1 && mu_istriggering==1'
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
              #selection = 'ismatched == 1 && trgmu_istriggering==1 && mu0_pt > {} && mu0_pt < {} && mu0_dxysig > {} && mu0_dxysig < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
              selection = 'ismatched == 1 && trgmu_istriggering==1 && trgmu_pt > {} && trgmu_pt < {} && abs(trgmu_dxysig) > {} && abs(trgmu_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
            elif muon_type == 'displaced':
              selection = 'ismatched == 1 && mu_istriggering==1 && mu_pt > {} && mu_pt < {} && abs(mu_dxysig) > {} && abs(mu_dxysig) < {}'.format(pt_min, pt_max, dxysig_min, dxysig_max)
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
        f = self.tools.getRootFile(signal_file.filename)
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
        selection = 'ismatched == 1 && trgmu_istriggering==1 && mu_istriggering==1'
        #selection = 'ismatched == 1 && mu0_istriggering==1 && mu_istriggering==1'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        tot = tot + entries
        hist2d.Fill(mu0_bins[0], mu_bins[0], entries/entries_incl)

        hist_name = 'hist01'
        selection = 'ismatched == 1 && trgmu_istriggering==1 && mu_istriggering==0'
        #selection = 'ismatched == 1 && mu0_istriggering==1 && mu_istriggering==0'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        tot = tot + entries
        hist2d.Fill(mu0_bins[0], mu_bins[1], entries/entries_incl)

        hist_name = 'hist10'
        selection = 'ismatched == 1 && trgmu_istriggering==0 && mu_istriggering==1'
        #selection = 'ismatched == 1 && mu0_istriggering==0 && mu_istriggering==1'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
        tot = tot + entries
        hist2d.Fill(mu0_bins[1], mu_bins[0], entries/entries_incl)

        hist_name = 'hist11'
        selection = 'ismatched == 1 && trgmu_istriggering==0 && mu_istriggering==0'
        #selection = 'ismatched == 1 && mu0_istriggering==0 && mu_istriggering==0'
        weight = '(1)'
        hist = self.tools.createHisto(tree, quantity, hist_name, weight, selection) 
        entries = hist.GetEntries()
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


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  do_studyHLTsyst = False
  do_plotOccupancy = True
  do_plotTriggeringTable = True

  signal_labels = []
  signal_label_m1 = 'central_V11_24Apr22_m1'
  signal_label_m1p5 = 'central_V11_24Apr22_m1p5'
  signal_label_m2 = 'central_V11_24Apr22_m2'
  signal_label_m3 = 'central_V11_24Apr22_m3'
  signal_label_m4p5 = 'central_V11_24Apr22_m4p5'
  signal_labels.append(signal_label_m1)
  signal_labels.append(signal_label_m1p5)
  signal_labels.append(signal_label_m2)
  signal_labels.append(signal_label_m3)
  signal_labels.append(signal_label_m4p5)

  categories = categories['inclusive']
  selection = selection['baseline_30Dec21'].flat

  if do_studyHLTsyst:
    weight_hlt = 'weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2' 
    Systematics(signal_labels=signal_labels, categories=categories, selection=selection, lumi=41.6).studyHLTsyst(weight_hlt=weight_hlt)

  if do_plotOccupancy:
    Systematics(signal_labels=signal_labels, categories=categories, selection=selection, lumi=41.6).plotOccupancy(muon_type='primary')
    Systematics(signal_labels=signal_labels, categories=categories, selection=selection, lumi=41.6).plotOccupancy(muon_type='displaced')

  if do_plotTriggeringTable:
    Systematics(signal_labels=signal_labels, categories=categories, selection=selection, lumi=41.6).plotTriggeringTable()



