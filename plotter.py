import os
import os.path
from os import path
import ROOT
from ROOT import gROOT, gStyle
from copy import copy
import math
import numpy as np
from itertools import product

from tools import Tools
from quantity import quantities_to_plot_small, quantities_to_plot_all, quantities_to_plot_selection, quantities_muonId_study_triggermuon, quantities_muonId_study_displacedmuon, quantities_muonId_study_triggermuon_small, quantities_muonId_study_displacedmuon_small, quantities_tag_and_probe, quantities_trackId
from samples import data_samples, data_samples_V02, data_samples_V03, data_samples_loose, data_samples_triggermuon_matching_check, data_samples_tag_and_probe, data_samples_loose_dsaonly, qcd_samples, qcd_samples_V03, qcd_samples_triggermuon_matching_check, signal_samples, signal_samples_loose, signal_samples_tag_and_probe, signal_samples_tag_and_probe_BToJPsiKstar, signal_samples_loose_dsaonly
from quantity import Quantity
from computeYields import ComputeYields

# temporary
import sys
sys.path.append('../../../../BHNLGen/CMSSW_10_2_15/src/HNLsGen/python/.')
from my_common import getVV


class Plotter(Tools):
  def __init__(self, quantity='', data_files='', qcd_files='', signal_files='', white_list=''):
    self.tools = Tools()
    self.quantity = quantity
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_files = signal_files
    self.white_list = white_list


  def getQCDMCLabel(self, label1, label2):
    '''
      based on MC label of the format ptXXtoYY 
    '''
    return label1[0:label1.find('to')] + label2[label2.find('to'):len(label2)]


  def getDataLabel(self, label):
    idx = 0
    data_label = ''
    while idx < len(label):
      period = label[label.find('Parking', idx)+10:label.find('Parking', idx)+11]
      run = label[label.find('Run2018', idx)+7:label.find('Run2018', idx)+8]
      dataset = run+period
      if idx == 0: data_label = dataset # or idx>len(label)-20:
      else: data_label += '+{}'.format(dataset)
      idx = idx + 21
    data_label = data_label.replace('A1+A2+A3+A4+A5+A6', 'A*')
    data_label = data_label.replace('B1+B2+B3+B4+B5+B6', 'B*')
    data_label = data_label.replace('C1+C2+C3+C4+C5', 'C*')
    data_label = data_label.replace('D1+D2+D3+D4+D5', 'D*')
    return data_label


  def getMaxRangeY(self, hist1, hist2, do_log=False, use_sig=False):
    margin = 0.15 if do_log==False else 0.5
    if use_sig: 
      max_hist1 = [hist.GetMaximum() for hist in hist1]
      the_max_hist1 = max(max_hist1)
      max_range = the_max_hist1+margin*the_max_hist1 if the_max_hist1>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    else:
      max_range = hist1.GetMaximum()+margin*hist1.GetMaximum() if hist1.GetMaximum()>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    return max_range


  def getLumiWeight(self, selection): # move
    '''
      weight = lumi_data / lumi_mc = N_data * sigma_mc / (N_mc * sigma_data) estimated as N_data / N_mc
    '''

    quantity_forweight = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=13000)
    n_obs_data = 0.
    n_err_data = 0.
    for data_file in self.data_files:
      f_data = self.tools.getRootFile(data_file.filename, with_ext=False) #  ROOT.TFile.Open(self.data_file.filename, 'READ')
      hist_data = self.tools.createHisto(f_data, 'signal_tree', quantity_forweight, branchname='flat', selection=selection)
      hist_data.Sumw2()
      n_obs_data += hist_data.GetBinContent(1)
      n_err_data += math.sqrt(n_obs_data) #hist_data.GetBinError(1)

    hist_mc_tot = self.tools.createWeightedHistoQCDMC(qcd_files=self.qcd_files, white_list=self.white_list, quantity=quantity_forweight, selection=selection)
    n_obs_mc = hist_mc_tot.GetBinContent(1)
    n_err_mc = math.sqrt(n_obs_mc) #hist_mc_tot.GetBinError(1)

    weight = float(n_obs_data) / float(n_obs_mc) if float(n_obs_mc)!= 0. else 0.

    #if n_obs_data != 0 and n_obs_mc != 0:
    #  err = weight* (n_err_data / n_obs_data + n_err_mc / n_obs_mc)
    #else: 
    #  err = 0

    return weight


  def plotHistFromTree(self, do_shape):
    data_file = data_samples[0]
    qcd_files = qcd_samples 

    filename_data = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano.root'
    filename_data = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano.root'
    filename_mc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root'
    f_data = ROOT.TFile.Open(filename_data, 'READ')
    f_mc = ROOT.TFile.Open(filename_mc, 'READ')
    
    hist_data = self.tools.createHisto(f_data, 'signal_tree', 'b_mass')
    hist_mc = self.tools.createHisto(f_mc, 'signal_tree', 'b_mass', weight=0) # 609

    hist_mc.SetFillColor(ROOT.kAzure-4)
    hist_mc.SetLineColor(2)

    #canv = ROOT.TCanvas('canv', 'canv', 900, 800)
    canv = self.tools.createTCanvas(name='canv', dimx=900, dimy=800)

    if do_shape: 
      int_data = hist_data.Integral()
      if int_data != 0: hist_data.Scale(1/int_data)
      int_mc = hist_mc.Integral()
      if int_mc != 0: hist_mc.Scale(1/int_mc)

    hist_mc_stack = ROOT.THStack('hist_mc_stack', '')
    hist_mc_stack.Add(hist_mc)

    hist_data.Draw()
    hist_mc_stack.Draw('histo same')
    #hist_mc.Draw('sames')

    canv.SaveAs('tmp.png')


  def plotDataMCComparison(self, selection='', title='', outdirlabel='', branchname='flat', treename='signal_tree', plot_data=False, plot_sig=False, plot_ratio=False, do_shape=True, do_luminorm=False, do_stack=True, do_log=False, add_overflow=False):
    if plot_data and self.data_files == '':
      raise RuntimeError('Please specify on which data sample to run')
    if plot_sig and self.signal_files == '':
      raise RuntimeError('Please specify on which signal sample to run')
    if not plot_data:
      plot_ratio = False

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, outdirlabel.replace('/', '_'), do_log, do_shape)
    #canv = ROOT.TCanvas(canv_name, canv_name, 1200, 1000)
    canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
    ROOT.SetOwnership(canv, False)
    #canv.SetGrid()
    canv.cd()

    # define the pads
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1) if plot_ratio else ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
    if plot_ratio: pad_up.SetBottomMargin(0.03)
    if do_log: pad_up.SetLogy()
    pad_up.Draw()
    canv.cd()
    if plot_ratio:
      pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
      pad_down.SetBottomMargin(0.25)
      pad_down.Draw()

    # prepare the legend
    if do_stack:
      #legend = self.tools.getRootTLegend(xmin=0.53, ymin=0.45, xmax=0.9, ymax=0.83, size=0.027)
      legend = self.tools.getRootTLegend(xmin=0.49, ymin=0.45, xmax=0.86, ymax=0.83, size=0.027)
    else:
      legend = self.tools.getRootTLegend(xmin=0.53, ymin=0.65, xmax=0.9, ymax=0.83, size=0.027)

    # write the top labels
    label_top_left = ROOT.TPaveText(0.1,0.92,0.5,0.93, "brNDC")
    label_top_left.SetBorderSize(0)
    label_top_left.SetFillColor(ROOT.kWhite)
    label_top_left.SetTextSize(0.042)
    label_top_left.SetTextAlign(11)
    label_top_left.SetTextFont(42)
    label_top_left.AddText('#textbf{Label} top left')
    #label_top_left.Draw()

    label_top_right = ROOT.TPaveText(0.5,0.92 if not plot_ratio else 0.95,0.91,0.93 if not plot_ratio else 0.96, "brNDC")
    label_top_right.SetBorderSize(0)
    label_top_right.SetFillColor(ROOT.kWhite)
    label_top_right.SetTextSize(0.042)
    label_top_right.SetTextAlign(31)
    label_top_right.SetTextFont(42)
    label_top_right.AddText(title)
    label_top_right.Draw()

    # data
    pad_up.cd()

    if plot_data:
      hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_data_tot.Sumw2()

      int_data_tot = 0.
      overflow_data_tot = 0
      error_overflow_data_tot = 0.

      hist_data_stack = ROOT.THStack('hist_data_stack', '')
      data_label = ''

      for idata, data_file in enumerate(self.data_files):
        f_data = self.tools.getRootFile(data_file.filename, with_ext=False) #  ROOT.TFile.Open(self.data_file.filename, 'READ')
        hist_data_name = 'hist_data_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        hist_data = self.tools.createHisto(f_data, treename, self.quantity, hist_name=hist_data_name, branchname=branchname, selection=selection)
        hist_data.Sumw2()
        if do_shape: 
          int_data_tot += hist_data.Integral()
        if add_overflow:
          overflow_data_tot += (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
          error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
          hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
          hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
          hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
          hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)

        hist_data_tot.Add(hist_data)
        if idata == 0: data_label = data_file.label[:20]
        else: data_label += '_{}'.format(data_file.label[:20])

        if do_shape and hist_data.Integral() != 0: 
          hist_data.Scale(1./hist_data.Integral())
          hist_data_stack.Add(hist_data)

      if do_shape and int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

      legend.AddEntry(hist_data_tot, 'data - {}'.format(self.getDataLabel(data_label) if len(self.data_files)>1 else data_file.label))

      ## set the style
      #hist_data_tot.SetLineWidth(0)
      hist_data_tot.SetMarkerStyle(20)
      #hist_data.SetTitle(self.title)

      #if not plot_ratio: 
      #  hist_data.GetXaxis().SetTitle(quantity.title)
      #  hist_data.GetXaxis().SetLabelSize(0.037)
      #  hist_data.GetXaxis().SetTitleSize(0.042)
      #  hist_data.GetXaxis().SetTitleOffset(1.1)
      #else:
      #  hist_data.GetXaxis().SetLabelSize(0.0)
      #  hist_data.GetXaxis().SetTitleSize(0.0)
      #hist_data.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
      #hist_data.GetYaxis().SetLabelSize(0.037)
      #hist_data.GetYaxis().SetTitleSize(0.042)
      #hist_data.GetYaxis().SetTitleOffset(1.1)
      #hist_data.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_data, hist_mc_stack, do_log))


    # signal
    signal_hists = []
    if plot_sig:
      for signal_file in self.signal_files:
        #f_signal = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+signal_file.filename, 'READ')
        f_signal = ROOT.TFile.Open(signal_file.filename, 'READ')
        hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        hist_signal = self.tools.createHisto(f_signal, treename, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection='ismatched==1' if selection=='' else 'ismatched==1 &&'+selection)
        hist_signal.Sumw2()
        #print '{} : {} entries'.format(signal_file.filename, int(hist_signal.Integral()))
        if do_shape: 
          int_signal = hist_signal.Integral()
          if int_signal != 0: hist_signal.Scale(1/int_signal)
        if add_overflow:
          overflow_signal = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
          error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
          hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
          hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
          hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
          hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)

        legend.AddEntry(hist_signal, 'signal - {}'.format(signal_file.label))

        hist_signal.SetLineWidth(3)
        hist_signal.SetLineColor(signal_file.colour)
        signal_hists.append(hist_signal)
        #hist_signal.Draw('same')

    # then, mc
    mc_hists = []

    hist_mc_tot = ROOT.TH1D('hist_mc_tot', 'hist_mc_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
    hist_mc_tot.Sumw2()

    int_mc_tot = 0.

    #colours = getColour()
    for ifile, qcd_file in enumerate(self.qcd_files):
      if qcd_file.label not in self.white_list: continue

      #f_mc = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+qcd_file.filename, 'READ')
      f_mc = ROOT.TFile.Open(qcd_file.filename, 'READ')
      
      #print qcd_file.label
      #weight_mc = self.tools.computeQCDMCWeight(self.tools.getTree(self, f_mc, 'signal_tree'), qcd_file.cross_section, qcd_file.filter_efficiency)
      weight_mc = self.tools.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
      #weight = '({}) * (weight_hlt)'.format(weight_mc)
      weight = '({})'.format(weight_mc)
      hist_mc_name = 'hist_mc_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
      hist_mc = self.tools.createHisto(f_mc, treename, self.quantity, hist_name=hist_mc_name, branchname=branchname, selection=selection, weight=weight) 
      hist_mc.Sumw2()

      #if qcd_file.label == 'V02_15to20':
      #  hist_mc = self.tools.createHisto(f_mc, 'signal_tree', 'b_mass', branchname, weight=1) # 609
      #else:
      #  hist_mc = self.tools.createHisto(f_mc, 'signal_tree', 'b_mass', branchname, weight=0.001) # 609

      hist_mc.SetFillColor(qcd_file.colour) #ROOT.kAzure-ifile)
      hist_mc.SetLineColor(1)
      
      if do_stack:
        legend.AddEntry(hist_mc, 'MC - {}'.format(qcd_file.label))

      if do_shape: 
        int_mc_tot += hist_mc.Integral()

      if add_overflow:
        overflow_mc = hist_mc.GetBinContent(hist_mc.GetNbinsX()) + hist_mc.GetBinContent(hist_mc.GetNbinsX()+1)
        error_overflow_mc = math.sqrt(math.pow(hist_mc.GetBinError(hist_mc.GetNbinsX()), 2) + math.pow(hist_mc.GetBinError(hist_mc.GetNbinsX()+1), 2)) 
        hist_mc.SetBinContent(hist_mc.GetNbinsX(), overflow_mc)
        hist_mc.SetBinError(hist_mc.GetNbinsX(), error_overflow_mc)
        hist_mc.SetBinContent(hist_mc.GetNbinsX()+1, 0)
        hist_mc.SetBinError(hist_mc.GetNbinsX()+1, 0)
  
      hist_mc_tot.Add(hist_mc)
      mc_hists.append(hist_mc)
  
    hist_mc_tot.SetFillColor(ROOT.kAzure-4)
    hist_mc_tot.SetLineColor(1)
    #print 'qcd mc: {} entries'.format(int(hist_mc_tot.Integral()))
  
    if not do_stack:
      legend.AddEntry(hist_mc_tot, 'MC - {}'.format(self.getQCDMCLabel(self.white_list[0], self.white_list[len(self.white_list)-1])))
      
    ## create stack histogram  
    hist_mc_stack = ROOT.THStack('hist_mc_stack', '')


    # compute the mc normalisation weight
    if do_luminorm: lumi_weight = self.getLumiWeight(selection=selection)
    #lumi_weight = Plotter()self.getLumiWeight(selection=selection)

    for hist_mc in mc_hists:
      if do_shape and int_mc_tot != 0: hist_mc.Scale(1/int_mc_tot)
      elif do_luminorm: hist_mc.Scale(lumi_weight)
      hist_mc_stack.Add(hist_mc)
    if do_shape and int_mc_tot!= 0: hist_mc_tot.Scale(1/int_mc_tot)
    elif do_luminorm: hist_mc_tot.Scale(lumi_weight)
    #hist_mc_tot.Scale(lumi_weight)

    #hist_data = self.tools.getRootXAxis(hist_data, label_size=0.0, title_size=0.0, offset=0)
    #hist_data = self.tools.getRootYAxis(hist_data, title='Entries' if not do_shape else 'Normalised to unity', label_size=0.037, title_size=0.042, offset=1.1, ymin=1e-9, ymax=self.getMaxRangeY(hist_data, hist_mc_stack, do_log))
    #canv.Update()

    #hist_data.GetXaxis().SetTitle(quantity.label)
    #hist_data.GetXaxis().SetLabelSize(0.0)
    #hist_data.GetXaxis().SetTitleSize(0.0)
    #hist_data.GetXaxis().SetTitleOffset(1.1)
    #hist_data.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
    #hist_data.GetYaxis().SetLabelSize(0.037)
    #hist_data.GetYaxis().SetTitleSize(0.042)
    #hist_data.GetYaxis().SetTitleOffset(1.1)
    #hist_data.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_data, hist_mc_stack, do_log))

    #if plot_data: hist = hist_data
    #elif not plot_data and not do_stack: hist = hist_mc_tot
    #else: 
    #  hist = hist_mc_stack

    hist_mc_tot.SetTitle('')
    if not plot_ratio: 
      hist_mc_tot.GetXaxis().SetTitle(quantity.title)
      hist_mc_tot.GetXaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
      hist_mc_tot.GetXaxis().SetTitleSize(0.042)
      hist_mc_tot.GetXaxis().SetTitleOffset(1.1)
    if plot_ratio:
      hist_mc_tot.GetXaxis().SetLabelSize(0.0)
      hist_mc_tot.GetXaxis().SetTitleSize(0.0)
    hist_mc_tot.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
    hist_mc_tot.GetYaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
    hist_mc_tot.GetYaxis().SetTitleSize(0.042)
    hist_mc_tot.GetYaxis().SetTitleOffset(1.3 if not plot_ratio else 1.1)
    if plot_data: hist_mc_tot.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_data_stack, hist_mc_stack, do_log))
    elif plot_sig: hist_mc_tot.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_mc_stack, do_log, use_sig=True))
    else: hist_mc_tot.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_mc_tot, hist_mc_stack, do_log))

    #ROOT.gStyle.SetPadLeftMargin(0.16) 
    ROOT.gStyle.SetOptStat(0)


    if plot_data: hist_data_tot.Draw()
    hist_mc_tot.Draw('histo')
    if do_stack:
      hist_mc_stack.Draw('histo same')
    else:
      hist_mc_tot.Draw('histo same')
    if plot_data: hist_data_tot.Draw('same') # making sure data points are always visible
    if plot_sig: 
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')

    # draw error bars
    hist_mc_tot_err = hist_mc_tot.Clone('hist_mc_tot_err')
    hist_mc_tot_err.SetLineWidth(0)
    hist_mc_tot_err.SetFillStyle(3244)
    hist_mc_tot_err.SetFillColor(ROOT.kGray+2)
    hist_mc_tot_err.Draw('E2 same')

    ## draw the legend
    legend.Draw('same')

    if do_luminorm:
      scale_text = ROOT.TPaveText(0.15, 0.83, 0.3, 0.88, "brNDC")
      scale_text.SetBorderSize(0)
      scale_text.SetFillColor(ROOT.kWhite)
      scale_text.SetTextSize(0.032)
      scale_text.SetTextAlign(31)
      scale_text.SetTextFont(42)
      scale_text.AddText('MC scaled by {}'.format(round(lumi_weight, 2)))
      #scale_text.Draw()

    # plot the ratio
    if plot_ratio:
      pad_down.cd()

      hist_ratio = self.tools.getRatioHistogram(hist_data_tot, hist_mc_tot)
      hist_ratio.Sumw2()

      for ibin in range(0, hist_ratio.GetNbinsX()+1):
        #print '{} {}'.format(hist_data.GetBinError(ibin), math.sqrt(hist_data.GetBinContent(ibin)))
        if hist_data_tot.GetBinContent(ibin) != 0 and hist_mc_tot.GetBinContent(ibin) != 0:
          #err = hist_ratio.GetBinContent(ibin) * (math.sqrt(hist_data.GetBinContent(ibin))/hist_data.GetBinContent(ibin) + math.sqrt(hist_mc_tot.GetBinContent(ibin))/hist_mc_tot.GetBinContent(ibin))
          #err = math.sqrt((math.sqrt(hist_data.GetBinContent(ibin))/hist_mc_tot.GetBinContent(ibin))**2 + (math.sqrt(hist_mc_tot.GetBinContent(ibin))*hist_data.GetBinContent(ibin)/(hist_mc_tot.GetBinContent(ibin))**2)**2)
          err = math.sqrt((hist_data_tot.GetBinError(ibin)/hist_mc_tot.GetBinContent(ibin))**2 + (hist_mc_tot.GetBinError(ibin)*hist_data_tot.GetBinContent(ibin)/(hist_mc_tot.GetBinContent(ibin))**2)**2)
        else: 
          err = 0
        if hist_ratio.GetBinContent(ibin) != 0: hist_ratio.SetBinError(ibin, err)

      hist_ratio.SetLineWidth(2)
      hist_ratio.SetMarkerStyle(20)
      hist_ratio.SetTitle('')
      hist_ratio.GetXaxis().SetTitle(self.quantity.title)

      #hist_ratio = self.tools.getRootXAxis(hist_ratio, title=self.quantity.title, label_size=0.1, title_size=0.13, offset=0.73)
      #val_min = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
      #val_max = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
      #hist_ratio = self.tools.getRootYAxis(hist_ratio, title='Data/MC', label_size=0.1, title_size=0.13, offset=0.345, ymin=val_min-0.15*val_min, ymax=val_max+0.15*val_max)

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

      hist_ratio.Draw('PE')

      ## draw line at ratio = 1
      line = ROOT.TLine(self.quantity.bin_min, 1, self.quantity.bin_max, 1)
      line.SetLineColor(4)
      line.SetLineWidth(2)
      line.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots/DataMCComparison', outdirlabel, do_shape, do_luminorm, do_stack, do_log)
    
    canv.SaveAs('{}/{}.png'.format(outputdir, self.quantity.label))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, self.quantity.label))


  def plotSignalBackgroundComparison(self, selection='', title='', outdirlabel='', branchname='nano', treename='Events', plot_ratio=False, do_shape=True, do_log=False):
    ROOT.gStyle.SetPadLeftMargin(0.12) 
    ROOT.gStyle.SetOptStat(0)

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, outdirlabel, do_log, do_shape)
    canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
    if do_log: canv.SetLogy()
    ROOT.SetOwnership(canv, False)
    #canv.SetGrid()
    canv.cd()

    # define the pads
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1) if plot_ratio else ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
    if plot_ratio: pad_up.SetBottomMargin(0.03)
    if do_log: pad_up.SetLogy()
    pad_up.Draw()
    canv.cd()
    if plot_ratio:
      pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
      pad_down.SetBottomMargin(0.25)
      pad_down.Draw()

    pad_up.cd()

    # prepare the legend
    legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.58, xmax=0.84, ymax=0.83, size=0.027)
    #legend = self.tools.getRootTLegend(xmin=0.53, ymin=0.65, xmax=0.9, ymax=0.83, size=0.027)

    label_top_right = ROOT.TPaveText(0.5,0.92,0.91,0.93, "brNDC")
    label_top_right.SetBorderSize(0)
    label_top_right.SetFillColor(ROOT.kWhite)
    label_top_right.SetTextSize(0.042)
    label_top_right.SetTextAlign(31)
    label_top_right.SetTextFont(42)
    label_top_right.AddText(title)
    label_top_right.Draw()

    # data
    hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
    hist_data_tot.Sumw2()
    int_data_tot = 0.
    for data_file in self.data_files:
      f_data = self.tools.getRootFile(data_file.filename, with_ext=False) #  ROOT.TFile.Open(self.data_file.filename, 'READ')
      hist_data_name = 'hist_data_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
      hist_data = self.tools.createHisto(f_data, treename, self.quantity, hist_name=hist_data_name, branchname=branchname, selection=selection)
      hist_data.Sumw2()
      if do_shape: 
        int_data_tot += hist_data.Integral()
      hist_data_tot.Add(hist_data)

    if int_data_tot != 0: hist_data_tot.Scale(1./int_data_tot)

    legend.AddEntry(hist_data_tot, 'data - {}'.format(self.getDataLabel(data_label) if len(self.data_files)>1 else data_file.label))

    ## set the style
    #hist_data.SetLineWidth(0)
    hist_data_tot.SetFillColor(ROOT.kBlue-3)
    hist_data_tot.SetFillStyle(3005)
    #hist_data.SetMarkerStyle(20)
    #hist_data.SetTitle(self.title)

    # signal
    signal_hists = []
    for signal_file in self.signal_files:
      #f_signal = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+signal_file.filename, 'READ')
      f_signal = ROOT.TFile.Open(signal_file.filename, 'READ')
      hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(self.quantity, outdirlabel, do_log, do_shape)
      #weight = 'weight_hlt'
      weight=-99
      hist_signal = self.tools.createHisto(f_signal, treename, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection='BToMuMuPi_isMatched==1 && BToMuMuPi_hnl_charge==0' if selection=='' else 'BToMuMuPi_isMatched==1 && BToMuMuPi_hnl_charge==0 &&'+selection)
      #hist_signal = self.tools.createHisto(f_signal, treename, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection=selection, weight=weight)
      hist_signal.Sumw2()
      if do_shape: 
        int_signal = hist_signal.Integral()
        if int_signal != 0: hist_signal.Scale(1/int_signal)
      legend.AddEntry(hist_signal, 'signal - {}'.format(signal_file.label))
      #legend.AddEntry(hist_signal, 'MC - {}'.format(signal_file.label))

      hist_signal.SetLineWidth(3)
      hist_signal.SetLineColor(signal_file.colour)
      signal_hists.append(hist_signal)
      #hist_signal.Draw('same')

    hist_data_tot.SetTitle('')
    hist_data_tot.GetXaxis().SetTitle(quantity.title)
    hist_data_tot.GetXaxis().SetLabelSize(0.033)
    hist_data_tot.GetXaxis().SetTitleSize(0.042)
    hist_data_tot.GetXaxis().SetTitleOffset(1.1)
    hist_data_tot.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
    hist_data_tot.GetYaxis().SetLabelSize(0.033)
    hist_data_tot.GetYaxis().SetTitleSize(0.042)
    hist_data_tot.GetYaxis().SetTitleOffset(1.3)
    #ymax = max(hist_data_tot.GetMaximum(), signal_hists[0].GetMaximum(), signal_hists[1].GetMaximum(), signal_hists[2].GetMaximum())
    ymax = max(hist_data_tot.GetMaximum(), signal_hists[0].GetMaximum())
    hist_data_tot.GetYaxis().SetRangeUser(1e-9, ymax+0.15*ymax)

    hist_data_tot.Draw('histo')
    for hist_sig in signal_hists:
      hist_sig.Draw('histo same')

    ## draw the legend
    legend.Draw('same')

    # plot the ratio
    if plot_ratio:
      pad_down.cd()

      hist_ratio = self.tools.getRatioHistogram(hist_data_tot, hist_signal)
      hist_ratio.Sumw2()

      for ibin in range(0, hist_ratio.GetNbinsX()+1):
        if hist_data_tot.GetBinContent(ibin) != 0 and hist_signal.GetBinContent(ibin) != 0:
          err = math.sqrt((hist_data_tot.GetBinError(ibin)/hist_signal.GetBinContent(ibin))**2 + (hist_signal.GetBinError(ibin)*hist_data_tot.GetBinContent(ibin)/(hist_signal.GetBinContent(ibin))**2)**2)
        else: 
          err = 0
        if hist_ratio.GetBinContent(ibin) != 0: hist_ratio.SetBinError(ibin, err)

      hist_ratio.SetLineWidth(2)
      hist_ratio.SetMarkerStyle(20)
      hist_ratio.SetTitle('')
      hist_ratio.GetXaxis().SetTitle(self.quantity.title)

      hist_ratio.GetXaxis().SetLabelSize(0.1)
      hist_ratio.GetXaxis().SetTitleSize(0.13)
      hist_ratio.GetXaxis().SetTitleOffset(0.73)
      hist_ratio.GetYaxis().SetTitle('Data/MC')
      hist_ratio.GetYaxis().SetLabelSize(0.1)
      hist_ratio.GetYaxis().SetTitleSize(0.13)
      hist_ratio.GetYaxis().SetTitleOffset(0.345)
      val_min = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
      val_max = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
      #hist_ratio.GetYaxis().SetRangeUser(val_min-0.15*val_min, val_max+0.15*val_max)
      hist_ratio.GetYaxis().SetRangeUser(0.5, 1.5)

      hist_ratio.Draw('PE')

      ## draw line at ratio = 1
      line = ROOT.TLine(self.quantity.bin_min, 1, self.quantity.bin_max, 1)
      line.SetLineColor(4)
      line.SetLineWidth(2)
      line.Draw('same')

    outputdir = self.tools.getOutDir('./myPlots/SignalBackgroundComparison', outdirlabel, do_shape, False, False, do_log)
    
    canv.SaveAs('{}/{}.png'.format(outputdir, self.quantity.label))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, self.quantity.label))


  def plotTwoSamples(file1, file2, branchname, tree1, tree2, selection1='', selection2='', legend1='legend1', legend2='legend2', do_printstat=False):
    f1 = ROOT.TFile.Open(file1, 'READ')
    f2 = ROOT.TFile.Open(file2, 'READ')

    #canv = ROOT.TCanvas('canv', 'canv', 900, 800)
    canv = self.tools.createTCanvas(name='canv', dimx=900, dimy=800)
    if self.do_log: canv.SetLogy()
    if not do_printstat: ROOT.gStyle.SetOptStat(0)
    
    hist1 = self.tools.createHisto(f1, tree1, self.quantity, hist_name=legend1, branchname=branchname, selection=selection1)
    hist2 = self.tools.createHisto(f2, tree2, self.quantity, hist_name=legend2, branchname=branchname, selection=selection2)
    
    if self.do_shape: 
      int1 = hist1.Integral()
      hist1.Scale(1/int1)
      int2 = hist2.Integral()
      hist2.Scale(1/int2)

    hist1.SetTitle(self.title)
    hist1.SetLineWidth(2)
    hist1.SetLineColor(ROOT.kOrange+1)
    hist1.GetXaxis().SetTitle(self.quantity.title)
    hist1.GetXaxis().SetLabelSize(0.037)
    hist1.GetXaxis().SetTitleSize(0.042)
    hist1.GetXaxis().SetTitleOffset(1.1)
    hist1.GetYaxis().SetTitle('Entries' if not self.do_shape else 'Normalised to unity')
    hist1.GetYaxis().SetLabelSize(0.037)
    hist1.GetYaxis().SetTitleSize(0.042)
    hist1.GetYaxis().SetTitleOffset(1.1)

    hist1.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist1, hist2, self.do_log))

    hist2.SetLineWidth(2)
    hist2.SetLineColor(ROOT.kBlue+2)

    hist1.Draw('histo')

    if do_printstat:
      ROOT.gPad.Update() 
      statsbox = ROOT.gPad.GetPrimitive("stats")
      y1 = statsbox.GetY1NDC()
      y2 = statsbox.GetY2NDC()
      newy1 = 2 * y1 - y2
      newy2 = y1
      statsbox.SetY1NDC(newy1)
      statsbox.SetY2NDC(newy2)

    hist2.Draw('histo sames')

    legend = self.tools.getRootTLegend(xmin=0.6, ymin=0.7, xmax=0.85, ymax=0.9, size=0.03)
    legend.AddEntry(hist2, legend2)
    legend.AddEntry(hist1, legend1)
    legend.Draw()

    outputdir = self.tools.getOutDir('./myPlots/comparison', self.outdirlabel, do_shape=self.do_shape, do_log=self.do_log)

    canv.SaveAs('{}/{}.png'.format(outputdir, self.quantity.label))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, self.quantity.label))


  def plotYields(self, lumi=0.774, selection='', title='', outdirlabel=''):

    canv = self.tools.createTCanvas('canv', dimx=900, dimy=800)
    canv.SetLogx()
    canv.SetLogy()
    canv.SetGrid()


    # gen yields
    gen_coupling_m1 = [5.37675839253e-06, 5.37675839253e-05, 0.000537675839253, 0.00537675839253, 0.0537675839253]
    #gen_yields_m1 = [0.233520553543, 21.1107914574, 2521.79127222, 198534.792249, 9787694.93131] # wrong weight
    #gen_yields_m1 = [0.222669464601, 21.7748154625, 2054.4123073, 164459.789617, 5551279.31665] # good weight, without trg mu eff
    gen_yields_m1 = [0.271587023481, 26.2785112766, 2480.63915671, 186917.31518, 6378867.68389] # good weight, with trg mu eff

    gen_coupling_m3 = [2.2126577747e-06, 2.2126577747e-05, 0.00022126577747, 0.0022126577747]
    #gen_yields_m3 = [0.0149191801294, 1.05017745728, 41.0813274363, 692.810973848]
    #gen_yields_m3 = [0.0199334286553, 1.14028241247, 23.582224495, 228.202498809]
    gen_yields_m3 = [0.0238935820319, 1.63795916047, 37.1400371471, 368.634025203]

    gen_coupling_m4p5 = [2.91378801607e-06, 2.91378801607e-05, 0.000291378801607, 0.00291378801607] 
    #gen_yields_m4p5 = [0.000314442951588, 0.0164840965898, 0.376472439332, 4.08626422153]
    #gen_yields_m4p5 = [0.000446445913129, 0.0127053455282, 0.160543306967, 1.68758530476]
    gen_yields_m4p5 = [0.000825665622275, 0.0250454472786, 0.322614030734, 3.41407673726]

    graph_gen_m1 = ROOT.TGraph()
    graph_gen_m3 = ROOT.TGraph()
    graph_gen_m4p5 = ROOT.TGraph()

    for pt in range(0, len(gen_coupling_m1)):
      point = graph_gen_m1.GetN()
      graph_gen_m1.SetPoint(point, gen_coupling_m1[pt], gen_yields_m1[pt])

    graph_gen_m1.SetMarkerStyle(22)
    graph_gen_m1.SetMarkerSize(2)
    graph_gen_m1.SetMarkerColor(ROOT.kOrange+0)
    graph_gen_m1.SetLineStyle(9)
    graph_gen_m1.SetLineWidth(2)
    graph_gen_m1.SetLineColor(ROOT.kOrange+0)

    for pt in range(0, len(gen_coupling_m3)):
      point = graph_gen_m3.GetN()
      graph_gen_m3.SetPoint(point, gen_coupling_m3[pt], gen_yields_m3[pt])

    graph_gen_m3.SetMarkerStyle(22)
    graph_gen_m3.SetMarkerSize(2)
    graph_gen_m3.SetMarkerColor(ROOT.kRed+1)
    graph_gen_m3.SetLineStyle(9)
    graph_gen_m3.SetLineWidth(2)
    graph_gen_m3.SetLineColor(ROOT.kRed+1)

    for pt in range(0, len(gen_coupling_m4p5)):
      point = graph_gen_m4p5.GetN()
      graph_gen_m4p5.SetPoint(point, gen_coupling_m4p5[pt], gen_yields_m4p5[pt])

    graph_gen_m4p5.SetMarkerStyle(22)
    graph_gen_m4p5.SetMarkerSize(2)
    graph_gen_m4p5.SetMarkerColor(ROOT.kRed+4)
    graph_gen_m4p5.SetLineStyle(9)
    graph_gen_m4p5.SetLineWidth(2)
    graph_gen_m4p5.SetLineColor(ROOT.kRed+4)

    from samples import signal_samples_limits_m1 as samples_m1
    from samples import signal_samples_limits_m3 as samples_m3
    from samples import signal_samples_limits_m4p5 as samples_m4p5

    graph_dummy = ROOT.TGraph() # move to TGraphAsymmErrors?
    graph_m1 = ROOT.TGraphAsymmErrors() #ROOT.TGraph() # move to TGraphAsymmErrors?
    graph_m3 = ROOT.TGraphAsymmErrors() #ROOT.TGraph() # move to TGraphAsymmErrors?
    graph_m4p5 = ROOT.TGraphAsymmErrors() #ROOT.TGraph() # move to TGraphAsymmErrors?

    graph_dummy.SetPoint(0, 1e-6, 1e-7)
    graph_dummy.SetPoint(1, 1e-1, 1e10)
    graph_dummy.SetMarkerStyle(0)
    graph_dummy.SetMarkerSize(0)
    graph_dummy.SetMarkerColor(0)
    #graph_dummy.SetTitle('Event {}'.format(ientry))
    graph_dummy.GetXaxis().SetTitle('|V^{2}|')
    graph_dummy.GetXaxis().SetLabelSize(0.037)
    graph_dummy.GetXaxis().SetTitleSize(0.042)
    graph_dummy.GetXaxis().SetTitleOffset(1.1)
    graph_dummy.GetYaxis().SetTitle('Signal yields')
    graph_dummy.GetYaxis().SetLabelSize(0.037)
    graph_dummy.GetYaxis().SetTitleSize(0.042)
    graph_dummy.GetYaxis().SetTitleOffset(1.1)

    print '\n mass 1'
    for signal_file in samples_m1:
      # get signal coupling
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau
      signal_v2 = getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

      # compute the signal yields
      signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=327.0e9) 

      # fill graph
      point = graph_m1.GetN()
      graph_m1.SetPoint(point, signal_v2, signal_yields)
      #graph_m1.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
      graph_m1.SetPointError(point, 0, 0, 0, 0)

    graph_m1.SetMarkerStyle(20)
    graph_m1.SetMarkerSize(2)
    graph_m1.SetMarkerColor(ROOT.kOrange+0)
    graph_m1.SetLineStyle(1)
    graph_m1.SetLineWidth(2)
    graph_m1.SetLineColor(ROOT.kOrange+0)


    print '\n mass 3'
    for signal_file in samples_m3:
      # get signal coupling
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau
      signal_v2 = getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

      # compute the signal yields
      signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=327.0e9) 

      # fill graph
      point = graph_m3.GetN()
      graph_m3.SetPoint(point, signal_v2, signal_yields)
      #graph_m3.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
      graph_m3.SetPointError(point, 0, 0, 0, 0)

    graph_m3.SetMarkerStyle(20)
    graph_m3.SetMarkerSize(2)
    graph_m3.SetMarkerColor(ROOT.kRed+1)
    graph_m3.SetLineStyle(1)
    graph_m3.SetLineWidth(2)
    graph_m3.SetLineColor(ROOT.kRed+1)

    print '\n mass 4.5'
    for signal_file in samples_m4p5:
      # get signal coupling
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau
      signal_v2 = getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

      # compute the signal yields
      signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=327.0e9) 

      # fill graph
      point = graph_m4p5.GetN()
      graph_m4p5.SetPoint(point, signal_v2, signal_yields)
      #graph_m4p5.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
      graph_m4p5.SetPointError(point, 0, 0, 0, 0)

    graph_m4p5.SetMarkerStyle(20)
    graph_m4p5.SetMarkerSize(2)
    graph_m4p5.SetMarkerColor(ROOT.kRed+4)
    graph_m4p5.SetLineStyle(1)
    graph_m4p5.SetLineWidth(2)
    graph_m4p5.SetLineColor(ROOT.kRed+4)

    graph_dummy.Draw('AP')  
    graph_m1.Draw('PL same')  
    graph_m3.Draw('PL same')  
    graph_m4p5.Draw('PL same')  
    graph_gen_m1.Draw('PL same')
    graph_gen_m3.Draw('PL same')
    graph_gen_m4p5.Draw('PL same')
      
    legend = self.tools.getRootTLegend(xmin=0.15, ymin=0.55, xmax=0.45, ymax=0.9, size=0.027)
    legend.AddEntry(graph_m1, 'm=1GeV, reco')
    legend.AddEntry(graph_gen_m1, 'm=1GeV, gen')
    legend.AddEntry(graph_m3, 'm=3GeV, reco')
    legend.AddEntry(graph_gen_m3, 'm=3GeV, gen')
    legend.AddEntry(graph_m4p5, 'm=4p5GeV, reco')
    legend.AddEntry(graph_gen_m4p5, 'm=4p5GeV, gen')
    legend.Draw()

    if not path.exists('./myPlots/yields'):
      os.system('mkdir -p ./myPlots/yields')

    canv.SaveAs('./myPlots/yields/signal_yields_gen_reco_corr.png')
    canv.SaveAs('./myPlots/yields/signal_yields_gen_reco_corr.pdf')



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  doDataMCComparison = False
  doSignalBackgroundComparison = True
  plotMCSR = False
  compareTwoDistributions = False
  plotYields = False

  plot_log = False
  plot_categories = False
  plot_SR = False
  plot_CR = True
  # add add_overflow


  if plotYields:
    from samples import signal_samples_limits_m4p5
    plotter = Plotter(signal_files=signal_samples_limits_m4p5)
    plotter.plotYields(lumi=41.6, selection='mu_isdsa!=1', title='', outdirlabel='')


  if doSignalBackgroundComparison:
    #baseline_selection = 'trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1'
    baseline_selection = 'BToMuMuPi_trg_mu_pt>7 && fabs(BToMuMuPi_trg_mu_eta)<1.5 && ProbeTracks_pt[BToMuMuPi_pi_idx]>0.7 && abs(ProbeTracks_eta[BToMuMuPi_pi_idx])<2 && abs(BToMuMuPi_pi_dz)>0.005 && abs(BToMuMuPi_pi_dxy)>0.001 && abs(BToMuMuPi_pi_dzS)>1.5 && abs(BToMuMuPi_pi_dxyS)>0.5 && BToMuMuPi_pi_DCASig>1 && Muon_pt[BToMuMuPi_sel_mu_idx]>2 && abs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 && abs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.001 && abs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.1 && BToMuMuPi_sv_prob>0.001 && BToMuMuPi_hnl_cos2D>0.95 && BToMuMuPi_mass<8 && BToMuMuPi_hnl_mass<6.3' 
    dirlabel = 'preselection_dsa_withpreselectionapplied'
    for quantity in quantities_to_plot_small:
      plotter = Plotter(quantity=quantity, data_files=data_samples_loose_dsaonly, signal_files=signal_samples_loose_dsaonly)
      #plotter.plotSignalBackgroundComparison(outdirlabel=dirlabel, branchname='nano', treename='Events', selection=baseline_selection, plot_ratio=False, do_shape=True, do_log=False)
    for quantity in quantities_trackId:
      plotter = Plotter(quantity=quantity, data_files=data_samples_loose_dsaonly, signal_files=signal_samples_loose_dsaonly)
      plotter.plotSignalBackgroundComparison(outdirlabel=dirlabel, branchname='nano', treename='Events', selection=baseline_selection, plot_ratio=False, do_shape=True, do_log=False)
    #dirlabel = 'test_JPsiToMuMu'
    #dirlabel = 'dataV06_tag_and_probe_v2_BToJPsiKstar_V0_sf_study15Sep21_A1_v1'
    #dirlabel = 'dataV07_BToJPsiKstar_V0_sf_study15Sep21_A1_v1'
    #for quantity in quantities_tag_and_probe:
    #  plotter = Plotter(quantity=quantity, data_files=data_samples_tag_and_probe, signal_files=signal_samples_tag_and_probe_BToJPsiKstar)
    #  plotter.plotSignalBackgroundComparison(outdirlabel=dirlabel, branchname='flat', treename='tree', selection=baseline_selection, plot_ratio=True, do_shape=True, do_log=False)


  #white_list_20to300 = ['QCD_pt20to30 (V04)', 'QCD_pt30to50 (V04)', 'QCD_pt50to80 (V04)', 'QCD_pt80to120 (V04)', 'QCD_pt80to120_ext (V04)', 'QCD_pt120to170 (V04)', 'QCD_pt120to170_ext (V04)', 'QCD_pt170to300 (V04)']
  white_list_20to300 = ['QCD_pt30to50 (V04)', 'QCD_pt50to80 (V04)', 'QCD_pt80to120 (V04)', 'QCD_pt120to170 (V04)', 'QCD_pt170to300 (V04)']
  white_list_20to30 = ['QCD_pt20to30 (V02)']

  for quantity in quantities_to_plot_small:
    small_plotter = Plotter(quantity=quantity, data_files=data_samples, qcd_files=qcd_samples, signal_files=signal_samples, white_list=white_list_20to300)
    #small_plotter.plotDataMCComparison(selection='hnl_charge!=0 && sv_lxysig>20', title='Control Region', outdirlabel='testing', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=True, do_log=False)

  if doDataMCComparison:

    #white_list_20to300 = ['QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)', 'QCD_pt170to300 (V02)']
    #white_list_15to300 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt170to300 (V02)']
    #white_list_15to300 = ['QCD_pt15to20 (V07_18Aug21)', 'QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
    white_list_15to300 = ['QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
    white_list_20to300 = ['QCD_pt20to30 (V05)', 'QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
    white_list_30to300 = ['QCD_pt30to50 (V05)', 'QCD_pt50to80 (V05)', 'QCD_pt80to120 (V05)', 'QCD_pt120to170 (V05)', 'QCD_pt170to300 (V05)']
    white_list_30to50 = ['QCD_pt30to50 (V05)']
    white_list_20to30 = ['QCD_pt20to30 (V05)']
    white_list_15to30 = ['QCD_pt15to20 (V05)', 'QCD_pt20to30 (V05)']
    #plotDataMCComparison(quantity, data_file=data_file, qcd_files=qcd_files, title='Charged #mu#pi CR', outdirlabel='dataV03_QCDV04_CR', selection='hnl_charge!=0', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #dirlabel = 'testing'
    #Plotter(quantity=quantity, title='Charged #mu#pi CR', outdirlabel=dirlabel, do_shape=True, do_stack=True, do_log=False).plotDataMCComparison(data_file=data_file, qcd_files=qcd_files, signal_files=signal_samples, branchname='flat', selection='hnl_charge!=0', white_list=white_list_20to300, plot_data=True, plot_qcdmc=True, plot_sig=False, plot_ratio=True)

    #baseline_selection = 'fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && hnl_cos2d>0.995 && sv_lxysig>20 && deltaeta_pi_fit_pi<0.015 && deltaphi_pi_fit_pi<0.03 && '
    #baseline_selection = 'fabs(mu_dzsig)>1 && fabs(mu_dxysig)>1.5 && sv_lxysig>20 && deltaeta_pi_fit_pi<0.015 && deltaphi_pi_fit_pi<0.03 && ' # if we want to categorise on cos2D, we dont want this cut to be too high
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && (trgmu_mu_mass<3.63 || trgmu_mu_mass>3.73) && '
    #baseline_selection = 'mu_isdsa !=1 && trgmu_looseid==1 && trgmu_softid==1 && '
    #baseline_selection = 'mu_isdsa !=1 && '
    #baseline_selection = 'mu_isdsa !=1 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && '
    #baseline_selection = 'mu_isdsa !=1 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6)) && '
    #baseline_selection = 'mu_isdsa !=1 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6 && '
    baseline_selection = 'mu_isdsa !=1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && '
    #baseline_selection = 'mu_isdsa != 1 && '

    if plot_SR:
      #dirlabel = 'dataV05_QCDV06_29Jun21'
      dirlabel = 'dataV06_tag_and_probe_v2_A1_nosf'
      for quantity in quantities_to_plot_small:
      #for quantity in quantities_muonId_study_displacedmuon_small:
        plotter = Plotter(quantity=quantity, data_files=data_samples, qcd_files=qcd_samples, signal_files=signal_samples, white_list=white_list_15to300)

        # signal region
        # inclusive
        plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0', title='Signal Region, inclusive', outdirlabel=dirlabel+'/SR/incl', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
        if plot_log:
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0', title='Signal Region, inclusive', outdirlabel=dirlabel+'/SR/incl', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)

        if plot_categories:
          # OS categories
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy<1 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}<1cm, OS', outdirlabel=dirlabel+'/SR/lxy0to1_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', title='Signal Region, (1<l_{xy}<5)cm, OS', outdirlabel=dirlabel+'/SR/lxy1to5_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', title='Signal Region, (5<l_{xy}<10)cm, OS', outdirlabel=dirlabel+'/SR/lxy5to10_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}>5cm, OS', outdirlabel=dirlabel+'/SR/lxygt5_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>10 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}>10cm, OS', outdirlabel=dirlabel+'/SR/lxygt10_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)


          if plot_log:
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy<1 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}<1cm, OS', outdirlabel=dirlabel+'/SR/lxy0to1_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', title='Signal Region, (1<l_{xy}<5)cm, OS', outdirlabel=dirlabel+'/SR/lxy1to5_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', title='Signal Region, (5<l_{xy}<10)cm, OS', outdirlabel=dirlabel+'/SR/lxy5to10_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}>5cm, OS', outdirlabel=dirlabel+'/SR/lxygt5_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>10 && trgmu_charge!=mu_charge', title='Signal Region, l_{xy}>10cm, OS', outdirlabel=dirlabel+'/SR/lxygt10_OS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)


          # SS categories
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy<1 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}<1cm, SS', outdirlabel=dirlabel+'/SR/lxy0to1_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', title='Signal Region, (1<l_{xy}<5)cm, SS', outdirlabel=dirlabel+'/SR/lxy1to5_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', title='Signal Region, (5<l_{xy}<10)cm, SS', outdirlabel=dirlabel+'/SR/lxy5to10_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}>5cm, SS', outdirlabel=dirlabel+'/SR/lxygt5_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>10 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}>10cm, SS', outdirlabel=dirlabel+'/SR/lxygt10_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=False, do_stack=False, do_log=False)

          if plot_log:
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy<1 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}<1cm, SS', outdirlabel=dirlabel+'/SR/lxy0to1_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=False, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', title='Signal Region, (1<l_{xy}<5)cm, SS', outdirlabel=dirlabel+'/SR/lxy1to5_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=False, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', title='Signal Region, (5<l_{xy}<10)cm, SS', outdirlabel=dirlabel+'/SR/lxy5to10_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=False, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge==0 && sv_lxy>5 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}>5cm, SS', outdirlabel=dirlabel+'/SR/lxygt5_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=False, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection='hnl_charge==0 && sv_lxy>10 && trgmu_charge==mu_charge', title='Signal Region, l_{xy}>10cm, SS', outdirlabel=dirlabel+'/SR/lxygt10_SS', branchname='flat', plot_data=False, plot_sig=True, plot_ratio=False, do_shape=True, do_stack=False, do_log=True)


    if plot_CR:
      #dirlabel = 'dataV04_QCDV05_v1'
      #dirlabel = 'test'
      dirlabel = 'dataV07_QCDV07_18Aug21_sf_study15Sep21_A1_v5_otherQCDrange_20_300'
      for quantity in quantities_to_plot_small:
        plotter = Plotter(quantity=quantity, data_files=data_samples, qcd_files=qcd_samples, signal_files=signal_samples, white_list=white_list_15to300)
        # control region
        # inclusive
        plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0', title='Control Region, inclusive', outdirlabel=dirlabel+'/CR/incl', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
        if plot_log:
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0', title='Control Region, inclusive', outdirlabel=dirlabel+'/CR/incl', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=True, do_log=True)
    
        if plot_categories:
          # OS categories
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy<1 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}<1cm, OS', outdirlabel=dirlabel+'/CR/lxy0to1_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', title='Control Region, (1<l_{xy}<5)cm, OS', outdirlabel=dirlabel+'/CR/lxy1to5_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', title='Control Region, (5<l_{xy}<10)cm, OS', outdirlabel=dirlabel+'/CR/lxy5to10_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}>5cm, OS', outdirlabel=dirlabel+'/CR/lxygt5_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>10 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}>10cm, OS', outdirlabel=dirlabel+'/CR/lxygt10_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)

          if plot_log:
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy<1 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}<1cm, OS', outdirlabel=dirlabel+'/CR/lxy0to1_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge!=mu_charge', title='Control Region, (1<l_{xy}<5)cm, OS', outdirlabel=dirlabel+'/CR/lxy1to5_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge!=mu_charge', title='Control Region, (5<l_{xy}<10)cm, OS', outdirlabel=dirlabel+'/CR/lxy5to10_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}>5cm, OS', outdirlabel=dirlabel+'/CR/lxygt5_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>10 && trgmu_charge!=mu_charge', title='Control Region, l_{xy}>10cm, OS', outdirlabel=dirlabel+'/CR/lxygt10_OS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)

          # SS categories
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy<1 && trgmu_charge==mu_charge', title='Control Region, l_{xy}<1cm, SS', outdirlabel=dirlabel+'/CR/lxy0to1_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', title='Control Region, (1<l_{xy}<5)cm, SS', outdirlabel=dirlabel+'/CR/lxy1to5_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', title='Control Region, (5<l_{xy}<10)cm, SS', outdirlabel=dirlabel+'/CR/lxy5to10_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && trgmu_charge==mu_charge', title='Control Region, l_{xy}>5cm, SS', outdirlabel=dirlabel+'/CR/lxygt5_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>10 && trgmu_charge==mu_charge', title='Control Region, l_{xy}>10cm, SS', outdirlabel=dirlabel+'/CR/lxygt10_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)

          if plot_log:
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy<1 && trgmu_charge==mu_charge', title='Control Region, l_{xy}<1cm, SS', outdirlabel=dirlabel+'/CR/lxy0to1_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>1 && sv_lxy<5 && trgmu_charge==mu_charge', title='Control Region, (1<l_{xy}<5)cm, SS', outdirlabel=dirlabel+'/CR/lxy1to5_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && sv_lxy<10 && trgmu_charge==mu_charge', title='Control Region, (5<l_{xy}<10)cm, SS', outdirlabel=dirlabel+'/CR/lxy5to10_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>5 && trgmu_charge==mu_charge', title='Control Region, l_{xy}>5cm, SS', outdirlabel=dirlabel+'/CR/lxygt5_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)
            plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0 && sv_lxy>10 && trgmu_charge==mu_charge', title='Control Region, l_{xy}>10cm, SS', outdirlabel=dirlabel+'/CR/lxygt10_SS', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=True)




  '''
    white_list_15to20 = ['QCD_pt15to20 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to20', selection='hnl_charge!=0', white_list=white_list_15to20, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to20', selection='hnl_charge!=0', white_list=white_list_15to20, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to30 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to30', selection='hnl_charge!=0', white_list=white_list_15to30, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to30', selection='hnl_charge!=0', white_list=white_list_15to30, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to50 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to50', selection='hnl_charge!=0', white_list=white_list_15to50, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to50', selection='hnl_charge!=0', white_list=white_list_15to50, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to80 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to80', selection='hnl_charge!=0', white_list=white_list_15to80, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to80', selection='hnl_charge!=0', white_list=white_list_15to80, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to120 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to120', selection='hnl_charge!=0', white_list=white_list_15to120, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to120', selection='hnl_charge!=0', white_list=white_list_15to120, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to170 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to170', selection='hnl_charge!=0', white_list=white_list_15to170, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to170', selection='hnl_charge!=0', white_list=white_list_15to170, do_shape=True, do_stack=True, do_log=True) 

    white_list_15to300 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)', 'QCD_pt170to300 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to300', selection='hnl_charge!=0', white_list=white_list_15to300, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt15to300', selection='hnl_charge!=0', white_list=white_list_15to300, do_shape=True, do_stack=True, do_log=True) 

    white_list_20to300 = ['QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)', 'QCD_pt170to300 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt20to300', selection='hnl_charge!=0', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt20to300', selection='hnl_charge!=0', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=True) 

    white_list_30to300 = ['QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)', 'QCD_pt170to300 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt30to300', selection='hnl_charge!=0', white_list=white_list_30to300, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt30to300', selection='hnl_charge!=0', white_list=white_list_30to300, do_shape=True, do_stack=True, do_log=True) 

    white_list_20to50 = ['QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)']
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt20to50', selection='hnl_charge!=0', white_list=white_list_20to50, do_shape=True, do_stack=True, do_log=False) 
    plotDataMCComparison(quantity, title='Data / MC comparison in the SS Control Region', outdirlabel='SS_CR_pt20to50', selection='hnl_charge!=0', white_list=white_list_20to50, do_shape=True, do_stack=True, do_log=True) 
    '''


if plotMCSR:
  for quantity in quantities_to_plot_small:
    #white_list_15to300 = ['QCD_pt15to20 (V02)', 'QCD_pt20to30 (V02)', 'QCD_pt30to50 (V02)', 'QCD_pt50to80 (V02)', 'QCD_pt80to120 (V02)', 'QCD_pt80to120_ext (V02)', 'QCD_pt120to170 (V02)', 'QCD_pt120to170_ext (V02)', 'QCD_pt170to300 (V02)']
    #plotMC(quantity, title='MC in the Signal Region', outdirlabel='SR_pt15to300_dispincl', selection='hnl_charge==0', white_list=white_list_15to300, do_shape=True, do_stack=True, do_log=False) 

    white_list_20to300 = ['QCD_pt20to30 (V04)', 'QCD_pt30to50 (V04)', 'QCD_pt50to80 (V04)', 'QCD_pt80to120 (V04)', 'QCD_pt120to170 (V04)', 'QCD_pt170to300 (V04)']
    Plotter(quantity=quantity, title='MC in the Signal Region', outdirlabel='SR_pt20to300_dispincl', do_shape=True, do_stack=True, do_log=False).plotMC(qcd_files=qcd_samples, selection='hnl_charge==0', white_list=white_list_20to300)
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region', outdirlabel='SR_pt20to300_dispincl', selection='hnl_charge==0', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    


    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, SV lxy<0.5cm', outdirlabel='SR_pt20to300_disp0to0p5', selection='hnl_charge==0 && sv_lxy<0.5', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (0.5 < SV lxy < 1)cm', outdirlabel='SR_pt20to300_disp0p5to1', selection='hnl_charge==0 && sv_lxy>0.5 && sv_lxy<1', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (1 < SV lxy < 3)cm', outdirlabel='SR_pt20to300_disp1to3', selection='hnl_charge==0 && sv_lxy>1 && sv_lxy<3', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (3 < SV lxy < 5)cm', outdirlabel='SR_pt20to300_disp3to5', selection='hnl_charge==0 && sv_lxy>3 && sv_lxy<5', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (3 < SV lxy < 15)cm', outdirlabel='SR_pt20to300_disp3to15', selection='hnl_charge==0 && sv_lxy>3 && sv_lxy<15', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (5<SV lxy<10)cm', outdirlabel='SR_pt20to300_disp5to10', selection='hnl_charge==0 && sv_lxy>5 && sv_lxy<10', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (10<SV lxy<20)cm', outdirlabel='SR_pt20to300_disp10to20', selection='hnl_charge==0 && sv_lxy>10 && sv_lxy<20', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (20<SV lxy<40)cm', outdirlabel='SR_pt20to300_disp20to40', selection='hnl_charge==0 && sv_lxy>20 && sv_lxy<40', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 
    #plotMC(quantity, qcd_files=qcd_samples, title='MC in the Signal Region, (40<SV lxy<100)cm', outdirlabel='SR_pt20to300_disp40to100', selection='hnl_charge==0 && sv_lxy>40 && sv_lxy<100', white_list=white_list_20to300, do_shape=True, do_stack=True, do_log=False) 


if compareTwoDistributions:
  for quantity in quantities_to_plot_small:
    file1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/old/bparknano_alltrgmus_looseselectionv2_matched.root'
    file2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_withlooseselection.root' 
    #plotTwoSamples(quantity, file1=file1, file2=file2, branchname='nano', tree1='Events', tree2='Events', legend1='signal (3GeV, 184mm)', legend2='background (data)', title='All candidates, no selection', outdirlabel='data_V02_B4_withlooseselection_signal_V20_emu_alltrgmus_looseselectionv2_matched', selection='BToMuMuPi_hnl_charge==0 || BToMuMuPi_hnl_charge!=0', do_shape=True, do_log=False)
    #plotTwoSamples(quantity, file1=file1, file2=file2, branchname='nano', tree1='Events', tree2='Events', legend1='signal (3GeV, 184mm)', legend2='background (data)', title='All candidates, no selection', outdirlabel='data_V02_B4_withlooseselection_signal_V20_emu_alltrgmus_looseselectionv2_matched', selection='BToMuMuPi_hnl_charge==0 || BToMuMuPi_hnl_charge!=0', do_shape=True, do_log=True)


  #file1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selected_alltrgmu_muonId.root'
  file1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_alltrgmu_muonId_iso.root'
  outdirlabel = 'V20_emu_looseselection_alltrgmu_muonId_triggermuon'
  for quantity in quantities_muonId_study_triggermuon:
    plotTwoSamples(quantity, file1=file1, file2=file1, branchname='nano', tree1='Events', tree2='Events', selection1='nBToMuMuPi>0 && Muon_pt[BToMuMuPi_trg_mu_idx]>7', selection2='nBToMuMuPi>0 && BToMuMuPi_isMatched==1 && Muon_pt[BToMuMuPi_trg_mu_idx]>7', legend1='all', legend2='matched', title='', outdirlabel=outdirlabel, do_shape=True, do_log=False, do_printstat=True)

  outdirlabel = 'V20_emu_looseselection_alltrgmu_muonId_displacedmuon'
  for quantity in quantities_muonId_study_displacedmuon_small:
    #selection = ' && Muon_inTimeMuon[BToMuMuPi_sel_mu_idx] && !((Muon_isGlobalMuon[BToMuMuPi_sel_mu_idx]==1 || Muon_isTrackerMuon[BToMuMuPi_sel_mu_idx]==1) && Muon_trackerHighPurityFlag[BToMuMuPi_sel_mu_idx]==1 && Muon_validHitFraction[BToMuMuPi_sel_mu_idx]>0.75)'
    #selection = ' && Muon_inTimeMuon[BToMuMuPi_sel_mu_idx] && ((Muon_isGlobalMuon[BToMuMuPi_sel_mu_idx]==1 || Muon_isTrackerMuon[BToMuMuPi_sel_mu_idx]==1) && Muon_trackerHighPurityFlag[BToMuMuPi_sel_mu_idx]==1 && Muon_validHitFraction[BToMuMuPi_sel_mu_idx]>0.75)'
    selection = ''
    #plotTwoSamples(quantity, file1=file1, file2=file1, branchname='nano', tree1='Events', tree2='Events', selection1='nBToMuMuPi>0 && Muon_pt[BToMuMuPi_trg_mu_idx]>7'+selection, selection2='nBToMuMuPi>0 && Muon_pt[BToMuMuPi_trg_mu_idx]>7 && BToMuMuPi_isMatched==1'+selection, legend1='all', legend2='matched', title='', outdirlabel=outdirlabel, do_shape=True, do_log=False, do_printstat=True)



