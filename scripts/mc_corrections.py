import os
import os.path
from os import path
import ROOT
from ROOT import gROOT, gStyle
import math
import numpy as np

from tools import Tools

import sys
sys.path.append('../objects')
from quantity import Quantity
from samples import data_samples, signal_samples
from baseline_selection import selection


class Plotter(Tools):
  def __init__(self, quantity='', data_files='', signal_files='', selection='', save_eos=False):
    self.tools = Tools()
    self.quantity = quantity
    self.data_files = data_files
    self.signal_files = signal_files
    self.selection = selection
    self.save_eos = save_eos


  def getQCDMCLabel(self, label1, label2, label3):
    '''
      based on MC label of the format ptXXtoYY 
    '''
    return label1[0:label1.find('to')] + label2[label2.find('to'):len(label2)] + label3[label3.find(' '):]


  def getDataLabel(self, label, version_label):
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
    if len(version_label) != 0:
      data_label += ' {}'.format(version_label)
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


  def createHisto(self, tree, quantity, hist_name='hist', branchname='flat', weight=-99, selection='', scale='', smear=''):
    ROOT.TH1.SetDefaultSumw2()
    hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)

    if 'signal' in hist_name:
      qte = '({})*({})*({})'.format(scale, smear, quantity.name_flat)
    else:
      qte = quantity.name_flat
      if '_uncorrected' in qte: qte = qte.replace('_uncorrected', '')
    #print qte

    if quantity.name_flat != 'probe_dxy_bs_uncorrected':
      selection += ' && probe_dxy_bs > 0.015' #TODO apply correction once it's derived
    else:
      selection += ' && {} > 0.015'.format(qte)
    if selection == '' and weight == -99: selection_string = selection
    elif selection != '' and weight == -99: selection_string = '({})'.format(selection)
    elif selection == ''  and weight != -99: selection_string = '({})'.format(weight)
    else: selection_string = '({sel}) * ({wght})'.format(sel=selection, wght=weight)
    #print selection_string

    tree.Project(hist_name, qte, selection_string)

    hist.SetDirectory(0)
    return hist


  def plot(self, 
      title='', 
      branchname='flat', 
      treename='tree', 
      add_weight_hlt=False, 
      add_weight_pu=False, 
      add_weight_muid=False, 
      weight_hlt='', 
      weight_puqcd='',
      weight_pusig='',
      weight_mu0id='', 
      weight_muid='', 
      plot_data=True, 
      plot_sig=True, 
      plot_ratio=True, 
      do_shape=True, 
      do_luminorm=False, 
      do_stack=True, 
      do_log=False, 
      add_overflow=True, 
      add_CMSlabel=True, 
      CMS_tag='Preliminary', 
      do_tdrstyle=True,
      scale='',
      smear=''
      ):

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, scale, smear, do_log)
    canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
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
      legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.45, xmax=0.84, ymax=0.83, size=0.027)
    else:
      if not do_tdrstyle:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.65, xmax=0.84, ymax=0.83, size=0.027)
      else:
        legend = self.tools.getRootTLegend(xmin=0.42, ymin=0.57, xmax=0.79, ymax=0.79, size=0.038)

    pad_up.cd()

    # data
    if plot_data:
      hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_data_tot.Sumw2()

      int_data_tot = 0.
      overflow_data_tot = 0
      error_overflow_data_tot = 0.

      hist_data_stack = ROOT.THStack('hist_data_stack', '')
      data_label = ''

      for idata, data_file in enumerate(self.data_files):
        f_data = self.tools.getRootFile(data_file.filename)
        tree_data = self.tools.getTree(f_data, treename)
        hist_data_name = 'hist_data_{}_{}_{}'.format(self.quantity.label, scale, smear)
        # apply sweight on data
        hist_data = self.createHisto(tree_data, self.quantity, hist_name=hist_data_name, branchname=branchname, selection=self.selection, weight='nbkg_sw')
        hist_data.Sumw2()
        if add_overflow:
          overflow_data_tot = (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
          error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
          hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
          hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
          hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
          hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)
        if do_shape: 
          int_data_tot += hist_data.Integral()

        hist_data_tot.Add(hist_data)
        if idata == 0: data_label = data_file.label[:20]
        else: data_label += '_{}'.format(data_file.label[:20])

        version_label = data_file.label[data_file.label.find('('):data_file.label.find(')')+1]

        if do_shape and hist_data.Integral() != 0: 
          hist_data.Scale(1./hist_data.Integral())
          hist_data_stack.Add(hist_data)

      if do_shape and int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

      if not do_tdrstyle:
        legend.AddEntry(hist_data_tot, 'data - {}'.format(self.getDataLabel(data_label, version_label) if len(self.data_files)>1 else data_file.label))
      else:
        legend.AddEntry(hist_data_tot, 'data-driven background')

      ## set the style
      hist_data_tot.SetMarkerStyle(20)
      #hist_data_tot.SetFillColor(ROOT.kBlue-3)
      #hist_data_tot.SetFillStyle(3005)

    # signal
    if plot_sig:
      signal_hists = []
      for signal_file in self.signal_files:
        f_sig = self.tools.getRootFile(signal_file.filename)
        tree_sig = self.tools.getTree(f_sig, treename)

        hist_signal_name = 'hist_signal_{}_{}_{}'.format(self.quantity, scale, smear)
        matching_selection = 'ismatched==1' if branchname == 'flat' else 'BToMuMuPi_isMatched==1'
        selection_signal = matching_selection if self.selection == '' else matching_selection + ' && ' + self.selection
        weight_sig = '(1)'
        if add_weight_hlt : weight_sig += ' * ({})'.format(weight_hlt)
        if add_weight_pu : weight_sig += ' * ({})'.format(weight_pusig)
        if add_weight_muid : weight_sig += ' * ({}) *({})'.format(weight_mu0id, weight_muid)

        hist_signal = self.createHisto(tree_sig, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection=selection_signal, weight=weight_sig, scale=scale, smear=smear)
        hist_signal.Sumw2()
        if add_overflow:
          overflow_signal_tot = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
          error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
          hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
          hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
          hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
          hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)
        if do_shape: 
          int_signal = hist_signal.Integral()
          if int_signal != 0: hist_signal.Scale(1/int_signal)

        legend.AddEntry(hist_signal, 'signal - {}'.format(signal_file.label))

        #hist_signal.SetLineWidth(3)
        #hist_signal.SetLineColor(signal_file.colour)
        hist_signal.SetFillColor(ROOT.kOrange+7)
        hist_signal.SetFillStyle(3005)
        signal_hists.append(hist_signal)

      
    frame = hist_data_tot.Clone('frame')

    frame.SetTitle('')
    if not plot_ratio: 
      frame.GetXaxis().SetTitle(quantity.title)
      frame.GetXaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
    if plot_ratio:
      frame.GetXaxis().SetLabelSize(0.0)
      frame.GetXaxis().SetTitleSize(0.0)
    frame.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
    frame.GetYaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
    frame.GetYaxis().SetTitleSize(0.042)
    frame.GetYaxis().SetTitleOffset(1.3 if not plot_ratio else 1.1)
    frame.GetYaxis().SetRangeUser(1e-4, self.getMaxRangeY(signal_hists, hist_data_tot, do_log, use_sig=True))

    #ROOT.gStyle.SetPadLeftMargin(0.16) 
    ROOT.gStyle.SetOptStat(0)

    # draw the distributions
    frame.Draw()
    
    hist_data_tot.Draw('same')
    for hist_sig in signal_hists:
      hist_sig.Draw('histo same')
    hist_data_tot.Draw('same')

    # draw the legend
    legend.Draw('same')

    # add labels
    if not do_tdrstyle:
      self.tools.printLatexBox(0.65, 0.86, title, size=0.04 if plot_ratio else 0.036)
    else:
      self.tools.printLatexBox(0.60, 0.84, title, size=0.04 if plot_ratio else 0.038)
    if add_CMSlabel: self.tools.printCMSTag(pad_up, CMS_tag, size=0.55 if plot_ratio else 0.43)
    if do_luminorm:
      scale_text = ROOT.TPaveText(0.15, 0.83, 0.3, 0.88, "brNDC")
      scale_text.SetBorderSize(0)
      scale_text.SetFillColor(ROOT.kWhite)
      scale_text.SetTextSize(0.032)
      scale_text.SetTextAlign(31)
      scale_text.SetTextFont(42)
      scale_text.AddText('MC scaled by {}'.format(round(lumi_weight/self.tools.getDataLumi(self.data_files), 2)))
      scale_text.Draw()

    # plot the ratio
    if plot_ratio:
      pad_down.cd()

      hist_ratio = self.tools.getRatioHistogram(hist_data_tot, signal_hists[0])
      hist_ratio.Sumw2()

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
      hist_ratio.GetYaxis().SetRangeUser(val_min-0.15*val_min, val_max+0.15*val_max)

      hist_ratio.Draw('PE')

      ## draw line at ratio = 1
      line = ROOT.TLine(self.quantity.bin_min, 1, self.quantity.bin_max, 1)
      line.SetLineColor(4)
      line.SetLineWidth(2)
      line.Draw('same')

    if not self.save_eos:
      outputdir = self.tools.getOutDir('./myPlots/mc_corrections', outdirlabel=quantity.label, do_log=do_log)
    else:
      outputdir = self.tools.getOutDir('/eos/home-a/anlyon/www/BHNL/mc_corrections', outdirlabel=quantity.label, do_log=do_log)
    
    canv.SaveAs('{}/{}_scale{}_smear{}.png'.format(outputdir, self.quantity.label, str(round(scale, 3)).replace('.', 'p'), smear[smear.rfind('_')+1:]))
    canv.SaveAs('{}/{}_scale{}_smear{}.pdf'.format(outputdir, self.quantity.label, str(round(scale, 3)).replace('.', 'p'), smear[smear.rfind('_')+1:]))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  dxy = Quantity(name_flat='probe_dxy_bs_uncorrected', label='probe_dxy_bs', title='probe #mu |d_{xy}| (BS) [cm]', nbins=90, bin_min=0, bin_max=0.5)
  dxy_sig = Quantity(name_flat='probe_dxy_sig_bs_uncorrected', label='probe_dxy_sig_bs', title='probe #mu |d_{xy}| significance (BS)', nbins=90, bin_min=0, bin_max=60)
  #quantities = [dxy, dxy_sig]
  quantities = [dxy_sig]

  data_files = data_samples['tag_and_probe']
  signal_files = signal_samples['tag_and_probe_BToJPsiKstar']

  baseline_selection = selection['tag_and_probe'].flat

  save_eos = True

  scales = np.arange(1., 1.3, 0.01)
  smears = [
            'smeared_corr_0p01',
            #'smeared_corr_0p02',
            'smeared_corr_0p03',
            #'smeared_corr_0p04',
            'smeared_corr_0p05', 
            #'smeared_corr_0p06',
            #'smeared_corr_0p07',
            #'smeared_corr_0p08',
            'smeared_corr_0p1',
            #'smeared_corr_0p15',
            'smeared_corr_0p2',
            'smeared_corr_0p25',
            'smeared_corr_0p3',
            'smeared_corr_0p35',
            'smeared_corr_0p4',
            'smeared_corr_0p45',
            'smeared_corr_0p5',
           ]
  #scales = [1.1, 1.15]
  #smears = ['smeared_corr_0p05', 'smeared_corr_0p1', 'smeared_corr_0p5']

  for quantity in quantities:
    plotter = Plotter(quantity=quantity, data_files=data_files, signal_files=signal_files, selection=baseline_selection, save_eos=save_eos)
    for scale in scales:
      for smear in smears:
        plotter.plot(scale=scale, smear=smear, do_log=True)
        plotter.plot(scale=scale, smear=smear, do_log=False)




