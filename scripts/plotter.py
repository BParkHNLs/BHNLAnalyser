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
from computeYields import ComputeYields

import sys
sys.path.append('../objects')
from quantity import Quantity, quantities
from samples import data_samples, qcd_samples, signal_samples
from categories import categories
from selection import selection


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to produce the main analysis plots', add_help=True)
  parser.add_argument('--outdirlabel'     , type=str, dest='outdirlabel'     , help='name of the outdir'                                            , default=None)
  parser.add_argument('--data_label'      , type=str, dest='data_label'      , help='which data samples to consider?'                               , default='V07_18Aug21')
  parser.add_argument('--qcd_label'       , type=str, dest='qcd_label'       , help='which qcd samples to consider?'                                , default='V07_18Aug21')
  parser.add_argument('--signal_label'    , type=str, dest='signal_label'    , help='which signal samples to consider?'                             , default='private')
  parser.add_argument('--quantities_label', type=str, dest='quantities_label', help='which quantity collection to consider?'                        , default='small')
  parser.add_argument('--selection_label' , type=str, dest='selection_label' , help='apply a baseline selection_label?'                             , default='standard')
  parser.add_argument('--categories_label', type=str, dest='categories_label', help='which phase space categorisation?'                             , default='standard')
  parser.add_argument('--sample_type'     , type=str, dest='sample_type'     , help='run the plotter on a nano or flat sample?'                     , default='flat', choices=['nano', 'flat'])
  parser.add_argument('--tree_name'       , type=str, dest='tree_name'       , help='name of the tree to analyse'                                   , default='signal_tree')
  parser.add_argument('--white_list '     , type=str, dest='white_list'      , help='pthat range to consider for qcd samples'                       , default='')
  parser.add_argument('--CMStag '         , type=str, dest='CMStag'          , help='CMS tag to be added if --add_CMSlabel'                         , default='Preliminary')
  parser.add_argument('--plot_CR'         ,           dest='plot_CR'         , help='plot QCDMC/data in the CR'                , action='store_true', default=False)
  parser.add_argument('--plot_SR'         ,           dest='plot_SR'         , help='plot QCDMC/signal in the SR'              , action='store_true', default=False)
  parser.add_argument('--plot_dataSig'    ,           dest='plot_dataSig'    , help='plot data (as bkg) and signal'            , action='store_true', default=False)
  parser.add_argument('--plot_ratio'      ,           dest='plot_ratio'      , help='plots the ratio'                          , action='store_true', default=False)
  parser.add_argument('--do_shape'        ,           dest='do_shape'        , help='normalise to unity'                       , action='store_true', default=False)
  parser.add_argument('--do_luminorm'     ,           dest='do_luminorm'     , help='normalise to the luminosity'              , action='store_true', default=False)
  parser.add_argument('--do_stack'        ,           dest='do_stack'        , help='create stack histtogram'                  , action='store_true', default=False)
  parser.add_argument('--do_log'          ,           dest='do_log'          , help='put the Y axis in log scale'              , action='store_true', default=False)
  parser.add_argument('--add_overflow'    ,           dest='add_overflow'    , help='add overflow bin'                         , action='store_true', default=False)
  parser.add_argument('--add_CMSlabel'    ,           dest='add_CMSlabel'    , help='add CMS label'                            , action='store_true', default=False)
  #parser.add_argument('--submit_batch', dest='submit_batch', help='submit on the batch?', action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):
  '''
    This function checks the validaty of the parser, and determines whether the parser was used by
    checking the outdirlabel against its default
  '''
  used_parser = (opt.outdirlabel != None)
  if used_parser and opt.plot_CR == False and opt.plot_SR == False and opt.plot_dataSig == False:
    raise RuntimeError('Please indicate the kind of plots you want to produce with --plot_CR and/or --plot_SR and/or --plot_dataSig')

  if used_parser and opt.do_shape == opt.do_luminorm:
    raise RuntimeError('Invalid arguments for --do_shape ({}) and --do_luminorm ({}) \nPlease only set exactly one option to True at a time'.format(opt.do_shape, opt.do_luminorm))

  return used_parser



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


  def plot(self, selection='', title='', outdirlabel='', plotdirlabel='', branchname='flat', treename='signal_tree', plot_data=False, plot_qcd=False, plot_sig=False, plot_ratio=False, do_shape=True, do_luminorm=False, do_stack=True, do_log=False, add_overflow=False, add_CMSlabel=True, CMS_tag='Preliminary'):

    # check the options
    if plot_data and self.data_files == '':
      raise RuntimeError('Please specify on which data sample to run')
    if plot_sig and self.signal_files == '':
      raise RuntimeError('Please specify on which signal sample to run')
    if not plot_data:
      plot_ratio = False

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, outdirlabel.replace('/', '_'), do_log, do_shape)
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
      legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.65, xmax=0.84, ymax=0.83, size=0.027)

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
        f_data = self.tools.getRootFile(data_file.filename, with_ext=False)
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
      if plot_data and plot_qcd:
        hist_data_tot.SetMarkerStyle(20)
      else:
        hist_data_tot.SetFillColor(ROOT.kBlue-3)
        hist_data_tot.SetFillStyle(3005)

    # signal
    if plot_sig:
      signal_hists = []
      for signal_file in self.signal_files:
        f_signal = self.tools.getRootFile(signal_file.filename, with_ext=False)
        hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        matching_selection = 'ismatched==1 && hnl_charge==0' if branchname == 'flat' else 'BToMuMuPi_isMatched==1 && BToMuMuPi_hnl_charge==0' # condition on charge added in the context of dsa study
        selection_signal = matching_selection if selection == '' else matching_selection + ' && ' + selection
        #TODO think of applying weight here
        hist_signal = self.tools.createHisto(f_signal, treename, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection=selection_signal)
        hist_signal.Sumw2()
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

    # qcd mc
    if plot_qcd:
      mc_hists = []
      hist_mc_tot = ROOT.TH1D('hist_mc_tot', 'hist_mc_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_mc_tot.Sumw2()
      int_mc_tot = 0.

      for ifile, qcd_file in enumerate(self.qcd_files):
        if qcd_file.label not in self.white_list: continue
        f_mc = self.tools.getRootFile(qcd_file.filename, with_ext=False)
        weight_mc = self.tools.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
        #TODO weight to be adapted here
        #weight = '({}) * (weight_hlt)'.format(weight_mc)
        weight = '({})'.format(weight_mc)
        hist_mc_name = 'hist_mc_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        hist_mc = self.tools.createHisto(f_mc, treename, self.quantity, hist_name=hist_mc_name, branchname=branchname, selection=selection, weight=weight) 
        hist_mc.Sumw2()
        hist_mc.SetFillColor(qcd_file.colour)
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
    
      hist_mc_tot.SetTitle('')
      hist_mc_tot.SetFillColor(ROOT.kAzure-4)
      hist_mc_tot.SetLineColor(1)
    
      if not do_stack:
        legend.AddEntry(hist_mc_tot, 'MC - {}'.format(self.getQCDMCLabel(self.white_list[0], self.white_list[len(self.white_list)-1])))
        
      ## create stack histogram  
      hist_mc_stack = ROOT.THStack('hist_mc_stack', '')

      ## compute the mc normalisation weight
      if do_luminorm: lumi_weight = self.tools.getLumiWeight(self.data_files, self.qcd_files, self.white_list, selection)

      for hist_mc in mc_hists:
        if do_shape and int_mc_tot != 0: hist_mc.Scale(1/int_mc_tot)
        elif do_luminorm: hist_mc.Scale(lumi_weight)
        hist_mc_stack.Add(hist_mc)
      if do_shape and int_mc_tot!= 0: hist_mc_tot.Scale(1/int_mc_tot)
      elif do_luminorm: hist_mc_tot.Scale(lumi_weight)


    if plot_qcd:
      frame = hist_mc_tot.Clone('frame')
    else:
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
    if plot_data and plot_qcd: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_data_stack, hist_mc_stack, do_log))
    elif plot_qcd and plot_sig: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_mc_stack, do_log, use_sig=True))
    elif plot_data and plot_sig: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_data_tot, do_log, use_sig=True))

    #ROOT.gStyle.SetPadLeftMargin(0.16) 
    ROOT.gStyle.SetOptStat(0)

    # draw the distributions
    frame.Draw()
    if plot_data and plot_qcd: hist_data_tot.Draw('same')
    if plot_data and not plot_qcd: hist_data_tot.Draw('histo same')
    if plot_qcd: 
      hist_mc_tot.Draw('histo same')
      if do_stack:
        hist_mc_stack.Draw('histo same')
      else:
        hist_mc_tot.Draw('histo same')
    if plot_data and plot_qcd: hist_data_tot.Draw('same') # making sure data points are always visible
    if plot_sig: 
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')

    # draw error bars
    if plot_qcd:
      hist_mc_tot_err = hist_mc_tot.Clone('hist_mc_tot_err')
      hist_mc_tot_err.SetLineWidth(0)
      hist_mc_tot_err.SetFillStyle(3244)
      hist_mc_tot_err.SetFillColor(ROOT.kGray+2)
      hist_mc_tot_err.Draw('E2 same')

    # draw the legend
    legend.Draw('same')

    # add labels
    self.tools.printLatexBox(0.65, 0.86, title, size=0.04 if plot_ratio else 0.036)
    if add_CMSlabel: self.tools.printCMSTag(pad_up, CMS_tag, size=0.55 if plot_ratio else 0.43)
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

    outputdir = self.tools.getOutDir('../outputs/{}/plots'.format(outdirlabel), plotdirlabel, do_shape, do_luminorm, do_stack, do_log)
    
    canv.SaveAs('{}/{}.png'.format(outputdir, self.quantity.label))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, self.quantity.label))



  def plotTwoSamples(file1, file2, branchname, tree1, tree2, selection1='', selection2='', legend1='legend1', legend2='legend2', do_printstat=False):
    f1 = ROOT.TFile.Open(file1, 'READ')
    f2 = ROOT.TFile.Open(file2, 'READ')

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
      signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

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
      signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

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
      signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

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


def printInfo(opt):
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '                   Running the plotter                          '
  print '\n'
  print ' data samples:    {}'.format(opt.data_label)
  print ' qcd samples:     {}'.format(opt.qcd_label)
  print ' signal samples:  {}'.format(opt.signal_label)
  print '\n'
  print ' categorisation:  {}'.format(opt.categories_label)
  print ' quantities:      {}'.format(opt.quantities_label)
  print ' selection:       {}'.format(opt.selection_label)
  print '\n'
  print ' outdir label:    {}'.format(opt.outdirlabel)
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '\n'

if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  opt = getOptions()
  used_parser = checkParser(opt)
  if used_parser:
    printInfo(opt)

    outdirlabel = opt.outdirlabel
    quantities = quantities[opt.quantities_label]
    categories = categories[opt.categories_label]
    data_files = data_samples[opt.data_label]
    qcd_files = qcd_samples[opt.qcd_label]
    signal_files = signal_samples[opt.signal_label]

    baseline_selection = selection[opt.selection_label].flat if opt.sample_type == 'flat' else selection[opt.selection_label].nano
    
    # to be modified
    white_list_15to300 = ['QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
    
    for category in categories:
      category_definition = category.definition_flat if opt.sample_type == 'flat' else category.definition_nano

      for quantity in quantities:
        plotter = Plotter(quantity=quantity, data_files=data_files, qcd_files=qcd_files, signal_files=signal_files, white_list=white_list_15to300)
        if opt.plot_CR:
          title = 'Control Region, {}'.format(category.title)
          plotdirlabel = 'CR/{}'.format(category.label)
          region_definition = 'hnl_charge != 0'
          plot_data = True
          plot_qcd = True
          plot_sig = False
          plot_ratio = True

          plotter.plot(selection = baseline_selection + ' && ' + region_definition + ' && ' + category_definition, 
                       title = title, 
                       outdirlabel = outdirlabel, 
                       plotdirlabel = plotdirlabel, 
                       branchname = opt.sample_type, 
                       treename = opt.tree_name,
                       plot_data = plot_data, 
                       plot_qcd = plot_qcd,
                       plot_sig = plot_sig, 
                       plot_ratio = plot_ratio, 
                       do_shape = opt.do_shape, 
                       do_luminorm = opt.do_luminorm, 
                       do_stack = opt.do_stack, 
                       do_log = opt.do_log,
                       add_overflow = opt.add_overflow,
                       add_CMSlabel = opt.add_CMSlabel,
                       CMS_tag = opt.CMStag
                       )


        if opt.plot_SR:
          title = 'Signal Region, {}'.format(category.title)
          plotdirlabel = 'SR/{}'.format(category.label)
          region_definition = 'hnl_charge == 0'
          plot_data = False
          plot_qcd = True
          plot_sig = True
          plot_ratio = False

          plotter.plot(selection = baseline_selection + ' && ' + region_definition + ' && ' + category_definition, 
                       title = title, 
                       outdirlabel = outdirlabel, 
                       plotdirlabel = plotdirlabel, 
                       branchname = opt.sample_type, 
                       treename = opt.tree_name,
                       plot_data = plot_data, 
                       plot_qcd = plot_qcd,
                       plot_sig = plot_sig, 
                       plot_ratio = plot_ratio, 
                       do_shape = opt.do_shape, 
                       do_luminorm = opt.do_luminorm, 
                       do_stack = opt.do_stack, 
                       do_log = opt.do_log,
                       add_overflow = opt.add_overflow,
                       add_CMSlabel = opt.add_CMSlabel,
                       CMS_tag = opt.CMStag
                       )
    
        if opt.plot_dataSig:
          title = 'Signal Region, {}'.format(category.title)
          plotdirlabel = 'DataSig/{}'.format(category.label)
          plot_data = True
          plot_qcd = False
          plot_sig = True
          plot_ratio = False

          plotter.plot(selection = baseline_selection + ' && ' + category_definition, 
                       title = title, 
                       outdirlabel = outdirlabel, 
                       plotdirlabel = plotdirlabel, 
                       branchname = opt.sample_type, 
                       treename = opt.tree_name,
                       plot_data = plot_data, 
                       plot_qcd = plot_qcd,
                       plot_sig = plot_sig, 
                       plot_ratio = plot_ratio, 
                       do_shape = opt.do_shape, 
                       do_luminorm = opt.do_luminorm, 
                       do_stack = opt.do_stack, 
                       do_log = opt.do_log,
                       add_overflow = opt.add_overflow,
                       add_CMSlabel = opt.add_CMSlabel,
                       CMS_tag = opt.CMStag
                       )
  else:

    doDataMCComparison = True
    doSignalBackgroundComparison = False
    compareTwoDistributions = False
    plotYields = False

    plot_log = False
    plot_categories = False
    plot_SR = False
    plot_CR = True
    # add add_overflow

    if doDataMCComparison:
      white_list_15to300 = ['QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
      #baseline_selection = 'mu_isdsa !=1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && '
      baseline_selection = ''
      if plot_CR:
        dirlabel = 'test'
        for quantity in quantities_to_plot_small:
          plotter = Plotter(quantity=quantity, data_files=data_samples, qcd_files=qcd_samples, signal_files=signal_samples, white_list=white_list_15to300)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0', title='Control Region, inclusive', outdirlabel=dirlabel+'/CR/incl', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)


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



