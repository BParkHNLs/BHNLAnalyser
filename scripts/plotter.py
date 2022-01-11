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
from compute_yields import ComputeYields

import sys
sys.path.append('../objects')
from quantity import Quantity, quantities
from samples import data_samples, qcd_samples, signal_samples
from categories import categories
from baseline_selection import selection
from qcd_white_list import white_list


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to produce the main analysis plots', add_help=True)
  parser.add_argument('--outdirlabel'     , type=str, dest='outdirlabel'     , help='name of the outdir'                                            , default=None)
  parser.add_argument('--subdirlabel'     , type=str, dest='subdirlabel'     , help='name of the subdir'                                            , default=None)
  parser.add_argument('--data_label'      , type=str, dest='data_label'      , help='which data samples to consider?'                               , default='V07_18Aug21')
  parser.add_argument('--qcd_label'       , type=str, dest='qcd_label'       , help='which qcd samples to consider?'                                , default='V07_18Aug21')
  parser.add_argument('--signal_label'    , type=str, dest='signal_label'    , help='which signal samples to consider?'                             , default='private')
  parser.add_argument('--quantities_label', type=str, dest='quantities_label', help='which quantity collection to consider?'                        , default='small')
  parser.add_argument('--selection_label' , type=str, dest='selection_label' , help='apply a baseline selection_label?'                             , default='standard')
  parser.add_argument('--categories_label', type=str, dest='categories_label', help='which phase space categorisation?'                             , default='standard')
  parser.add_argument('--category_label'  , type=str, dest='category_label'  , help='label of a given category within this list'                    , default=None)
  parser.add_argument('--sample_type'     , type=str, dest='sample_type'     , help='run the plotter on a nano or flat sample?'                     , default='flat', choices=['nano', 'flat'])
  parser.add_argument('--tree_name'       , type=str, dest='tree_name'       , help='name of the tree to analyse'                                   , default='signal_tree')
  parser.add_argument('--qcd_white_list ' , type=str, dest='qcd_white_list'  , help='pthat range to consider for qcd samples'                       , default='20to300')
  parser.add_argument('--CMStag '         , type=str, dest='CMStag'          , help='CMS tag to be added if --add_CMSlabel'                         , default='Preliminary')
  parser.add_argument('--weight_hlt'      , type=str, dest='weight_hlt'      , help='name of the hlt weight branch'                                 , default='weight_hlt_A1')
  parser.add_argument('--weight_puqcd'    , type=str, dest='weight_puqcd'    , help='name of the pu qcd weight branch'                                  , default='weight_pu_A')
  parser.add_argument('--add_weight_hlt'  ,           dest='add_weight_hlt'  , help='add hlt weight'                           , action='store_true', default=False)
  parser.add_argument('--add_weight_pu'   ,           dest='add_weight_pu'   , help='add pile-up weight'                       , action='store_true', default=False)
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


def printInfo(opt):
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '                   Running the plotter                          '
  print '\n'
  print ' data samples:    {}'.format(opt.data_label)
  print ' qcd samples:     {}'.format(opt.qcd_label)
  print ' signal samples:  {}'.format(opt.signal_label)
  print '\n'
  print ' categorisation:  {}'.format(opt.categories_label)
  print ' category:        {}'.format(opt.category_label)
  print ' quantities:      {}'.format(opt.quantities_label)
  print ' selection:       {}'.format(opt.selection_label)
  print '\n'
  print ' outdir label:    {}'.format(opt.outdirlabel)
  print ' subdir label:    {}'.format(opt.subdirlabel)
  print '-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
  print '\n'



class Plotter(Tools):
  def __init__(self, quantity='', data_files='', qcd_files='', signal_files='', white_list=''):
    self.tools = Tools()
    self.quantity = quantity
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_files = signal_files
    self.white_list = white_list


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


  def plot(self, selection='', title='', outdirloc='', outdirlabel='', subdirlabel='', plotdirlabel='', branchname='flat', treename='signal_tree', add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd='', plot_data=False, plot_qcd=False, plot_sig=False, plot_ratio=False, do_shape=True, do_luminorm=False, do_stack=True, do_log=False, add_overflow=False, add_CMSlabel=True, CMS_tag='Preliminary'):

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
        f_data = self.tools.getRootFile(data_file.filename)
        tree_data = self.getTree(f_data, treename)
        hist_data_name = 'hist_data_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        hist_data = self.tools.createHisto(tree_data, self.quantity, hist_name=hist_data_name, branchname=branchname, selection=selection)
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

      legend.AddEntry(hist_data_tot, 'data - {}'.format(self.getDataLabel(data_label, version_label) if len(self.data_files)>1 else data_file.label))

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
        f_sig = self.tools.getRootFile(signal_file.filename)
        tree_sig = self.getTree(f_sig, treename)

        hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        matching_selection = 'ismatched==1' if branchname == 'flat' else 'BToMuMuPi_isMatched==1'
        selection_signal = matching_selection if selection == '' else matching_selection + ' && ' + selection
        weight_sig = '(1)'
        if add_weight_hlt : weight_sig += ' * ({})'.format(weight_hlt)
        #if add_weight_pu : weight_sig += ' * (weight_pu_qcd)' #TODO modify pileup weight

        hist_signal = self.tools.createHisto(tree_sig, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection=selection_signal, weight=weight_sig)
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

        hist_signal.SetLineWidth(3)
        hist_signal.SetLineColor(signal_file.colour)
        signal_hists.append(hist_signal)

    # qcd mc
    if plot_qcd:
      qcd_hists = []
      hist_qcd_tot = ROOT.TH1D('hist_qcd_tot', 'hist_qcd_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_qcd_tot.Sumw2()
      int_qcd_tot = 0.

      #lumi_tot = self.tools.getQCDMCLumi(self.qcd_files, self.white_list)
      for ifile, qcd_file in enumerate(self.qcd_files):
        qcd_file_pthatrange = self.tools.getPthatRange(qcd_file.label)
        if qcd_file_pthatrange not in self.white_list: continue

        f_qcd = self.tools.getRootFile(qcd_file.filename)
        tree_qcd = self.getTree(f_qcd, treename)
        tree_run = self.getTree(f_qcd, 'run_tree')

        weight_qcd = self.tools.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
        #weight_qcd = 1. / (self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency) * lumi_tot)
        #weight_qcd = '(1)'
        weight_qcd = '({}*4.644)'.format(weight_qcd)
        #weight_qcd = '({})'.format(weight_qcd)
        print 'weight_qcd ',weight_qcd
        if add_weight_hlt : weight_qcd += ' * ({})'.format(weight_hlt)
        if add_weight_pu : weight_qcd += ' * ({}) '.format(weight_puqcd)

        hist_qcd_name = 'hist_qcd_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), do_log, do_shape)
        hist_qcd = self.tools.createHisto(tree_qcd, self.quantity, hist_name=hist_qcd_name, branchname=branchname, selection=selection, weight=weight_qcd) 
        hist_qcd.Sumw2()
        hist_qcd.SetFillColor(qcd_file.colour)
        hist_qcd.SetLineColor(1)
        
        if do_stack:
          legend.AddEntry(hist_qcd, 'MC - {}'.format(qcd_file.label))
        if add_overflow:
          overflow_qcd = hist_qcd.GetBinContent(hist_qcd.GetNbinsX()) + hist_qcd.GetBinContent(hist_qcd.GetNbinsX()+1)
          error_overflow_qcd = math.sqrt(math.pow(hist_qcd.GetBinError(hist_qcd.GetNbinsX()), 2) + math.pow(hist_qcd.GetBinError(hist_qcd.GetNbinsX()+1), 2)) 
          hist_qcd.SetBinContent(hist_qcd.GetNbinsX(), overflow_qcd)
          hist_qcd.SetBinError(hist_qcd.GetNbinsX(), error_overflow_qcd)
          hist_qcd.SetBinContent(hist_qcd.GetNbinsX()+1, 0)
          hist_qcd.SetBinError(hist_qcd.GetNbinsX()+1, 0)
        if do_shape: 
          int_qcd_tot += hist_qcd.Integral()
    
        hist_qcd_tot.Add(hist_qcd)
        qcd_hists.append(hist_qcd)
    
      hist_qcd_tot.SetTitle('')
      hist_qcd_tot.SetFillColor(ROOT.kAzure-4)
      hist_qcd_tot.SetLineColor(1)
    
      if not do_stack:
        legend.AddEntry(hist_qcd_tot, 'MC - {}'.format(self.getQCDMCLabel(self.white_list[0], self.white_list[len(self.white_list)-1], qcd_file.label)))
        
      ## create stack histogram  
      hist_qcd_stack = ROOT.THStack('hist_qcd_stack', '')

      ## compute the mc normalisation weight
      #if do_luminorm: lumi_weight = self.tools.getLumiWeight(self.data_files, self.qcd_files, self.white_list, selection, add_weight_hlt, add_weight_pu, weight_hlt, weight_puqcd)
      if do_luminorm:
        lumi_qcd = self.tools.getQCDMCLumi(self.qcd_files, self.white_list) #*112915./3884564. #9e-04 #0.018 #self.tools.getQCDMCLumi(self.qcd_files, self.white_list) 
        print 'lumi_qcd ',lumi_qcd
        lumi_data = self.tools.getDataLumi(self.data_files)
        print 'lumi_data ',lumi_data
        lumi_weight = 1#lumi_data #1#lumi_data / lumi_qcd
        print 'lumi_weight ',lumi_weight


      for hist_qcd in qcd_hists:
        if do_shape and int_qcd_tot != 0: hist_qcd.Scale(1/int_qcd_tot)
        elif do_luminorm: hist_qcd.Scale(lumi_weight)
        hist_qcd_stack.Add(hist_qcd)
      if do_shape and int_qcd_tot!= 0: hist_qcd_tot.Scale(1/int_qcd_tot)
      elif do_luminorm: hist_qcd_tot.Scale(lumi_weight)

    if plot_qcd:
      frame = hist_qcd_tot.Clone('frame')
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
    if plot_data and plot_qcd: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(hist_data_tot, hist_qcd_tot, do_log))
    elif plot_qcd and plot_sig: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_qcd_stack, do_log, use_sig=True))
    elif plot_data and plot_sig: frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_data_tot, do_log, use_sig=True))

    #ROOT.gStyle.SetPadLeftMargin(0.16) 
    ROOT.gStyle.SetOptStat(0)

    # draw the distributions
    frame.Draw()
    # tmp
    f1 = ROOT.TFile.Open('/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/data/pileup/profiles/pileup_data_2018.root', 'READ')
    hist_data_pu = f1.Get('pileup_2018A')
    hist_data_pu.Scale(1./hist_data_pu.Integral());
    hist_data_pu.SetMarkerStyle(41)
    hist_data_pu.SetMarkerSize(2)
    hist_data_pu.SetMarkerColor(ROOT.kBlue+2)
    #hist_data_pu.Draw('same')
    
    
    if plot_data and plot_qcd: hist_data_tot.Draw('same')
    if plot_data and not plot_qcd: hist_data_tot.Draw('histo same')
    if plot_qcd: 
      hist_qcd_tot.Draw('histo same')
      if do_stack:
        hist_qcd_stack.Draw('histo same')
      else:
        hist_qcd_tot.Draw('histo same')
    if plot_data and plot_qcd: hist_data_tot.Draw('same') # making sure data points are always visible
    if plot_sig: 
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')
    #hist_data_pu.Draw('same')

    # draw error bars
    if plot_qcd:
      hist_qcd_tot_err = hist_qcd_tot.Clone('hist_qcd_tot_err')
      hist_qcd_tot_err.SetLineWidth(0)
      hist_qcd_tot_err.SetFillStyle(3244)
      hist_qcd_tot_err.SetFillColor(ROOT.kGray+2)
      hist_qcd_tot_err.Draw('E2 same')

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

      hist_ratio = self.tools.getRatioHistogram(hist_data_tot, hist_qcd_tot)
      #hist_ratio = self.tools.getRatioHistogram(hist_data_tot, signal_hists[0])
      hist_ratio.Sumw2()

      for ibin in range(0, hist_ratio.GetNbinsX()+1):
        if hist_data_tot.GetBinContent(ibin) != 0 and hist_qcd_tot.GetBinContent(ibin) != 0:
          #err = hist_ratio.GetBinContent(ibin) * (math.sqrt(hist_data.GetBinContent(ibin))/hist_data.GetBinContent(ibin) + math.sqrt(hist_qcd_tot.GetBinContent(ibin))/hist_qcd_tot.GetBinContent(ibin))
          #err = math.sqrt((math.sqrt(hist_data.GetBinContent(ibin))/hist_qcd_tot.GetBinContent(ibin))**2 + (math.sqrt(hist_qcd_tot.GetBinContent(ibin))*hist_data.GetBinContent(ibin)/(hist_qcd_tot.GetBinContent(ibin))**2)**2)
          err = math.sqrt((hist_data_tot.GetBinError(ibin)/hist_qcd_tot.GetBinContent(ibin))**2 + (hist_qcd_tot.GetBinError(ibin)*hist_data_tot.GetBinContent(ibin)/(hist_qcd_tot.GetBinContent(ibin))**2)**2)
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

    outputdir = self.tools.getOutDir('{}/{}/plots/{}'.format(outdirloc, outdirlabel, subdirlabel), plotdirlabel, do_shape, do_luminorm, do_stack, do_log, add_overflow)
    
    canv.SaveAs('{}/{}.png'.format(outputdir, self.quantity.label))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, self.quantity.label))


  def plotTwoSamples(file1, file2, branchname, tree1, tree2, selection1='', selection2='', legend1='legend1', legend2='legend2', do_printstat=False):
    canv = self.tools.createTCanvas(name='canv', dimx=900, dimy=800)
    if self.do_log: canv.SetLogy()
    if not do_printstat: ROOT.gStyle.SetOptStat(0)

    hist1 = self.tools.createHisto(file1, tree1, self.quantity, hist_name=legend1, branchname=branchname, selection=selection1)
    hist2 = self.tools.createHisto(file2, tree2, self.quantity, hist_name=legend2, branchname=branchname, selection=selection2)
    
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



  def getBinnedPrefitMass(self, selection='', branchname='flat', treename='signal_tree', categories='', n_bins=50, add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd=''):

    # this will be to avoid memory issues with python
    signal_hists_all = []
    qcd_hists_all = []
    qcd_stack_hists_all = [] 
    frames_all = [] 

    # get signal histograms
    for category in categories:
      signal_hists = []
      for signal_file in self.signal_files:
        f_sig = self.tools.getRootFile(signal_file.filename)
        tree_sig = self.getTree(f_sig, treename)
        
        mass = signal_file.mass
        sigma = signal_file.resolution
    
        quantity_yields = Quantity(name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=n_bins, bin_min=mass-2*sigma, bin_max=mass+2*sigma) #TODO modify 

        hist_signal_name = 'hist_signal_{}_{}'.format(quantity_yields, outdirlabel.replace('/', '_'))
        matching_selection = 'ismatched==1' if branchname == 'flat' else 'BToMuMuPi_isMatched==1'
        selection_signal = matching_selection if selection == '' else matching_selection + ' && ' + selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
        weight_sig = '(1)'
        if add_weight_hlt : weight_sig += ' * ({})'.format(weight_hlt)
        #if add_weight_pu : weight_sig += ' * (weight_pu_qcd)' #TODO modify pileup weight

        hist_signal = self.tools.createHisto(tree_sig, quantity_yields, hist_name=hist_signal_name, branchname=branchname, selection=selection_signal, weight=weight_sig)
        hist_signal.Sumw2()
      
        # normalise signal to actual yields
        lumi_target = 41.6 #TODO
        sigma_B = 472.8e9 #TODO 
        #print 'signal selection ',selection_signal
        signal_yields = ComputeYields(signal_file=signal_file, selection=selection_signal).computeSignalYields(lumi=lumi_target, sigma_B=sigma_B, add_weight_hlt=False, weight_hlt='', isBc=False)[0]
        int_signal = hist_signal.Integral()
        hist_signal.Scale(signal_yields/int_signal)

        hist_signal.SetLineWidth(3)
        hist_signal.SetLineColor(signal_file.colour)
        signal_hists.append(hist_signal)

      signal_hists_all.append(signal_hists)

      # get qcd mc histograms
      qcd_hists = []
      hist_qcd_tot = ROOT.TH1D('hist_qcd_tot', 'hist_qcd_tot', quantity_yields.nbins, quantity_yields.bin_min, quantity_yields.bin_max)
      hist_qcd_tot.Sumw2()
      int_qcd_tot = 0.

      for ifile, qcd_file in enumerate(self.qcd_files):
        qcd_file_pthatrange = self.tools.getPthatRange(qcd_file.label)
        if qcd_file_pthatrange not in self.white_list: continue

        f_qcd = self.tools.getRootFile(qcd_file.filename)
        tree_qcd = self.getTree(f_qcd, treename)
        tree_run = self.getTree(f_qcd, 'run_tree')

        weight_qcd = self.tools.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
        weight_qcd = '({})'.format(weight_qcd)
        if add_weight_hlt : weight_qcd += ' * ({})'.format(weight_hlt)
        if add_weight_pu : weight_qcd += ' * ({}) '.format(weight_puqcd)

        hist_qcd_name = 'hist_qcd_{}_{}'.format(quantity_yields, outdirlabel.replace('/', '_'))
        selection_qcd = selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
        #print 'qcd selection ',selection_qcd
        hist_qcd = self.tools.createHisto(tree_qcd, quantity_yields, hist_name=hist_qcd_name, branchname=branchname, selection=selection_qcd, weight=weight_qcd) 
        hist_qcd.Sumw2()
        hist_qcd.SetFillColor(qcd_file.colour)
        hist_qcd.SetLineColor(1)
        
        int_qcd_tot += hist_qcd.Integral()
    
        hist_qcd_tot.Add(hist_qcd)
        qcd_hists.append(hist_qcd)
    
      hist_qcd_tot.SetTitle('')
      hist_qcd_tot.SetFillColor(ROOT.kAzure-4)
      hist_qcd_tot.SetLineColor(1)
    
      ## create stack histogram  
      hist_qcd_stack = ROOT.THStack('hist_qcd_stack', '')

      ## compute the mc normalisation weight
      #TODO use function
      lumi_true = 0
      for data_file in self.data_files:
        lumi_true += data_file.lumi
      print 'lumi true ',lumi_true

      from qcd_white_list import white_list
      from ABCD_regions import ABCD_regions
      background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=selection_qcd, white_list=white_list['20to300']).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=sigma, ABCD_regions=ABCD_regions['cos2d_svprob_0p996'], sigma_mult_window=2)[0] #TODO modify
      background_yields *= 41.6/lumi_true
      #print '{} bkg {}'.format(category.label, background_yields)

      # normalise to estimated yields
      for hist_qcd in qcd_hists:
        hist_qcd.Scale(background_yields/int_qcd_tot)
        hist_qcd_stack.Add(hist_qcd)
      hist_qcd_tot.Scale(background_yields/int_qcd_tot)

      qcd_hists_all.append(hist_qcd_tot)
      qcd_stack_hists_all.append(hist_qcd_stack)

      frame = ROOT.TH1D('frame', 'frame', n_bins, mass-2.5*sigma, mass+2.5*sigma) #TODO modify?
      frame.SetTitle('')
      #frame.GetYaxis().SetRangeUser(0, 10000)
      frame.GetXaxis().SetTitle(quantity.title)
      frame.GetXaxis().SetLabelSize(0.033)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetTitle('Yields')
      frame.GetYaxis().SetLabelSize(0.033)
      frame.GetYaxis().SetTitleSize(0.042)
      frame.GetYaxis().SetTitleOffset(1.3)
      frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_qcd_stack, do_log=True, use_sig=True))
      frames_all.append(frame)

    return signal_hists_all, qcd_hists_all, qcd_stack_hists_all


  def plotBinnedPrefitMass(self, selection='', title='', outdirloc='', outdirlabel='', subdirlabel='', plotdirlabel='', branchname='flat', treename='signal_tree', n_bins=25, add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd='', do_stack=True, do_log=False, add_CMSlabel=True, CMS_tag='Preliminary'):

    #ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
    #ROOT.gROOT.ProcessLine('.L /t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNLAnalyser/scripts/example_style.C')
    #ROOT.gROOT.ProcessLine('set_tdr_style()')
    #ROOT.gROOT.ProcessLine('set_root_style()')
    #ROOT.gROOT.ProcessLine('set_mult_pad_style()')
    #ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
    #ROOT.gROOT.ProcessLine('setTDRStyle()')

    signal_hists_all, qcd_hists_all, qcd_stack_hists_all = self.getBinnedPrefitMass(selection, branchname, treename, categories, n_bins, add_weight_hlt, add_weight_pu, weight_hlt, weight_puqcd) 

    ## plot each category separately 
    for icat, category in enumerate(categories):
      # create the canvas
      canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, outdirlabel.replace('/', '_'), do_log, category.label)
      canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
      canv.cd()

      # define the pad
      pad_up = ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
      if do_log: pad_up.SetLogy()
      pad_up.Draw()
      canv.cd()

      # prepare the legend
      if do_stack:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.45, xmax=0.84, ymax=0.83, size=0.027)
      else:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.65, xmax=0.84, ymax=0.83, size=0.027)

      # get histograms
      signal_hists = signal_hists_all[icat]
      hist_qcd_tot = qcd_hists_all[icat]
      hist_qcd_stack = qcd_hists_all[icat]

      #for hist_signal in signal_hists:
      #  legend.AddEntry(hist_signal, 'signal - {}'.format(signal_file.label))
      #legend.AddEntry(hist_qcd_tot, 'MC - {}'.format(self.getQCDMCLabel(self.white_list[0], self.white_list[len(self.white_list)-1], qcd_file.label)))

      pad_up.cd()

      #ROOT.gStyle.SetPadLeftMargin(0.16) 
      ROOT.gStyle.SetOptStat(0)

      mass = 1
      sigma = 0.00853 
      frame = ROOT.TH1D('frame', 'frame', 50, mass-2.5*sigma, mass+2.5*sigma) #TODO modify?
      frame.SetTitle('')
      #frame.GetYaxis().SetRangeUser(0, 10000)
      frame.GetXaxis().SetTitle(quantity.title)
      frame.GetXaxis().SetLabelSize(0.033)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetTitle('Yields')
      frame.GetYaxis().SetLabelSize(0.033)
      frame.GetYaxis().SetTitleSize(0.042)
      frame.GetYaxis().SetTitleOffset(1.3)
      frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_qcd_stack, do_log=True, use_sig=True))
      #frame.GetYaxis().SetRangeUser(1e-9, 1e6)

      # draw the distributions
      frame.Draw()
      hist_qcd_tot.Draw('histo same')
      if do_stack:
        hist_qcd_stack.Draw('histo same')
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')

      # draw error bars
      hist_qcd_tot_err = hist_qcd_tot.Clone('hist_qcd_tot_err')
      hist_qcd_tot_err.SetLineWidth(0)
      hist_qcd_tot_err.SetFillStyle(3244)
      hist_qcd_tot_err.SetFillColor(ROOT.kGray+2)
      hist_qcd_tot_err.Draw('E2 same')

      # draw the legend
      legend.Draw('same')

      # add labels
      self.tools.printLatexBox(0.65, 0.86, category.title, size=0.036)
      if add_CMSlabel: self.tools.printCMSTag(pad_up, CMS_tag, size=0.43)

      plotname = 'hnl_mass_{}'.format(category.label)
      outputdir = self.tools.getOutDir('{}/{}/datacards/{}'.format(outdirloc, outdirlabel, subdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
      
      canv.SaveAs('{}/{}.png'.format(outputdir, plotname))
      canv.SaveAs('{}/{}.pdf'.format(outputdir, plotname))


    # plot categories on same canvas 
    # create the canvas
    canv_name = 'canv_mult_{}_{}_{}'.format(self.quantity.label, outdirlabel.replace('/', '_'), do_log)
    canv_mult = self.tools.createTCanvas(name=canv_name, dimx=1500, dimy=1000)
    #canv_mult.Divide((len(categories)-1)/2, 2)
    canv_mult.Divide(1, 2)
    
    for icat, category in enumerate(categories):
      print category.label
      if 'incl' in category.label: continue

      pad_min = 0.
      pad_max = 1.
      pad_size = float((pad_max-pad_min)/(len(categories)-1)*2)
      if icat <= len(categories)/2:
        canv_mult.cd(1)
        pad_up = ROOT.TPad("pad_up","pad_up",pad_min+(icat-1)*pad_size,0,pad_min+icat*pad_size,1)
        pad_up.SetLeftMargin(0.15)
        pad_up.SetRightMargin(0.09)
        pad_up.SetTopMargin(0.15)
        pad_up.SetBottomMargin(0.03)
      else: 
        canv_mult.cd(2)
        pad_up = ROOT.TPad("pad_up","pad_up",pad_min+(icat-1-(len(categories)-1)/2)*pad_size,0,pad_min+(icat-(len(categories)-1)/2)*pad_size,1)
        pad_up.SetLeftMargin(0.15)
        pad_up.SetRightMargin(0.09)
        pad_up.SetTopMargin(0.03)
        pad_up.SetBottomMargin(0.15)
  
      #pad_up = ROOT.TPad("pad_up","pad_up",pad_min+icat*pad_size,0,pad_min+(icat+1)*pad_size,1)
      pad_up.SetLogy()
      pad_up.Draw()
      pad_up.cd()

      # prepare the legend
      if do_stack:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.45, xmax=0.84, ymax=0.83, size=0.027)
      else:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.65, xmax=0.84, ymax=0.83, size=0.027)

      # get histograms
      signal_hists = signal_hists_all[icat]
      hist_qcd_tot = qcd_hists_all[icat]
      hist_qcd_stack = qcd_hists_all[icat]

      #pad_up.cd()

      #ROOT.gStyle.SetPadLeftMargin(0.16) 
      ROOT.gStyle.SetOptStat(0)

      #mass = 1
      #sigma = 0.00853 
      #frame = ROOT.TH1D('frame'+category.label, 'frame'+category.label, 30, mass-2.5*sigma, mass+2.5*sigma) #TODO modify?
      #frame.SetTitle('')
      ##frame.GetYaxis().SetRangeUser(0, 10000)
      #frame.GetXaxis().SetTitle(quantity.title)
      #frame.GetXaxis().SetLabelSize(0.033)
      #frame.GetXaxis().SetTitleSize(0.042)
      #frame.GetXaxis().SetTitleOffset(1.1)
      #frame.GetYaxis().SetTitle('Yields')
      #frame.GetYaxis().SetLabelSize(0.033)
      #frame.GetYaxis().SetTitleSize(0.042)
      #frame.GetYaxis().SetTitleOffset(1.3)
      ##frame.GetYaxis().SetRangeUser(1e-9, self.getMaxRangeY(signal_hists, hist_qcd_stack, do_log=True, use_sig=True))
      #frame.GetYaxis().SetRangeUser(1e-9, 1e6)

      hist_qcd_tot.SetTitle('')
      hist_qcd_tot.GetXaxis().SetTitle(quantity.title if icat==len(categories)-1 else '')
      hist_qcd_tot.GetXaxis().SetLabelSize(0.06 if icat>(len(categories)-1)/2 else 0)
      hist_qcd_tot.GetXaxis().SetTitleSize(0.06 if icat>(len(categories)-1)/2 else 0)
      hist_qcd_tot.GetXaxis().SetTickLength(0.03)
      hist_qcd_tot.GetXaxis().SetNdivisions(405)
      hist_qcd_tot.GetXaxis().SetTitleOffset(1.1)
      hist_qcd_tot.GetYaxis().SetTitle('Yields' if icat==1 else '')
      hist_qcd_tot.GetYaxis().SetLabelSize(0.06 if (icat==1 or icat-1==(len(categories)-1)/2) else 0)
      hist_qcd_tot.GetYaxis().SetTitleSize(0.06 if (icat==1 or icat-1==(len(categories)-1)/2) else 0)
      ##hist_qcd_tot.GetYaxis().SetLabelSize(0.05)
      ##hist_qcd_tot.GetYaxis().SetTitleSize(0.05)
      hist_qcd_tot.GetYaxis().SetTickLength(0.03)
      hist_qcd_tot.GetYaxis().SetNdivisions(405)
      hist_qcd_tot.GetYaxis().SetTitleOffset(1.3)
      hist_qcd_tot.GetYaxis().SetRangeUser(1e-9, 1e6)

      # draw the distributions
      #frame.Draw()
      #hist_qcd_tot.Draw('histo same')
      hist_qcd_tot.Draw('histo')
      if do_stack:
        hist_qcd_stack.Draw('histo same')
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')

      # draw error bars
      hist_qcd_tot_err = hist_qcd_tot.Clone('hist_qcd_tot_err')
      hist_qcd_tot_err.SetLineWidth(0)
      hist_qcd_tot_err.SetFillStyle(3244)
      hist_qcd_tot_err.SetFillColor(ROOT.kGray+2)
      hist_qcd_tot_err.Draw('E2 same')

      # draw the legend
      legend.Draw('same')

      # add labels
      self.tools.printLatexBox(0.65, 0.86, category.title, size=0.036)
      if add_CMSlabel and icat==1: self.tools.printCMSTag(pad_up, CMS_tag, size=0.6)

    plotname = 'hnl_mass_mult'
    outputdir = self.tools.getOutDir('{}/{}/datacards/{}'.format(outdirloc, outdirlabel, subdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
    
    canv_mult.cd()
    canv_mult.SaveAs('{}/{}.png'.format(outputdir, plotname))
    canv_mult.SaveAs('{}/{}.pdf'.format(outputdir, plotname))




if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  opt = getOptions()
  used_parser = checkParser(opt)
  if used_parser:
    '''
      Automatic handling 
    '''
    printInfo(opt)

    outdirlabel = opt.outdirlabel
    subdirlabel = opt.subdirlabel
    quantities = quantities[opt.quantities_label]
    categories = categories[opt.categories_label]
    data_files = data_samples[opt.data_label]
    qcd_files = qcd_samples[opt.qcd_label]
    signal_files = signal_samples[opt.signal_label]

    baseline_selection = selection[opt.selection_label].flat if opt.sample_type == 'flat' else selection[opt.selection_label].nano
    
    white_list = white_list[opt.qcd_white_list]
    
    for category in categories:
      if opt.category_label != None and category.label != opt.category_label: continue # needed for category parallelisation on the batch
      category_definition = category.definition_flat if opt.sample_type == 'flat' else category.definition_nano
      category_cutbased_selection = category.cutbased_selection

      for quantity in quantities:
        plotter = Plotter(quantity=quantity, data_files=data_files, qcd_files=qcd_files, signal_files=signal_files, white_list=white_list)
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
                       outdirloc = './outputs',
                       outdirlabel = outdirlabel, 
                       subdirlabel = subdirlabel, 
                       plotdirlabel = plotdirlabel, 
                       branchname = opt.sample_type, 
                       treename = opt.tree_name,
                       add_weight_hlt = 1 if opt.add_weight_hlt else 0,
                       add_weight_pu = 1 if opt.add_weight_pu else 0,
                       weight_hlt = opt.weight_hlt,
                       weight_puqcd = opt.weight_puqcd,
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
          do_shape = False #True
          do_luminorm = True #False

          ##plotter.plot(selection = baseline_selection + ' && ' + region_definition + ' && ' + category_definition, 
          #             title = title, 
          #             outdirloc = './outputs',
          #             outdirlabel = outdirlabel, 
          #             subdirlabel = subdirlabel, 
          #             plotdirlabel = plotdirlabel, 
          #             branchname = opt.sample_type, 
          #             treename = opt.tree_name,
          #             add_weight_hlt = 1 if opt.add_weight_hlt else 0,
          #             add_weight_pu = 1 if opt.add_weight_pu else 0,
          #             weight_hlt = opt.weight_hlt,
          #             weight_puqcd = opt.weight_puqcd,
          #             plot_data = plot_data, 
          #             plot_qcd = plot_qcd,
          #             plot_sig = plot_sig, 
          #             plot_ratio = plot_ratio, 
          #             do_shape = do_shape, 
          #             do_luminorm = do_luminorm, 
          #             do_stack = opt.do_stack, 
          #             do_log = opt.do_log,
          #             add_overflow = opt.add_overflow,
          #             add_CMSlabel = opt.add_CMSlabel,
          #             CMS_tag = opt.CMStag
          #             )

          plotter.plotBinnedPrefitMass(selection = baseline_selection + ' && ' + region_definition + ' && ' + category_definition + ' && ' + category_cutbased_selection, 
                       title = title, 
                       outdirloc = './outputs',
                       outdirlabel = outdirlabel, 
                       subdirlabel = subdirlabel, 
                       plotdirlabel = plotdirlabel, 
                       branchname = opt.sample_type, 
                       treename = opt.tree_name,
                       n_bins = 25,
                       add_weight_hlt = 1 if opt.add_weight_hlt else 0,
                       add_weight_pu = 1 if opt.add_weight_pu else 0,
                       weight_hlt = opt.weight_hlt,
                       weight_puqcd = opt.weight_puqcd,
                       do_stack = opt.do_stack, 
                       #do_log = opt.do_log,
                       do_log = True, #opt.do_log,
                       add_CMSlabel = opt.add_CMSlabel,
                       CMS_tag = opt.CMStag
                       )
    
        if opt.plot_dataSig:
          title = 'Signal Region, {}'.format(category.title)
          plotdirlabel = 'DataSig/{}'.format(category.label)
          plot_data = True
          plot_qcd = False
          plot_sig = True
          plot_ratio = opt.plot_ratio

          plotter.plot(selection = baseline_selection + ' && ' + category_definition, 
                       title = title, 
                       outdirloc = './outputs',
                       outdirlabel = outdirlabel, 
                       subdirlabel = subdirlabel, 
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
    '''
       ## Interactive board ##

    For an interactive usage of the script, use the code below
    '''

    doDataMCComparison = True
    doSignalBackgroundComparison = False
    compareTwoDistributions = False

    plot_log = False
    plot_categories = False
    plot_SR = False
    plot_CR = True
    # add add_overflow

    if doDataMCComparison:
      white_list_15to300 = white_list['15to300']
      #white_list_15to300 = ['QCD_pt20to30 (V07_18Aug21)', 'QCD_pt30to50 (V07_18Aug21)', 'QCD_pt50to80 (V07_18Aug21)', 'QCD_pt80to120 (V07_18Aug21)', 'QCD_pt120to170 (V07_18Aug21)', 'QCD_pt170to300 (V07_18Aug21)']
      #baseline_selection = 'mu_isdsa !=1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && '
      baseline_selection = ''
      if plot_CR:
        dirlabel = 'test'
        for quantity in quantities_to_plot_small:
          plotter = Plotter(quantity=quantity, data_files=data_samples, qcd_files=qcd_samples, signal_files=signal_samples, white_list=white_list_15to300)
          plotter.plotDataMCComparison(selection=baseline_selection+'hnl_charge!=0', title='Control Region, inclusive', outdirlabel=dirlabel+'/CR/incl', branchname='flat', plot_data=True, plot_sig=False, plot_ratio=True, do_shape=False, do_luminorm=True, do_stack=False, do_log=False)


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



