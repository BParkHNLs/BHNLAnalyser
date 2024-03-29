import ROOT
import numpy as np
import os
from os import path
from tools import Tools
import sys
import math
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples, qcd_samples, data_samples
from categories import Category
from baseline_selection import selection
from qcd_white_list import white_list
from quantity import Quantity as Qte


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to study the selection', add_help=True)
  parser.add_argument('--dirlabel'             , type=str, dest='dirlabel'            , help='name of the outdir'  , default=None)
  parser.add_argument('--mass'                 , type=str, dest='mass'                , help='mass'                , default=None)
  parser.add_argument('--category'             , type=str, dest='category'            , help='category'            , default=None)
  parser.add_argument('--quantity'             , type=str, dest='quantity'            , help='quantity'            , default=None)
  parser.add_argument('--action'               , type=str, dest='action'              , help='action'              , default=None)
  parser.add_argument('--cut_b_mass'           , type=str, dest='cut_b_mass'          , help='cut_b_mass'          , default=None)
  parser.add_argument('--cut_pi_pt'            , type=str, dest='cut_pi_pt'           , help='cut_pi_pt'           , default=None)
  parser.add_argument('--cut_sv_lxysig'        , type=str, dest='cut_sv_lxysig'       , help='cut_sv_lxysig'       , default=None)
  parser.add_argument('--cut_mu_dxysig'        , type=str, dest='cut_mu_dxysig'       , help='cut_mu_dxysig'       , default=None)
  parser.add_argument('--cut_pi_dxysig'        , type=str, dest='cut_pi_dxysig'       , help='cut_pi_dxysig'       , default=None)
  parser.add_argument('--cut_max_mu_pi_dxysig' , type=str, dest='cut_max_mu_pi_dxysig', help='cut_max_mu_pi_dxysig', default=None)
  parser.add_argument('--cut_min_mu_pi_dxysig' , type=str, dest='cut_min_mu_pi_dxysig', help='cut_min_mu_pi_dxysig', default=None)
  parser.add_argument('--cut_hnl_cos2d'        , type=str, dest='cut_hnl_cos2d'       , help='cut_hnl_cos2d'       , default=None)
  return parser.parse_args()


class Quantity(object):
  def __init__(self, name_flat='', label='', title='', logic='', units='', binMin=0., binMax=0.):
    self.name_flat = name_flat
    self.label = label
    self.title = title
    self.logic = logic
    self.units = units
    self.binMin = binMin
    self.binMax = binMax


class SelectedQuantity(object):
  '''
  this class will allow to impose a selection when studying 
  the impact of the cut on a given parameter
  '''
  def __init__(self, quantity, chosen_cut):
    self.quantity   = quantity
    self.chosen_cut = chosen_cut

#TODO apply weights on signal, move dxysig definition to BS (DCASig for pions), add selection on the b mass difference with mB

class Selection(object):
  def __init__(self, signal_label='', signal_files='', data_file='', qcd_files='', baseline_selection='', quantity='', category=None, use_data=True, dirlabel='', preexisting_selection=None, proposed_cut=None, white_list=''):
    self.tools                 = Tools()
    self.signal_label          = signal_label
    self.signal_files          = signal_files
    self.data_file             = data_file
    self.qcd_files             = qcd_files
    self.baseline_selection    = baseline_selection
    self.quantity              = quantity
    self.category              = category
    self.use_data              = use_data
    self.dirlabel              = dirlabel
    self.preexisting_selection = preexisting_selection
    self.proposed_cut          = proposed_cut
    self.white_list            = white_list

    #self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2' 
    self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2'
    print 'weight_hlt: {}'.format(self.weight_hlt)
    self.weight_pusig = 'weight_pu_sig_D'
    self.weight_mu0id = 'weight_mu0_softid'
    self.weight_muid = 'weight_mu_looseid'

    self.resolution_p0 = 0.0002747
    self.resolution_p1 = 0.008302


  def createOutDir(self, outputdir):
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))


  def getLabel(self, str_):
    new_str_ = str_
    if 'abs(' in str_: new_str_ = new_str_.replace('abs(', '')
    if ')'in str_: new_str_ = new_str_.replace(')', '')
    return new_str_


  def getSelectionString(self):
    '''
    function to write the already-applied cuts in a string
    '''
    preselection_str = []
    for item, _ in enumerate(self.preexisting_selection):
      name_variable = self.preexisting_selection[item].quantity.name_flat
      preselection_str.append('{}{}{}'.format(name_variable,self.preexisting_selection[item].quantity.logic,self.preexisting_selection[item].chosen_cut))
    return ' && '.join(preselection_str)


  def getEff(self, entries_selected, entries_initial):
    eff = entries_selected / entries_initial if entries_initial!=0 else 0.
    return eff


  def getScanGraph(self, npoints=20):

    '''
    plot the signal efficiency and background rejection as a function of the cut applied on the quantity of interest
    possibility to draw a line at the chosen cut
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, npoints) 

    canv = self.tools.createTCanvas('canv', 900, 800)
    canv.SetGrid()
    
    gs_sig = []
    gs_bkg = []

    for signal_file in self.signal_files:
      print '\nsignal_file ',signal_file.ctau
      # define the mass window
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)
      print '{} {}'.format(signal_file.mass-2*signal_file.resolution, signal_file.mass+2*signal_file.resolution)

      # signal efficiency
      g_sig = ROOT.TGraph()

      selection_sig_ini = 'ismatched == 1 && {} &&  {}'.format(self.baseline_selection, self.category.definition_flat) 
      if self.preexisting_selection !=None: selection_sig_ini += ' && {}'.format(self.getSelectionString())
      #print 'selection_sig_ini ',selection_sig_ini
      f_sig = self.tools.getRootFile(signal_file.filename)
      tree_sig = self.tools.getTree(f_sig, 'signal_tree')
      initial_sig_entries = self.tools.createHisto(tree_sig, self.hnl_mass, selection=selection_sig_ini).Integral()
      #print 'initial_sig_entries {} '.format(initial_sig_entries)

      for idx, cut in enumerate(points):
        selection_sig_sel = selection_sig_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        selected_sig_entries = self.tools.createHisto(tree_sig, self.hnl_mass, selection=selection_sig_sel).Integral()
        #print 'selection_sig_sel {} {}'.format(selection_sig_sel, selected_sig_entries)
        g_sig.SetPoint(idx, cut, self.getEff(selected_sig_entries, initial_sig_entries))
        g_sig.SetLineWidth(2)
        g_sig.SetLineStyle(9)
        g_sig.SetLineColor(signal_file.colour)
        g_sig.SetMarkerColor(signal_file.colour)
        g_sig.SetMarkerStyle(20)
        print 'efficiency signal selection ',selection_sig_sel
        print 'efficiency cut {} ini signal {} selected signal {}'.format(cut, initial_sig_entries, selected_sig_entries)

      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetRangeUser(0, 1)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('APL')
        else:
          g_sig.Draw('PL')

      # background rejection
      g_bkg = ROOT.TGraph()
      selection_bkg_ini = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
      if self.preexisting_selection !=None: selection_bkg_ini += ' && {}'.format(self.getSelectionString())
      #print 'selection_bkg_ini ',selection_bkg_ini
      if not self.use_data:
        initial_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_ini).Integral() # no weight for the moment
      else:
        f_data = self.tools.getRootFile(self.data_file)
        tree_data = self.tools.getTree(f_data, 'signal_tree')
        hist_data_name = 'hist_data_{}'.format(self.quantity)
        hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=selection_bkg_ini)
        initial_bkg_entries = hist_data.Integral()
      #print 'initial_bkg_entries ',initial_bkg_entries

      for idx, cut in enumerate(points):
        selection_bkg_sel = selection_bkg_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        #print 'selection_bkg_sel ',selection_bkg_sel
        if not self.use_data:
          selected_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_sel).Integral() # no weight for the moment
        else:
          f_data = self.tools.getRootFile(self.data_file)
          tree_data = self.tools.getTree(f_data, 'signal_tree')
          hist_data_name = 'hist_data_sel_{}'.format(self.quantity)
          hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=selection_bkg_sel)
          selected_bkg_entries = hist_data.Integral()
        g_bkg.SetPoint(idx, cut, 1-self.getEff(selected_bkg_entries, initial_bkg_entries))
        g_bkg.SetLineWidth(2)
        g_bkg.SetLineColor(signal_file.colour)
        g_bkg.SetMarkerColor(signal_file.colour)
        g_bkg.SetMarkerSize(2)
        g_bkg.SetMarkerStyle(45)
        print 'efficiency background selection ',selection_bkg_sel
        print 'efficiency cut {} ini background {} selected background {}'.format(cut, initial_bkg_entries, selected_bkg_entries)

      gs_bkg.append(g_bkg)

      for g_bkg in gs_bkg:
        g_bkg.Draw('PL')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, 0, self.proposed_cut, 1)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      line.Draw('same')
    
    self.tools.printLatexBox(0.73, 0.925, self.category.title, size=0.041)
          
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    for ifile, signal_file in enumerate(self.signal_files):
      legend.AddEntry(gs_sig[ifile], 'sig efficiency ({}GeV, {}mm)'.format(signal_file.mass, signal_file.ctau))
      legend.AddEntry(gs_bkg[ifile], 'bkg rejection around {}GeV'.format(signal_file.mass))
    legend.Draw()

    canv.cd()
    self.createOutDir('myPlots/selection/{}'.format(self.dirlabel))
    canv.SaveAs('myPlots/selection/{}/scan_{}_{}.png'.format(self.dirlabel, self.getLabel(self.quantity.label), self.category.label))


  def getSignificanceScan(self, npoints=20):
    '''
    plot the significance as a function of the cut applied on the quantity of interest
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, npoints) 
      
    canv = self.tools.createTCanvas('canv', 900, 800)
    canv.SetGrid()

    canv_log = self.tools.createTCanvas('canv_log', 900, 800)
    canv_log.SetGrid()
    canv_log.SetLogy()
    
    gs_sig = []

    significance_max = -1
    significance_min = -1

    for signal_file in self.signal_files:
      print '\nsignal_file ',signal_file.ctau
      g_sig = ROOT.TGraph()
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)

      for idx, cut in enumerate(points):
        # compute signal yields
        signal_selection = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
        if self.preexisting_selection !=None: signal_selection += ' && {}'.format(self.getSelectionString())
        signal_selection += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)

        lumi_D1 = 5.302 # restrict to this lumi, as background taken from unblinded D1 data
        #signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi_D1, sigma_B=472.8e9) 
        signal_yields = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mui0d=self.weight_mu0id, weight_muid=self.weight_muid) 
        #print '{} {} {}'.format(signal_file.mass, signal_file.ctau, int(signal_yields))

        # compute background yields
        #background_selection = self.category.definition_flat
        #if self.preexisting_selection !=None: background_selection += ' && {}'.format(self.getSelectionString())
        #background_selection += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)

        #data_files = data_samples['V09_06Nov21']
        #qcd_files = qcd_samples['V09_06Nov21']
        #from ABCD_regions import ABCD_regions 
        #ABCD_regions = ABCD_regions['cos2d_svprob_0p996']

        #background_yields = ComputeYields(data_files=data_files, qcd_files=qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromABCDData(mass=signal_file.mass, resolution=signal_file.resolution, ABCD_regions=ABCD_regions)[0]

        #lumi_true = 0.
        #for data_file in data_samples['V09_06Nov21']:
        #  lumi_true += data_file.lumi

        #background_yields = background_yields * 41.6/lumi_true
        # get background yields from unblinded data
        f_data = self.tools.getRootFile(self.data_file)
        tree_data = self.tools.getTree(f_data, 'signal_tree')
        background_selection = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
        if self.preexisting_selection !=None: background_selection += ' && {}'.format(self.getSelectionString())
        background_selection += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        hist_data_name = 'hist_data_{}'.format(self.quantity)
        hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=background_selection)
        background_yields = hist_data.Integral()

        # compute significance
        if background_yields != 0:
          significance = signal_yields / math.sqrt(background_yields)
          print 'significance signal selection ',signal_selection
          print 'significance background selection ',background_selection
          print 'significance cut {} signal yields {} background yields {} significance {}'.format(cut, signal_yields, background_yields, significance)
          #print '{} {}'.format(signal_file.ctau, significance)
        else:
          significance = -99.

        if significance_max == -1 or significance > significance_max:
          significance_max = significance

        if significance_min == -1 or significance < significance_min:
          significance_min = significance

        g_sig.SetPoint(idx, cut, significance)
        g_sig.SetLineWidth(2)
        g_sig.SetLineStyle(9)
        g_sig.SetLineColor(signal_file.colour)
        g_sig.SetMarkerColor(signal_file.colour)
        g_sig.SetMarkerStyle(20)

      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetRangeUser(significance_min-0.15*significance_min, significance_max+0.15*significance_max)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          canv.cd()
          g_sig.Draw('APL')
          canv_log.cd()
          g_sig.Draw('APL')
        else:
          canv.cd()
          g_sig.Draw('PL')
          canv_log.cd()
          g_sig.Draw('PL')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, significance_min-0.15*significance_min, self.proposed_cut, significance_max+0.15*significance_max)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      canv.cd()
      line.Draw('same')
      canv_log.cd()
      line.Draw('same')
    
    self.tools.printLatexBox(0.73, 0.925, self.category.title, size=0.041)
          
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    for ifile, signal_file in enumerate(self.signal_files):
      legend.AddEntry(gs_sig[ifile], 'significance ({}GeV, {}mm)'.format(signal_file.mass, signal_file.ctau))
    canv.cd()
    legend.Draw()
    canv_log.cd()
    legend.Draw()

    self.createOutDir('myPlots/selection/{}'.format(self.dirlabel))

    canv.cd()
    canv.SaveAs('myPlots/selection/{}/significance_scan_{}_{}.png'.format(self.dirlabel, self.getLabel(self.quantity.label), self.category.label))

    canv_log.cd()
    canv_log.SaveAs('myPlots/selection/{}/significance_scan_{}_{}_log.png'.format(self.dirlabel, self.getLabel(self.quantity.label), self.category.label))


  def getSignificanceDiffScan(self, npoints=20):
    '''
    plot the significance gain/loss as a function of the cut applied on the quantity of interest
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, npoints) 
      
    canv = self.tools.createTCanvas('canv', 900, 800)
    #canv.SetGrid()

    #canv_log = self.tools.createTCanvas('canv_log', 900, 800)
    #canv_log.SetGrid()
    #canv_log.SetLogy()

    ROOT.gStyle.SetPadLeftMargin(0.8)
    ROOT.gStyle.SetPadRightMargin(0.8)
    
    gs_sig = []

    significance_diff_max = -1
    significance_diff_min = -1

    for signal_file in self.signal_files:
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau

      g_sig = ROOT.TGraph()
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)
      #print '{} {} {} {}'.format(signal_mass, signal_ctau, signal_file.mass-2*signal_file.resolution, signal_file.mass+2*signal_file.resolution)

      lumi_D1 = 5.302 # when computing the relative difference in significance, this value does not matter

      # compute the initial significance
      signal_selection_ini = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
      if self.preexisting_selection !=None: signal_selection_ini += ' && {}'.format(self.getSelectionString())
      #signal_yields_ini, err_signal_yields_ini = ComputeYields(signal_file=signal_file, selection=signal_selection_ini).computeSignalYields(lumi=lumi_D1, sigma_B=472.8e9) 
      signal_yields_ini = ComputeYields(signal_label=signal_label, selection=signal_selection_ini).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid)[0] 
      #print 'signal_yields_ini ',signal_yields_ini

      f_data = self.tools.getRootFile(self.data_file)
      tree_data = self.tools.getTree(f_data, 'signal_tree')
      background_selection_ini = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
      if self.preexisting_selection !=None: background_selection_ini += ' && {}'.format(self.getSelectionString())
      hist_data_name = 'hist_data_ini_{}'.format(self.quantity)
      hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=background_selection_ini)
      background_yields_ini = hist_data.Integral()
      #print 'background_yields_ini ',background_yields_ini

      if background_yields_ini != 0:
        #if 2*((signal_yields_ini + background_yields_ini)*math.log(1 + signal_yields_ini/background_yields_ini) - signal_yields_ini) > 0:
        #  significance_ini = math.sqrt(2*((signal_yields_ini + background_yields_ini)*math.log(1 + signal_yields_ini/background_yields_ini) - signal_yields_ini))
        #else:
        #  significance_ini = signal_yields_ini / math.sqrt(background_yields_ini)
        significance_ini = signal_yields_ini / math.sqrt(background_yields_ini)
        #print 'ini {} {} {} {}'.format(signal_yields_ini, background_yields_ini, significance_ini_1, significance_ini)
      else:
        significance_ini = -99.
      #print 'ini {} / sqrt({}) = {}'.format(signal_yields_ini, background_yields_ini, significance_ini)

      for idx, cut in enumerate(points):
        print cut
        # compute signal yields
        signal_selection_sel = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
        if self.preexisting_selection !=None: signal_selection_sel += ' && {}'.format(self.getSelectionString())
        signal_selection_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        #print signal_selection_sel

        #signal_yields_sel, err_signal_yields_sel = ComputeYields(signal_file=signal_file, selection=signal_selection_sel).computeSignalYields(lumi=lumi_D1, sigma_B=472.8e9) 
        signal_yields_sel = ComputeYields(signal_label=signal_label, selection=signal_selection_sel).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid)[0]
        #print 'signal_yields_sel ',signal_yields_sel

        # get background yields from unblinded data
        f_data = self.tools.getRootFile(self.data_file)
        tree_data = self.tools.getTree(f_data, 'signal_tree')
        background_selection_sel = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
        if self.preexisting_selection !=None: background_selection_sel += ' && {}'.format(self.getSelectionString())
        background_selection_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        #print background_selection_sel
        hist_data_name = 'hist_data_sel_{}_{}'.format(self.quantity, cut)
        hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=background_selection_sel)
        background_yields_sel = hist_data.Integral()
        #print 'background_yields_sel ',background_yields_sel

        # compute significance
        if background_yields_sel != 0:
          significance_sel = signal_yields_sel / math.sqrt(background_yields_sel)
          #if 2*((signal_yields_sel + background_yields_sel)*math.log(1 + signal_yields_sel/background_yields_sel) - signal_yields_sel) > 0:
          #  significance_sel = math.sqrt(2*((signal_yields_sel + background_yields_sel)*math.log(1 + signal_yields_sel/background_yields_sel) - signal_yields_sel))
          #else:
          #  significance_sel = signal_yields_sel / math.sqrt(background_yields_sel)
          #print 'sel {} {} {} {}'.format(signal_yields_sel, background_yields_sel, significance_sel_1, significance_sel)
        else:
          significance_sel = -99.
        #print '{} / sqrt({}) = {}'.format(signal_yields_sel, background_yields_sel, significance_sel)

        # compute relative difference in significance
        significance_diff = ((significance_sel / significance_ini) - 1) * 100

        if significance_diff_max == -1 or significance_diff > significance_diff_max:
          significance_diff_max = significance_diff

        if significance_diff_min == -1 or significance_diff < significance_diff_min:
          significance_diff_min = significance_diff

        if significance_diff_min < -50: significance_diff_min = -50

        g_sig.SetPoint(idx, cut, significance_diff)
        g_sig.SetLineWidth(2)
        g_sig.SetLineStyle(9)
        g_sig.SetLineColor(signal_file.colour)
        g_sig.SetMarkerColor(signal_file.colour)
        g_sig.SetMarkerStyle(20)

      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetTitle('Relative difference in significance [%]')
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetTitleSize(0.045)
      gs_sig[0].GetYaxis().SetTitleOffset(1.15)
      gs_sig[0].GetYaxis().SetRangeUser(significance_diff_min-0.15*significance_diff_min, significance_diff_max+0.15*significance_diff_max)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          #canv.cd()
          g_sig.Draw('APL')
          #canv_log.cd()
          #g_sig.Draw('APL')
        else:
          #canv.cd()
          g_sig.Draw('PL')
          #canv_log.cd()
          #g_sig.Draw('PL')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, significance_diff_min-0.15*significance_diff_min, self.proposed_cut, significance_diff_max+0.15*significance_diff_max)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      #canv.cd()
      line.Draw('same')
      #canv_log.cd()
      #line.Draw('same')
    
    self.tools.printLatexBox(0.67, 0.25, self.category.title, size=0.041)
          
    #legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.25, xmax=0.82, ymax=0.4, size=0.04)
    for ifile, signal_file in enumerate(self.signal_files):
      legend.AddEntry(gs_sig[ifile], '{}GeV, {}mm'.format(signal_file.mass, signal_file.ctau))
    #canv.cd()
    legend.Draw()
    #canv_log.cd()
    #legend.Draw()

    self.createOutDir('myPlots/selection/{}'.format(self.dirlabel))

    canv.cd()
    canv.SaveAs('myPlots/selection/{}/significance_diff_scan_{}_{}.png'.format(self.dirlabel, self.getLabel(self.quantity.label), self.category.label))

    #canv_log.cd()
    #canv_log.SaveAs('myPlots/selection/{}/significance_diff_scan_{}_{}_log.png'.format(self.dirlabel, self.getLabel(self.quantity.label), self.category.label))



  def printCutflowLine(self, printNum=False):
    cutflow_line = '{qte} {log} {cut} & '.format(
        qte = self.quantity.label, 
        log = self.quantity.logic, 
        cut = self.proposed_cut, 
    )

    for ifile, signal_file in enumerate(self.signal_files):
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)
      f_sig = self.tools.getRootFile(signal_file.filename)
      tree_sig = self.tools.getTree(f_sig, 'signal_tree')

      selection_sig_ini = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
      if self.preexisting_selection !=None: selection_sig_ini += ' && {}'.format(self.getSelectionString())
      initial_sig_entries = self.tools.createHisto(tree_sig, self.hnl_mass, selection=selection_sig_ini).Integral()

      selection_sig_sel = selection_sig_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, self.proposed_cut)
      selected_sig_entries = self.tools.createHisto(tree_sig, self.hnl_mass, selection=selection_sig_sel).Integral()

      selection_bkg_ini = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
      if self.preexisting_selection !=None: selection_bkg_ini += ' && {}'.format(self.getSelectionString())
      if not self.use_data:
        initial_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_ini).Integral() # no weight for the moment
      else:
        f_data = self.tools.getRootFile(self.data_file)
        tree_data = self.tools.getTree(f_data, 'signal_tree')
        hist_data_name = 'hist_data_{}'.format(self.quantity)
        hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=selection_bkg_ini)
        initial_bkg_entries = hist_data.Integral()

      selection_bkg_sel = selection_bkg_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, self.proposed_cut)
      if not self.use_data:
        selected_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_sel).Integral() # no weight for the moment
      else:
        f_data = self.tools.getRootFile(self.data_file)
        tree_data = self.tools.getTree(f_data, 'signal_tree')
        hist_data_name = 'hist_data_sel_{}'.format(self.quantity)
        hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=selection_bkg_sel)
        selected_bkg_entries = hist_data.Integral()

      if printNum: print 'sig m{} ctau{} {} {}'.format(signal_file.mass, signal_file.ctau, initial_sig_entries, selected_sig_entries)
      if printNum: print 'bkg m{} {} {}'.format(signal_file.mass, initial_bkg_entries, selected_bkg_entries)

      if ifile == 0:
        #cutflow_line += ' -{sig_per}\% & -{bkg_per}\% '.format(
        cutflow_line += ' -{sig_per}\% '.format(
          sig_per = round((1 - selected_sig_entries / initial_sig_entries)*100, 2) if initial_sig_entries!=0. else '-', 
          #bkg_per = round((1 - selected_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
        )
      else:
        #cutflow_line += ' & -{sig_per}\% & -{bkg_per}\% '.format(
        cutflow_line += ' & -{sig_per}\% '.format(
          sig_per = round((1 - selected_sig_entries / initial_sig_entries)*100, 2) if initial_sig_entries!=0. else '-', 
          #bkg_per = round((1 - selected_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
        )

    cutflow_line += ' & -{bkg_per}\% '.format(
      bkg_per = round((1 - selected_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
    )

    cutflow_line+= '\\\ '

    print cutflow_line


  def do_print_significance(self):
    '''
      Get the significance before and after the selection
    '''

    for signal_file in self.signal_files:
      signal_mass = signal_file.mass
      signal_ctau = signal_file.ctau
      resolution = self.resolution_p0 + self.resolution_p1 * signal_mass

      g_sig = ROOT.TGraph()
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_mass-2*resolution, bin_max=signal_mass+2*resolution)

      lumi_D1 = 5.302 # when computing the relative difference in significance, this value does not matter

      # compute significance before selection
      signal_selection_ini = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
      #print signal_selection_ini
      #signal_yields_ini, err_signal_yields_ini = ComputeYields(signal_file=signal_file, selection=signal_selection_ini).computeSignalYields(lumi=lumi_D1, sigma_B=472.8e9) 
      signal_yields_ini = ComputeYields(signal_label=signal_label, selection=signal_selection_ini).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid)[0] 

      f_data = self.tools.getRootFile(self.data_file)
      tree_data = self.tools.getTree(f_data, 'signal_tree')
      background_selection_ini = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
      #print background_selection_ini
      hist_data_name = 'hist_data_ini_{}'.format(self.quantity)
      hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=background_selection_ini)
      background_yields_ini = hist_data.Integral()

      if background_yields_ini != 0:
        significance_ini = signal_yields_ini / math.sqrt(background_yields_ini)
        #print '{} {}'.format(signal_yields_ini, background_yields_ini)
        #significance_ini = math.sqrt(2*((signal_yields_ini + background_yields_ini)*math.log(1 + signal_yields_ini/background_yields_ini) - signal_yields_ini))
      else:
        significance_ini = -99.

      # compute significance after selection
      signal_selection_sel = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, self.category.definition_flat) 
      if self.preexisting_selection !=None: signal_selection_sel += ' && {}'.format(self.getSelectionString())
      #print signal_selection_sel

      #signal_yields_sel, err_signal_yields_sel = ComputeYields(signal_file=signal_file, selection=signal_selection_sel).computeSignalYields(lumi=lumi_D1, sigma_B=472.8e9) 
      signal_yields_sel = ComputeYields(signal_label=signal_label, selection=signal_selection_sel).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid)[0] 

      # get background yields from unblinded data
      f_data = self.tools.getRootFile(self.data_file)
      tree_data = self.tools.getTree(f_data, 'signal_tree')
      background_selection_sel = '{} && {}'.format(self.baseline_selection, self.category.definition_flat)
      if self.preexisting_selection !=None: background_selection_sel += ' && {}'.format(self.getSelectionString())
      #print background_selection_sel
      #background_selection_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
      hist_data_name = 'hist_data_sel_{}'.format(self.quantity)
      hist_data = self.tools.createHisto(tree_data, self.hnl_mass, hist_name=hist_data_name, branchname='flat', selection=background_selection_sel)
      background_yields_sel = hist_data.Integral()

      # compute significance
      if background_yields_sel != 0:
        significance_sel = signal_yields_sel / math.sqrt(background_yields_sel)
        #print '{} {}'.format(signal_yields_sel, background_yields_sel)
        #significance_sel = math.sqrt(2*((signal_yields_sel + background_yields_sel)*math.log(1 + signal_yields_sel/background_yields_sel) - signal_yields_sel))
        #print 'sel {} {} {} {}'.format(signal_yields_sel, background_yields_sel, significance_sel_1, significance_sel)
      else:
        significance_sel = -99.

      # compute relative difference in significance
      significance_diff = ((significance_sel / significance_ini) - 1) * 100

      print '{} & {} & {} & {} & {}\% \\\ '.format(signal_file.mass, signal_file.ctau, '{:.2e}'.format(significance_ini), '{:.2e}'.format(significance_sel), round(significance_diff))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  signal_label_m1 = 'V12_08Aug22_m1'
  #signal_label_m1 = 'V12_08Aug22_benchmark'
  signal_label_m3 = 'V12_08Aug22_m3'
  signal_label_m4p5 = 'V12_08Aug22_m4p5'

  signal_m1_ctau1000 = signal_samples[signal_label_m1][0]
  signal_m1_ctau100 = signal_samples[signal_label_m1][1]
  signal_m1_ctau10 = signal_samples[signal_label_m1][2]

  signal_m1p5_ctau1000 = signal_samples['V10_30Dec21_m1p5'][0]
  signal_m2_ctau1000 = signal_samples['V10_30Dec21_m2'][0]

  signal_m3_ctau1000 = signal_samples[signal_label_m3][0]
  signal_m3_ctau100 = signal_samples[signal_label_m3][1]
  #signal_m3_ctau10 = signal_samples[signal_label_m3][2]
  signal_m3_ctau1 = signal_samples[signal_label_m3][2]

  signal_m4p5_ctau100 = signal_samples[signal_label_m4p5][0]
  signal_m4p5_ctau10 = signal_samples[signal_label_m4p5][1]
  signal_m4p5_ctau1 = signal_samples[signal_label_m4p5][2]
  signal_m4p5_ctau0p1 = signal_samples[signal_label_m4p5][3]

  #data_file = data_samples['V10_30Dec21'][0].filename # we only consider D1 which is unblinded
  #data_file = data_samples['V11_24Apr22_small'][0].filename # we only consider D1 which is unblinded
  data_file = data_samples['V12_08Aug22'][0].filename # we only consider D1 which is unblinded
  #data_file = data_samples['V12_08Aug22_small'][0].filename # we only consider D1 which is unblinded
  #data_file = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21_10Chunks.root'
  #data_file = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_30Dec21.root'
  #data_file = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_5files.root'
  #data_file = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22.root'

  baseline_selection = selection['baseline_08Aug22'].flat

  qcd_files = qcd_samples['V10_30Dec21']

  white_list = white_list['20to300']

  hnl_charge = Quantity('hnl_charge', 'hnl_charge', '#mu#pi charge', '==', '', -2, 2)
  #hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.9993, 1) 
  pi_pt = Quantity('pi_pt', 'pi_pt', 'displaced pion pT', '>', 'GeV', 0, 8)
  mu_dxy = Quantity('abs(mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 0.3)
  pi_dxy = Quantity('abs(pi_dxy)', 'pi_dxy', 'displaced pion |IP| on xy', '>', 'cm', 0, 1)
  max_mu_pi_dxysig = Quantity('max(abs(pi_dxysig), abs(mu_dxysig))', 'max_mu_pi_dxysig', 'max(muon, pion) |IP| significance on xy', '>', 'cm', 0, 80)
  min_mu_pi_dxysig = Quantity('min(abs(pi_dcasig), abs(mu_dxysig_bs))', 'min_mu_pi_dxysig', 'min(muon, pion) |IP| significance on xy (BS)', '>', 'cm', 0, 200)
  mu_dxysig = Quantity('abs(mu_dxysig)', 'mu_dxysig', 'displaced muon |IP| significance on xy', '>', 'cm', 0, 50)
  pi_dxysig = Quantity('abs(pi_dxysig)', 'pi_dxysig', 'displaced pion |IP| significance on xy', '>', 'cm', 0, 100)
  pi_dcasig = Quantity('abs(pi_dcasig)', 'pi_dcasig', 'pion DCA significance', '>', 'cm', 0, 50)
  cos_theta_star_pion = Quantity('fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.5, 1) 
  mu0_softid = Quantity('mu0_softid', 'mu0_softid', 'trigger muon soft ID', '==', '', 0, 1)
  mu_looseid = Quantity('mu_looseid', 'mu_looseid', 'muon loose ID', '==', '', 0, 1)
  #mu_customid = Quantity('mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))', 'mu_customid', 'muon custom ID', '==', '', 0, 1)
  mu_customid = Quantity('mu_customisedid', 'mu_customisedid', 'muon custom ID', '==', '', 0, 1)
  mu_whnlid = Quantity('mu_whnlid', 'mu_whnlid', 'muon whnl ID', '==', '', 0, 1)
  b_mass = Quantity('abs(b_mass-5.3)', 'b_mass', '#mu#mu#pi mass', '<', 'GeV', 0, 1.5) 
  sv_lxy_sig = Quantity('sv_lxysig', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 500)
  mu_numberofpixellayers = Quantity('mu_numberofpixellayers', 'mu_numberofpixellayers', 'mu_numberofpixellayers', '<', '', 0, 6)
  mu_numberoftrackerlayers = Quantity('mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', '<', '', 0, 18)
  deltar_trgmu_pi = Quantity('abs(deltar_trgmu_pi)', 'abs(deltar_trgmu_pi)', '|#DeltaR(trg#mu, #pi)|', '<', '', 0.3, 1.) 
  deltar_trgmu_hnl = Quantity('abs(deltar_trgmu_hnl)', 'abs(deltar_trgmu_hnl)', '|#DeltaR(trg#mu, hnl)|', '<', '', 0., 0.5) 
  mu0_pfiso03_rel = Quantity('mu0_pfiso03_rel', 'mu0_pfiso03_rel', 'primary muon relative PF iso03', '<', '', 0.3, 3)
  mu0_pfiso04_rel = Quantity('mu0_pfiso04_rel', 'mu0_pfiso04_rel', 'primary muon relative PF iso04', '<', '', 0.3, 3)
  mu_pfiso03_rel = Quantity('mu_pfiso03_rel', 'mu_pfiso03_rel', 'displaced muon relative PF iso03', '<', '', 0., 0.5)
  mu_pfiso04_rel = Quantity('mu_pfiso04_rel', 'mu_pfiso04_rel', 'displaced muon relative PF iso04', '<', '', 0.3, 3)
  
  #signal_files = [signal_m1_ctau1000, signal_m3_ctau100, signal_m4p5_ctau1]
  #signal_files = [signal_m1_ctau1000]
  signal_files = [signal_samples['V12_08Aug22_benchmark'][0]]
  signal_label = 'V12_08Aug22_benchmark'
  #signal_label = signal_label_m1
  #signal_files = [signal_m1_ctau1000, signal_m1_ctau100, signal_m1_ctau10]

  printSelectionString = True

  opt = getOptions()

  dirlabel = opt.dirlabel
  mass = opt.mass
  category = opt.category
  quantity = opt.quantity
  action = opt.action

  if mass == '1': 
    signal_label = signal_label_m1
    signal_files = [signal_m1_ctau1000, signal_m1_ctau100, signal_m1_ctau10]
  if mass == '3': 
    signal_label = signal_label_m3
    signal_files = [signal_m3_ctau1000, signal_m3_ctau100, signal_m3_ctau1]
  if mass == '4p5': 
    signal_label = signal_label_m4p5
    signal_files = [signal_m4p5_ctau100, signal_m4p5_ctau10, signal_m4p5_ctau1]

  do_lxy0to1_OS = False
  do_lxy0to1_SS = False
  do_lxy1to5_OS = False
  do_lxy1to5_SS = False
  do_lxygt5_OS = False
  do_lxygt5_SS = False

  if category == 'lxy0to1_OS': do_lxy0to1_OS = True
  if category == 'lxy0to1_SS': do_lxy0to1_SS = True
  if category == 'lxy1to5_OS': do_lxy1to5_OS = True
  if category == 'lxy1to5_SS': do_lxy1to5_SS = True
  if category == 'lxygt5_OS': do_lxygt5_OS = True
  if category == 'lxygt5_SS': do_lxygt5_SS = True

  if category == 'lxy0to1_OS' or category == 'lxy0to1_SS': hnl_cos2D = Quantity('1-hnl_cos2d', 'hnl_cos2D', '(1 - SV cos2D)', '<', '', 0, 3e-3) 
  if category == 'lxy1to5_OS' or category == 'lxy1to5_SS': hnl_cos2D = Quantity('1-hnl_cos2d', 'hnl_cos2D', '(1 - SV cos2D)', '<', '', 0, 3e-4) 
  if category == 'lxygt5_OS' or category == 'lxygt5_SS': hnl_cos2D = Quantity('1-hnl_cos2d', 'hnl_cos2D', '(1 - SV cos2D)', '<', '', 0, 3e-5) 

  do_b_mass = False
  do_pi_pt = False
  do_sv_lxysig = False
  do_mu_dxysig = False
  do_pi_dxysig = False
  do_max_mu_pi_dxysig = False
  do_min_mu_pi_dxysig = False
  do_hnl_cos2d = False

  if quantity == 'b_mass': do_b_mass = True
  if quantity == 'pi_pt': do_pi_pt = True
  if quantity == 'sv_lxysig': do_sv_lxysig = True
  if quantity == 'mu_dxysig': do_mu_dxysig = True
  if quantity == 'pi_dxysig': do_pi_dxysig = True
  if quantity == 'max_mu_pi_dxysig': do_max_mu_pi_dxysig = True
  if quantity == 'min_mu_pi_dxysig': do_min_mu_pi_dxysig = True
  if quantity == 'hnl_cos2d': do_hnl_cos2d = True

  do_scan_efficiency = False
  do_scan_significance = False
  do_scan_significance_diff = False
  do_print_cutflow = False
  do_print_significance = False

  if action == 'scan_efficiency': do_scan_efficiency = True
  if action == 'scan_significance': do_scan_significance = True
  if action == 'scan_significance_diff': do_scan_significance_diff = True
  if action == 'print_cutflow': do_print_cutflow = True
  if action == 'print_significance': do_print_significance = True 

  ### CATEGORY inclusive ###

  category_incl = Category(
      label = 'incl',
      title = 'Inclusive',
      definition_flat = 'hnl_pt>0',
  )

  selection_sequential_incl = [] 

  #cut_mu0_pfiso04_rel = 0.5
  ##selection_sequential_incl.append(SelectedQuantity(mu0_pfiso04_rel, cut_mu0_pfiso04_rel))
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_pfiso04_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu0_pfiso04_rel, preexisting_selection=selection_sequential_incl)
  #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_pfiso04_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu0_pfiso04_rel) 
  ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
  #if do_print_cutflow: selection.printCutflowLine(printNum=False)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
  #selection.getSignificanceDiffScan(npoints=10)
  ##selection_sequential_incl.append(SelectedQuantity(mu0_pfiso04_rel, cut_mu0_pfiso04_rel))

  #selection.do_print_significance()

  #cut_mu0_pfiso03_rel = 0.5
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_pfiso03_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu0_pfiso03_rel, preexisting_selection=selection_sequential_incl)
  #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_pfiso03_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu0_pfiso03_rel) 
  ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
  #if do_print_cutflow: selection.printCutflowLine(printNum=False)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
  #selection.getSignificanceDiffScan(npoints=10)
  ##selection_sequential_incl.append(SelectedQuantity(mu0_pfiso03_rel, cut_mu0_pfiso03_rel))

  #cut_mu_pfiso04_rel = 0.5
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_pfiso04_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu_pfiso04_rel, preexisting_selection=selection_sequential_incl)
  #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_pfiso04_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu_pfiso04_rel) 
  ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
  #if do_print_cutflow: selection.printCutflowLine(printNum=False)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
  #selection.getSignificanceDiffScan(npoints=10)
  ##selection_sequential_incl.append(SelectedQuantity(mu_pfiso04_rel, cut_mu_pfiso04_rel))

  cut_mu_pfiso03_rel = 0.5
  #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_pfiso03_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu_pfiso03_rel, preexisting_selection=selection_sequential_incl)
  selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_pfiso03_rel, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu_pfiso03_rel) 
  #if do_scan_efficiency: selection.getScanGraph(npoints=2)
  if do_print_cutflow: selection.printCutflowLine(printNum=False)
  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
  selection.getSignificanceDiffScan(npoints=10)
  #selection_sequential_incl.append(SelectedQuantity(mu_pfiso03_rel, cut_mu_pfiso03_rel))

  cut_mu0_softid = 1
  selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_incl, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
  #if do_scan_efficiency: selection.getScanGraph(npoints=2)
  if do_print_cutflow: selection.printCutflowLine(printNum=False)
  selection_sequential_incl.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_incl, dirlabel=dirlabel, preexisting_selection=selection_sequential_incl, proposed_cut=cut_mu_looseid)
  #if do_scan_efficiency: selection.getScanGraph(npoints=2)
  if do_print_cutflow: selection.printCutflowLine()
  #selection_sequential_incl.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_whnlid = 1
  selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_incl, dirlabel=dirlabel, preexisting_selection=selection_sequential_incl, proposed_cut=cut_mu_whnlid)
  #if do_scan_efficiency: selection.getScanGraph(npoints=2)
  if do_print_cutflow: selection.printCutflowLine()
  ##selection_sequential_incl.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_mu_customid = 1
  selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_incl, dirlabel=dirlabel, preexisting_selection=selection_sequential_incl, proposed_cut=cut_mu_customid)
  #if do_scan_efficiency: selection.getScanGraph(npoints=2)
  if do_print_cutflow: selection.printCutflowLine()
  ##selection_sequential_incl.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  ### CATEGORY lxygt5_OS ###

  if do_lxygt5_OS:
    category_lxygt5_OS = Category(
        label = 'lxygt5_OS',
        title = 'l_{xy}>5cm, OS',
        definition_flat = 'sv_lxy>5 && mu0_charge!=mu_charge',
    )

    selection_sequential_lxygt5_OS = [] 

    #signal_files = [signal_m1_ctau1000, signal_m1p5_ctau1000, signal_m2_ctau1000, signal_m3_ctau100]

    print 'category lxygt5_OS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxygt5_OS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxygt5_OS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxygt5_OS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    ##cut_deltar_trgmu_pi = 1.5
    ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=deltar_trgmu_pi, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_deltar_trgmu_pi)
    ##if do_scan_efficiency: selection.getScanGraph(npoints=30)
    ##selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(deltar_trgmu_pi, cut_deltar_trgmu_pi))

    ##cut_pi_dcasig = 20
    ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dcasig, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_pi_dcasig)
    ##if do_scan_efficiency: selection.getScanGraph(npoints=30)
    ###selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxygt5_OS.append(SelectedQuantity(pi_dcasig, cut_pi_dcasig))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()

  #### CATEGORY lxygt5_SS ###

  if do_lxygt5_SS:
    category_lxygt5_SS = Category(
        label = 'lxygt5_SS',
        title = 'l_{xy}>5cm, SS',
        definition_flat = 'sv_lxy>5 && mu0_charge==mu_charge',
    )

    selection_sequential_lxygt5_SS = [] 

    print 'category lxygt5_SS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxygt5_SS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxygt5_SS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
      if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
      if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxygt5_SS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxygt5_SS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxygt5_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()


  ### CATEGORY lxy1to5_OS ###

  if do_lxy1to5_OS:
    category_lxy1to5_OS = Category(
        label = 'lxy1to5_OS',
        title = '(1<l_{xy}<=5)cm, OS',
        definition_flat = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge!=mu_charge',
    )

    selection_sequential_lxy1to5_OS = [] 

    print 'category lxy1to5_OS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxy1to5_OS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxy1to5_OS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to5_SS ###

  if do_lxy1to5_SS:
    category_lxy1to5_SS = Category(
        label = 'lxy1to5_SS',
        title = '(1<l_{xy}<=5)cm, SS',
        definition_flat = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge==mu_charge',
    )

    selection_sequential_lxy1to5_SS = [] 

    print 'category lxy1to5_SS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxy1to5_SS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxy1to5_SS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy1to5_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()


  ### CATEGORY lxy0to1_OS ###

  if do_lxy0to1_OS:
    category_lxy0to1_OS = Category(
        label = 'lxy0to1_OS',
        title = '(0<l_{xy}<=1)cm, OS',
        definition_flat = 'sv_lxy<=1 && mu0_charge!=mu_charge',
    )

    selection_sequential_lxy0to1_OS = [] 

    print 'category lxy0to1_OS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxy0to1_OS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=3)
    #if do_scan_significance: selection.getSignificanceScan(npoints=3)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    #hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.99, 1) 
    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()


### CATEGORY lxy0to1_SS ###

  if do_lxy0to1_SS:
    category_lxy0to1_SS = Category(
        label = 'lxy0to1_SS',
        title = '(0<l_{xy}<=1)cm, SS',
        definition_flat = 'sv_lxy<=1 && mu0_charge==mu_charge',
    )

    selection_sequential_lxy0to1_SS = [] 

    print 'category lxy0to1_SS'

    cut_mu0_softid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu0_softid, category=category_lxy0to1_SS, dirlabel=dirlabel, proposed_cut=cut_mu0_softid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_looseid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

    cut_mu_whnlid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_whnlid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_mu_customid = 1
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_customid)
    #if do_scan_efficiency: selection.getScanGraph(npoints=2)
    if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

    cut_hnl_charge = 0
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_hnl_charge)
    #if do_scan_efficiency: selection.getScanGraph(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

    cut_b_mass = float(opt.cut_b_mass)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_b_mass)
    if do_b_mass:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(b_mass, cut_b_mass))

    cut_pi_pt = float(opt.cut_pi_pt)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_pi_pt)
    if do_pi_pt:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

    cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_sv_lxy_sig)
    if do_sv_lxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

    #cut_mu_dxysig = float(opt.cut_mu_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_dxysig)
    #if do_mu_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine()
    ##selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

    #cut_pi_dxysig = float(opt.cut_pi_dxysig)
    #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_pi_dxysig)
    #if do_pi_dxysig:
    #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
    #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
    #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    #if do_print_cutflow: selection.printCutflowLine(printNum=False)
    ##selection_sequential_lxy0to1_SS.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

    cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_min_mu_pi_dxysig)
    if do_min_mu_pi_dxysig:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine(printNum=False)
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

    #hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.99, 1) 
    cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
    selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_hnl_cos2D)
    if do_hnl_cos2d:
      if do_scan_efficiency: selection.getScanGraph(npoints=30)
      if do_scan_significance: selection.getSignificanceScan(npoints=30)
      if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
    if do_print_cutflow: selection.printCutflowLine()
    selection_sequential_lxy0to1_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    if do_print_significance: selection.do_print_significance()

    if printSelectionString: print selection.getSelectionString()



  ##cut_mu_numberofpixellayers = 2 
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_numberofpixellayers, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_numberofpixellayers)
  ##if do_scan_efficiency: selection.getScanGraph(npoints=7)
  ##if do_print_cutflow: selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(mu_numberofpixellayers, cut_mu_numberofpixellayers))

  ##cut_mu_numberoftrackerlayers = 13 
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_numberoftrackerlayers, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_numberoftrackerlayers)
  ##if do_scan_efficiency: selection.getScanGraph(npoints=19)
  ##if do_print_cutflow: selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(mu_numberoftrackerlayers, cut_mu_numberoftrackerlayers))

  ##cut_sv_lxy_sig = 30
  ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_sv_lxy_sig)
  ##if do_scan_efficiency: selection.getScanGraph(npoints=30)
  ##if do_print_cutflow: selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

