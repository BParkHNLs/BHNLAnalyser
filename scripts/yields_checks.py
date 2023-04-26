import os
import sys
import numpy as np
from os import path
import math
from tools import Tools
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples, data_samples, qcd_samples
from categories import categories
from ABCD_regions import ABCD_regions
from baseline_selection import selection
from qcd_white_list import white_list
from ctau_points import ctau_points

import ROOT


class YieldsChecks(Tools):
  def __init__(self, data_files='', qcd_files='', signal_files='', white_list='', categories='', selection='', ABCD_regions=''):
    self.tools = Tools()
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_files = signal_files
    self.white_list = white_list
    self.categories = categories
    self.selection = selection
    self.ABCD_regions = ABCD_regions


  def plotSigYields(self, lumi=0.774, selection='', label='', outdirlabel='', category=None):
    #canv = self.tools.createTCanvas('canv', dimx=900, dimy=800)
    canv = self.tools.createTCanvas('canv', dimx=1200, dimy=800)
    canv.SetLogx()
    canv.SetLogy()
    canv.SetGrid()

    # gen yields
    gen_coupling_m1 = [5.37675839253e-06, 5.37675839253e-05, 0.000537675839253, 0.00537675839253, 0.0537675839253] #[5.37675839253e-06, 5.37675839253e-05, 0.000537675839253, 0.00537675839253, 0.0537675839253]
    #gen_yields_m1 = [0.233520553543, 21.1107914574, 2521.79127222, 198534.792249, 9787694.93131] # wrong weight
    #gen_yields_m1 = [0.222669464601, 21.7748154625, 2054.4123073, 164459.789617, 5551279.31665] # good weight, without trg mu eff
    gen_yields_m1 = [0.271587023481, 26.2785112766, 2480.63915671, 186917.31518, 6378867.68389] #[0.271587023481, 26.2785112766, 2480.63915671, 186917.31518, 6378867.68389] # good weight, with trg mu eff

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
      graph_gen_m1.SetPoint(point, gen_coupling_m1[pt], gen_yields_m1[pt]*472.8/327)

    graph_gen_m1.SetMarkerStyle(22)
    graph_gen_m1.SetMarkerSize(2)
    graph_gen_m1.SetMarkerColor(ROOT.kOrange+0)
    graph_gen_m1.SetLineStyle(9)
    graph_gen_m1.SetLineWidth(2)
    graph_gen_m1.SetLineColor(ROOT.kOrange+0)

    for pt in range(0, len(gen_coupling_m3)):
      point = graph_gen_m3.GetN()
      graph_gen_m3.SetPoint(point, gen_coupling_m3[pt], gen_yields_m3[pt]*472.8/327)

    graph_gen_m3.SetMarkerStyle(22)
    graph_gen_m3.SetMarkerSize(2)
    graph_gen_m3.SetMarkerColor(ROOT.kRed+1)
    graph_gen_m3.SetLineStyle(9)
    graph_gen_m3.SetLineWidth(2)
    graph_gen_m3.SetLineColor(ROOT.kRed+1)

    for pt in range(0, len(gen_coupling_m4p5)):
      point = graph_gen_m4p5.GetN()
      graph_gen_m4p5.SetPoint(point, gen_coupling_m4p5[pt], gen_yields_m4p5[pt]*472.8/327)

    graph_gen_m4p5.SetMarkerStyle(22)
    graph_gen_m4p5.SetMarkerSize(2)
    graph_gen_m4p5.SetMarkerColor(ROOT.kRed+4)
    graph_gen_m4p5.SetLineStyle(9)
    graph_gen_m4p5.SetLineWidth(2)
    graph_gen_m4p5.SetLineColor(ROOT.kRed+4)

    #signal_label_m1 = 'V10_30Dec21_m1_large'
    #signal_label_m1 = 'V10_30Dec21_m1'
    signal_label_m1 = 'V12_08Aug22_m1'
    #signal_label_m1 = 'V42_08Aug22_m1p02'
    #signal_label_m3 = 'V10_30Dec21_m3'
    signal_label_m3 = 'V12_08Aug22_m3'
    #signal_label_m3 = 'V42_08Aug22_m3p0'
    #signal_label_m4p5 = 'V10_30Dec21_m4p5_large'
    #signal_label_m3 = 'V11_24Apr22_m3'
    #signal_label_m4p5 = 'central_V11_24Apr22_m4p5'
    #signal_label_m4p5 = 'V10_30Dec21_m4p5'
    signal_label_m4p5 = 'V12_08Aug22_m4p5'
    #signal_label_m4p5 = 'V42_08Aug22_m4p5'
    #signal_label_m4p5 = 'V11_24Apr22_m4p5'
    #signal_label_m1 = 'V00_looseselection_m1'
    #signal_label_m3 = 'V00_looseselection_m3'
    #signal_label_m4p5 = 'V00_looseselection_m4p5'
    samples_m1 = signal_samples[signal_label_m1]
    samples_m3 = signal_samples[signal_label_m3]
    samples_m4p5 = signal_samples[signal_label_m4p5]

    #isMixed = True

    #ctau_points = ctau_points['baseline']
    strategy = 'inclusive'

    resolution_p0 = 0.0002747 
    resolution_p1 = 0.008302 
    window_size = 4

    #ctau_points_m1 = [10., 15., 20., 30., 40., 50., 70., 100., 150., 200., 300., 400., 500., 700., 1000., 1500., 2000., 3000., 4000., 5000., 7000., 10000., 15000., 20000., 30000., 40000., 50000., 70000.] 
    #ctau_points_m3 = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1., 1.5, 2., 3., 4., 5., 7., 10., 15., 20., 30., 40., 50., 70., 100., 150., 200., 300., 400., 500., 700., 1000.] 
    #ctau_points_m3 = [0.5, 0.7 ,1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0]
    #ctau_points_m3 = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1., 1.5, 2., 3., 4., 5., 7., 10., 15., 20., 30., 40., 50., 70., 100., 150.] 
    ctau_points_m3 = [1., 10., 100., 1000.] 
    #ctau_points_m4p5 = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1., 1.5, 2., 3., 4., 5., 7., 10., 15., 20., 30., 40., 50., 70., 100.] 
    #ctau_points_m4p5 = [0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07]
    ctau_points_m4p5 = [0.1, 1., 10., 100.] 

    #is_bc = False if '_Bc' not in category.label else True
    is_bc = False
    addBc = False

    line_style = 1 if not is_bc else 9

    #selection = selection + ' && sv_lxysig>=150 && mu0_charge!=mu_charge && b_mass>5.7' #FIXME remove
    selection = selection + ' && ' + category.definition_flat
    weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable'
    weight_pusig = 'weight_pu_sig_D'
    weight_mu0id = 'weight_mu0_softid'
    weight_muid = 'weight_mu_looseid'
    add_weight_hlt = True
    add_weight_pu = True
    add_weight_muid = True

    graph_dummy = ROOT.TGraph()
    graph_m1 = ROOT.TGraphAsymmErrors()
    graph_m3 = ROOT.TGraphAsymmErrors()
    graph_m4p5 = ROOT.TGraphAsymmErrors()

    #graph_dummy.SetPoint(0, 1e-6, 1e-4)
    graph_dummy.SetPoint(0, 1e-6, 1e-3)
    #graph_dummy.SetPoint(0, 1e-6, 1e-0)
    #graph_dummy.SetPoint(1, 1e-1, 1e10)
    graph_dummy.SetPoint(1, 1e-2, 1e4)
    #graph_dummy.SetPoint(1, 1e-1, 1e4)
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

    #print '\n mass 1'
    #for ctau_point in ctau_points_m1:
    #  # get signal coupling
    #  signal_mass = 1.0
    #  signal_ctau = ctau_point
    #  signal_resolution = resolution_p0 + resolution_p1 * signal_mass
    #  print '\n',signal_ctau
    #  signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

    #  # compute the signal yields
    #  #signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
    #  signal_selection = 'hnl_mass>{} && hnl_mass<{} && {}'.format(signal_mass-window_size*signal_resolution, signal_mass+window_size*signal_resolution, selection)
    #  signal_yields, err = ComputeYields(signal_label=signal_label_m1, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=False, strategy=strategy) 
    #  print signal_yields
    #  if addBc:
    #    signal_yields += ComputeYields(signal_file=signal_file, signal_label=signal_label_m1, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9, is_bc=True)[0]
    #    print signal_yields

    #  # fill graph
    #  point = graph_m1.GetN()
    #  graph_m1.SetPoint(point, signal_v2, signal_yields)
    #  graph_m1.SetPointError(point, 0, 0, err, err)
    #  #graph_m1.SetPointError(point, 0, 0, 0, 0)

    #graph_m1.SetMarkerStyle(20)
    #graph_m1.SetMarkerSize(2)
    #graph_m1.SetMarkerColor(ROOT.kOrange+0)
    ##graph_m1.SetMarkerColor(ROOT.kBlue)
    #graph_m1.SetLineStyle(1)
    #graph_m1.SetLineWidth(2)
    #graph_m1.SetLineColor(ROOT.kOrange+0)

    print '\n mass 3'
    for ctau_point in ctau_points_m3:
      # get signal coupling
      signal_mass = 3.
      signal_ctau = ctau_point
      signal_resolution = resolution_p0 + resolution_p1 * signal_mass
      print '\n',signal_ctau
      signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

      # compute the signal yields
      signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
      #signal_selection = 'hnl_mass>{} && hnl_mass<{} && {}'.format(signal_mass-window_size*signal_resolution, signal_mass+window_size*signal_resolution, selection)
      #signal_selection += ' && sv_lxysig>150 && mu0_charge!=mu_charge' 
      #signal_selection += ' && sv_lxysig>0 && sv_lxysig<=50 && mu0_charge!=mu_charge' 
      #print signal_selection
      signal_yields, err = ComputeYields(signal_label=signal_label_m3, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=is_bc, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_pusig=weight_pusig, weight_mu0id=weight_mu0id, weight_muid=weight_muid, strategy=strategy) 
      print '{} {}%'.format(signal_yields, round(err/signal_yields*100, 1))
      if addBc:
        signal_yields += ComputeYields(signal_label=signal_label_m3, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=True, strategy=strategy)[0]
        print signal_yields

      # fill graph
      point = graph_m3.GetN()
      graph_m3.SetPoint(point, signal_v2, signal_yields)
      graph_m3.SetPointError(point, 0, 0, err, err)
      #graph_m3.SetPointError(point, 0, 0, 0, 0)

    graph_m3.SetMarkerStyle(20)
    graph_m3.SetMarkerSize(2)
    graph_m3.SetMarkerColor(ROOT.kRed+1)
    #graph_m3.SetMarkerColor(ROOT.kBlue)
    graph_m3.SetLineStyle(line_style)
    graph_m3.SetLineWidth(2)
    graph_m3.SetLineColor(ROOT.kRed+1)

    print '\n mass 4.5'
    for ctau_point in ctau_points_m4p5:
      # get signal coupling
      signal_mass = 4.5
      signal_ctau = ctau_point
      signal_resolution = resolution_p0 + resolution_p1 * signal_mass
      print '\n',signal_ctau
      signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

      # compute the signal yields
      signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
      #signal_selection = 'hnl_mass>{} && hnl_mass<{} && {}'.format(signal_mass-window_size*signal_resolution, signal_mass+window_size*signal_resolution, selection)
      print '\n non-Bc'
      signal_yields, err = ComputeYields(signal_label=signal_label_m4p5, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=is_bc, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_pusig=weight_pusig, weight_mu0id=weight_mu0id, weight_muid=weight_muid, strategy=strategy) 
      print signal_yields
      if addBc:
        print '\n Bc'
        signal_yields_bc = ComputeYields(signal_label=signal_label_m4p5, selection=signal_selection).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=True, strategy=strategy)[0]
        print '{} {} {}'.format(signal_yields, signal_yields_bc, (signal_yields+signal_yields_bc)/signal_yields)
        signal_yields += signal_yields_bc
        print signal_yields

      # fill graph
      point = graph_m4p5.GetN()
      graph_m4p5.SetPoint(point, signal_v2, signal_yields)
      graph_m4p5.SetPointError(point, 0, 0, err, err)
      #graph_m4p5.SetPointError(point, 0, 0, 0, 0)

    graph_m4p5.SetMarkerStyle(20)
    graph_m4p5.SetMarkerSize(2)
    graph_m4p5.SetMarkerColor(ROOT.kRed+4)
    #graph_m4p5.SetMarkerColor(ROOT.kBlue)
    graph_m4p5.SetLineStyle(line_style)
    graph_m4p5.SetLineWidth(2)
    graph_m4p5.SetLineColor(ROOT.kRed+4)

    graph_dummy.Draw('AP')  
    #graph_m1.Draw('PL same')  
    graph_m3.Draw('PL same')  
    graph_m4p5.Draw('PL same')  
    #graph_gen_m1.Draw('PL same')
    #graph_gen_m3.Draw('PL same')
    #graph_gen_m4p5.Draw('PL same')
      
    legend = self.tools.getRootTLegend(xmin=0.15, ymin=0.55, xmax=0.45, ymax=0.9, size=0.033)
    #legend = self.tools.getRootTLegend(xmin=0.15, ymin=0.35, xmax=0.45, ymax=0.55, size=0.033)
    #legend.AddEntry(graph_m1, 'm=1GeV, reco')
    #legend.AddEntry(graph_m1, 'm=1GeV')
    #legend.AddEntry(graph_gen_m1, 'm=1GeV, gen')
    #legend.AddEntry(graph_m3, 'm=3GeV, reco')
    legend.AddEntry(graph_m3, 'm=3GeV')
    #legend.AddEntry(graph_gen_m3, 'm=3GeV, gen')
    #legend.AddEntry(graph_m4p5, 'm=4p5GeV, reco')
    legend.AddEntry(graph_m4p5, 'm=4p5GeV')
    #legend.AddEntry(graph_m4p5, 'generated')
    #legend.AddEntry(graph_gen_m4p5, 'm=4p5GeV, gen')
    legend.Draw()

    if not path.exists('./myPlots/yields'):
      os.system('mkdir -p ./myPlots/yields')

    label = label + '_' + category.label
    canv.SaveAs('./myPlots/yields/signal_yields_{}.png'.format(label))
    canv.SaveAs('./myPlots/yields/signal_yields_{}.pdf'.format(label))


  def getMCWeight(self, signal_files, signal_mass, signal_ctau, lumi, sigma_B=472.8e9, lhe_efficiency=0.08244, add_weight_hlt=False, add_weight_pu=False, add_weight_muid=False, weight_hlt='', weight_pusig='', weight_mu0id='', weight_muid='', is_bc=False):
    weight_sig_list = ['gen_hnl_ct']
    weight_ctau = self.tools.getCtauWeight(signal_files=signal_files, ctau=signal_ctau, is_bc=is_bc)
    weight_signal = self.tools.getSignalWeight(signal_files=signal_files, mass=signal_mass, ctau=signal_ctau, sigma_B=sigma_B, lumi=lumi, lhe_efficiency=lhe_efficiency, is_bc=is_bc)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
    if add_weight_hlt: 
      weight_sig += ' * ({})'.format(weight_hlt)
      weight_sig_list.append(weight_hlt)
    if add_weight_pu: 
      weight_sig += ' * ({})'.format(weight_pusig)
      weight_sig_list.append(weight_pusig)
    if add_weight_muid: 
      weight_sig += ' * ({}) * ({})'.format(weight_mu0id, weight_muid)
      weight_sig_list.append(weight_mu0id)
      weight_sig_list.append(weight_muid)

    return weight_sig, weight_sig_list


  def writeYieldsTable(self):
    baseline_selection = selection['baseline_08Aug22'].flat

    masses = [3., 4.5, 5.5]
    #signal_labels = ['V12_08Aug22_m1', 'V12_08Aug22_m3', 'V12_08Aug22_m4p5']
    signal_labels = ['V12_08Aug22_genmatching']
    ctau_points = [0.1, 1, 10, 100, 1000]
    strategy = 'inclusive'

    lumi = 41.6

    weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2'
    weight_pusig = 'weight_pu_sig_D'
    weight_mu0id = 'weight_mu0_softid'
    weight_muid = 'weight_mu_looseid'
    add_weight_hlt = True
    add_weight_pu = True
    add_weight_muid = True

    resolution_p0 = 6.98338e-04 
    resolution_p1 = 7.78382e-03 
    window_size = 4

    is_bc = True

    for imass, mass in enumerate(masses):
      for ctau_point in ctau_points:
        #print '\n'
        # get signal coupling
        signal_mass = mass
        signal_ctau = ctau_point
        signal_resolution = resolution_p0 + resolution_p1 * signal_mass
        signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)

        # compute the signal yields
        signal_selection_1 = 'ismatched==1 && {}'.format(baseline_selection)
        signal_selection_2 = 'hnl_mass>{} && hnl_mass<{} && {}'.format(signal_mass-window_size*signal_resolution, signal_mass+window_size*signal_resolution, baseline_selection)
        print signal_selection_2

        signal_yields_1, err_1 = ComputeYields(signal_label=signal_labels[0], selection=signal_selection_1).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=is_bc, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_pusig=weight_pusig, weight_mu0id=weight_mu0id, weight_muid=weight_muid, strategy=strategy) 
        signal_yields_2, err_2 = ComputeYields(signal_label=signal_labels[0], selection=signal_selection_2).computeSignalYields(mass=signal_mass, ctau=signal_ctau, lumi=lumi, sigma_B=472.8e9, is_bc=is_bc, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_pusig=weight_pusig, weight_mu0id=weight_mu0id, weight_muid=weight_muid, strategy=strategy) 
        #if addBc:
        #  signal_yields += ComputeYields(signal_file=signal_file, signal_label=signal_label_m1, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9, isBc=True)[0]
        #  print signal_yields
        print '{} GeV & {} mm & {} & {} & {} & {} \\\ '.format(signal_mass, signal_ctau, self.tools.getCouplingLabel(signal_v2), '{:.2e}'.format(signal_yields_1), '{:.2e}'.format(signal_yields_2), round(signal_yields_2/signal_yields_1, 2))


  def plotSigBkgYields(self, lumi=0.774, selection='', label='', doABCD=False, doABCDHybrid=False, doTF=False):
    for category in self.categories:
      print category.label
      canv = self.tools.createTCanvas('canv'+category.label, dimx=1200, dimy=1000)
      canv.SetLogx()
      canv.SetLogy()
      canv.SetGrid()

      graph_dummy = ROOT.TGraph()

      #graph_dummy.SetPoint(0, 1e-6, 1e-7)
      graph_dummy.SetPoint(0, 1e0, 1e-5)
      #graph_dummy.SetPoint(1, 1e-1, 1e10)
      #graph_dummy.SetPoint(1, 1, 1e10)
      graph_dummy.SetPoint(1, 1, 1e4)
      graph_dummy.SetMarkerStyle(0)
      graph_dummy.SetMarkerSize(0)
      graph_dummy.SetMarkerColor(0)
      graph_dummy.GetXaxis().SetTitle('|V^{2}|')
      graph_dummy.GetXaxis().SetLabelSize(0.037)
      graph_dummy.GetXaxis().SetTitleSize(0.042)
      graph_dummy.GetXaxis().SetTitleOffset(1.1)
      graph_dummy.GetYaxis().SetTitle('Yields')
      graph_dummy.GetYaxis().SetLabelSize(0.037)
      graph_dummy.GetYaxis().SetTitleSize(0.042)
      graph_dummy.GetYaxis().SetTitleOffset(1.1)

      graph_dummy.Draw('AP')  

      legend = self.tools.getRootTLegend(xmin=0.15, ymin=0.55, xmax=0.45, ymax=0.9, size=0.027)
      #legend = self.tools.getRootTLegend(xmin=0.25, ymin=0.2, xmax=0.75, ymax=0.35, size=0.033)

      samples_m1 = signal_samples['central_V09_06Nov21_m1']
      samples_m3 = signal_samples['central_V09_06Nov21_m3']
      samples_m4p5 = signal_samples['central_V09_06Nov21_m4p5']

      #  print signal_files
      #colour = [ROOT.kOrange+0, ROOT.kRed+1, ROOT.kRed+4]
      #resolution = [1.59e-03, 5.58e-03, 0.0368] 
      #for imass, mass in enumerate(['1', '3', '4.5']):
      #for ifile, signal_files_ in enumerate(self.signal_files):

      graph_m1 = ROOT.TGraphAsymmErrors()
      #for signal_file in signal_files_:
      for signal_file in samples_m1:
      #for signal_file in self.signal_files:
        #print signal_file.filename
        # get signal coupling
        signal_mass = signal_file.mass
        #if float(signal_mass) != float(mass): continue
        signal_ctau = signal_file.ctau
        signal_resolution = signal_file.resolution
        signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
        signal_colour = signal_file.colour
        print '{}Gev {}mm'.format(signal_mass, signal_ctau)

        # compute the signal yields
        signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
        signal_selection += ' && ' + category.definition_flat + '&& ' + category.cutbased_selection 
        #print signal_selection
        signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9)[0] 

        # fill graph
        point = graph_m1.GetN()
        graph_m1.SetPoint(point, signal_v2, signal_yields)
        #graph_m1.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
        graph_m1.SetPointError(point, 0, 0, 0, 0)

      graph_m1.SetMarkerStyle(20)
      graph_m1.SetMarkerSize(2)
      graph_m1.SetMarkerColor(signal_colour)
      graph_m1.SetLineStyle(1)
      graph_m1.SetLineWidth(2)
      graph_m1.SetLineColor(signal_colour)

      graph_m1.Draw('PL same')  
      legend.AddEntry(graph_m1, 'signal, m={}GeV'.format(signal_mass))

      # background yields
      background_selection = selection + ' && ' + category.definition_flat + '&& ' + category.cutbased_selection
      if doABCD: background_yields, background_err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, background_err =  ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, background_err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=signal_mass, resolution=signal_resolution)

      # for the moment, taking 30% uncertainty
      background_err = 0.3 * background_yields

      lumi_true = 0.
      for data_file in self.data_files:
        lumi_true += data_file.lumi

      background_yields = background_yields * lumi/lumi_true
      background_err = background_err * lumi/lumi_true

      bkg_box_m1 = ROOT.TBox(1e-6, background_yields-background_err, 1, background_yields+background_err)
      bkg_box_m1.SetFillColor(signal_colour)
      bkg_box_m1.SetFillStyle(3005)
      bkg_box_m1.Draw('same')

      bkg_line_m1 = ROOT.TLine(1e-6, background_yields, 1, background_yields)
      bkg_line_m1.SetLineColor(signal_colour)
      bkg_line_m1.SetLineWidth(3)
      bkg_line_m1.SetLineStyle(9)

      bkg_line_m1.Draw('same')
      legend.AddEntry(bkg_line_m1, 'background in 2#sigma window around m={}GeV'.format(signal_mass))


      graph_m3 = ROOT.TGraphAsymmErrors()
      #for signal_file in signal_files_:
      for signal_file in samples_m3:
      #for signal_file in self.signal_files:
        #print signal_file.filename
        # get signal coupling
        signal_mass = signal_file.mass
        #if float(signal_mass) != float(mass): continue
        signal_ctau = signal_file.ctau
        signal_resolution = signal_file.resolution
        signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
        signal_colour = signal_file.colour
        print '{}Gev {}mm'.format(signal_mass, signal_ctau)

        # compute the signal yields
        signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
        signal_selection += ' && ' + category.definition_flat + '&& ' + category.cutbased_selection 
        #print signal_selection
        signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9)[0] 

        # fill graph
        point = graph_m3.GetN()
        graph_m3.SetPoint(point, signal_v2, signal_yields)
        #graph_m3.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
        graph_m3.SetPointError(point, 0, 0, 0, 0)

      graph_m3.SetMarkerStyle(20)
      graph_m3.SetMarkerSize(2)
      graph_m3.SetMarkerColor(signal_colour)
      graph_m3.SetLineStyle(1)
      graph_m3.SetLineWidth(2)
      graph_m3.SetLineColor(signal_colour)

      graph_m3.Draw('PL same')  
      legend.AddEntry(graph_m3, 'signal, m={}GeV'.format(signal_mass))

      # background yields
      background_selection = selection + ' && ' + category.definition_flat + '&& ' + category.cutbased_selection
      if doABCD: background_yields, background_err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, background_err =  ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, background_err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=signal_mass, resolution=signal_resolution)

      lumi_true = 0.
      for data_file in self.data_files:
        lumi_true += data_file.lumi

      # for the moment, taking 30% uncertainty
      background_err = 0.3 * background_yields

      background_yields = background_yields * lumi/lumi_true
      background_err = background_err * lumi/lumi_true

      bkg_box_m3 = ROOT.TBox(1e-6, background_yields-background_err, 1, background_yields+background_err)
      bkg_box_m3.SetFillColor(signal_colour)
      bkg_box_m3.SetFillStyle(3005)
      bkg_box_m3.Draw('same')

      bkg_line_m3 = ROOT.TLine(1e-6, background_yields, 1, background_yields)
      bkg_line_m3.SetLineColor(signal_colour)
      bkg_line_m3.SetLineWidth(3)
      bkg_line_m3.SetLineStyle(9)

      bkg_line_m3.Draw('same')
      legend.AddEntry(bkg_line_m3, 'background in 2#sigma window around m={}GeV'.format(signal_mass))


      graph_m4p5 = ROOT.TGraphAsymmErrors()
      #for signal_file in signal_files_:
      for signal_file in samples_m4p5:
      #for signal_file in self.signal_files:
        #print signal_file.filename
        # get signal coupling
        signal_mass = signal_file.mass
        #if float(signal_mass) != float(mass): continue
        signal_ctau = signal_file.ctau
        signal_resolution = signal_file.resolution
        signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
        signal_colour = signal_file.colour
        print '{}Gev {}mm'.format(signal_mass, signal_ctau)

        # compute the signal yields
        signal_selection = 'ismatched==1' if selection=='' else 'ismatched==1 && {}'.format(selection)
        signal_selection += ' && ' + category.definition_flat + '&& ' + category.cutbased_selection 
        #print signal_selection
        signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9)[0] 

        # fill graph
        point = graph_m4p5.GetN()
        graph_m4p5.SetPoint(point, signal_v2, signal_yields)
        #graph_m4p5.SetPointError(point, 0, 0, err_signal_yields, err_signal_yields)
        graph_m4p5.SetPointError(point, 0, 0, 0, 0)

      graph_m4p5.SetMarkerStyle(20)
      graph_m4p5.SetMarkerSize(2)
      graph_m4p5.SetMarkerColor(signal_colour)
      graph_m4p5.SetLineStyle(1)
      graph_m4p5.SetLineWidth(2)
      graph_m4p5.SetLineColor(signal_colour)

      graph_m4p5.Draw('PL same')  
      legend.AddEntry(graph_m4p5, 'signal, m={}GeV'.format(signal_mass))

      # background yields
      background_selection = selection + ' && ' + category.definition_flat + '&& ' + category.cutbased_selection
      if doABCD: background_yields, background_err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, background_err =  ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=signal_mass, resolution=signal_resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, background_err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=signal_mass, resolution=signal_resolution)

      lumi_true = 0.
      for data_file in self.data_files:
        lumi_true += data_file.lumi

      # for the moment, taking 30% uncertainty
      background_err = 0.3 * background_yields

      background_yields = background_yields * lumi/lumi_true
      background_err = background_err * lumi/lumi_true

      bkg_box_m4p5 = ROOT.TBox(1e-6, background_yields-background_err, 1, background_yields+background_err)
      bkg_box_m4p5.SetFillColor(signal_colour)
      bkg_box_m4p5.SetFillStyle(3005)
      bkg_box_m4p5.Draw('same')

      bkg_line_m4p5 = ROOT.TLine(1e-6, background_yields, 1, background_yields)
      bkg_line_m4p5.SetLineColor(signal_colour)
      bkg_line_m4p5.SetLineWidth(3)
      bkg_line_m4p5.SetLineStyle(9)

      bkg_line_m4p5.Draw('same')
      legend.AddEntry(bkg_line_m4p5, 'background in 2#sigma window around m={}GeV'.format(signal_mass))

      legend.Draw()
      self.tools.printLatexBox(0.69, 0.86, category.title, size=0.04)

      if not path.exists('./myPlots/yields'):
        os.system('mkdir -p ./myPlots/yields')

      canv.cd()
      canv.SaveAs('./myPlots/yields/signalbkg_yields_{}_{}.png'.format(label, category.label))
      canv.SaveAs('./myPlots/yields/signalbkg_yields_{}_{}.pdf'.format(label, category.label))


  def plotBkgYields(self, label='', lumi=41.6, doABCD=False, doABCDHybrid=False, doTF=False):
    hist_m1 = ROOT.TH1D('hist_m1', 'hist_m1', len(self.categories)-1, 0, len(self.categories)-1)
    hist_m3 = ROOT.TH1D('hist_m3', 'hist_m3', len(self.categories)-1, 0, len(self.categories)-1)
    hist_m4p5 = ROOT.TH1D('hist_m4p5', 'hist_m4p5', len(self.categories)-1, 0, len(self.categories)-1)
    canv = self.tools.createTCanvas('canv', dimx=1500, dimy=1000)
    canv.SetGrid()
    canv_log = self.tools.createTCanvas('canv_log', dimx=1500, dimy=1000)
    canv_log.SetLogy()
    canv_log.SetGrid()

    ROOT.gStyle.SetPadRightMargin(0.16)
    ROOT.gStyle.SetPadLeftMargin(0.16)

    lumi_true = 0.
    for data_file in self.data_files:
      lumi_true += data_file.lumi
    print lumi_true

    for icat, category in enumerate(self.categories):
      if 'incl' in category.label: continue

      background_selection = self.selection + ' && ' + category.definition_flat + '&& ' + category.cutbased_selection

      signal_file_m1 = signal_samples['central_V09_06Nov21_benchmark'][2]
      mass = signal_file_m1.mass
      resolution = signal_file_m1.resolution
      if doABCD: background_yields, err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, err = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=mass, resolution=resolution)
      background_yields = background_yields * lumi/lumi_true
      err = err * lumi/lumi_true
      if background_yields == 0.: 
        background_yields = 1e-1
      print 'mass 1 {} {} +- {}'.format(category.label, int(background_yields), int(err))
      hist_m1.Fill(category.label, 0)
      hist_m1.SetBinContent(icat, background_yields)
      hist_m1.SetBinError(icat, err)

      signal_file_m3 = signal_samples['central_V09_06Nov21_benchmark'][1]
      mass = signal_file_m3.mass
      resolution = signal_file_m3.resolution
      if doABCD: background_yields, err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, err = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=mass, resolution=resolution)
      background_yields = background_yields * lumi/lumi_true
      err = err * lumi/lumi_true
      if background_yields == 0.: 
        background_yields = 1e-1
      print 'mass 3 {} {} +- {}'.format(category.label, int(background_yields), int(err))
      hist_m3.Fill(category.label, 0)
      hist_m3.SetBinContent(icat, background_yields)
      hist_m3.SetBinError(icat, err)

      signal_file_m4p5 = signal_samples['central_V09_06Nov21_benchmark'][0]
      mass = signal_file_m4p5.mass
      resolution = signal_file_m4p5.resolution
      background_yields = 1.
      err = 1.
      if doABCD: background_yields, err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doABCDHybrid: background_yields, err = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      elif doTF: background_yields, err, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, white_list=self.white_list, selection=background_selection).computeBkgYieldsFromMC(mass=mass, resolution=resolution)
      background_yields = background_yields * lumi/lumi_true
      err = err * lumi/lumi_true
      if background_yields == 0.: 
        background_yields = 1e-1
      print 'mass 4p5 {} {} +- {}'.format(category.label, int(background_yields), int(err))
      hist_m4p5.Fill(category.label, 0)
      hist_m4p5.SetBinContent(icat, background_yields)
      hist_m4p5.SetBinError(icat, err)

    hist_frame = ROOT.TH1D('hist_frame', 'hist_frame', len(self.categories)-1, 0, len(self.categories)-1)
    for icat, category in enumerate(self.categories):
      if 'incl' in category.label: continue
      hist_frame.Fill(category.label, 0)
    hist_frame.SetTitle('')
    hist_frame.GetYaxis().SetTitle('Background yields')
    hist_frame.GetYaxis().SetLabelSize(0.029)
    range_min = min(hist_m1.GetMinimum(), hist_m3.GetMinimum(), hist_m4p5.GetMinimum())
    if range_min == 0.: range_min = 1e-1
    range_max = max(hist_m1.GetMaximum(), hist_m3.GetMaximum(), hist_m4p5.GetMaximum())
    range_max += 0.6*range_max
    #hist_frame.GetYaxis().SetRangeUser(range_min, range_max)
    hist_frame.GetYaxis().SetRangeUser(9e-02, 1e06)

    hist_m1.SetMarkerStyle(20)
    hist_m1.SetMarkerSize(2)
    hist_m1.SetMarkerColor(ROOT.kOrange+0)
    hist_m1.SetLineWidth(2)
    hist_m1.SetLineColor(ROOT.kOrange+0)

    hist_m3.SetMarkerStyle(20)
    hist_m3.SetMarkerSize(2)
    hist_m3.SetMarkerColor(ROOT.kRed+1)
    hist_m3.SetLineWidth(2)
    hist_m3.SetLineColor(ROOT.kRed+1)
    #hist_m3.GetYaxis().SetRangeUser(1e-01, 1e06)

    hist_m4p5.SetMarkerStyle(20)
    hist_m4p5.SetMarkerSize(2)
    hist_m4p5.SetMarkerColor(ROOT.kRed+4)
    hist_m4p5.SetLineWidth(2)
    hist_m4p5.SetLineColor(ROOT.kRed+4)

    leg = self.tools.getRootTLegend(xmin=0.54, ymin=0.65, xmax=0.89, ymax=0.89, size=0.027)
    leg.SetFillColorAlpha(0, 0)
    leg.AddEntry(hist_m1, '2#sigma window around 1 GeV')
    leg.AddEntry(hist_m3, '2#sigma window around 3 GeV')
    leg.AddEntry(hist_m4p5, '2#sigma window around 4.5 GeV')

    ROOT.gStyle.SetOptStat(0)

    canv.cd()
    hist_frame.Draw()
    hist_m3.Draw('PE0 same')
    hist_m1.Draw('PE0 same')
    hist_m4p5.Draw('PE0 same')
    leg.Draw()

    if not path.exists('myPlots/yields'):
      os.system('mkdir -p myPlots/yields')

    #label = 'bkg_yields_V08_29Sep21_partial_scaledtolumi_sel13Oct21_nocatsel_ABCDcos2dsvprob_dsa_nonmatched_allmasses'

    canv.SaveAs('myPlots/yields/bkg_yields_{}.png'.format(label))

    canv_log.cd()
    hist_frame.Draw()
    hist_m3.Draw('PE0 same')
    hist_m1.Draw('PE0 same')
    hist_m4p5.Draw('PE0 same')
    leg.Draw()

    canv_log.SaveAs('myPlots/yields/bkg_yields_{}_log.png'.format(label))


  def getSignificance(self, lumi=41.6):

    significances = []  

    # compute the significances
    for category in self.categories:
      print category.label

      masses = [1.0, 1.5, 2.0, 3.0, 4.5]
      resolutions = [0.00853, 0.0126, 0.0164, 0.0237, 0.0377]

      for imass, mass in enumerate(masses):
        print 'mass {} {}'.format(mass, resolutions[imass])
        background_selection = self.selection + ' && ' + category.definition_flat# + '&& ' + category.cutbased_selection
        #background_yields, err, _, _, _ = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolutions[imass], ABCD_regions=self.ABCD_regions)

        from qcd_white_list import white_list
        self.white_list = white_list['20to300']
        background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=background_selection, white_list=self.white_list).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolutions[imass], ABCD_regions=self.ABCD_regions)[0]

        lumi_true = 0.
        for data_file in self.data_files:
          lumi_true += data_file.lumi

        # scale yields to lumi
        background_yields = background_yields * lumi/lumi_true

        for signal_file in self.signal_files:
          if signal_file.mass != mass: continue
          # compute the signal yields
          signal_selection = 'ismatched==1 && ' + self.selection + ' && ' + category.definition_flat

          signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9) 
          #print '{} {} {}'.format(signal_file.mass, signal_file.ctau, int(signal_yields))

          significance = signal_yields / math.sqrt(background_yields)
          significances.append(Significance(category=category.label, mass=signal_file.mass, ctau=signal_file.ctau, significance=significance))

    # make the ranking
    significances.sort(key = lambda x: float(x.significance), reverse=True)
    
    ## per category
    print '\nRanking per category'
    for category in self.categories:
      print '\nCategory {}'.format(category.label)
      for item in significances:
        if item.category==category.label: print '{} {} {}'.format(item.mass, item.ctau, round(item.significance, 2))

    ## per signal point
    print '\nRanking per signal point'
    for signal_file in self.signal_files:
      print '\nSignal point m{} ctau{}'.format(signal_file.mass, signal_file.ctau)
      for item in significances:
        if item.mass==signal_file.mass and item.ctau==signal_file.ctau: print '{} {}'.format(item.category, round(item.significance, 2))
        if item.mass==signal_file.mass and item.ctau==signal_file.ctau: print '{} & {} \\\ '.format(item.category, round(item.significance, 2))



      
      
          

class Significance(object):
  def __init__(self, category, mass, ctau, significance):
    self.category = category
    self.mass = mass
    self.ctau = ctau 
    self.significance = significance



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  plotSigYields = False
  plotSigBkgYields = False
  plotBkgYields = False
  getSignificance = False
  printYieldsTable = True

  if printYieldsTable:
    plotter = YieldsChecks()
    plotter.writeYieldsTable()

  if plotSigYields:
    data_files = data_samples['V10_30Dec21']
    qcd_files = qcd_samples['V10_30Dec21']
    white_list = white_list['20to300']
    #categories = categories['categories_0_50_150_Bc']
    categories = categories['V12_08Aug22_permass']
    ABCD = ABCD_regions['cos2d_svprob_0p996']
    #doABCD = False
    #doABCDHybrid = False
    #doTF = True

    #selection = selection['study_Nov21'].flat #'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge==mu_charge'
    #selection = 'mu_isdsa==1 && sv_lxy>20 && pi_pt>1.1 && mu_ismatchedtoslimmedmuon==0'
    #selection = selection['baseline_30Dec21'].flat
    selection = selection['baseline_08Aug22'].flat

    #plotter = YieldsChecks()
    plotter = YieldsChecks(data_files=data_files, qcd_files=qcd_files, white_list=white_list, categories=categories, selection=selection, ABCD_regions=ABCD)

    #label = 'central_reweighting_allmasses_strategyctaumerged_updatedsigma'
    #label = 'V08_29Sep21_dsa_study_ABCD_m1_nonmatched'
    label = 'test'
    for category in categories:
      plotter.plotSigYields(lumi=41.6, selection=selection, label=label, outdirlabel='', category=category)


  if plotSigBkgYields:
    data_files = data_samples['V09_06Nov21']
    qcd_files = qcd_samples['V09_06Nov21']
    samples_m1 = signal_samples['central_V09_06Nov21_m1']
    samples_m3 = signal_samples['central_V09_06Nov21_m3']
    samples_m4p5 = signal_samples['central_V09_06Nov21_m4p5']
    signal_files = [samples_m1, samples_m3]
    white_list = white_list['20to300']
    categories = categories['3cat_0_1_5_significance']
    #categories = categories['inclusive']
    ABCD = ABCD_regions['cos2d_svprob_0p996']
    doABCD = False
    doABCDHybrid = True
    doTF = False
    selection = selection['study_Nov21'].flat
    label = 'V09_06Nov21_ABCDHybrid_uncert0p3'

    plotter = YieldsChecks(data_files=data_files, qcd_files=qcd_files, signal_files=signal_files, white_list=white_list, categories=categories, selection=selection, ABCD_regions=ABCD)
    #for category in categories:
    #  label += '_{}'.format(category.label)
    #plotter.plotSigBkgYields(lumi=41.6, selection=selection, label=label, category=category, doABCD=doABCD, doABCDHybrid=doABCDHybrid, doTF=doTF)
    plotter.plotSigBkgYields(lumi=41.6, selection=selection, label=label, doABCD=doABCD, doABCDHybrid=doABCDHybrid, doTF=doTF)


  if plotBkgYields:
    data_files = data_samples['V09_06Nov21']
    qcd_files = qcd_samples['V09_06Nov21']
    white_list = white_list['50to300']
    categories = categories['3cat_0_1_5_significance']
    selection = selection['study_Nov21'].flat
    ABCD = ABCD_regions['cos2d_svprob_0p996']
    doABCD = True
    doABCDHybrid = False
    doTF = False

    label = 'test' #V09_06Nov21_fulllumi_selstudyNov21_cos2dsvprob0p996_catselsignificance_20sigmawindow_qcd50to300_fullA_TF'
    #label = 'V09_06Nov21_fulllumi_selstudyNov21_cos2dsvprob0p996_fullA_nocatsel'
    plotter = YieldsChecks(data_files=data_files, qcd_files=qcd_files, white_list=white_list, categories=categories, selection=selection, ABCD_regions=ABCD)
    plotter.plotBkgYields(label=label, doABCD=doABCD, doABCDHybrid=doABCDHybrid, doTF=doTF)

  if getSignificance:
    signal_files = signal_samples['central_V09_06Nov21_m1'] + signal_samples['central_V09_06Nov21_m1p5'] + signal_samples['central_V09_06Nov21_m2'] + signal_samples['central_V09_06Nov21_m3'] + signal_samples['central_V09_06Nov21_m4p5']
    data_files = data_samples['V09_06Nov21']
    qcd_files = qcd_samples['V09_06Nov21']
    categories = categories['standard']
    selection = selection['study_Nov21'].flat + ' && hnl_charge==0'
    ABCD = ABCD_regions['cos2d_svprob_0p996']

    plotter = YieldsChecks(data_files=data_files, qcd_files=qcd_files, signal_files=signal_files, categories=categories, selection=selection, ABCD_regions=ABCD)
    plotter.getSignificance(lumi=41.6)



