import os
import sys
import numpy as np
from os import path
from tools import Tools
from compute_yields import ComputeYields
sys.path.append('../objects')
from samples import signal_samples, data_samples, qcd_samples
from categories import categories
from ABCD_regions import ABCD_regions
from baseline_selection import selection

import ROOT

class YieldsChecks(Tools):
  def __init__(self, data_file='', qcd_files='', categories='', selection='', ABCD_regions=''):
    self.tools = Tools()
    self.data_file = data_file
    self.qcd_files = qcd_files
    self.categories = categories
    self.selection = selection
    self.ABCD_regions = ABCD_regions


  def plotSigYields(self, lumi=0.774, selection='', label='', outdirlabel=''):
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

    #samples_m1 = signal_samples['limits_m1_29Sep21']
    #samples_m3 = signal_samples['limits_m3_29Sep21']
    #samples_m4p5 = signal_samples['limits_m4p5_29Sep21']

    samples_m1 = signal_samples['central_V09_06Nov21_m1_large']
    samples_m3 = signal_samples['central_V09_06Nov21_m3_large']
    samples_m4p5 = signal_samples['central_V09_06Nov21_m4p5_large']

    graph_dummy = ROOT.TGraph()
    graph_m1 = ROOT.TGraphAsymmErrors()
    graph_m3 = ROOT.TGraphAsymmErrors()
    graph_m4p5 = ROOT.TGraphAsymmErrors()

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
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9) 

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
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9) 

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
      signal_yields, err_signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=lumi, sigma_B=472.8e9) 

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
    legend.AddEntry(graph_m1, 'm=1GeV')
    legend.AddEntry(graph_gen_m1, 'm=1GeV, gen')
    legend.AddEntry(graph_m3, 'm=3GeV, reco')
    #legend.AddEntry(graph_m3, 'm=3GeV')
    legend.AddEntry(graph_gen_m3, 'm=3GeV, gen')
    legend.AddEntry(graph_m4p5, 'm=4p5GeV, reco')
    #legend.AddEntry(graph_m4p5, 'm=4p5GeV')
    legend.AddEntry(graph_gen_m4p5, 'm=4p5GeV, gen')
    legend.Draw()

    if not path.exists('./myPlots/yields'):
      os.system('mkdir -p ./myPlots/yields')

    canv.SaveAs('./myPlots/yields/signal_yields_{}.png'.format(label))
    canv.SaveAs('./myPlots/yields/signal_yields_{}.pdf'.format(label))


  def plotBkgYields(self, label):
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

    for icat, category in enumerate(self.categories):
      if 'incl' in category.label: continue

      background_selection = self.selection + ' && ' + category.definition_flat# + '&& ' + category.cutbased_selection

      signal_file_m1 = signal_samples['central_V09_06Nov21_benchmark'][2]
      mass = signal_file_m1.mass
      resolution = signal_file_m1.resolution
      #background_yields, err, _, _, _ = ComputeYields(data_file=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields, err = ComputeYields(data_files=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields = background_yields * 41.6/(0.774) #* 0.108)
      err = err * 41.6/(0.774) #* 0.108)
      if background_yields == 0.: 
        background_yields = 1e-1
      print 'mass 1 {} {} +- {}'.format(category.label, int(background_yields), int(err))
      hist_m1.Fill(category.label, 0)
      hist_m1.SetBinContent(icat, background_yields)
      hist_m1.SetBinError(icat, err)

      signal_file_m3 = signal_samples['central_V09_06Nov21_benchmark'][1]
      mass = signal_file_m3.mass
      resolution = signal_file_m3.resolution
      #background_yields, err, _, _, _ = ComputeYields(data_file=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields, err = ComputeYields(data_files=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields = background_yields * 41.6/(0.774) # * 0.108)
      err = err * 41.6/(0.774) # * 0.108)
      if background_yields == 0.: 
        background_yields = 1e-1
      print 'mass 3 {} {} +- {}'.format(category.label, int(background_yields), int(err))
      hist_m3.Fill(category.label, 0)
      hist_m3.SetBinContent(icat, background_yields)
      hist_m3.SetBinError(icat, err)

      signal_file_m4p5 = signal_samples['central_V09_06Nov21_benchmark'][0]
      mass = signal_file_m4p5.mass
      resolution = signal_file_m4p5.resolution
      #background_yields, err, _, _, _ = ComputeYields(data_file=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDData(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields, err = ComputeYields(data_files=self.data_file, qcd_files=self.qcd_files, selection=background_selection).computeBkgYieldsFromABCDHybrid(mass=mass, resolution=resolution, ABCD_regions=self.ABCD_regions)
      background_yields = background_yields * 41.6/(0.774) # * 0.108)
      err = err * 41.6/(0.774) #* 0.108)
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
    hist_frame.GetYaxis().SetRangeUser(range_min, range_max)

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
    #hist_m3.GetYaxis().SetRangeUser(0, 50000)

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



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  plotSigYields = False
  plotBkgYields = True

  if plotSigYields:
    plotter = YieldsChecks()
    selection = selection['study_Nov21'].flat #'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge==mu_charge'
    label = 'central_reweighting_allmasses_strategy3_updatedsigma'
    plotter.plotSigYields(lumi=41.6, selection=selection, label=label, outdirlabel='')

  if plotBkgYields:
    data_file = data_samples['V09_06Nov21']
    qcd_files = qcd_samples['V09_06Nov21']
    categories = categories['study_Nov21']
    selection = selection['study_Nov21'].flat
    ABCD = ABCD_regions['bmass_hnlcharge']

    #label = 'V09_06Nov21_scaledtolumi_sel13Oct21_ABCDcos2dsvprob_withcatsel_allmasses'
    label = 'test'
    plotter = YieldsChecks(data_file=data_file, categories=categories, selection=selection, ABCD_regions=ABCD)
    plotter.plotBkgYields(label=label)





