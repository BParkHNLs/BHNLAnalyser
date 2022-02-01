import os
from os import path
import ROOT

from tools import Tools
from compute_yields import ComputeYields

import sys
sys.path.append('../objects')
from quantity import Quantity, quantities
from samples import data_samples, qcd_samples, signal_samples
from categories import categories
from baseline_selection import selection
from qcd_white_list import white_list


class PrefitPlotter(Tools):
  def __init__(self, data_files='', qcd_files='', signal_files='', white_list='', window_size=2, nbins=25):
    self.tools = Tools()
    self.data_files = data_files
    self.qcd_files = qcd_files
    self.signal_files = signal_files
    self.white_list = white_list
    self.window_size = window_size
    self.nbins = nbins


  def getMaxRangeY(self, hist1, hist2, do_log=False, use_sig=False):
    margin = 0.15 if do_log==False else 0.5
    if use_sig: 
      max_hist1 = [hist.GetMaximum() for hist in hist1]
      the_max_hist1 = max(max_hist1)
      max_range = the_max_hist1+margin*the_max_hist1 if the_max_hist1>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    else:
      max_range = hist1.GetMaximum()+margin*hist1.GetMaximum() if hist1.GetMaximum()>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    return max_range


  def getBinnedPrefitMass(self, selection='', use_qcd=False, use_sigyield_fromdatacard=False, lumi_target=41.6, outdirlabel='', subdirlabel='', homedir='', treename='signal_tree', categories='', add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd=''):

    # this will be to avoid memory issues with python
    signal_hists_all = []
    bkg_hists_all = []
    bkg_stack_hists_all = []
    bkg_err_hists_all = []
    frames_all = [] 

    # get signal histograms
    for category in categories:
      print category.label
      signal_hists = []
      for signal_file in self.signal_files:
        f_sig = self.tools.getRootFile(signal_file.filename)
        tree_sig = self.getTree(f_sig, treename)
        
        signal_mass = signal_file.mass
        sigma = signal_file.resolution
        signal_ctau = signal_file.ctau
        signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
        signal_coupling = self.tools.getCouplingLabel(signal_v2)

        # prepare histogram
        quantity_yields = Quantity(name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=self.nbins, bin_min=signal_mass-self.window_size*sigma, bin_max=signal_mass+self.window_size*sigma) 
        hist_signal_name = 'hist_signal_{}_{}'.format(quantity_yields.label, outdirlabel.replace('/', '_'))
        matching_selection = 'ismatched==1'
        selection_signal = matching_selection + ' && ' + selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
        weight_ctau = self.tools.getCtauWeight(signal_file)
        weight_signal = self.tools.getSignalWeight(signal_file=signal_file, sigma_B=472.8e9, lumi=lumi_target)
        if add_weight_hlt : weight_sig += ' * ({})'.format(weight_hlt)
        #if add_weight_pu : weight_sig += ' * (weight_pu_qcd)' #TODO modify pileup weight
    
        if use_sigyield_fromdatacard:
          weight_sig = '({})'.format(weight_ctau)
          hist_signal = self.tools.createHisto(tree_sig, quantity_yields, hist_name=hist_signal_name, branchname='flat', selection=selection_signal, weight=weight_sig)

          # get yields
          file_label = 'bhnl_cat_{}_m_{}_v2_{}'.format(category.label, signal_mass, signal_coupling).replace('-', 'm').replace('.', 'p')
          summary_filename = '{}/outputs/{}/datacards/{}/summary_yields_{}.txt'.format(homedir, outdirlabel, subdirlabel, file_label)
          try:
            summary_file = open(summary_filename)
            lines = summary_file.readlines()
            for line in lines:
              if not 'sig' in line: continue
              signal_yields = float(line[line.find('sig')+4:len(line)])
          except:
            sigma_B = 472.8e9 #TODO 
            signal_yields = ComputeYields(signal_file=signal_file, selection=selection_signal).computeSignalYields(lumi=lumi_target, sigma_B=sigma_B, add_weight_hlt=False, weight_hlt='', isBc=False)[0]

          # restrict to the 2sigma window
          int_signal = hist_signal.Integral(hist_signal.FindBin(signal_mass-2*sigma), hist_signal.FindBin(signal_mass+2*sigma))
          
          hist_signal.Scale(signal_yields/int_signal) 

        else:
          weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
          hist_signal = self.tools.createHisto(tree_sig, quantity_yields, hist_name=hist_signal_name, branchname='flat', selection=selection_signal, weight=weight_sig)

        hist_signal.SetLineWidth(3)
        hist_signal.SetLineColor(signal_file.colour)
        signal_hists.append(hist_signal)

      signal_hists_all.append(signal_hists)

      if use_qcd:
        # get qcd mc histograms
        qcd_hists = []
        hist_qcd_name = 'hist_qcd_tot_{}'.format(category.label)
        hist_qcd_tot = ROOT.TH1D(hist_qcd_name, hist_qcd_name, quantity_yields.nbins, quantity_yields.bin_min, quantity_yields.bin_max)
        hist_qcd_tot.Sumw2()
        hist_qcd_tot.SetDirectory(ROOT.gROOT)
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

          hist_qcd_name = 'hist_qcd_{}_{}'.format(quantity_yields.label, outdirlabel.replace('/', '_'))
          selection_qcd = selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
          hist_qcd = self.tools.createHisto(tree_qcd, quantity_yields, hist_name=hist_qcd_name, branchname='flat', selection=selection_qcd, weight=weight_qcd) 
          hist_qcd.Sumw2()
          hist_qcd.SetFillColor(qcd_file.colour)
          hist_qcd.SetLineColor(1)
          
          int_qcd_tot += hist_qcd.Integral(hist_qcd.FindBin(signal_mass-2*sigma), hist_qcd.FindBin(signal_mass+2*sigma))
      
          hist_qcd_tot.Add(hist_qcd)
          qcd_hists.append(hist_qcd)
      
        hist_qcd_tot.SetTitle('')
        hist_qcd_tot.SetFillColor(ROOT.kAzure-4)
        hist_qcd_tot.SetLineColor(1)
      
        ## create stack histogram  
        hist_qcd_stack = ROOT.THStack('hist_qcd_stack', '')

        ## get the mc normalisation weight
        try:
          summary_file = open(summary_filename)
          lines = summary_file.readlines()
          for line in lines:
            if not 'bkg' in line: continue
            background_yields = float(line[line.find('bkg')+4:len(line)])
        except:
          from ABCD_regions import ABCD_regions
          from qcd_white_list import white_list
          background_yields = ComputeYields(data_files=self.data_files, qcd_files=self.qcd_files, selection=selection_qcd, white_list=white_list['20to300']).computeBkgYieldsFromABCDHybrid(mass=signal_mass, resolution=sigma, ABCD_regions=ABCD_regions['cos2d_svprob_0p996'], sigma_mult_window=2)[0] #TODO modify
          lumi_true = self.tools.getDataLumi(self.data_files)
          background_yields *= lumi_target/lumi_true

        ## normalise to estimated yields
        for hist_qcd in qcd_hists:
          if int_qcd_tot!=0: hist_qcd.Scale(background_yields/int_qcd_tot)
          hist_qcd_stack.Add(hist_qcd)
        if int_qcd_tot!=0: hist_qcd_tot.Scale(background_yields/int_qcd_tot)

        # draw error bars
        hist_qcd_tot_err = hist_qcd_tot.Clone('hist_qcd_tot_err')
        hist_qcd_tot_err.SetDirectory(ROOT.gROOT)
        hist_qcd_tot_err.SetLineWidth(0)
        hist_qcd_tot_err.SetFillStyle(3244)
        hist_qcd_tot_err.SetFillColor(ROOT.kGray+2)

        bkg_hists_all.append(hist_qcd_tot)
        bkg_stack_hists_all.append(hist_qcd_stack)
        bkg_err_hists_all.append(hist_qcd_tot_err)

      else:
        hist_data_name = 'hist_data_tot_{}'.format(category.label)
        hist_data_tot = ROOT.TH1D(hist_data_name, hist_data_name, quantity_yields.nbins, quantity_yields.bin_min, quantity_yields.bin_max)
        hist_data_tot.Sumw2()
        hist_data_stack = ROOT.THStack('hist_data_stack', '')
        hist_data_tot.SetDirectory(ROOT.gROOT)

        for idata, data_file in enumerate(self.data_files):
          f_data = self.tools.getRootFile(data_file.filename)
          tree_data = self.getTree(f_data, treename)
          hist_data_name = 'hist_data_{}_{}_{}'.format(quantity_yields.label, outdirlabel.replace('/', '_'), do_log)
          selection_data = selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection
          hist_data = self.tools.createHisto(tree_data, quantity_yields, hist_name=hist_data_name, branchname='flat', selection=selection_data)
          hist_data.Sumw2()

          hist_data_stack.Add(hist_data)
          hist_data_tot.Add(hist_data)

        hist_data_tot.SetTitle('')
        hist_data_tot.SetFillColor(ROOT.kBlue-3)
        hist_data_tot.SetFillStyle(3005)
        hist_data_tot.SetLineColor(1)

        # draw error bars
        hist_data_tot_err = hist_data_tot.Clone('hist_data_tot_err')
        hist_data_tot_err.SetDirectory(ROOT.gROOT)
        hist_data_tot_err.SetLineWidth(0)
        hist_data_tot_err.SetFillStyle(3244)
        hist_data_tot_err.SetFillColor(ROOT.kGray+2)

        bkg_hists_all.append(hist_data_tot)
        bkg_stack_hists_all.append(hist_data_stack)
        bkg_err_hists_all.append(hist_data_tot_err)

    return signal_hists_all, bkg_hists_all, bkg_stack_hists_all, bkg_err_hists_all


  def plotBinnedPrefitMass(self, selection='', use_qcd=False, use_sigyield_fromdatacard=False, lumi_target=41.6, outdirloc='', outdirlabel='', subdirlabel='', homedir='', treename='signal_tree', categories='', add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd='', do_stack=True, do_log=False, add_CMSlabel=True, CMS_tag='Preliminary'):

    signal_hists_all, bkg_hists_all, bkg_stack_hists_all, bkg_err_hists_all = self.getBinnedPrefitMass(selection, use_qcd, use_sigyield_fromdatacard, lumi_target, outdirlabel, subdirlabel, homedir, treename, categories, add_weight_hlt, add_weight_pu, weight_hlt, weight_puqcd) 

    ## plot each category separately 
    for icat, category in enumerate(categories):
      # create the canvas
      canv_name = 'canv_{}_{}_{}'.format(outdirlabel.replace('/', '_'), do_log, category.label)
      canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
      canv.cd()

      # define the pad
      pad = ROOT.TPad("pad","pad",0.02,0,1,1)
      if do_log: pad.SetLogy()
      pad.Draw()
      canv.cd()

      # get histograms
      signal_hists = signal_hists_all[icat]
      hist_bkg_tot = bkg_hists_all[icat]
      hist_bkg_stack = bkg_hists_all[icat]
      hist_bkg_tot_err = bkg_err_hists_all[icat]

      pad.cd()

      ROOT.gStyle.SetOptStat(0)

      hist_bkg_tot.SetTitle('')
      hist_bkg_tot.GetXaxis().SetTitle('#mu#pi invariant mass [GeV]')
      hist_bkg_tot.GetXaxis().SetLabelSize(0.033)
      hist_bkg_tot.GetXaxis().SetTitleSize(0.042)
      hist_bkg_tot.GetXaxis().SetTitleOffset(1.1)
      hist_bkg_tot.GetYaxis().SetTitle('Yields')
      hist_bkg_tot.GetYaxis().SetLabelSize(0.033)
      hist_bkg_tot.GetYaxis().SetTitleSize(0.042)
      hist_bkg_tot.GetYaxis().SetTitleOffset(1.3)
      hist_bkg_tot.GetYaxis().SetRangeUser(0 if not do_log else 1e-9, self.getMaxRangeY(signal_hists, hist_bkg_stack, do_log=True, use_sig=True) if not do_log else 5e6)

      # draw the distributions
      hist_bkg_tot.Draw('histo')
      if do_stack:
        hist_bkg_stack.Draw('histo same')
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')
      hist_bkg_tot_err.Draw('E2 same')

      # draw the legend
      legend = self.tools.getRootTLegend(xmin=0.15, ymin=0.6, xmax=0.4, ymax=0.8, size=0.033, do_alpha=False)
      for isig, signal_hist in enumerate(signal_hists):
        legend.AddEntry(signal_hist, 'signal - ({}GeV, {}mm)'.format(self.signal_files[isig].mass, self.signal_files[isig].ctau))
      legend.AddEntry(hist_bkg_tot, 'QCD' if use_qcd else 'Data')
      legend.Draw('same')

      # add labels
      self.tools.printLatexBox(0.65, 0.86, category.title, size=0.036)
      if add_CMSlabel: self.tools.printCMSTag(pad, CMS_tag, size=0.43)
      if add_CMSlabel: self.tools.printLumiTag(pad, lumi_target, size=0.43, offset=0.57)

      plotname = 'hnl_mass_{}_window{}sigma_nbins{}_lumi{}_m{}'.format(category.label, self.window_size, self.nbins if not use_qcd else str(self.nbins)+'_qcd', lumi_target, self.signal_files[0].mass).replace('.', 'p')
      #outputdir = self.tools.getOutDir('{}/{}/datacards/{}'.format(outdirloc, outdirlabel, subdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
      outputdir = self.tools.getOutDir('myPlots/{}'.format(outdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
      
      canv.SaveAs('{}/{}.png'.format(outputdir, plotname))
      canv.SaveAs('{}/{}.pdf'.format(outputdir, plotname))


    ## plot categories on same canvas 
    # create the canvas
    canv_name = 'canv_mult_{}_{}'.format(outdirlabel.replace('/', '_'), do_log)
    canv_mult = self.tools.createTCanvas(name=canv_name, dimx=1500, dimy=1000)
    canv_mult.Divide(1, 2)
    
    for icat, category in enumerate(categories):
      if 'incl' in category.label: continue

      pad_min = 0.
      pad_max = 1.
      pad_size = float((pad_max-pad_min)/(len(categories)-1)*2)
      if icat <= len(categories)/2:
        canv_mult.cd(1)
        pad = ROOT.TPad("pad","pad",pad_min+(icat-1)*pad_size,0,pad_min+icat*pad_size,1)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.07)
        pad.SetTopMargin(0.15)
        pad.SetBottomMargin(0.03)
      else: 
        canv_mult.cd(2)
        pad = ROOT.TPad("pad","pad",pad_min+(icat-1-(len(categories)-1)/2)*pad_size,0,pad_min+(icat-(len(categories)-1)/2)*pad_size,1)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.07)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.15)
  
      if do_log: pad.SetLogy()
      pad.Draw()
      pad.cd()

      # get histograms
      signal_hists = signal_hists_all[icat]
      hist_bkg_tot = bkg_hists_all[icat]
      hist_bkg_stack = bkg_hists_all[icat]
      hist_bkg_tot_err = bkg_err_hists_all[icat]

      ROOT.gStyle.SetOptStat(0)

      hist_bkg_tot.SetTitle('')
      hist_bkg_tot.GetXaxis().SetTitle('#mu#pi invariant mass [GeV]' if icat==len(categories)-1 else '')
      hist_bkg_tot.GetXaxis().SetLabelSize(0.06 if icat>(len(categories)-1)/2 else 0)
      hist_bkg_tot.GetXaxis().SetTitleSize(0.07 if icat>(len(categories)-1)/2 else 0)
      hist_bkg_tot.GetXaxis().SetTickLength(0.03)
      hist_bkg_tot.GetXaxis().SetNdivisions(405)
      hist_bkg_tot.GetXaxis().SetTitleOffset(1.1)
      hist_bkg_tot.GetYaxis().SetTitle('Yields' if icat==1 else '')
      hist_bkg_tot.GetYaxis().SetLabelSize(0.06 if (icat==1 or icat-1==(len(categories)-1)/2) else 0)
      hist_bkg_tot.GetYaxis().SetTitleSize(0.07 if (icat==1 or icat-1==(len(categories)-1)/2) else 0)
      hist_bkg_tot.GetYaxis().SetTickLength(0.03)
      hist_bkg_tot.GetYaxis().SetNdivisions(405)
      hist_bkg_tot.GetYaxis().SetTitleOffset(1.3)
      hist_bkg_tot.GetYaxis().SetRangeUser(1e-9, 5e6)

      # draw the distributions
      hist_bkg_tot.Draw('histo')
      if do_stack:
        hist_bkg_stack.Draw('histo same')
      for hist_sig in signal_hists:
        hist_sig.Draw('histo same')
      hist_bkg_tot_err.Draw('E2 same')

      # draw the legend
      if icat==1:
        legend = self.tools.getRootTLegend(xmin=0.26, ymin=0.06, xmax=0.84, ymax=0.37, size=0.04, do_alpha=False)
        for isig, signal_hist in enumerate(signal_hists):
          legend.AddEntry(signal_hist, 'signal - ({}GeV, {}mm)'.format(self.signal_files[isig].mass, self.signal_files[isig].ctau))
        legend.AddEntry(hist_bkg_tot, 'QCD' if use_qcd else 'Data')
        legend.Draw('same')

      # draw 2sigma window
      #line_right = ROOT.TLine(self.signal_files[0].mass-2*self.signal_files[0].resolution, 1e-9, self.signal_files[0].mass-2*self.signal_files[0].resolution, 5e6)
      #line_right.SetLineColor(1)
      #line_right.SetLineWidth(3)
      #line_right.SetLineStyle(9)
      #line_right.Draw('same')

      #line_left = ROOT.TLine(self.signal_files[0].mass+2*self.signal_files[0].resolution, 1e-9, self.signal_files[0].mass+2*self.signal_files[0].resolution, 5e6)
      #line_left.SetLineColor(1)
      #line_left.SetLineWidth(3)
      #line_left.SetLineStyle(9)
      #line_left.Draw('same')

      # add labels
      offset = 1.4 if icat<=(len(categories)-1)/2 else 2.8
      self.tools.printLatexBox(0.73, 1-offset*pad.GetTopMargin(), category.title, size=0.056)
      if add_CMSlabel and icat==1: self.tools.printCMSTag(pad, CMS_tag, size=0.45, offset=0.15)
      if add_CMSlabel and icat==(len(categories)-1)/2: self.tools.printLumiTag(pad, lumi_target, size=0.45, offset=0.34)

    plotname = 'hnl_mass_mult_window{}sigma_nbins{}_lumi{}_m{}'.format(self.window_size, self.nbins if not use_qcd else str(self.nbins)+'_qcd', lumi_target, self.signal_files[0].mass).replace('.', 'p')
    #outputdir = self.tools.getOutDir('{}/{}/datacards/{}'.format(outdirloc, outdirlabel, subdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
    outputdir = self.tools.getOutDir('myPlots/{}'.format(outdirlabel), 'binned_prefit', do_shape=False, do_luminorm=False, do_stack=do_stack, do_log=do_log, add_overflow=False)
    
    canv_mult.cd()
    canv_mult.SaveAs('{}/{}.png'.format(outputdir, plotname))
    canv_mult.SaveAs('{}/{}.pdf'.format(outputdir, plotname))


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  
  homedir = '/t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNLAnalyser'
  data_files = data_samples['V10_30Dec21']
  qcd_files = qcd_samples['V10_30Dec21']
  signal_files = signal_samples['V10_30Dec21_m3']
  categories = categories['3cat_0_1_5_significance']
  baseline_selection = selection['baseline_30Dec21'].flat
  white_list = white_list['20to300']

  add_weight_hlt = False
  weight_hlt = ''

  add_weight_pu = True
  weight_puqcd = 'weight_pu_qcd_D'

  nbins = 30
  window_size = 4
  lumi_target = 5.302 

  outdirlabel = 'V10_30Dec21'
  #subdirlabel = 'test_shape_analysis_v3_2'
  subdirlabel = 'test_shape_analysis_lumiD1_v2'

  do_log = True

  use_qcd = False
  use_sigyield_fromdatacard = False

  plotter = PrefitPlotter(data_files=data_files, qcd_files=qcd_files, signal_files=signal_files, white_list=white_list, window_size=window_size, nbins=nbins)

  region_definition = 'hnl_charge == 0'
  plotter.plotBinnedPrefitMass(selection = baseline_selection + ' && ' + region_definition, 
               use_qcd = use_qcd,
               use_sigyield_fromdatacard = use_sigyield_fromdatacard,
               lumi_target = lumi_target,
               outdirloc = './outputs',
               outdirlabel = outdirlabel, 
               subdirlabel = subdirlabel, 
               homedir = homedir,
               treename = 'signal_tree',
               categories = categories,
               add_weight_hlt = add_weight_hlt,
               add_weight_pu = add_weight_pu,
               weight_hlt = weight_hlt, 
               weight_puqcd = weight_puqcd,
               do_stack = False, 
               do_log = do_log,
               add_CMSlabel = True,
               CMS_tag = 'Preliminary'
               )
