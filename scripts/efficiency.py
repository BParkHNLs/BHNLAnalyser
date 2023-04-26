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
from samples import qcd_samples
from qcd_white_list import white_list as wt


class EfficiencyAnalyser(Tools):
  def __init__(self, filename='', sample_type='flat', process_type='mumupi', add_dsa='yes', matchings=None, cuts=None, displacement_bins=None, pt_bins=None, title='', outdirlabel='default'):
    self.tools = Tools()
    self.filename = filename
    self.sample_type = sample_type
    self.process_type = process_type
    self.add_dsa = add_dsa
    self.matchings = matchings
    self.displacement_bins = displacement_bins
    self.pt_bins = pt_bins
    self.cuts = cuts
    self.title = title
    self.outdirlabel = outdirlabel
    self.outputdir = self.tools.getOutDir('./myPlots/efficiency', self.outdirlabel)

    self.treename = ''
    if self.sample_type == 'nano': self.treename = 'Events'
    else:
      if self.process_type == 'mumupi': self.treename = 'signal_tree'
      else: self.treename = 'hnl_tree'

    if self.sample_type not in ['flat', 'nano']:
      raise RuntimeError("Unknown sample type. Please choose among ['flat', 'nano']")

    if self.process_type not in ['mumupi', 'mupi']:
      raise RuntimeError("Unknown process type. Please choose among ['mumupi', 'mupi']")

    if self.add_dsa not in ['yes', 'no']:
      raise RuntimeError("Unknown add-dsa decision. Please choose among ['yes', 'no']")

    if self.matchings != None:
      for matching in self.matchings:
        if matching not in ['candidate', 'trigger_muon', 'muon', 'pion']:
          raise RuntimeError("Unknown matching strategy. Please choose among ['candidate', 'trigger_muon', 'muon', 'pion']")

    if path.exists('{}/numbers.txt'.format(self.outputdir)):
      os.system('rm {}/numbers.txt'.format(self.outputdir))


  def getEfficiencyFromNano(self, displacement_bins=None, pt_bins=None):
    f = self.tools.getRootFile(self.filename, with_ext=False) 
    tree = self.tools.getTree(f, self.treename)

    if displacement_bins != None and pt_bins != None:
      n_num = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings)))  #[[0] * len(displacement_bins)] * len(pt_bins)
      n_deno = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings))) #[[0] * len(displacement_bins)] * len(pt_bins)
      efficiency = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings))) #[[0] * len(displacement_bins)] * len(pt_bins)
      error = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings))) #[[0] * len(displacement_bins)] * len(pt_bins)
    else:
      n_num = np.zeros(len(self.cuts))
      n_deno = np.zeros(len(self.cuts))
      efficiency = np.zeros(len(self.cuts))
      error = np.zeros(len(self.cuts))

    for entry in tree:
      # only consider events with at least one candidate
      #if entry.nBToMuMuPi == 0: continue

      #print '\n'
      # retrieve candidate matching information
      trgmu_ismatched = 0
      mu_ismatched = 0
      pi_ismatched = 0
      cand_ismatched = 0
      for icand in range(0, entry.nBToMuMuPi):
        if entry.BToMuMuPi_trg_mu_isMatched[icand] == 1: trgmu_ismatched = 1
        if entry.BToMuMuPi_sel_mu_isMatched[icand] == 1: mu_ismatched = 1
        if entry.BToMuMuPi_pi_isMatched[icand] == 1: pi_ismatched = 1
        if entry.BToMuMuPi_isMatched[icand] == 1: cand_ismatched = 1
        #if entry.BToMuMuPi_sel_mu_isMatched[icand] == 1 and entry.BToMuMuPi_pi_isMatched[icand] == 1 and entry.BToMuMuPi_mupi_mass_reco_gen_reldiff[icand] < 0.5: 
        #  mupi_ismatched = 1

      matching_cond = {}
      matching_cond['candidate'] = cand_ismatched==1 
      matching_cond['trigger_muon'] = trgmu_ismatched==1 
      matching_cond['muon'] = mu_ismatched==1 
      matching_cond['pion'] = pi_ismatched==1 

      hnl_idx = -1
      mother_idx = -1
      trgmu_idx = -1
      mu_idx = -1
      pi_idx = -1

      # search for hnl and its mother
      for igen in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[igen]) == 9900015: 
          hnl_idx = igen 
          mother_idx = entry.GenPart_genPartIdxMother[hnl_idx]
          break

      # search for trigger muon and hnl daughters
      for igen in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[igen]) == 13 and entry.GenPart_genPartIdxMother[igen] == mother_idx: trgmu_idx = igen
        if abs(entry.GenPart_pdgId[igen]) == 13 and entry.GenPart_genPartIdxMother[igen] == hnl_idx: mu_idx = igen
        if abs(entry.GenPart_pdgId[igen]) == 211 and entry.GenPart_genPartIdxMother[igen] == hnl_idx: pi_idx = igen

      # remove e-channel
      if mu_idx == -1: continue

      for imatch, matching in enumerate(self.matchings):
        if displacement_bins != None and pt_bins != None:
          # fetch n_deno and n_num with acceptance cuts
          #FIXME make sure to use latest deno cuts
          #if entry.GenPart_pt[mu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
          #   and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
          #   and entry.GenPart_pt[trgmu_idx] > 7.:
          if (entry.GenPart_pt[trgmu_idx] > 7. or entry.GenPart_pt[mu_idx] > 7) and entry.GenPart_pt[mu_idx] > 1.5 and entry.GenPart_pt[trgmu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
             and abs(entry.GenPart_eta[trgmu_idx]) < 2.5 and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5:
          #if entry.GenPart_pt[mu_idx] > 3.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
          #   and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
          #   and entry.GenPart_pt[trgmu_idx] > 9.5:
          #if entry.GenPart_pt[mu_idx] > 0:
               displacement = math.sqrt(pow(entry.GenPart_vx[mu_idx] - entry.GenPart_vx[trgmu_idx], 2) + pow(entry.GenPart_vy[mu_idx] - entry.GenPart_vy[trgmu_idx], 2))
               #if matching == 'candidate': obj_pt = entry.GenPart_pt[mother_idx]
               if matching == 'candidate': obj_pt = entry.GenPart_pt[hnl_idx]
               elif matching == 'mupi_candidate': obj_pt = entry.GenPart_pt[hnl_idx]
               elif matching == 'trigger_muon': obj_pt = entry.GenPart_pt[trgmu_idx]
               elif matching == 'muon': obj_pt = entry.GenPart_pt[mu_idx]
               elif matching == 'pion': obj_pt = entry.GenPart_pt[pi_idx]

               for ibin_disp, displacement_bin in enumerate(displacement_bins):
                 for ibin_pt, pt_bin in enumerate(pt_bins):
                   bin_min_disp, bin_max_disp = displacement_bin
                   bin_min_pt, bin_max_pt = pt_bin

                   if displacement > bin_min_disp and displacement < bin_max_disp and obj_pt > bin_min_pt and obj_pt < bin_max_pt: 
                     n_deno[ibin_disp][ibin_pt][imatch] = n_deno[ibin_disp][ibin_pt][imatch] + 1
                     #print '{} {}'.format(matching, matching_cond[matching])
                     if matching_cond[matching]: n_num[ibin_disp][ibin_pt][imatch] = n_num[ibin_disp][ibin_pt][imatch] + 1
        else:
          for icut, cut in enumerate(self.cuts):
            # fetch n_deno and n_num with acceptance cuts
            #if entry.GenPart_pt[mu_idx] > float(cut) and abs(entry.GenPart_eta[mu_idx]) < 2.5: # \
            if entry.GenPart_pt[mu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
               and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
               and entry.GenPart_pt[trgmu_idx] > float(cut):
            #if (entry.GenPart_pt[trgmu_idx] > 7. or entry.GenPart_pt[mu_idx] > 7) and entry.GenPart_pt[mu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
            #   and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5:
                 n_deno[icut] = n_deno[icut] + 1
                 if matching == 'candidate' and cand_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'trigger_muon' and trgmu_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'muon' and mu_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'pion' and pi_ismatched == 1: n_num[icut] = n_num[icut] + 1

    # compute efficiency
    f = open('{}/numbers.txt'.format(self.outputdir), 'a')
    f.write('\n')
    if displacement_bins != None and pt_bins != None:
      for imatch, matching in enumerate(self.matchings):
        for ibin_disp, displacement_bin in enumerate(displacement_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            efficiency[ibin_disp][ibin_pt][imatch] = float(n_num[ibin_disp][ibin_pt][imatch]) / float(n_deno[ibin_disp][ibin_pt][imatch]) if float(n_deno[ibin_disp][ibin_pt][imatch]) != 0. else 0.
            if n_num[ibin_disp][ibin_pt][imatch] == 0: n_num[ibin_disp][ibin_pt][imatch] = 1e-11
            if float(n_num[ibin_disp][ibin_pt][imatch]) != 0 and float(n_deno[ibin_disp][ibin_pt][imatch]) != 0:
              error[ibin_disp][ibin_pt][imatch] = efficiency[ibin_disp][ibin_pt][imatch] * ( math.sqrt(float(n_num[ibin_disp][ibin_pt][imatch]))/float(n_num[ibin_disp][ibin_pt][imatch])  + math.sqrt(float(n_deno[ibin_disp][ibin_pt][imatch]))/float(n_deno[ibin_disp][ibin_pt][imatch]) )     
            else:
              error[ibin_disp][ibin_pt][imatch] = 0
            # for aesthetics
            if efficiency[ibin_disp][ibin_pt][imatch] == 0.: efficiency[ibin_disp][ibin_pt][imatch] = 1e-9

      #for imatch, matching in enumerate(self.matchings):
      #  print '\n'
      #  for ibin_disp, displacement_bin in enumerate(displacement_bins):
      #    for ibin_pt, pt_bin in enumerate(pt_bins):
      #      print '{} {} {} {} {} {}+-{}'.format(matching, displacement_bin, pt_bin, n_deno[ibin_disp][ibin_pt][imatch], n_num[ibin_disp][ibin_pt][imatch], efficiency[ibin_disp][ibin_pt][imatch], error[ibin_disp][ibin_pt][imatch])
      for imatch, matching in enumerate(self.matchings):
        for ibin_disp, displacement_bin in enumerate(displacement_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            f.write('\n{} {} {} {} {} {}+-{}'.format(matching, displacement_bin, pt_bin, n_deno[ibin_disp][ibin_pt][imatch], n_num[ibin_disp][ibin_pt][imatch], efficiency[ibin_disp][ibin_pt][imatch], error[ibin_disp][ibin_pt][imatch]))
    else:
      for icut, cut in enumerate(self.cuts):
        efficiency[icut] = float(n_num[icut]) / float(n_deno[icut]) if float(n_deno[icut]) != 0. else 0. 
        error[icut] = efficiency[icut] * ( math.sqrt(n_num[icut])/n_num[icut]  + math.sqrt(n_deno[icut])/n_deno[icut] )     
    f.close()

      #for icut, cut in enumerate(self.cuts):
      #  print '{} {} {} {}'.format(cut, n_deno[icut], n_num[icut], efficiency[icut])

    return efficiency, error


  def getEfficiency(self, displacement_bins=None, pt_bins=None):
    f = self.tools.getRootFile(self.filename, with_ext=False) 
    tree = self.tools.getTree(f, self.treename)

    n_num = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings)))
    n_deno = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings)))
    efficiency = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings)))
    error = np.zeros((len(displacement_bins), len(pt_bins), len(self.matchings)))

    for entry in tree:
      # retrieve candidate matching information
      trgmu_ismatched = entry.trgmu_ismatched if self.process_type == 'mumupi' else 0
      mu_ismatched = entry.mu_ismatched
      pi_ismatched = entry.pi_ismatched
      cand_ismatched = entry.ismatched # either mumupi candidate or mupi candidate, depending on the sample (BToMuMuPi vs HNLToMuPi)

      matching_cond = {}
      if self.add_dsa == 'yes':
        matching_cond['candidate'] = cand_ismatched==1 
        matching_cond['trigger_muon'] = trgmu_ismatched==1 
        matching_cond['muon'] = mu_ismatched==1 
        matching_cond['pion'] = pi_ismatched==1 
      else:
        matching_cond['candidate'] = cand_ismatched==1 and entry.mu_isdsa==0 and entry.hnl_charge==0 
        matching_cond['trigger_muon'] = trgmu_ismatched==1 and entry.mu_isdsa==0 and entry.hnl_charge==0
        matching_cond['muon'] = mu_ismatched==1 and entry.mu_isdsa==0 and entry.hnl_charge==0
        matching_cond['pion'] = pi_ismatched==1 and entry.mu_isdsa==0 and entry.hnl_charge==0
      
      #hnl_idx = -1
      #mother_idx = -1
      #trgmu_idx = -1
      #mu_idx = -1
      #pi_idx = -1
      
      for imatch, matching in enumerate(self.matchings):
        # fetch n_deno and n_num with acceptance cuts
        #FIXME update the following cuts
        if entry.gen_mu_pt > 1. and abs(entry.gen_mu_eta) < 2.5 \
           and entry.gen_pi_pt > 1. and abs(entry.gen_pi_eta) < 2.5:
             #and entry.gen_trgmu_pt > 7.:
             displacement = entry.gen_hnl_lxy
             if matching == 'candidate': obj_pt = entry.gen_hnl_pt if self.process_type == 'mupi' else entry.gen_b_pt
             elif matching == 'trigger_muon': obj_pt = entry.gen_trgmu_pt
             elif matching == 'muon': obj_pt = entry.gen_mu_pt
             elif matching == 'pion': obj_pt = entry.gen_pi_pt

             for ibin_disp, displacement_bin in enumerate(displacement_bins):
               for ibin_pt, pt_bin in enumerate(pt_bins):
                 bin_min_disp, bin_max_disp = displacement_bin
                 bin_min_pt, bin_max_pt = pt_bin

                 if displacement > bin_min_disp and displacement < bin_max_disp and obj_pt > bin_min_pt and obj_pt < bin_max_pt: 
                   n_deno[ibin_disp][ibin_pt][imatch] = n_deno[ibin_disp][ibin_pt][imatch] + 1
                   if matching_cond[matching]: n_num[ibin_disp][ibin_pt][imatch] = n_num[ibin_disp][ibin_pt][imatch] + 1
                   
    # compute efficiency
    f = open('{}/numbers.txt'.format(self.outputdir), 'a')
    f.write('\n')
    if displacement_bins != None and pt_bins != None:
      for imatch, matching in enumerate(self.matchings):
        for ibin_disp, displacement_bin in enumerate(displacement_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            efficiency[ibin_disp][ibin_pt][imatch] = float(n_num[ibin_disp][ibin_pt][imatch]) / float(n_deno[ibin_disp][ibin_pt][imatch]) if float(n_deno[ibin_disp][ibin_pt][imatch]) != 0. else 0.
            if n_num[ibin_disp][ibin_pt][imatch] == 0: n_num[ibin_disp][ibin_pt][imatch] = 1e-11
            if float(n_num[ibin_disp][ibin_pt][imatch]) != 0 and float(n_deno[ibin_disp][ibin_pt][imatch]) != 0:
              error[ibin_disp][ibin_pt][imatch] = efficiency[ibin_disp][ibin_pt][imatch] * ( math.sqrt(float(n_num[ibin_disp][ibin_pt][imatch]))/float(n_num[ibin_disp][ibin_pt][imatch])  + math.sqrt(float(n_deno[ibin_disp][ibin_pt][imatch]))/float(n_deno[ibin_disp][ibin_pt][imatch]) )     
            else:
              error[ibin_disp][ibin_pt][imatch] = 0
            # for aesthetics
            if efficiency[ibin_disp][ibin_pt][imatch] == 0.: efficiency[ibin_disp][ibin_pt][imatch] = 1e-9

      for imatch, matching in enumerate(self.matchings):
        for ibin_disp, displacement_bin in enumerate(displacement_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            f.write('\n{} {} {} {} {} {}+-{}'.format(matching, displacement_bin, pt_bin, n_deno[ibin_disp][ibin_pt][imatch], n_num[ibin_disp][ibin_pt][imatch], efficiency[ibin_disp][ibin_pt][imatch], error[ibin_disp][ibin_pt][imatch]))
    else:
      for icut, cut in enumerate(self.cuts):
        efficiency[icut] = float(n_num[icut]) / float(n_deno[icut]) if float(n_deno[icut]) != 0. else 0. 
        error[icut] = efficiency[icut] * ( math.sqrt(n_num[icut])/n_num[icut]  + math.sqrt(n_deno[icut])/n_deno[icut] )     
    f.close()

    return efficiency, error


  def plot2DEfficiency(self):
    gStyle.SetPadRightMargin(0.16)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetOptStat(0)
    gStyle.SetPaintTextFormat(".2f")

    efficiency = -99.
    error = -99.
    if self.sample_type == 'flat':
      efficiency, error = self.getEfficiency(displacement_bins=self.displacement_bins, pt_bins=self.pt_bins)
    else:
      efficiency, error = self.getEfficiencyFromNano(displacement_bins=self.displacement_bins, pt_bins=self.pt_bins)

    for imatch, matching in enumerate(self.matchings):
      canv_name = '2d_canv_{}_{}'.format(matching, self.outdirlabel)
      canv = self.tools.createTCanvas(canv_name, 900, 800)

      bin_min_disp, a = self.displacement_bins[0]
      b ,bin_max_disp = self.displacement_bins[len(self.displacement_bins)-1]
      step_disp = (bin_max_disp-bin_min_disp) / (len(self.displacement_bins)-1)

      hist_name = 'hist_{}_{}'.format(matching, self.outdirlabel)
      hist = ROOT.TH2D(hist_name, hist_name, len(self.displacement_bins), 1, 100, len(self.pt_bins), 1, 100)
      for ibin_disp, displacement_bin in enumerate(self.displacement_bins):
        for ibin_pt, pt_bin in enumerate(self.pt_bins):
          hist.Fill(str(displacement_bin), str(pt_bin), efficiency[ibin_disp][ibin_pt][imatch])

      hist.Draw('colz')
      hist.SetMarkerSize(2)
      hist.Draw('text' +'same')

      hist.SetTitle(self.title)
      hist.GetXaxis().SetTitle('l_{xy} [cm]')
      hist.GetXaxis().SetLabelSize(0.037)
      hist.GetXaxis().SetTitleSize(0.042)
      hist.GetXaxis().SetTitleOffset(1.1)
      if matching == 'candidate': ylabel = '#mu#pi pT [GeV]' #if self.process_type == 'mupi' else '#mu#mu#pi pT [GeV]'
      elif matching == 'trigger_muon': ylabel = 'trigger #mu pT [GeV]'
      elif matching == 'muon': ylabel = '#mu pT [GeV]'
      elif matching == 'pion': ylabel = '#pi pT [GeV]'
      hist.GetYaxis().SetTitle(ylabel)
      hist.GetYaxis().SetLabelSize(0.037)
      hist.GetYaxis().SetTitleSize(0.042)
      hist.GetYaxis().SetTitleOffset(1.8)
      hist.GetZaxis().SetTitle('Efficiency')
      hist.GetZaxis().SetLabelSize(0.037)
      hist.GetZaxis().SetTitleSize(0.042)
      hist.GetZaxis().SetTitleOffset(1.2)
      hist.GetZaxis().SetRangeUser(-1e-9, 1)

      ROOT.gPad.Modified()
      ROOT.gPad.Update()

      canv.cd()
      canv.SaveAs('{}/2d_matching_{}.png'.format(self.outputdir, matching))
      canv.SaveAs('{}/2d_matching_{}.pdf'.format(self.outputdir, matching))

    
  def plot1DEfficiency(self, binning='displacement'):
    if binning not in ['displacement', 'pt']:
      raise RuntimeError("Please choose in which quantity to bin ('displacement', 'bin')")

    if binning == 'displacement':
      bins = self.displacement_bins
      displacement_bins = bins
      pt_bins = [(0, 1000)]
    elif binning == 'pt':
      bins = self.pt_bins
      displacement_bins = [(0, 1000)]
      pt_bins = bins

    #efficiency, error = self.getEfficiency(displacement_bins=displacement_bins, pt_bins=pt_bins)
    if self.sample_type == 'flat':
      efficiency, error = self.getEfficiency(displacement_bins=displacement_bins, pt_bins=pt_bins)
    else:
      efficiency, error = self.getEfficiencyFromNano(displacement_bins=displacement_bins, pt_bins=pt_bins)

    for imatch, matching in enumerate(self.matchings):
      canv_name = 'canv_{}_{}_{}'.format(matching, binning, self.outdirlabel)
      canv = self.tools.createTCanvas(canv_name, 900, 800)

      pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
      pad.Draw()
      pad.SetGrid()
      pad.cd()

      graph = ROOT.TGraphAsymmErrors()
      for ibin, bin_ in enumerate(bins):
        bin_min, bin_max = bin_
        x = (float(bin_min) + float(bin_max))/2.
        y = efficiency[ibin][0][imatch] if binning == 'displacement' else efficiency[0][ibin][imatch]
        err = error[ibin][0][imatch] if binning == 'displacement' else error[0][ibin][imatch]
        point = graph.GetN()
        graph.SetPoint(point, x, y)
        graph.SetPointError(point, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err, err)
        
      graph.Draw('AP')  

      graph.SetTitle(self.title)
      graph.SetLineColor(ROOT.kBlue+2 if self.add_dsa == 'no' else ROOT.kRed)
      graph.SetMarkerStyle(20)
      graph.SetMarkerColor(ROOT.kBlue+2 if self.add_dsa == 'no' else ROOT.kRed)
      if binning == 'displacement': xlabel = 'l_{xy} [cm]'
      elif binning == 'pt' and matching == 'candidate': xlabel = '#mu#pi pT [GeV]' if self.process_type == 'mupi' else '#mu#mu#pi pT [GeV]'
      elif binning == 'pt' and matching == 'trigger_muon': xlabel = 'trigger #mu pT [GeV]'
      elif binning == 'pt' and matching == 'muon': xlabel = '#mu pT [GeV]'
      elif binning == 'pt' and matching == 'pion': xlabel = '#pi pT [GeV]'
      graph.GetXaxis().SetTitle(xlabel)
      graph.GetXaxis().SetLabelSize(0.037)
      graph.GetXaxis().SetTitleSize(0.042)
      graph.GetXaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetTitle('Efficiency')
      graph.GetYaxis().SetLabelSize(0.037)
      graph.GetYaxis().SetTitleSize(0.042)
      graph.GetYaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetRangeUser(0, 1.)

      #leg = self.tools.getRootTLegend(xmin=0.5, ymin=0.3, xmax=0.7, ymax=0.5, size=0.035)
      #leg.AddEntry(graph, 'slimmed muons')
      #leg.Draw()

      label = ROOT.TPaveText(0.6,0.7,0.85,0.8,"brNDC")
      text_process = 'HNL#rightarrow #mu#pi' if self.process_type == 'mupi' else 'B#rightarrow#mu#mu#pi'
      text = 'inclusive pT' if binning == 'displacement' else 'inclusive displacement'
      label.AddText(text_process)
      label.AddText(text)
      label.SetBorderSize(0)
      label.SetFillColorAlpha(0, 0)
      label.SetTextSize(0.035)
      label.SetTextFont(42)
      label.SetTextAlign(11)
      label.Draw()

      canv.cd()
      canv.SaveAs('{}/1d_{}_matching_{}.png'.format(self.outputdir, binning, matching))
      canv.SaveAs('{}/1d_{}_matching_{}.pdf'.format(self.outputdir, binning, matching))


  def plotAcceptanceScan(self, quantity='trigger_muon'):
    efficiency, error = self.getEfficiency(displacement_bins=None, pt_bins=None)

    for imatch, matching in enumerate(self.matchings):
      canv_name = 'canv_{}_{}'.format(matching, self.outdirlabel)
      canv = self.tools.createTCanvas(canv_name, 900, 800)

      pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
      pad.Draw()
      pad.cd()

      graph = ROOT.TGraphAsymmErrors()
      for icut, cut_ in enumerate(self.cuts):
        x = cut_ 
        y = efficiency[icut]
        err = error[icut]
        point = graph.GetN()
        graph.SetPoint(point, x, y)
        #graph.SetPointError(point, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err, err)
        graph.SetPointError(point, 0., 0., 0., 0.)
        
      graph.Draw('AP')  

      graph.SetTitle(self.title)
      graph.SetLineColor(ROOT.kBlue+2)
      graph.SetMarkerStyle(20)
      graph.SetMarkerColor(ROOT.kBlue+2)
      graph.GetXaxis().SetTitle('Cut on {} pT [GeV]'.format(quantity))
      graph.GetXaxis().SetLabelSize(0.037)
      graph.GetXaxis().SetTitleSize(0.042)
      graph.GetXaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetTitle('Efficiency')
      graph.GetYaxis().SetLabelSize(0.037)
      graph.GetYaxis().SetTitleSize(0.042)
      graph.GetYaxis().SetTitleOffset(1.1)

      label = ROOT.TPaveText(0.15,0.1,0.45,0.2,"brNDC")
      text = 'inclusive pT and displacement'
      label.AddText(text)
      label.SetBorderSize(0)
      label.SetFillColor(ROOT.kWhite)
      label.SetTextSize(0.035)
      label.SetTextFont(42)
      label.SetTextAlign(11)
      label.Draw()

      canv.cd()
      canv.SaveAs('{}/acceptance_cut_{}_pt_matching_{}.png'.format(self.outputdir, quantity, matching))
      canv.SaveAs('{}/acceptance_cut_{}_pt_matching_{}.pdf'.format(self.outputdir, quantity, matching))


  def plotMuonIDEfficiency(self, muon_type, muon_id, binning):
    if muon_type not in ['primary_muon', 'displaced_muon']:
      raise RuntimeError('Muon type "{}" not valid.'.format(muon_type))

    if muon_id not in ['soft', 'loose']:
      raise RuntimeError('Muon ID "{}" not valid.'.format(muon_id))

    if binning not in ['pt', 'displacement']:
      raise RuntimeError('Binning "{}" not valid'.format(binning))

    filename_m1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root' 
    filename_m3 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root'
    filename_m4p5 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root'

    f_m1 = self.tools.getRootFile(filename_m1, with_ext=False) 
    f_m3 = self.tools.getRootFile(filename_m3, with_ext=False) 
    f_m4p5 = self.tools.getRootFile(filename_m4p5, with_ext=False) 

    tree_m1 = self.tools.getTree(f_m1, 'signal_tree')
    tree_m3 = self.tools.getTree(f_m3, 'signal_tree')
    tree_m4p5 = self.tools.getTree(f_m4p5, 'signal_tree')

    # for fake rate
    qcd_files = qcd_samples['V11_24Apr22_sources']
    white_list = wt['20to300']

    if binning == 'pt':
      #bins = [(1, 7), (7, 8), (8, 9), (9, 10), (10, 12), (12, 20), (20, 30), (30, 100)] 
      #bins = [(1, 7), (7, 10), (10, 20), (20, 30)] 
      #bins = [(1, 7), (7, 10), (10, 15), (15, 20), (20, 30)] 
      #bins = [(1, 7), (7, 10), (10, 15), (15, 20), (20, 25), (25, 30)] 
      #bins = [(1, 7), (7, 10), (10, 15), (15, 30), (30, 100)] 
      #bins = [(1, 3), (3, 6), (6, 10), (10, 15), (15, 20)] 
      bins = [(1, 4), (4, 7), (7, 10), (10, 15), (15, 20), (20, 30)] 
    elif binning == 'displacement':
      #bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 30), (30, 50), (50, 100)]
      bins = [(0, 3), (3, 5), (5, 10), (10, 20), (20, 100)]
      #bins = [(0, 3), (3, 5)]

    efficiency_m1 = np.zeros(len(bins))
    efficiency_m3 = np.zeros(len(bins))
    efficiency_m4p5 = np.zeros(len(bins))
    fake_rate = np.zeros(len(bins))

    err_m1 = np.zeros(len(bins))
    err_m3 = np.zeros(len(bins))
    err_m4p5 = np.zeros(len(bins))
    err_fake = np.zeros(len(bins))
      
    hnl_mass = Quantity(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=0, bin_max=1000)

    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_

      if binning == 'pt': qte = '{part}_pt'
      elif binning == 'displacement': qte = 'sv_lxy'

      if muon_type == 'primary_muon': part = 'mu0'
      elif muon_type == 'displaced_muon': part = 'mu'

      if muon_id == 'soft': muid = 'softid'
      elif muon_id == 'loose': muid = 'looseid'

      if ibin != len(bins)-1:
        #selection_num = '{part}_{qte}>= {mn} && {part}_{qte}<{mx} && {part}_{muid}==1'.format(qte=qte, part=part, muid=muid, mn=bin_min, mx=bin_max)
        #selection_deno = '{part}_{qte}>= {mn} && {part}_{qte}<{mx}'.format(qte=qte, part=part, mn=bin_min, mx=bin_max)
        selection_num = '{qte}>= {mn} && {qte}<{mx} && {part}_{muid}==1'.format(qte=qte, part=part, muid=muid, mn=bin_min, mx=bin_max)
        selection_deno = '{qte}>= {mn} && {qte}<{mx}'.format(qte=qte, part=part, mn=bin_min, mx=bin_max)
      else: # add overflow
        #selection_num = '{part}_{qte}>= {mn} && {part}_{muid}==1'.format(qte=qte, part=part, muid=muid, mn=bin_min)
        #selection_deno = '{part}_{qte}>= {mn}'.format(qte=qte, part=part, mn=bin_min)
        selection_num = '{qte}>= {mn} && {part}_{muid}==1'.format(qte=qte, part=part, muid=muid, mn=bin_min)
        selection_deno = '{qte}>= {mn}'.format(qte=qte, part=part, mn=bin_min)

      #weight = '(weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2)'
      weight = '(1.)'

      hist_num_m1 = self.tools.createHisto(tree_m1, hnl_mass, hist_name='hist_m1', branchname='flat', selection=selection_num, weight=weight)
      hist_deno_m1 = self.tools.createHisto(tree_m1, hnl_mass, hist_name='hist_m1', branchname='flat', selection=selection_deno, weight=weight)
      hist_num_m3 = self.tools.createHisto(tree_m3, hnl_mass, hist_name='hist_m3', branchname='flat', selection=selection_num, weight=weight)
      hist_deno_m3 = self.tools.createHisto(tree_m3, hnl_mass, hist_name='hist_m3', branchname='flat', selection=selection_deno, weight=weight)
      hist_num_m4p5 = self.tools.createHisto(tree_m4p5, hnl_mass, hist_name='hist_m4p5', branchname='flat', selection=selection_num, weight=weight)
      hist_deno_m4p5 = self.tools.createHisto(tree_m4p5, hnl_mass, hist_name='hist_m4p5', branchname='flat', selection=selection_deno, weight=weight)
        
      err_num_m1 = ROOT.double(0.)
      err_deno_m1 = ROOT.double(0.)
      num_m1 = hist_num_m1.IntegralAndError(0, 1000, err_num_m1)
      deno_m1 = hist_deno_m1.IntegralAndError(0, 1000, err_deno_m1)

      err_num_m3 = ROOT.double(0.)
      err_deno_m3 = ROOT.double(0.)
      num_m3 = hist_num_m3.IntegralAndError(0, 1000, err_num_m3)
      deno_m3 = hist_deno_m3.IntegralAndError(0, 1000, err_deno_m3)

      err_num_m4p5 = ROOT.double(0.)
      err_deno_m4p5 = ROOT.double(0.)
      num_m4p5 = hist_num_m4p5.IntegralAndError(0, 1000, err_num_m4p5)
      deno_m4p5 = hist_deno_m4p5.IntegralAndError(0, 1000, err_deno_m4p5)

      efficiency_m1[ibin] = num_m1 / deno_m1 if deno_m1 != 0 else -1
      efficiency_m3[ibin] = num_m3 / deno_m3 if deno_m3 != 0 else -1
      #print '{} {}'.format(num_m4p5, deno_m4p5)
      efficiency_m4p5[ibin] = num_m4p5 / deno_m4p5 if deno_m4p5 != 0 else -1

      err_m1[ibin] = num_m1 / deno_m1 * math.sqrt(math.pow(err_num_m1 / num_m1, 2) + math.pow(err_deno_m1 / deno_m1, 2)) if deno_m1 != 0 else 0 
      err_m3[ibin] = num_m3 / deno_m3 * math.sqrt(math.pow(err_num_m3 / num_m3, 2) + math.pow(err_deno_m3 / deno_m3, 2)) if deno_m3 != 0 else 0
      err_m4p5[ibin] = num_m4p5 / deno_m4p5 * math.sqrt(math.pow(err_num_m4p5 / num_m4p5, 2) + math.pow(err_deno_m4p5 / deno_m4p5, 2)) if deno_m4p5 != 0 else 0

      # fake rate
      if muon_type == 'primary_muon': part2 = 'trgmu'
      elif muon_type == 'displaced_muon': part2 = 'mu'
      if ibin != len(bins)-1:
        #selection_fake_num = '{part}_{qte}>= {mn} && {part}_{qte}<{mx} && {part2}_{muid}==1 && {part}_isfake==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min, mx=bin_max)
        #selection_fake_deno = '{part}_{qte}>= {mn} && {part}_{qte}<{mx} && {part2}_{muid}==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min, mx=bin_max)
        selection_fake_num = '{qte}>= {mn} && {qte}<{mx} && {part2}_{muid}==1 && {part}_isfake==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min, mx=bin_max)
        selection_fake_deno = '{qte}>= {mn} && {qte}<{mx} && {part2}_{muid}==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min, mx=bin_max)
      else: # add overflow
        #selection_fake_num = '{part}_{qte}>= {mn} && {part2}_{muid}==1 && {part}_isfake==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min)
        #selection_fake_deno = '{part}_{qte}>= {mn} && {part2}_{muid}==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min)
        selection_fake_num = '{qte}>= {mn} && {part2}_{muid}==1 && {part}_isfake==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min)
        selection_fake_deno = '{qte}>= {mn} && {part2}_{muid}==1'.format(qte=qte, part=part, part2=part2, muid=muid, mn=bin_min)

      hist_fake_num = self.tools.createWeightedHistoQCDMC(qcd_files, white_list=white_list, quantity=hnl_mass, selection=selection_fake_num, treename='sources')
      hist_fake_deno = self.tools.createWeightedHistoQCDMC(qcd_files, white_list=white_list, quantity=hnl_mass, selection=selection_fake_deno, treename='sources')

      err_fake_num = ROOT.double(0.)
      err_fake_deno = ROOT.double(0.)
      fake_num = hist_fake_num.IntegralAndError(0, 1000, err_fake_num)
      fake_deno = hist_fake_deno.IntegralAndError(0, 1000, err_fake_deno)

      fake_rate[ibin] = fake_num / fake_deno
      err_fake[ibin] = fake_num / fake_deno * math.sqrt(math.pow(err_fake_num / fake_num, 2) + math.pow(err_fake_deno / fake_deno, 2)) 

    print efficiency_m1
    print efficiency_m3
    print efficiency_m4p5

    canv = self.tools.createTCanvas('canv', 900, 800)

    pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
    pad.Draw()
    pad.SetGrid()
    pad.cd()

    leg = self.tools.getRootTLegend(xmin=0.15, ymin=0.4, xmax=0.4, ymax=0.6, size=0.035)

    graph_m1 = ROOT.TGraphAsymmErrors()
    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_
      x_m1 = (float(bin_min) + float(bin_max))/2.
      y_m1 = efficiency_m1[ibin]
      #err = 0.

      point_m1 = graph_m1.GetN()
      graph_m1.SetPoint(point_m1, x_m1, y_m1)
      graph_m1.SetPointError(point_m1, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err_m1[ibin], err_m1[ibin])

    graph_m1.SetLineColor(ROOT.kOrange+0)
    graph_m1.SetMarkerStyle(20)
    graph_m1.SetMarkerColor(ROOT.kOrange+0)
    leg.AddEntry(graph_m1, 'mass 1 GeV, ctau 1000 mm')

    graph_m3 = ROOT.TGraphAsymmErrors()
    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_
      x_m3 = (float(bin_min) + float(bin_max))/2.
      y_m3 = efficiency_m3[ibin]
      #err = 0.

      point_m3 = graph_m3.GetN()
      graph_m3.SetPoint(point_m3, x_m3, y_m3)
      graph_m3.SetPointError(point_m3, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err_m3[ibin], err_m3[ibin])

    graph_m3.SetLineColor(ROOT.kRed+1)
    graph_m3.SetMarkerStyle(20)
    graph_m3.SetMarkerColor(ROOT.kRed+1)
    leg.AddEntry(graph_m3, 'mass 3 GeV, ctau 100 mm')

    graph_m4p5 = ROOT.TGraphAsymmErrors()
    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_
      x_m4p5 = (float(bin_min) + float(bin_max))/2.
      y_m4p5 = efficiency_m4p5[ibin]
      #err = 0.

      point_m4p5 = graph_m4p5.GetN()
      graph_m4p5.SetPoint(point_m4p5, x_m4p5, y_m4p5)
      graph_m4p5.SetPointError(point_m4p5, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err_m4p5[ibin], err_m4p5[ibin])

    graph_m4p5.SetLineColor(ROOT.kRed+4)
    graph_m4p5.SetMarkerStyle(20)
    graph_m4p5.SetMarkerColor(ROOT.kRed+4)
    leg.AddEntry(graph_m4p5, 'mass 4.5 GeV, ctau 1 mm')

    graph_fake = ROOT.TGraphAsymmErrors()
    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_
      x_fake = (float(bin_min) + float(bin_max))/2.
      y_fake = fake_rate[ibin]
      #err = 0.

      point_fake = graph_fake.GetN()
      graph_fake.SetPoint(point_fake, x_fake, y_fake)
      graph_fake.SetPointError(point_fake, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err_fake[ibin], err_fake[ibin])

    graph_fake.SetLineColor(ROOT.kBlue+2)
    graph_fake.SetMarkerStyle(45)
    graph_fake.SetMarkerSize(2)
    graph_fake.SetMarkerColor(ROOT.kBlue+2)
    leg.AddEntry(graph_fake, 'fake rate')

    graph_m1.Draw('AP')  
    graph_m3.Draw('P same')  
    graph_m4p5.Draw('P same')  
    graph_fake.Draw('P same')
    leg.Draw()

    if muon_type == 'primary_muon': label1 = 'Primary muon'
    elif muon_type == 'displaced_muon': label1 = 'Displaced muon'

    graph_m1.SetTitle('{} ({} ID)'.format(label1, muon_id))
    if binning == 'pt':
      if muon_type == 'primary_muon': xlabel = 'primary muon p_{T} [GeV]'
      elif muon_type == 'displaced_muon': xlabel = 'displaced muon p_{T} [GeV]'
    elif binning == 'displacement':
      xlabel = 'SV l_{xy} [cm]'
    graph_m1.GetXaxis().SetTitle(xlabel)
    graph_m1.GetXaxis().SetLabelSize(0.037)
    graph_m1.GetXaxis().SetTitleSize(0.042)
    graph_m1.GetXaxis().SetTitleOffset(1.1)
    graph_m1.GetYaxis().SetTitle('Efficiency')
    graph_m1.GetYaxis().SetLabelSize(0.037)
    graph_m1.GetYaxis().SetTitleSize(0.042)
    graph_m1.GetYaxis().SetTitleOffset(1.1)
    graph_m1.GetYaxis().SetRangeUser(0., 1.1)

    canv.cd()
    #canv.SaveAs('{}/efficiency_{}_{}id.png'.format(self.outputdir, muon_type, muon_id))
    #canv.SaveAs('{}/efficiency_{}_{}id.pdf'.format(self.outputdir, muon_type, muon_id))
    canv.SaveAs('{}/efficiency_{}_{}id_{}.png'.format(self.outputdir, muon_type, muon_id, binning))
    canv.SaveAs('{}/efficiency_{}_{}id_{}.pdf'.format(self.outputdir, muon_type, muon_id, binning))




if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  do_studyRecoEfficiency = False
  do_studyMuonIDEfficiency = True

  if do_studyRecoEfficiency:
    #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_unresolved_motherpdgid_unfittedmass.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_unresolved_fittedmass_looseselection_originalmatching.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looselection_updatedmatching_fullgen.root' # used for the july study
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looselection_updatedmatching_fullgen.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looseselection_mupi_v1.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looseselection_mupi_minigenmatching_v2.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_efficiency_study.root' # used for study on efficiency gain with dsa (Nov 21)
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V36/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looseselection.root' # for reconstruction efficiency without trigger matching efficiency
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V36/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looseselection_fullgen.root'

    displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 30), (30, 50), (50, 100)]
    #displacement_bins = [(0, 3), (3, 5), (5, 10), (10, 20), (20, 100)]
    pt_bins = [(1, 7), (7, 10), (10, 15), (15, 30), (30, 100)] 
    #pt_bins = [(1, 7), (7, 8), (8, 9), (9, 10), (10, 12), (12, 20), (20, 30), (30, 100)] 

    cuts_mu_pt = [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    cuts_trigger_mu_pt = [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]

    #outdirlabel = 'signalV20_test_unresolved_fittedmass_looseselection_originalmatching_withoutdsa'
    #outdirlabel = 'dsa_study_Nov21_dsa'
    outdirlabel = 'test_V36_fullgen_v1'

    title = ''

    matchings = ['candidate', 'trigger_muon', 'muon', 'pion']
    #matchings = ['candidate', 'trigger_muon']

    sample_type = 'nano'
    process_type = 'mumupi'
    add_dsa = 'no'

    analyser = EfficiencyAnalyser(filename=filename, sample_type=sample_type, process_type=process_type, add_dsa=add_dsa, matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel)

    analyser.plot2DEfficiency()
    analyser.plot1DEfficiency(binning='displacement')
    analyser.plot1DEfficiency(binning='pt')


  if do_studyMuonIDEfficiency:
    outdirlabel = 'muon_ID'
    analyser = EfficiencyAnalyser(
        outdirlabel=outdirlabel,
        )

    muon_type = 'primary_muon'
    muon_id = 'soft'
    binning = 'pt'

    #analyser.plotMuonIDEfficiency(
    #    muon_type = muon_type, 
    #    muon_id = muon_id, 
    #    binning = binning,
    #    )

    muon_type = 'displaced_muon'
    muon_id = 'loose'
    binning = 'displacement' #'pt'

    analyser.plotMuonIDEfficiency(
        muon_type = muon_type, 
        muon_id = muon_id, 
        binning = binning,
        )


