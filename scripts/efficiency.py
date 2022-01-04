import os
import os.path
from os import path
import ROOT
from ROOT import gROOT, gStyle
import math
import numpy as np

from tools import Tools


class EfficiencyAnalyser(Tools):
  def __init__(self, filename, sample_type='flat', process_type='mumupi', add_dsa='yes', matchings=None, cuts=None, displacement_bins=None, pt_bins=None, title='', outdirlabel='default'):
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

    for matching in self.matchings:
      if matching not in ['candidate', 'trigger_muon', 'muon', 'pion']:
        raise RuntimeError("Unknown matching strategy. Please choose among ['candidate', 'trigger_muon', 'muon', 'pion']")

    if path.exists('{}/numbers.txt'.format(self.outputdir)):
      os.system('rm {}/numbers.txt'.format(self.outputdir))


  def getEfficiencyFromNano(self, displacement_bins=None, pt_bins=None):
    f = PlottingTools.getRootFile(self, self.filename, with_ext=False) 
    tree = PlottingTools.getTree(self, f, self.treename)

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
        #if entry.BToMuMuPi_trg_mu_isMatched[icand] == 1 and entry.Muon_isDSAMuon[entry.BToMuMuPi_sel_mu_idx[icand]] != 1: trgmu_ismatched = 1
        if entry.BToMuMuPi_trg_mu_isMatched[icand] == 1: trgmu_ismatched = 1
        #if entry.BToMuMuPi_sel_mu_isMatched[icand] == 1 and entry.Muon_isDSAMuon[entry.BToMuMuPi_sel_mu_idx[icand]] != 1: mu_ismatched = 1
        if entry.BToMuMuPi_sel_mu_isMatched[icand] == 1: mu_ismatched = 1
        #if entry.BToMuMuPi_pi_isMatched[icand] == 1 and entry.Muon_isDSAMuon[entry.BToMuMuPi_sel_mu_idx[icand]] != 1: pi_ismatched = 1
        if entry.BToMuMuPi_pi_isMatched[icand] == 1: pi_ismatched = 1
        #if entry.BToMuMuPi_isMatched[icand] == 1 and entry.Muon_isDSAMuon[entry.BToMuMuPi_sel_mu_idx[icand]] != 1: 
        if entry.BToMuMuPi_isMatched[icand] == 1: 
          cand_ismatched = 1
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
          if entry.GenPart_pt[mu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
             and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
             and entry.GenPart_pt[trgmu_idx] > 7.:
          #if entry.GenPart_pt[mu_idx] > 3.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
          #   and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
          #   and entry.GenPart_pt[trgmu_idx] > 9.5:
          #if entry.GenPart_pt[mu_idx] > 0:
               displacement = math.sqrt(pow(entry.GenPart_vx[mu_idx] - entry.GenPart_vx[trgmu_idx], 2) + pow(entry.GenPart_vy[mu_idx] - entry.GenPart_vy[trgmu_idx], 2))
               if matching == 'candidate': obj_pt = entry.GenPart_pt[mother_idx]
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
                 n_deno[icut] = n_deno[icut] + 1
                 if matching == 'candidate' and cand_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'trigger_muon' and trgmu_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'muon' and mu_ismatched == 1: n_num[icut] = n_num[icut] + 1
                 elif matching == 'pion' and pi_ismatched == 1: n_num[icut] = n_num[icut] + 1

    # compute efficiency
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
        print '\n'
        for ibin_disp, displacement_bin in enumerate(displacement_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            print '{} {} {} {} {} {}+-{}'.format(matching, displacement_bin, pt_bin, n_deno[ibin_disp][ibin_pt][imatch], n_num[ibin_disp][ibin_pt][imatch], efficiency[ibin_disp][ibin_pt][imatch], error[ibin_disp][ibin_pt][imatch])
    else:
      for icut, cut in enumerate(self.cuts):
        efficiency[icut] = float(n_num[icut]) / float(n_deno[icut]) if float(n_deno[icut]) != 0. else 0. 
        error[icut] = efficiency[icut] * ( math.sqrt(n_num[icut])/n_num[icut]  + math.sqrt(n_deno[icut])/n_deno[icut] )     

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
      
      hnl_idx = -1
      mother_idx = -1
      trgmu_idx = -1
      mu_idx = -1
      pi_idx = -1
      
      for imatch, matching in enumerate(self.matchings):
        # fetch n_deno and n_num with acceptance cuts
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
      if matching == 'candidate': ylabel = '#mu#pi pT [GeV]' if self.process_type == 'mupi' else '#mu#mu#pi pT [GeV]'
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

    efficiency, error = self.getEfficiency(displacement_bins=displacement_bins, pt_bins=pt_bins)

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
      graph.GetYaxis().SetRangeUser(0, 0.5)

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



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano.root'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_unresolved_motherpdgid_unfittedmass.root'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_test/mass3.0_ctau184.256851021/nanoFiles/merged/flat_bparknano_unresolved_fittedmass_looseselection_originalmatching.root'
  filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looselection_updatedmatching_fullgen.root'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looseselection_mupi_v1.root'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looseselection_mupi_minigenmatching_v2.root'
  #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_efficiency_study.root' # used for study on efficiency gain with dsa (Nov 21)

  displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 30), (30, 50), (50, 100)]
  #displacement_bins = [(0, 3), (3, 5), (5, 10), (10, 20), (20, 100)]
  pt_bins = [(1, 7), (7, 10), (10, 15), (15, 30), (30, 100)] 

  cuts_mu_pt = [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
  cuts_trigger_mu_pt = [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]

  #outdirlabel = 'signalV20_test_unresolved_fittedmass_looseselection_originalmatching_withoutdsa'
  #outdirlabel = 'dsa_study_Nov21_dsa'
  outdirlabel = 'test'

  title = ''

  #matchings = ['candidate', 'trigger_muon', 'muon', 'pion']
  matchings = ['candidate', 'trigger_muon']

  sample_type = 'flat'
  process_type = 'mumupi'
  add_dsa = 'no'

  analyser = EfficiencyAnalyser(filename=filename, sample_type=sample_type, process_type=process_type, add_dsa=add_dsa, matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel)

  analyser.plot2DEfficiency()
  analyser.plot1DEfficiency(binning='displacement')
  analyser.plot1DEfficiency(binning='pt')

