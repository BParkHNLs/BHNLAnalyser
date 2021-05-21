import os
import os.path
from os import path
import ROOT
from ROOT import gROOT, gStyle
import math
import numpy as np

from common import PlottingTools


class EfficiencyAnalyser(PlottingTools):
  def __init__(self, filename, treename='Events', matchings=None, cuts=None, displacement_bins=None, pt_bins=None, title='', outdirlabel='default'):
    self.filename = filename
    self.treename = treename
    self.matchings = matchings
    self.displacement_bins = displacement_bins
    self.pt_bins = pt_bins
    #self.bins = bins # for 1D plot bins can either be in displacement or in pt
    self.cuts = cuts
    self.title = title
    self.outdirlabel = outdirlabel
    self.outputdir = PlottingTools.getOutDir(self, './myPlots/efficiency', self.outdirlabel)


  def getEfficiency(self, displacement_bins=None, pt_bins=None):
    f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+self.filename, 'READ')
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
      if entry.nBToMuMuPi == 0: continue
        
      #print '\n'
      # retrieve candidate matching information
      trgmu_ismatched = 0
      mu_ismatched = 0
      pi_ismatched = 0
      cand_ismatched = 0
      matched_idx = -1
      for icand in range(0, entry.nBToMuMuPi):
        if entry.BToMuMuPi_trg_mu_isMatched[icand] == 1: trgmu_ismatched = 1
        if entry.BToMuMuPi_sel_mu_isMatched[icand] == 1: mu_ismatched = 1
        if entry.BToMuMuPi_pi_isMatched[icand] == 1: pi_ismatched = 1
        if entry.BToMuMuPi_isMatched[icand] == 1: 
          matched_idx = icand
          cand_ismatched = 1

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
        if entry.GenPart_pdgId[igen] == 9900015: 
          hnl_idx = igen 
          mother_idx = entry.GenPart_genPartIdxMother[hnl_idx]

      # search for trigger muon and hnl daughters
      for igen in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[igen]) == 13 and entry.GenPart_genPartIdxMother[igen] == mother_idx: trgmu_idx = igen
        if abs(entry.GenPart_pdgId[igen]) == 13 and entry.GenPart_pdgId[entry.GenPart_genPartIdxMother[igen]] == 9900015: mu_idx = igen
        if abs(entry.GenPart_pdgId[igen]) == 211 and entry.GenPart_pdgId[entry.GenPart_genPartIdxMother[igen]] == 9900015: pi_idx = igen

      # remove e-channel
      if mu_idx == -1: continue

      for imatch, matching in enumerate(self.matchings):
        if displacement_bins != None and pt_bins != None:
          # fetch n_deno and n_num with acceptance cuts
          if entry.GenPart_pt[mu_idx] > 1.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
             and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
             and entry.GenPart_pt[trgmu_idx] > 7.5:
          #if entry.GenPart_pt[mu_idx] > 3.5 and abs(entry.GenPart_eta[mu_idx]) < 2.5 \
          #   and entry.GenPart_pt[pi_idx] > 0.7 and abs(entry.GenPart_eta[pi_idx]) < 2.5 \
          #   and entry.GenPart_pt[trgmu_idx] > 9.5:
          #if entry.GenPart_pt[mu_idx] > 0:
               displacement = math.sqrt(pow(entry.GenPart_vx[mu_idx] - entry.GenPart_vx[trgmu_idx], 2) + pow(entry.GenPart_vy[mu_idx] - entry.GenPart_vy[trgmu_idx], 2))
               if matching == 'candidate': obj_pt = entry.GenPart_pt[hnl_idx]
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


  def plot2DEfficiency(self):
    gStyle.SetPadRightMargin(0.16)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetOptStat(0)
    gStyle.SetPaintTextFormat(".2f")

    efficiency, error = self.getEfficiency(displacement_bins=self.displacement_bins, pt_bins=self.pt_bins)

    for imatch, matching in enumerate(self.matchings):
      if matching not in ['candidate', 'trigger_muon', 'muon', 'pion']:
        raise RuntimeError("Unknown matching strategy. Please choose among ['candidate', 'trigger_muon', 'muon', 'pion']")

      canv = ROOT.TCanvas('2d_canv_{}_{}'.format(matching, self.outdirlabel), '2d_canv_{}'.format(matching, self.outdirlabel), 900, 800)
      #ROOT.SetOwnership(canv, False)

      #pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
      #pad.Draw()
      #pad.cd()

      error_list = []

      bin_min_disp, a = self.displacement_bins[0]
      b ,bin_max_disp = self.displacement_bins[len(self.displacement_bins)-1]
      step_disp = (bin_max_disp-bin_min_disp) / (len(self.displacement_bins)-1)

      errors = []
      hist_name = 'hist_{}_{}'.format(matching, self.outdirlabel)
      hist = ROOT.TH2D(hist_name, hist_name, len(self.displacement_bins), 1, 100, len(self.pt_bins), 1, 100)
      for ibin_disp, displacement_bin in enumerate(self.displacement_bins):
        for ibin_pt, pt_bin in enumerate(self.pt_bins):
          hist.Fill(str(displacement_bin), str(pt_bin), efficiency[ibin_disp][ibin_pt][imatch])

          # print the error
          #bin_min_disp, bin_max_disp = displacement_bin
          bin_min_pt, bin_max_pt = pt_bin
          #x1 = (bin_max_disp + bin_min_disp)/2 - (bin_max_disp - bin_min_disp)*0.25
          #x1 = (bin_max_disp + bin_min_disp)/2.
          x1 = ibin_disp/float((len(self.displacement_bins)-1)) - ibin_disp/float((len(self.displacement_bins)-1))*0.05
          x2 = ibin_disp/float((len(self.displacement_bins)-1)) + ibin_disp/float((len(self.displacement_bins)-1))*0.05
          y1 = ibin_pt/float((len(self.pt_bins)-1)) - ibin_pt/float((len(self.pt_bins)-1))*0.05
          y2 = ibin_pt/float((len(self.pt_bins)-1)) + ibin_pt/float((len(self.pt_bins)-1))*0.05
          #x1 = bin_min_disp + ibin_disp*step_disp
          #print '{} {} {} {}'.format(ibin_disp/float((len(self.displacement_bins)-1)), ibin_disp/float((len(self.displacement_bins)-1))*0.05, x1, x2)
          #x2 = (bin_max_disp + bin_min_disp)/2 + (bin_max_disp - bin_min_disp)*0.25
          #y1 = (bin_max_pt + bin_min_pt)/2 - (bin_max_pt - bin_min_pt)*0.25
          ##print '{} {} {}'.format(bin_min_pt, bin_max_pt, y1)
          #y1 = (bin_max_pt + bin_min_pt)/2.
          #y2 = (bin_max_pt + bin_min_pt)/2 + (bin_max_pt - bin_min_pt)*0.25
          error_print = ROOT.TPaveText(x1, y1, x2, y2, "brNDC")
          ##error_print = ROOT.TLatex(x1, y1, '#pm {}'.format(round(float(error[ibin_disp][ibin_pt][imatch], 2)))) 
          #error_print = ROOT.TLatex(x1, y1, 'X')
          #errors.append(error_print)
          #error_print.AddText('#pm {}'.format(round(error[ibin_disp][ibin_pt][imatch], 2)))
          error_print.AddText('X')
          error_print.SetFillColorAlpha(0, 0)
          error_list.append(error_print)

      #gStyle.SetPadRightMargin(0.16)
      #gStyle.SetPadLeftMargin(0.16)

      #hist.SetBarOffset(0.3)

      hist.Draw('colz')
      hist.SetMarkerSize(2)
      #hist.Draw('text,error' +'same')
      hist.Draw('text' +'same')

      #for error in errors:
      #  error.Draw('same')

      for error in error_list:
        #error.Draw('same')
        error.SetBorderSize(0)
        error.SetTextSize(0.015)
        error.SetTextFont(62)
        error.SetTextAlign(11)

      hist.SetTitle(self.title)
      #hist.SetTitleSize(0.8)
      hist.GetXaxis().SetTitle('l_{xy}(trgmu, mu) [cm]')
      hist.GetXaxis().SetLabelSize(0.037)
      hist.GetXaxis().SetTitleSize(0.042)
      hist.GetXaxis().SetTitleOffset(1.1)
      if matching == 'candidate': ylabel = '#mu#pi pT [GeV]'
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

      #gStyle.SetOptStat(0)
      #ROOT.gStyle.SetPadTopMargin(0.05)
      #ROOT.gStyle.SetPadBottomMargin(0.13)
      #ROOT.gStyle.SetPadLeftMargin(0.05)
      #ROOT.gStyle.SetPadRightMargin(0.70)
      #gStyle.SetPaintTextFormat(".2f")
        
      #ROOT.gStyle = ROOT.gROOT.GetGlobal( "gStyle", 1 )

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
      if matching not in ['candidate', 'trigger_muon', 'muon', 'pion']:
        raise RuntimeError("Unknown matching strategy. Please choose among ['candidate', 'trigger_muon', 'muon', 'pion']")

      canv = ROOT.TCanvas('canv_{}_{}_{}'.format(matching, binning, self.outdirlabel), 'canv_{}_{}_{}'.format(matching, binning, self.outdirlabel), 900, 800)
      #ROOT.SetOwnership(canv, False)
      pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
      pad.Draw()
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
      #graph.SetLineWidth(1.5)
      graph.SetLineColor(ROOT.kBlue+2)
      #graph.SetMarkerSize(2)
      graph.SetMarkerStyle(20)
      graph.SetMarkerColor(ROOT.kBlue+2)
      if binning == 'displacement': xlabel = 'l_{xy}(trgmu, mu) [cm]'
      elif binning == 'pt' and matching == 'candidate': xlabel = '#mu#pi pT [GeV]'
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

      label = ROOT.TPaveText(0.5,0.75,0.75,0.8,"brNDC")
      text = 'inclusive pT' if binning == 'displacement' else 'inclusive displacement'
      label.AddText(text)
      label.SetBorderSize(0)
      label.SetFillColor(ROOT.kWhite)
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
      if matching not in ['candidate', 'trigger_muon', 'muon', 'pion']:
        raise RuntimeError("Unknown matching strategy. Please choose among ['candidate', 'trigger_muon', 'muon', 'pion']")

      canv = ROOT.TCanvas('canv_{}_{}'.format(matching, self.outdirlabel), 'canv_{}'.format(matching, self.outdirlabel), 900, 800)
      #ROOT.SetOwnership(canv, False)
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
      #graph.SetLineWidth(1.5)
      graph.SetLineColor(ROOT.kBlue+2)
      #graph.SetMarkerSize(2)
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

  #plotHistFromTree(do_shape=True)

  #filename1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looseselection_muononly_alltrgmu.root'
  filename1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looseselection_muononly_stdtrgmu.root'
  #filename1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/Chunk0_n500/bparknano_looseselection_muononly_stdtrgmu_nj1.root'
  #filename1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/Chunk0_n500/merged/bparknano_looseselection_muononly_stdtrgmu.root'
  #filename1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selected_muononly_stdtrgmu.root'
  #filename1 = '/scratch/anlyon/bparknano_looseselection_muononly_alltrgmu.root'
  #displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 100)]

  displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 30), (30, 50), (50, 100)]
  #displacement_bins = [(0, 1), (1, 5), (5, 10), (10, 15), (15, 30), (30, 50), (50, 100)]
  #displacement_bins = [(0, 1), (1, 5), (5, 15), (15, 100)]
  pt_bins = [(1, 7), (7, 10), (10, 15), (15, 30), (30, 100)] 
  #pt_bins = [(1, 5), (5, 10), (10, 20), (20, 100)] 

  cuts_mu_pt = [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
  cuts_trigger_mu_pt = [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]
  #plotAcceptanceScan(filename=filename1, treename='Events', cuts=cuts_mu_pt, matching='candidate', quantity='trigger_muon', title='Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel='signalV25_looseselection_alltrgmu')

  #outdirlabel = 'signalV25_looseselection_alltrgmu_tightacceptancecuts_mupt3p5_trgmupt_9p5'
  #outdirlabel = 'signalV25_looseselection_stdtrgmu_styled'
  outdirlabel = 'testing'

  #title = 'Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons'
  title = ''

  matchings = ['candidate', 'trigger_muon', 'muon', 'pion']
  #matchings = ['candidate', 'trigger_muon']
  #matchings = ['trigger_muon']

  EfficiencyAnalyser(filename=filename1, treename='Events', matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel).plot2DEfficiency()
  #EfficiencyAnalyser(filename=filename1, treename='Events', matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel).plot1DEfficiency(binning='displacement')
  #EfficiencyAnalyser(filename=filename1, treename='Events', matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel).plot1DEfficiency(binning='pt')

  #EfficiencyAnalyser(filename=filename1, treename='Events', matchings=matchings, displacement_bins=displacement_bins, pt_bins=pt_bins, title='Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel).plot1DEfficiency(binning='pt')
  #EfficiencyAnalyser(filename=filename1, treename='Events', matchings=matchings, cuts=cuts_trigger_mu_pt, title='Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel).plotAcceptanceScan(quantity='trigger_muon')

  #plot1DEfficiency(filename=filename1, treename='Events', bins=displacement_bins, matching='trigger_muon', binning='displacement', title='Trigger muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename1, treename='Events', bins=displacement_bins, matching='muon', binning='displacement', title='Muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename1, treename='Events', bins=displacement_bins, matching='pion', binning='displacement', title='Pion matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)

  #plot1DEfficiency(filename=filename1, treename='Events', bins=pt_bins, matching='candidate', binning='pt', title='Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename1, treename='Events', bins=pt_bins, matching='trigger_muon', binning='pt', title='Trigger muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename1, treename='Events', bins=pt_bins, matching='muon', binning='pt', title='Muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename1, treename='Events', bins=pt_bins, matching='pion', binning='pt', title='Pion matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)

  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='candidate', title='Candidate matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='trigger_muon', title='Trigger muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='muon', title='Muon matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='pion', title='Pion matching, signal (3GeV,2000mm), loose selection, all trigger muons', outdirlabel=outdirlabel)


  #filename2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_selected_muononly_alltrgmu.root'
  #outdirlabel = 'signalV25_selected_alltrgmu_tightacceptancecuts_mupt3p5_trgmupt_9p5'

  #plot1DEfficiency(filename=filename2, treename='Events', bins=displacement_bins, matching='candidate', binning='displacement', title='Candidate matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=displacement_bins, matching='trigger_muon', binning='displacement', title='Trigger muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=displacement_bins, matching='muon', binning='displacement', title='Muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=displacement_bins, matching='pion', binning='displacement', title='Pion matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)

  #plot1DEfficiency(filename=filename2, treename='Events', bins=pt_bins, matching='candidate', binning='pt', title='Candidate matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=pt_bins, matching='trigger_muon', binning='pt', title='Trigger muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=pt_bins, matching='muon', binning='pt', title='Muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot1DEfficiency(filename=filename2, treename='Events', bins=pt_bins, matching='pion', binning='pt', title='Pion matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)

  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='candidate', title='Candidate matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='trigger_muon', title='Trigger muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='muon', title='Muon matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='pion', title='Pion matching, signal (3GeV,2000mm), selected, all trigger muons', outdirlabel=outdirlabel)




  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='trigger_muon', title='Trigger muon matching, signal (3GeV,184cm), selected, all trigger muons', outdirlabel='signalV20_emu_selected_alltrgmu')
  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='muon', title='Muon matching, signal (3GeV,184cm), selected, all trigger muons', outdirlabel='signalV20_emu_selected_alltrgmu')
  #plot2DEfficiency(filename=filename1, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='pion', title='Pion matching, signal (3GeV,184cm), selected, all trigger muons', outdirlabel='signalV20_emu_selected_alltrgmu')


  #filename2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_muononly_alltrgmu.root'
  #displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 100)]
  #pt_bins = [(1, 5), (5, 10), (10, 15), (15, 20), (20, 100)] #(20, 30), (30, 100)]
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='candidate', title='Candidate matching, signal (3GeV,184cm), looseselection, all trigger muons', outdirlabel='signalV20_emu_looseselection_alltrgmu')
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='trigger_muon', title='Trigger muon matching, signal (3GeV,184cm), looseselection, all trigger muons', outdirlabel='signalV20_emu_looseselection_alltrgmu')
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='muon', title='Muon matching, signal (3GeV,184cm), loose selection, all trigger muons', outdirlabel='signalV20_emu_looseselection_alltrgmu')
  #plot2DEfficiency(filename=filename2, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='pion', title='Pion matching, signal (3GeV,184cm), looseselection, all trigger muons', outdirlabel='signalV20_emu_looseselection_alltrgmu')

  #filename3 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_muononly_stdtrgmu.root'
  #displacement_bins = [(0, 1), (1, 3), (3, 5), (5, 10), (10, 15), (15, 100)]
  #pt_bins = [(1, 5), (5, 10), (10, 15), (15, 20), (20, 100)] #(20, 30), (30, 100)]
  #plot2DEfficiency(filename=filename3, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='candidate', title='Candidate matching, signal (3GeV,184cm), looseselection, standard trigger muons', outdirlabel='signalV20_emu_looseselection_stdtrgmu')
  #plot2DEfficiency(filename=filename3, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='trigger_muon', title='Trigger muon matching, signal (3GeV,184cm), looseselection, standard trigger muons', outdirlabel='signalV20_emu_looseselection_stdtrgmu')
  #plot2DEfficiency(filename=filename3, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='muon', title='Muon matching, signal (3GeV,184cm), loose selection, standard trigger muons', outdirlabel='signalV20_emu_looseselection_stdtrgmu')
  #plot2DEfficiency(filename=filename3, treename='Events', displacement_bins=displacement_bins, pt_bins=pt_bins, matching='pion', title='Pion matching, signal (3GeV,184cm), looseselection, standard trigger muons', outdirlabel='signalV20_emu_looseselection_stdtrgmu')


