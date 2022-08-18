import os
from os import path
import ROOT
from tools import Tools

import sys
sys.path.append('../objects')
from quantity import Quantity, quantities
from samples import qcd_samples
from qcd_white_list import white_list
from baseline_selection import selection
from categories import categories


class BackgroundSource(object):
  def __init__(self, name, label, colour):
    self.name = name
    self.label = label
    self.colour = colour


class BackgroundPlotter(object):
  def __init__(self, quantity, qcd_files, white_list='', categories='', baseline_selection='', add_weight_hlt='', add_weight_muid='',  add_weight_pu='', weight_hlt='', weight_puqcd='', weight_mu0id='', weight_muid='', CMS_tag=''):
    self.tools = Tools()
    self.quantity = quantity
    self.qcd_files = qcd_files
    self.white_list = white_list
    self.categories = categories
    self.baseline_selection = baseline_selection
    self.add_weight_hlt = add_weight_hlt
    self.add_weight_pu = add_weight_pu
    self.add_weight_muid = add_weight_muid
    self.weight_hlt = weight_hlt
    self.weight_puqcd = weight_puqcd
    self.weight_mu0id = weight_mu0id
    self.weight_muid = weight_muid
    self.CMS_tag = CMS_tag


  def plotStackHistogram(self, sources='', title='', plot_label=''):

    for category in self.categories:

      # will contain the stack of the different sources
      hist_stack = ROOT.THStack('hist_stack_{}'.format(category.label), '')

      # is used to check on the normalisation
      hist_tot_control_name = 'hist_tot_control_{}'.format(category.label)
      hist_tot_control = ROOT.TH1D(hist_tot_control_name, hist_tot_control_name, self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_tot_control.SetDirectory(0)

      # get the legend
      legend = self.tools.getRootTLegend(xmin=0.53, ymin=0.6, xmax=0.87, ymax=0.87, size=0.027)

      # for each source, collect all the pthat ranges
      hist_tot_all = []
      for isource, source in enumerate(sources):
        hist_tot_name = 'hist_tot_{}'.format(category.label)
        hist_tot = ROOT.TH1D(hist_tot_name, hist_tot_name, self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
        hist_tot.Sumw2()
        hist_tot.SetDirectory(0)

        hist_qcd_tot = []

        for qcd_file in self.qcd_files:
          qcd_file_pthatrange = self.tools.getPthatRange(qcd_file.label)
          if qcd_file_pthatrange not in self.white_list: continue

          f_qcd = self.tools.getRootFile(qcd_file.filename)
          tree = self.tools.getTree(f_qcd, 'sources')
          tree_run = self.tools.getTree(f_qcd, 'run_tree')

          weight_qcd = self.tools.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
          weight_qcd = '({})'.format(weight_qcd)
          if self.add_weight_hlt : weight_qcd += ' * ({})'.format(self.weight_hlt)
          if self.add_weight_pu : weight_qcd += ' * ({}) '.format(self.weight_puqcd)
          if self.add_weight_muid : weight_qcd += ' * ({}) *({})'.format(self.weight_mu0id, self.weight_mu0id)

          hist_name = 'hist_{}_{}'.format(qcd_file_pthatrange, category.label)
          selection_qcd = self.baseline_selection + ' && ' + category.definition_flat + ' && ' + category.cutbased_selection + ' && ' + source.name
          hist = self.tools.createHisto(tree, self.quantity, hist_name=hist_name, branchname='flat', selection=selection_qcd, weight=weight_qcd)

          hist_qcd_tot.append(hist)

          hist_tot.Add(hist)
          hist_tot.SetFillColor(source.colour)
          hist_tot.SetLineWidth(0)

        hist_tot_all.append(hist_tot)
          
        legend.AddEntry(hist_tot, source.label)

      int_tot = 0
      for hist_tot in hist_tot_all:
        int_tot += hist_tot.Integral()

      for hist_tot in hist_tot_all:
        if int_tot != 0.: hist_tot.Scale(1./int_tot)

        # add source to stack histogram
        hist_stack.Add(hist_tot)
        hist_tot_control.Add(hist_tot)
  
      print hist_tot_control.Integral()

      frame = hist_tot.Clone('frame')
      frame.SetTitle('')
      frame.GetXaxis().SetTitle(quantity.title)
      frame.GetYaxis().SetLabelSize(0.033)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetTitle('Normalised to unity')
      frame.GetYaxis().SetLabelSize(0.033)
      frame.GetYaxis().SetTitleSize(0.042)
      frame.GetYaxis().SetTitleOffset(1.1)
      frame.GetYaxis().SetRangeUser(0., hist_stack.GetMaximum() + 0.15*hist_stack.GetMaximum())

      canv = self.tools.createTCanvas('canv_{}_{}'.format(category.label, plot_label), 800, 700)
      canv.cd()

      pad = ROOT.TPad("pad","pad",0,0,1,1)
      pad.Draw()
      pad.cd()

      ROOT.gStyle.SetOptStat(0)

      frame.Draw()
      hist_stack.Draw('histo same')
      hist_tot_control.Draw('histo same')

      legend.Draw()

      self.tools.printCMSTag(pad, self.CMS_tag, size=0.43)
      self.tools.printLatexBox(0.55, 0.56, title, size=0.033, pos='left')
      self.tools.printLatexBox(0.55, 0.5, category.title, size=0.033, pos='left')

      if not path.exists('./myPlots/background_sources'):
        os.system('mkdir -p ./myPlots/background_sources')

      canv.cd()
      canv.SaveAs('myPlots/background_sources/{}_{}.png'.format(plot_label, category.label))
      canv.SaveAs('myPlots/background_sources/{}_{}.pdf'.format(plot_label, category.label))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  qcd_files = qcd_samples['V11_24Apr22_sources']
  white_list = white_list['20to300']
  baseline_selection = selection['baseline_sources'].flat
  categories = categories['inclusive']

  CMS_tag = 'Preliminary'
  
  add_weight_hlt = False
  add_weight_pu = False
  add_weight_muid = False
  weight_hlt = ''
  weight_puqcd = ''
  weight_mu0id = ''
  weight_muid = ''

  mu0_nonmatched = BackgroundSource('mu0_ismatched!=1', 'non-genmatched', ROOT.kGray)
  mu0_matched = BackgroundSource('mu0_ismatched==1 && mu0_isfake==0', 'genmatched', ROOT.kMagenta-10)
  mu0_fake = BackgroundSource('mu0_ismatched==1 && mu0_isfake==1', 'genmatched + fake', ROOT.kAzure-9)
  mu0_directBdecay = BackgroundSource('mu0_directBdecay==1 && mu0_ismatched==1', 'direct B decay', ROOT.kRed-9)
  mu0_directBdecay_fullmatched = BackgroundSource('mu0_directBdecay==1 && cand_isfullmatched==1', 'direct B decay', ROOT.kRed-9)
  mu0_BtoDdecay = BackgroundSource('mu0_BtoDdecay==1 && mu0_ismatched==1', 'B->D decay', ROOT.kRed-4)
  mu0_BtoDdecay_fullmatched = BackgroundSource('mu0_BtoDdecay==1 && cand_isfullmatched==1', 'B->D decay', ROOT.kRed-4)
  mu0_BtoDtoKdecay = BackgroundSource('mu0_BtoDtoKdecay==1 && mu0_ismatched==1', 'B->D->K decay', ROOT.kViolet+1)
  mu0_Btotaudecay = BackgroundSource('mu0_Btotaudecay==1 && mu0_ismatched==1', 'B->tau decay', ROOT.kViolet-9)
  mu0_BtoJPsidecay = BackgroundSource('mu0_BtoJPsidecay==1 && mu0_ismatched==1', 'B->JPsi decay', ROOT.kPink+1)
  mu0_BtoDexciteddecay = BackgroundSource('mu0_BtoDexciteddecay==1 && mu0_ismatched==1', 'excited D decay chain', ROOT.kPink-9)
  mu0_other = BackgroundSource('mu0_other==1 && mu0_ismatched==1', 'other sources from B', ROOT.kRed+2)
  mu0_nonBdecay = BackgroundSource('mu0_nonBdecay==1 && mu0_ismatched==1', 'non B decay', ROOT.kRed+4)
  mu0_nonBdecay_fullmatched = BackgroundSource('mu0_nonBdecay==1 && cand_isfullmatched==1', 'non B decay', ROOT.kRed+4)
  mu0_cascade = BackgroundSource('mu0_directBdecay==0 && mu0_BtoDdecay==0 && mu0_nonBdecay==0 && mu_ismatched==1', 'B cascade decay', ROOT.kRed+2)
  mu0_cascade_fullmatched = BackgroundSource('mu0_directBdecay==0 && mu0_BtoDdecay==0 && mu0_nonBdecay==0 && cand_isfullmatched==1', 'B cascade decay', ROOT.kRed+2)

  mu_nonmatched = BackgroundSource('mu_ismatched!=1', 'non-genmatched', ROOT.kGray)
  mu_matched = BackgroundSource('mu_ismatched==1 && mu_isfake==0', 'genmatched', ROOT.kMagenta-10)
  mu_fake = BackgroundSource('mu_ismatched==1 && mu_isfake==1', 'genmatched + fake', ROOT.kAzure-9)
  mu_directBdecay = BackgroundSource('mu_directBdecay==1 && mu_ismatched==1', 'direct B decay', ROOT.kRed-9)
  mu_directBdecay_fullmatched = BackgroundSource('mu_directBdecay==1 && cand_isfullmatched==1', 'direct B decay', ROOT.kRed-9)
  mu_BtoDdecay = BackgroundSource('mu_BtoDdecay==1 && mu_ismatched==1', 'B->D decay', ROOT.kRed-4)
  mu_BtoDdecay_fullmatched = BackgroundSource('mu_BtoDdecay==1 && cand_isfullmatched==1', 'B->D decay', ROOT.kRed-4)
  mu_BtoDtoKdecay = BackgroundSource('mu_BtoDtoKdecay==1 && mu_ismatched==1', 'B->D->K decay', ROOT.kViolet+1)
  mu_Btotaudecay = BackgroundSource('mu_Btotaudecay==1 && mu_ismatched==1', 'B->tau decay', ROOT.kViolet-9)
  mu_BtoJPsidecay = BackgroundSource('mu_BtoJPsidecay==1 && mu_ismatched==1', 'B->JPsi decay', ROOT.kPink+1)
  mu_BtoDexciteddecay = BackgroundSource('mu_BtoDexciteddecay==1 && mu_ismatched==1', 'excited D decay chain', ROOT.kPink-9)
  mu_other = BackgroundSource('mu_other==1 && mu_ismatched==1', 'other sources from B', ROOT.kRed+2)
  mu_nonBdecay = BackgroundSource('mu_nonBdecay==1 && mu_ismatched==1', 'non B decay', ROOT.kRed+4)
  mu_nonBdecay_fullmatched = BackgroundSource('mu_nonBdecay==1 && cand_isfullmatched==1', 'non B decay', ROOT.kRed+4)
  mu_cascade = BackgroundSource('mu_directBdecay==0 && mu_BtoDdecay==0 && mu_nonBdecay==0 && mu_ismatched==1', 'B cascade decay', ROOT.kRed+2)
  mu_cascade_fullmatched = BackgroundSource('mu_directBdecay==0 && mu_BtoDdecay==0 && mu_nonBdecay==0 && cand_isfullmatched==1', 'B cascade decay', ROOT.kRed+2)

  pi_nonmatched = BackgroundSource('pi_ismatched!=1', 'non-genmatched', ROOT.kGray)
  pi_matched = BackgroundSource('pi_ismatched==1 && pi_isfake==0', 'genmatched', ROOT.kMagenta-10)
  pi_fake = BackgroundSource('pi_ismatched==1 && pi_isfake==1', 'genmatched + fake', ROOT.kAzure-9)
  pi_directBdecay = BackgroundSource('pi_directBdecay==1 && pi_ismatched==1', 'direct B decay', ROOT.kRed-9)
  pi_directBdecay_fullmatched = BackgroundSource('pi_directBdecay==1 && cand_isfullmatched==1', 'direct B decay', ROOT.kRed-9)
  pi_BtoDdecay = BackgroundSource('pi_BtoDdecay==1 && pi_ismatched==1', 'B->D decay', ROOT.kRed-4)
  pi_BtoDdecay_fullmatched = BackgroundSource('pi_BtoDdecay==1 && cand_isfullmatched==1', 'B->D decay', ROOT.kRed-4)
  pi_BtoDtoKdecay = BackgroundSource('pi_BtoDtoKdecay==1 && pi_ismatched==1', 'B->D->K decay', ROOT.kViolet+1)
  pi_Btotaudecay = BackgroundSource('pi_Btotaudecay==1 && pi_ismatched==1', 'B->tau decay', ROOT.kViolet-9)
  pi_BtoJPsidecay = BackgroundSource('pi_BtoJPsidecay==1 && pi_ismatched==1', 'B->JPsi decay', ROOT.kPink+1)
  pi_BtoDexciteddecay = BackgroundSource('pi_BtoDexciteddecay==1 && pi_ismatched==1', 'excited D decay chain', ROOT.kPink-9)
  pi_other = BackgroundSource('pi_other==1 && pi_ismatched==1', 'other sources from B', ROOT.kRed+2)
  pi_nonBdecay = BackgroundSource('pi_nonBdecay==1 && pi_ismatched==1', 'non B decay', ROOT.kRed+4)
  pi_nonBdecay_fullmatched = BackgroundSource('pi_nonBdecay==1 && cand_isfullmatched==1', 'non B decay', ROOT.kRed+4)
  pi_cascade = BackgroundSource('pi_directBdecay==0 && pi_BtoDdecay==0 && pi_nonBdecay==0 && pi_ismatched==1', 'B cascade decay', ROOT.kRed+2)
  pi_cascade_fullmatched = BackgroundSource('pi_directBdecay==0 && pi_BtoDdecay==0 && pi_nonBdecay==0 && cand_isfullmatched==1', 'B cascade decay', ROOT.kRed+2)

  cand_nonmatched = BackgroundSource('cand_ismatched!=1', 'non-genmatched', ROOT.kGray)
  cand_mu0mu_samemother = BackgroundSource('cand_isfullmatched==1 && cand_mu0mu_samemother==1 && cand_mumupi_samemother==0', 'genmatched + mu0mu same mother', ROOT.kRed-10)
  cand_mu0pi_samemother = BackgroundSource('cand_isfullmatched==1 && cand_mu0pi_samemother==1 && cand_mumupi_samemother==0', 'genmatched + mu0pi same mother', ROOT.kRed-9)
  cand_mupi_samemother = BackgroundSource('cand_isfullmatched==1 && cand_mupi_samemother==1 && cand_mumupi_samemother==0', 'genmatched + mupi same mother', ROOT.kRed+0)
  #cand_mumupi_samemother = BackgroundSource('cand_ismatched && cand_mumupi_samemother', 'genmatched + mumupi same mother', ROOT.kMagenta-10)
  cand_differentmothers = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==0 && cand_mupi_samemother==0 && cand_mu0pi_samemother==0 && cand_mu0mu_samemother==0', 'genmatched + different mothers', ROOT.kRed+4)

  cand_partial_mu0mu_samemother = BackgroundSource('cand_ispartialmatched==1 && cand_mu0mu_samemother==1 && cand_mumupi_samemother==0', 'partial genmatched + mu0mu same mother', ROOT.kGreen-10)
  cand_partial_mu0pi_samemother = BackgroundSource('cand_ispartialmatched==1 && cand_mu0pi_samemother==1 && cand_mumupi_samemother==0', 'partial genmatched + mu0pi same mother', ROOT.kGreen-6)
  cand_partial_mupi_samemother = BackgroundSource('cand_ispartialmatched==1 && cand_mupi_samemother==1 && cand_mumupi_samemother==0', 'partial genmatched + mupi same mother', ROOT.kSpring-6)
  #cand_mumupi_samemother = BackgroundSource('cand_ismatched && cand_mumupi_samemother', 'genmatched + mumupi same mother', ROOT.kMagenta-10)
  cand_partial_differentmothers = BackgroundSource('cand_ispartialmatched==1 && cand_mumupi_samemother==0 && cand_mupi_samemother==0 && cand_mu0pi_samemother==0 && cand_mu0mu_samemother==0', 'partial genmatched + different mothers', ROOT.kGreen+4)

  cand_mumupi_samemother = BackgroundSource('cand_ismatched==1 && cand_mumupi_samemother==1', 'B cascade decay', ROOT.kRed-10)
  #cand_differentmothers = BackgroundSource('cand_ismatched==1 && cand_mumupi_samemother==0', 'combinatorial', ROOT.kGreen-10)

  #cand_unmatched = BackgroundSource('cand_ismatched==0 && cand_ispartially_matched==0', 'non-genmatched', ROOT.kGray)
  cand_unmatched = BackgroundSource('cand_isnotmatched==1', 'non-genmatched', ROOT.kGray)
  cand_partiallymatched = BackgroundSource('cand_ispartialmatched==1', 'partially genmatched', ROOT.kCyan-10)
  cand_combinatorial = BackgroundSource('cand_ispartialmatched==1 || (cand_isfullmatched==1 && cand_mumupi_samemother==0)', 'combinatorial', ROOT.kGreen-10)
  cand_cascade = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1', 'B cascade decay', ROOT.kRed-10)

  cand_directBdecay = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_directBdecay==1', 'direct B decay', ROOT.kRed-9)
  cand_BtoDdecay = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_BtoDdecay==1', 'B->D decay', ROOT.kRed-4)
  cand_BtoJPsidecay = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_BtoJPsidecay==1', 'B->JPsi decay', ROOT.kPink+1)
  #cand_other = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_directBdecay==0 && mu0_BtoDdecay==0 && mu0_BtoJPsidecay==0 && mu0_nonBdecay==0', 'B cascade decay', ROOT.kOrange+0)
  cand_other = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_directBdecay==0 && mu0_BtoDdecay==0 && mu0_nonBdecay==0', 'B cascade decay', ROOT.kRed+2)
  cand_nonBdecay = BackgroundSource('cand_isfullmatched==1 && cand_mumupi_samemother==1 && mu0_nonBdecay==1', 'non B decay', ROOT.kRed+4)

  #cand_differentmothers = BackgroundSource('cand_ismatched==1 && cand_samemother!=1', 'genmatched + different mothers', ROOT.kRed-10)
  #cand_samemother = BackgroundSource('cand_ismatched==1 && cand_samemother==1', 'genmatched + same mother', ROOT.kGreen-10)

  mu0_matching = [
    mu0_nonmatched,
    mu0_matched,
    mu0_fake,
    ]

  mu_matching = [
    mu_nonmatched,
    mu_matched,
    mu_fake,
    ]

  pi_matching = [
    pi_nonmatched,
    pi_matched,
    pi_fake,
    ]

  mu0_sources = [
    mu0_directBdecay,
    mu0_BtoDdecay,
    mu0_BtoDtoKdecay,
    mu0_Btotaudecay,
    mu0_BtoJPsidecay,
    mu0_BtoDexciteddecay,
    mu0_other,
    mu0_nonBdecay,
    mu0_nonmatched,
    ]

  mu0_sources_cascade_fullmatched = [
    mu0_directBdecay_fullmatched,
    mu0_BtoDdecay_fullmatched,
    mu0_cascade_fullmatched,
    mu0_nonBdecay_fullmatched,
    ]

  mu_sources = [
    mu_directBdecay,
    mu_BtoDdecay,
    mu_BtoDtoKdecay,
    mu_Btotaudecay,
    mu_BtoJPsidecay,
    mu_BtoDexciteddecay,
    mu_other,
    mu_nonBdecay,
    mu_nonmatched,
    ]

  mu_sources_cascade_fullmatched = [
    mu_directBdecay_fullmatched,
    mu_BtoDdecay_fullmatched,
    mu_cascade_fullmatched,
    mu_nonBdecay_fullmatched,
    ]

  pi_sources = [
    pi_directBdecay,
    pi_BtoDdecay,
    pi_BtoDtoKdecay,
    pi_Btotaudecay,
    pi_BtoJPsidecay,
    pi_BtoDexciteddecay,
    pi_other,
    pi_nonBdecay,
    pi_nonmatched,
    ]

  pi_sources_cascade_fullmatched = [
    pi_directBdecay_fullmatched,
    pi_BtoDdecay_fullmatched,
    pi_cascade_fullmatched,
    pi_nonBdecay_fullmatched,
    ]

  pi_sources_cascade = [
    #pi_nonmatched,
    pi_directBdecay,
    pi_BtoDdecay,
    pi_cascade,
    pi_nonBdecay,
    ]

  cand_sources_general = [
    cand_cascade,
    cand_combinatorial,
    #cand_partiallymatched,
    cand_unmatched,
    #cand_mumupi_samemother,
    #cand_differentmothers,
    #cand_nonmatched,
    #cand_differentmothers,
    #cand_mu0mu_samemother,
    #cand_mu0pi_samemother,
    #cand_mupi_samemother,
    #cand_mumupi_samemother,
    #cand_nonmatched,
    ]

  cand_sources_cascade = [
    cand_directBdecay,
    cand_BtoDdecay,
    #cand_BtoJPsidecay,
    cand_other,
    cand_nonBdecay,
    #cand_combinatorial,
    #cand_unmatched,
    ]

  cand_sources_combinatorial = [
    #cand_cascade,
    #cand_partiallymatched,
    cand_partial_mu0mu_samemother,
    cand_partial_mu0pi_samemother,
    cand_partial_mupi_samemother,
    cand_partial_differentmothers,
    cand_mu0mu_samemother,
    cand_mu0pi_samemother,
    cand_mupi_samemother,
    cand_differentmothers,
    #cand_unmatched,
    ]

  hnl_mass = Quantity(name_flat='hnl_mass', label='hnl_mass', title='#mu#mu#pi invariant mass [GeV]', nbins=80, bin_min=0, bin_max=5.4)
  sv_lxy = Quantity(name_flat='sv_lxy', label='sv_lxy', title='SV l_{xy} [cm]', nbins=80, bin_min=0, bin_max=5)
  mu0_pt = Quantity(name_flat='mu0_pt', label='mu0_pt', title='primary #mu p_{T} [GeV]', nbins=80, bin_min=0, bin_max=30)
  mu0_eta = Quantity(name_flat='mu0_eta', label='mu0_eta', title='primary #mu #eta', nbins=80, bin_min=-2.1, bin_max=2.1)
  mu0_dxysig = Quantity(name_flat='fabs(mu0_dxysig)', label='mu0_dxysig', title='primary #mu dxysig', nbins=80, bin_min=0, bin_max=150)
  mu_pt = Quantity(name_flat='mu_pt', label='mu_pt', title='displaced #mu p_{T} [GeV]', nbins=80, bin_min=0, bin_max=30)
  mu_eta = Quantity(name_flat='mu_eta', label='mu_eta', title='displaced #mu #eta', nbins=80, bin_min=-2.1, bin_max=2.1)
  mu_dxysig = Quantity(name_flat='fabs(mu_dxysig)', label='mu_dxysig', title='displaced #mu dxysig', nbins=80, bin_min=0, bin_max=150)
  pi_pt = Quantity(name_flat='pi_pt', label='pi_pt', title='displaced #pi p_{T} [GeV]', nbins=80, bin_min=0, bin_max=15)
  pi_eta = Quantity(name_flat='pi_eta', label='pi_eta', title='displaced #pi #eta', nbins=80, bin_min=-2.1, bin_max=2.1)
  pi_dxysig = Quantity(name_flat='fabs(pi_dxysig)', label='pi_dxysig', title='displaced #pi dxysig', nbins=80, bin_min=0, bin_max=150)

  quantities = [
    hnl_mass,
    sv_lxy,
    #mu0_pt,
    #mu0_eta,
    #mu0_dxysig,
    #mu_pt,
    #mu_eta,
    #mu_dxysig,
    #pi_pt,
    #pi_eta,
    #pi_dxysig,
    ]

  for quantity in quantities:

    plotter = BackgroundPlotter(
        quantity = quantity,
        qcd_files = qcd_files,
        white_list = white_list,
        categories = categories,
        baseline_selection = baseline_selection,
        add_weight_hlt = add_weight_hlt,
        add_weight_pu = add_weight_pu,
        weight_hlt = weight_hlt,
        weight_puqcd = weight_puqcd,
        CMS_tag = CMS_tag,
        )

    plotter.plotStackHistogram(sources=mu0_matching, title='Primary muon', plot_label='matching_{}_mu0'.format(quantity.label))
    plotter.plotStackHistogram(sources=mu_matching, title='Displaced muon', plot_label='matching_{}_mu'.format(quantity.label))
    plotter.plotStackHistogram(sources=pi_matching, title='Displaced pion', plot_label='matching_{}_pi'.format(quantity.label))

    plotter.plotStackHistogram(sources=mu0_sources, title='Primary muon', plot_label='sources_{}_mu0'.format(quantity.label))
    plotter.plotStackHistogram(sources=mu0_sources_cascade_fullmatched, title='Primary muon', plot_label='sources_{}_mu0_cascade_fullmatched'.format(quantity.label))
    plotter.plotStackHistogram(sources=mu_sources, title='Displaced muon', plot_label='sources_{}_mu'.format(quantity.label))
    plotter.plotStackHistogram(sources=mu_sources_cascade_fullmatched, title='Displaced muon', plot_label='sources_{}_mu_cascade_fullmatched'.format(quantity.label))
    plotter.plotStackHistogram(sources=pi_sources, title='Displaced pion', plot_label='sources_{}_pi'.format(quantity.label))
    plotter.plotStackHistogram(sources=pi_sources_cascade_fullmatched, title='Displaced pion', plot_label='sources_{}_pi_cascade_fullmatched'.format(quantity.label))
    plotter.plotStackHistogram(sources=pi_sources_cascade, title='Displaced pion', plot_label='sources_{}_pi_cascade'.format(quantity.label))
    plotter.plotStackHistogram(sources=cand_sources_general, title='#mu#mu#pi candidate', plot_label='sources_{}_cand_general'.format(quantity.label))
    plotter.plotStackHistogram(sources=cand_sources_cascade, title='#mu#mu#pi candidate', plot_label='sources_{}_cand_cascade'.format(quantity.label))
    plotter.plotStackHistogram(sources=cand_sources_combinatorial, title='#mu#mu#pi candidate', plot_label='sources_{}_cand_combinatorial'.format(quantity.label))

