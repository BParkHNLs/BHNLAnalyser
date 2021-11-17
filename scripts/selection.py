import ROOT
import numpy as np
import os
from os import path
from tools import Tools
import sys
sys.path.append('../objects')
from samples import signal_samples, qcd_samples
from categories import Category
from baseline_selection import selection
from qcd_white_list import white_list
from quantity import Quantity as Qte


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


class Selection(object):
  def __init__(self, signal_files, qcd_files, quantity, category=None, dirlabel='', preexisting_selection=None, proposed_cut=None, white_list=''):
    self.tools                 = Tools()
    self.signal_files          = signal_files
    self.qcd_files             = qcd_files
    self.quantity              = quantity
    self.category              = category
    self.dirlabel              = dirlabel
    self.preexisting_selection = preexisting_selection
    self.proposed_cut          = proposed_cut
    self.white_list            = white_list



  def createOutDir(self, outputdir):
    if not path.exists(outputdir):
      os.system('mkdir {}'.format(outputdir))


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
      # define the mass window
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)
      print '{} {}'.format(signal_file.mass-2*signal_file.resolution, signal_file.mass+2*signal_file.resolution)

      # signal efficiency
      g_sig = ROOT.TGraph()

      #selection_sig_ini = 'ismatched == 1 && trgmu_softid==1 && mu_looseid==1 && {}'.format(self.category.definition_flat) 
      selection_sig_ini = 'ismatched == 1 && {}'.format(self.category.definition_flat) 
      if self.preexisting_selection !=None: selection_sig_ini += ' && {}'.format(self.getSelectionString())
      #print 'selection_sig_ini ',selection_sig_ini
      f_sig = self.tools.getRootFile(signal_file.filename, with_ext=False) 
      initial_sig_entries = self.tools.createHisto(f_sig, 'signal_tree', self.hnl_mass, selection=selection_sig_ini).Integral()
      #print 'initial_sig_entries {} '.format(initial_sig_entries)

      for idx, cut in enumerate(points):
        #selection_sig_sel = 'ismatched == 1 && trgmu_softid==1 && mu_looseid==1 && {}'.format(self.category.definition_flat) 
        selection_sig_sel = 'ismatched == 1 && {}'.format(self.category.definition_flat) 
        if self.preexisting_selection !=None: selection_sig_sel += ' && {}'.format(self.getSelectionString())
        selection_sig_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        selected_sig_entries = self.tools.createHisto(f_sig, 'signal_tree', self.hnl_mass, selection=selection_sig_sel).Integral()
        #print 'selection_sig_sel {} {}'.format(selection_sig_sel, selected_sig_entries)
        g_sig.SetPoint(idx, cut, self.getEff(selected_sig_entries, initial_sig_entries))
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
      gs_sig[0].GetYaxis().SetRangeUser(0, 1)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('APL')
        else:
          g_sig.Draw('PL')

      # background rejection
      g_bkg = ROOT.TGraph()
      #selection_bkg_ini = 'trgmu_softid==1 && mu_looseid==1 && {}'.format(self.category.definition_flat)
      selection_bkg_ini = '{}'.format(self.category.definition_flat)
      if self.preexisting_selection !=None: selection_bkg_ini += ' && {}'.format(self.getSelectionString())
      #print 'selection_bkg_ini ',selection_bkg_ini
      initial_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_ini).Integral() # no weight for the moment
      #print 'initial_bkg_entries ',initial_bkg_entries

      for idx, cut in enumerate(points):
        #selection_bkg_sel = 'trgmu_softid==1 && mu_looseid==1 && {}'.format(self.category.definition_flat)
        selection_bkg_sel = '{}'.format(self.category.definition_flat)
        if self.preexisting_selection !=None: selection_bkg_sel += ' && {}'.format(self.getSelectionString())
        selection_bkg_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        #print 'selection_bkg_sel ',selection_bkg_sel
        selected_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_sel).Integral() # no weight for the moment
        g_bkg.SetPoint(idx, cut, 1-self.getEff(selected_bkg_entries, initial_bkg_entries))
        g_bkg.SetLineWidth(2)
        g_bkg.SetLineColor(signal_file.colour)
        g_bkg.SetMarkerColor(signal_file.colour)
        g_bkg.SetMarkerSize(2)
        g_bkg.SetMarkerStyle(45)

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


  def printCutflowLine(self, printNum=False):
    cutflow_line = '{qte} {log} {cut} & '.format(
        qte = self.quantity.title, 
        log = self.quantity.logic, 
        cut = self.proposed_cut, 
    )

    for signal_file in self.signal_files:
      self.hnl_mass = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_file.mass-2*signal_file.resolution, bin_max=signal_file.mass+2*signal_file.resolution)

      f_sig = self.tools.getRootFile(signal_file.filename, with_ext=False) 
      #selection_sig_ini = 'ismatched == 1 && trgmu_softid==1 && mu_looseid==1 &&  {}'.format(self.category.definition_flat) 
      selection_sig_ini = 'ismatched == 1 && {}'.format(self.category.definition_flat) 
      if self.preexisting_selection !=None: selection_sig_ini += ' && {}'.format(self.getSelectionString())
      initial_sig_entries = self.tools.createHisto(f_sig, 'signal_tree', self.hnl_mass, selection=selection_sig_ini).Integral()

      selection_sig_sel = selection_sig_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, self.proposed_cut)
      selected_sig_entries = self.tools.createHisto(f_sig, 'signal_tree', self.hnl_mass, selection=selection_sig_sel).Integral()

      #selection_bkg_ini = 'trgmu_softid==1 && mu_looseid==1 && {}'.format(self.category.definition_flat)
      selection_bkg_ini = '{}'.format(self.category.definition_flat)
      if self.preexisting_selection !=None: selection_bkg_ini += ' && {}'.format(self.getSelectionString())
      initial_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_ini).Integral() # no weight for the moment

      selection_bkg_sel = selection_bkg_ini + ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, self.proposed_cut)
      selected_bkg_entries = self.tools.createWeightedHistoQCDMC(self.qcd_files, white_list=self.white_list, quantity=self.hnl_mass, selection=selection_bkg_sel).Integral() # no weight for the moment

      if printNum: print 'sig m{} ctau{} {} {}'.format(signal_file.mass, signal_file.ctau, initial_sig_entries, selected_sig_entries)
      if printNum: print 'bkg m{} {} {}'.format(signal_file.mass, initial_bkg_entries, selected_bkg_entries)

      cutflow_line += ' -{sig_per}\% & -{bkg_per}\% '.format(
        sig_per = round((1 - selected_sig_entries / initial_sig_entries)*100, 2) if initial_sig_entries!=0. else '-', 
        bkg_per = round((1 - selected_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
      )

    cutflow_line+= '\\'

    print cutflow_line


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  signal_files = signal_samples['central_V09_06Nov21_benchmark']
  qcd_files = qcd_samples['V09_06Nov21']

  white_list = white_list['20to300']

  hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.9999, 1) 
  pi_pt = Quantity('pi_pt', 'pi_pt', 'displaced pion pT', '>', 'GeV', 0, 6)
  mu_dxy = Quantity('abs(mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 0.3)
  pi_dxy = Quantity('abs(pi_dxy)', 'pi_dxy', 'displaced pion |IP| on xy', '>', 'cm', 0, 1)
  cos_theta_star_pion = Quantity('fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.5, 1) 
  trgmu_softid = Quantity('trgmu_softid', 'trgmu_softid', 'trigger muon soft ID', '==', '', 0, 1)
  mu_looseid = Quantity('mu_looseid', 'mu_looseid', 'muon loose ID', '==', '', 0, 1)
  mu_customid = Quantity('mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))', 'mu_customid', 'muon custom ID', '==', '', 0, 1)
  b_mass = Quantity('b_mass', 'b_mass', '#mu#mu#pi mass', '>', 'GeV', 1, 7) 
  sv_lxy_sig = Quantity('sv_lxysig', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 100)
  mu_numberofpixellayers = Quantity('mu_numberofpixellayers', 'mu_numberofpixellayers', 'mu_numberofpixellayers', '<', '', 0, 6)
  mu_numberoftrackerlayers = Quantity('mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', '<', '', 0, 18)
  

  printSelectionString = True

  dirlabel = 'benchmark'

  ### CATEGORY lxygt20_OS ###

  category_lxygt20_OS = Category(
      label = 'lxygt20_OS',
      title = 'l_{xy}>20cm, OS',
      definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxygt20_OS = [] 

  print 'category lxygt20_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt20_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt20_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99998 # this can be tighten
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.8
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.25
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.5
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt20_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  #cut_cos_theta_star_pion = 0.7
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=cos_theta_star_pion, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_OS, proposed_cut=cut_cos_theta_star_pion)
  ##selection.getScanGraph(npoints=30)
  ##selection.printCutflowLine()
  ##selection_sequential_lxygt20_OS.append(SelectedQuantity(cos_theta_star_pion, cut_cos_theta_star_pion))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxygt20_SS ###

  category_lxygt20_SS = Category(
      label = 'lxygt20_SS',
      title = 'l_{xy}>20cm, SS',
      definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxygt20_SS = [] 

  print 'category lxygt20_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt20_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt20_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  cos_theta_star_pion = Quantity('fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.5, 1) 
  selection_sequential_lxygt20_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.8
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt20_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  ##mu_dxy = Quantity('abs(mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 1)
  cut_mu_dxy = 0.25
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt20_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  ##cos_theta_star_pion = Quantity('fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.3, 1) 
  ##cut_cos_theta_star_pion = 0.7
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=cos_theta_star_pion, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_cos_theta_star_pion)
  ##selection.getScanGraph(npoints=30)
  ##selection.printCutflowLine()
  ##selection_sequential_lxygt20_SS.append(SelectedQuantity(cos_theta_star_pion, cut_cos_theta_star_pion))

  cut_pi_dxy = 0.5
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt20_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt20_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxygt10_OS ###

  category_lxygt10_OS = Category(
      label = 'lxygt10_OS',
      title = 'l_{xy}>10cm, OS',
      definition_flat = 'sv_lxy>10 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxygt10_OS = [] 

  print 'category lxygt10_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt10_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt10_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.3
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt10_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxygt10_SS ###

  category_lxygt10_SS = Category(
      label = 'lxygt10_SS',
      title = 'l_{xy}>10cm, SS',
      definition_flat = 'sv_lxy>10 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxygt10_SS = [] 

  print 'category lxygt10_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt10_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt10_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.3
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt10_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt10_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt10_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxygt5_OS ###

  category_lxygt5_OS = Category(
      label = 'lxygt5_OS',
      title = 'l_{xy}>5cm, OS',
      definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxygt5_OS = [] 

  print 'category lxygt5_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt5_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.3
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt5_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxygt5_SS ###

  category_lxygt5_SS = Category(
      label = 'lxygt5_SS',
      title = 'l_{xy}>5cm, SS',
      definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxygt5_SS = [] 

  print 'category lxygt5_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxygt5_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.3
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxygt5_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxygt5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxygt5_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxygt5_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxy3to20_OS ###

  category_lxy3to20_OS = Category(
      label = 'lxy3to20_OS',
      title = '(3<l_{xy}<=20)cm, OS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy3to20_OS = [] 

  print 'category lxy3to20_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to20_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to20_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99994
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.0
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.07
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxy3to20_SS ###

  category_lxy3to20_SS = Category(
      label = 'lxy3to20_SS',
      title = '(3<l_{xy}<=20)cm, SS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy3to20_SS = [] 

  print 'category lxy3to20_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to20_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to20_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99994
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1.0
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.07
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.2
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to20_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to20_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxy3to10_OS ###

  category_lxy3to10_OS = Category(
      label = 'lxy3to10_OS',
      title = '(3<l_{xy}<=10)cm, OS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy3to10_OS = [] 

  print 'category lxy3to10_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to10_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to10_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99991
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.06
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.17
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

  ### CATEGORY lxy3to10_SS ###

  category_lxy3to10_SS = Category(
      label = 'lxy3to10_SS',
      title = '(3<l_{xy}<=10)cm, SS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy3to10_SS = [] 

  print 'category lxy3to10_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to10_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to10_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.99991
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.06
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.17
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to10_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to10_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy3to5_OS ###

  category_lxy3to5_OS = Category(
      label = 'lxy3to5_OS',
      title = '(3<l_{xy}<=5)cm, OS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy3to5_OS = [] 

  print 'category lxy3to5_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to5_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.9990, 1) 
  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy3to5_SS ###

  category_lxy3to5_SS = Category(
      label = 'lxy3to5_SS',
      title = '(3<l_{xy}<=5)cm, SS',
      definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy3to5_SS = [] 

  print 'category lxy3to5_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy3to5_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy3to5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.9990, 1) 
  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy3to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy3to5_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy3to5_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to20_OS ###

  category_lxy1to20_OS = Category(
      label = 'lxy1to20_OS',
      title = '(1<l_{xy}<=20)cm, OS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy1to20_OS = [] 

  print 'category lxy1to20_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to20_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to20_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to20_SS ###

  category_lxy1to20_SS = Category(
      label = 'lxy1to20_SS',
      title = '(1<l_{xy}<=20)cm, SS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy1to20_SS = [] 

  print 'category lxy1to20_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to20_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to20_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to20_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to20_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to10_OS ###

  category_lxy1to10_OS = Category(
      label = 'lxy1to10_OS',
      title = '(1<l_{xy}<=10)cm, OS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy1to10_OS = [] 

  print 'category lxy1to10_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to10_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to10_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to10_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to10_SS ###

  category_lxy1to10_SS = Category(
      label = 'lxy1to10_SS',
      title = '(1<l_{xy}<=10)cm, SS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy1to10_SS = [] 

  print 'category lxy1to10_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to10_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to10_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.04
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to10_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to10_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to10_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to5_OS ###

  category_lxy1to5_OS = Category(
      label = 'lxy1to5_OS',
      title = '(1<l_{xy}<=5)cm, OS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy1to5_OS = [] 

  print 'category lxy1to5_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to5_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to5_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to5_SS ###

  category_lxy1to5_SS = Category(
      label = 'lxy1to5_SS',
      title = '(1<l_{xy}<=5)cm, SS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy1to5_SS = [] 

  print 'category lxy1to5_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to5_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9997
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to5_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to5_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to5_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to3_OS ###

  category_lxy1to3_OS = Category(
      label = 'lxy1to3_OS',
      title = '(1<l_{xy}<=3)cm, OS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy1to3_OS = [] 

  print 'category lxy1to3_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to3_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to3_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9994
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.02
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to3_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy1to3_SS ###

  category_lxy1to3_SS = Category(
      label = 'lxy1to3_SS',
      title = '(1<l_{xy}<=3)cm, SS',
      definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy1to3_SS = [] 

  print 'category lxy1to3_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy1to3_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy1to3_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  cut_hnl_cos2D = 0.9994
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_mu_dxy = 0.02
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_mu_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  cut_pi_dxy = 0.03
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_pi_dxy)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine(printNum=False)
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 70
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy1to3_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy1to3_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy1to3_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy0to1_OS ###

  category_lxy0to1_OS = Category(
      label = 'lxy0to1_OS',
      title = '(0<l_{xy}<=1)cm, OS',
      definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
  )

  selection_sequential_lxy0to1_OS = [] 

  print 'category lxy0to1_OS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy0to1_OS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_OS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.99, 1) 
  cut_hnl_cos2D = 0.998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_OS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_OS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  ##cut_mu_dxy = 0.02
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_dxy, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_mu_dxy)
  ##selection.getScanGraph(npoints=30)
  ##selection.printCutflowLine()
  ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(mu_dxy, cut_mu_dxy))

  ##cut_pi_dxy = 0.03
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_dxy, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_pi_dxy)
  ##selection.getScanGraph(npoints=30)
  ##selection.printCutflowLine(printNum=False)
  ##selection_sequential_lxy0to1_OS.append(SelectedQuantity(pi_dxy, cut_pi_dxy))

  cut_sv_lxy_sig = 25
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy0to1_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_OS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_OS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()

### CATEGORY lxy0to1_SS ###

  category_lxy0to1_SS = Category(
      label = 'lxy0to1_SS',
      title = '(0<l_{xy}<=1)cm, SS',
      definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
  )

  selection_sequential_lxy0to1_SS = [] 

  print 'category lxy0to1_SS'

  cut_trgmu_softid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=trgmu_softid, category=category_lxy0to1_SS, dirlabel=dirlabel, proposed_cut=cut_trgmu_softid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_SS.append(SelectedQuantity(trgmu_softid, cut_trgmu_softid))

  cut_mu_looseid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_looseid, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_looseid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

  cut_mu_customid = 1
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_customid, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_mu_customid)
  #selection.getScanGraph(npoints=2)
  #selection.printCutflowLine()
  ##selection_sequential_lxy0to1_SS.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

  hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.99, 1) 
  cut_hnl_cos2D = 0.998
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=hnl_cos2D, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_hnl_cos2D)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_SS.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_pi_pt = 0.9
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=pi_pt, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_pi_pt)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_SS.append(SelectedQuantity(pi_pt, cut_pi_pt))

  cut_sv_lxy_sig = 25
  selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxy0to1_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential_lxy0to1_SS, proposed_cut=cut_sv_lxy_sig)
  #selection.getScanGraph(npoints=30)
  #selection.printCutflowLine()
  selection_sequential_lxy0to1_SS.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  if printSelectionString: print selection.getSelectionString()



  ##cut_mu_numberofpixellayers = 2 
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_numberofpixellayers, category=category_lxygt20_SS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_numberofpixellayers)
  ##selection.getScanGraph(npoints=7)
  ##selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(mu_numberofpixellayers, cut_mu_numberofpixellayers))

  ##cut_mu_numberoftrackerlayers = 13 
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=mu_numberoftrackerlayers, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_numberoftrackerlayers)
  ##selection.getScanGraph(npoints=19)
  ##selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(mu_numberoftrackerlayers, cut_mu_numberoftrackerlayers))

  ##cut_sv_lxy_sig = 30
  ##selection = Selection(signal_files=signal_files, qcd_files=qcd_files, white_list=white_list, quantity=sv_lxy_sig, category=category_lxygt20_OS, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_sv_lxy_sig)
  ##selection.getScanGraph(npoints=30)
  ##selection.printCutflowLine()
  ##selection_sequential.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

