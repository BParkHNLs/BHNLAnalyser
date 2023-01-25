import numpy as np
import ROOT
import os
import os.path
from os import path
from tools import Tools



class Quantity(object):
  def __init__(self, name_nano='', name_flat='', label='', title='', logic='', units='', binMin=0., binMax=0.):
    self.name_nano = name_nano
    self.name_flat = name_flat
    self.label = label
    self.title = title
    self.logic = logic
    self.units = units
    self.binMin = binMin
    self.binMax = binMax
                              

class Selection(Tools):
  def __init__(self, files, quantity, preexisting_selection=None, npoints=50, sample_type='nano', write_cut_analysis=False, proposed_cut=None):
    self.tools                 = Tools()
    self.files                 = files
    self.quantity              = quantity
    self.preexisting_selection = preexisting_selection
    self.npoints               = npoints
    self.sample_type           = sample_type
    self.write_cut_analysis    = write_cut_analysis
    self.proposed_cut          = proposed_cut


    # efficiency will be assessed based on the mass distributions
    self.hnl_mass = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass', 'hnl_mass', '', binMin=0, binMax=10000) # with very loose selection some candidates have very large mupi invariant mass
    
    # baseline selection
    #self.baseline_selection = ' && '.join(['{}==0'.format('hnl_charge' if str_=='sig' else 'b_hnl_charge')])
    #                                  , 'b_mass<6.35']) 


  def createOutDir(self, outputdir):
    if not path.exists(outputdir):
      os.system('mkdir {}'.format(outputdir))


  def getLabel(self, str_):
    new_str_ = str_
    if 'abs(' in str_: new_str_ = new_str_.replace('abs(', '')
    if ')'in str_: new_str_ = new_str_.replace(')', '')
    return new_str_


  def getPreselectionString(self, str_):
    '''
    function to write the already-applied cuts in a string
    str_ is to make the distinction between nano and matched samples,
    where variables are named differently
    '''
    preselection_str = []
    for item, _ in enumerate(self.preexisting_selection):
      name_variable = self.preexisting_selection[item].quantity.name_nano if self.sample_type == 'nano' else self.preexisting_selection[item].quantity.name_flat
      preselection_str.append('{}{}{}'.format(name_variable,self.preexisting_selection[item].quantity.logic,self.preexisting_selection[item].chosen_cut))
    return ' && '.join(preselection_str)


  def createHisto(self, file_, str_, with_preexisting_selection, with_extra_selection, cut=0):
    '''
    str_ makes the difference between signal and background
    '''

    filename = file_.sample_name 
    tree_name = 'Events' if self.sample_type == 'nano' else 'signal_tree'

    cut_variable = self.quantity.name_nano if self.sample_type == 'nano' else self.quantity.name_flat

    if self.sample_type == 'nano':
      baseline_selection = 'BToMuMuPi_isMatched==1' if str_=='sig' else 'BToMuMuPi_isMatched>-99'
    else:
      baseline_selection = 'ismatched==1' if str_=='sig' else 'ismatched>-99'
    
    c = ROOT.TCanvas()
    f = ROOT.TFile.Open(filename, 'READ')
    tree = f.Get(tree_name)
    
    # elements to fill the histogram
    distr_todraw = self.hnl_mass.name_nano if self.sample_type == 'nano' else self.hnl_mass.name_flat

    preselection = baseline_selection 
    if with_preexisting_selection and self.preexisting_selection != None: 
      preselection += ' && ' + self.getPreselectionString(str_)

    if with_extra_selection:
      preselection += ' && {}{}{}'.format(cut_variable, self.quantity.logic, cut)

    #print preselection
    
    hist = ROOT.TH1D('hist', 'hist', 1500, self.hnl_mass.binMin, self.hnl_mass.binMax)
    tree.Draw('{}>>hist'.format(distr_todraw), preselection)

    hist.SetDirectory(0)
    return hist


  def getEff(self, entries_selected, entries_initial):
    eff = entries_selected / entries_initial if entries_initial!=0 else 0.
    return eff


  def getScanGraph(self):

    '''
    plot the signal efficiency and background rejection as a function of the cut applied on the quantity of interest
    possibility to draw a line at the chosen cut
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints) 

    canv = self.tools.createTCanvas('canv', 900, 800)
    canv.SetGrid()
    
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
    pad_up.SetBottomMargin(0.1)
    pad_up.Draw()
    canv.cd()
    pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
    pad_down.SetBottomMargin(0.15)
    pad_down.Draw()
    canv.cd()
    pad_leg = ROOT.TPad("pad_leg","pad_leg",0.5,0,1,0.25)
    pad_leg.SetBottomMargin(0.15)

    if self.write_cut_analysis:
      pad_down.cd()

    gs_sig = []
    gs_bkg = []

    # signal efficiency
    for ifile, file_ in enumerate(self.files):
      if file_.process != 'signal': continue
      g_sig = ROOT.TGraph()

      initial_sig_entries = self.createHisto(file_, 'sig', True, False).GetEntries()

      for idx, cut in enumerate(points):
        selected_sig_entries = self.createHisto(file_, 'sig', True, True, cut).GetEntries()
        g_sig.SetPoint(idx, cut, self.getEff(selected_sig_entries, initial_sig_entries))
        g_sig.SetLineWidth(0)
        g_sig.SetMarkerColor(ROOT.kOrange+1)
        g_sig.SetMarkerStyle(20+ifile)
  
      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetRangeUser(0, 1)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('AP')
        else:
          g_sig.Draw('P')

    # draw background rejection
    for ifile, file_ in enumerate(self.files):
      if file_.process != 'background': continue
      g_bkg = ROOT.TGraph()
      initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()

      for idx, cut in enumerate(points):
        selected_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, cut).GetEntries()
        g_bkg.SetPoint(idx, cut, 1-self.getEff(selected_bkg_entries, initial_bkg_entries))
        g_bkg.SetLineWidth(0)
        g_bkg.SetMarkerColor(ROOT.kBlue+2)
        g_bkg.SetMarkerStyle(20)

      g_bkg.Draw('P')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, 0, self.proposed_cut, 1)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      line.Draw('same')
          
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    for ifile, file_ in enumerate(self.files):
      if file_.process == 'signal':
        # only consistent if the background sample is fed before the signal samples
        legend.AddEntry(gs_sig[ifile-1], 'sig efficiency ({}GeV, {}mm)'.format(file_.signal_mass, file_.signal_ctau))
      else:
        legend.AddEntry(g_bkg, 'bkg rejection')
    legend.Draw()

    # write cutflow information
    if self.write_cut_analysis:
      pad_down.cd()
      proposed_sig_entries = self.createHisto(self.files[0], 'sig', True, True, self.proposed_cut).GetEntries()
      proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()

      box1 = self.tools.getTextBox("NDC", 0.05, 0.7, 0.4, 0.98, 'Additional proposed cut: {}{}{}'.format(self.quantity.label, self.quantity.logic, self.proposed_cut), ROOT.kRed)
      box2 = self.tools.getTextBox("brNDC", 0.02, 0.3, 0.08, 0.4, '{}GeV'.format(self.files[0].signal_mass), ROOT.kBlack)
      box3 = self.tools.getTextBox("brNDC", 0.08, 0.45, 0.4, 0.65, 'N_matched ini: {}'.format(initial_sig_entries), ROOT.kOrange+1)
      box4 = self.tools.getTextBox("brNDC", 0.08, 0.25, 0.4, 0.55, 'N_matched new: {}'.format(proposed_sig_entries), ROOT.kOrange+1)
      box5 = self.tools.getTextBox("brNDC", 0.35, 0.35, 0.5, 0.65, '==> -{}%'.format(round((1 - proposed_sig_entries / initial_sig_entries)*100, 2)), ROOT.kOrange+1)
      box6 = self.tools.getTextBox("brNDC", 0.08, 0.15, 0.4, 0.3, 'N_nano ini: {}'.format(initial_bkg_entries), ROOT.kBlue+2)
      box7 = self.tools.getTextBox("brNDC", 0.08, 0., 0.4, 0.15, 'N_nano new: {}'.format(proposed_bkg_entries), ROOT.kBlue+2)
      box8 = self.tools.getTextBox("brNDC", 0.35, 0., 0.5, 0.3, '==> -{}%'.format(round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2)), ROOT.kBlue+2)

      box1.Draw('same')
      box2.Draw('same')
      box3.Draw('same')
      box4.Draw('same')
      box5.Draw('same')
      box6.Draw('same')
      box7.Draw('same')
      box8.Draw('same')

    canv.cd()
    self.createOutDir('myPlots/preselection')
    canv.SaveAs('myPlots/preselection/scan_{}.png'.format(self.getLabel(self.quantity.label)))


  def getROCGraph(self):

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints)

    canv = self.tools.createTCanvas('canv', 'canv', 1200, 1000)
    g = ROOT.TGraph2D()
    g.SetTitle('Cut on {} ({})'.format(self.quantity.title, self.quantity.logic))

    # do it for first file only
    file_ = self.files[0]
    
    initial_sig_entries = self.createHisto(file_, 'sig', True, False).GetEntries()
    initial_bkg_entries = self.createHisto(file_, 'bkg', True, False).GetEntries()

    for idx, cut in enumerate(points):
      selected_sig_entries = self.createHisto(file_, 'sig', True, True, cut).GetEntries()
      selected_bkg_entries = self.createHisto(file_, 'bkg', True, True, cut).GetEntries()
      g.SetPoint(idx, 1-self.getEff(selected_bkg_entries, initial_bkg_entries), self.getEff(selected_sig_entries, initial_sig_entries), cut)
    
    g.GetXaxis().SetTitle('background rejection')
    g.GetXaxis().SetLabelSize(0.038)
    g.GetXaxis().SetTitleSize(0.04)
    g.GetXaxis().SetTitleOffset(1.4)
    g.GetYaxis().SetTitle('signal efficiency')
    g.GetYaxis().SetLabelSize(0.038)
    g.GetYaxis().SetTitleSize(0.04)
    g.GetYaxis().SetTitleOffset(1.4)
    
    
    #ROOT.gStyle.SetPadRightMargin(5.5) 
    #g.GetZaxis().SetTitle('Cut on {}'.format(self.quantity.title))
    g.GetZaxis().SetTitle('Cut')
    #g.GetZaxis().SetLabelSize(0.038)
    g.GetZaxis().SetTitleSize(0.14)
    #g.GetZaxis().SetTitleOffset(1.2)

    ROOT.gStyle.SetPalette(53)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.7)
    g.Draw("PCOLZ")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ROOT.gPad.GetView().TopView()

    #ROOT.gPad.Update()
    #ROOT.gPad.Modified()
    ROOT.gStyle.SetPadRightMargin(1.5) 
   
    # draw diagonal line
    line = ROOT.TLine()
    line.SetLineColor(1)
    line.SetLineWidth(3)
    line.SetLineStyle(9)
    line.DrawLineNDC(0.1, 0.9, 0.9, 0.1)
    

    self.createOutDir('myPlots/preselection')
    canv.SaveAs('myPlots/preselection/roc_{}.png'.format(self.getLabel(self.quantity.label)))


  def printCutflowLine(self):
    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()
    initial_sig_entries = self.createHisto(self.files[1], 'sig', True, False).GetEntries()
  
    proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()
    proposed_sig_entries = self.createHisto(self.files[1], 'sig', True, True, self.proposed_cut).GetEntries()

    #print 'sig1 {} {}'.format(int(initial_sig_entries), int(initial_bkg_entries))
    #print 'sig1 {} {}'.format(int(proposed_sig_entries), int(proposed_bkg_entries))

    if len(self.files)==1:
      cutflow_line = '{qte} {log} {cut} & -{sig_per}\% & -{bkg_per}\% \\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
      )

    elif len(self.files)==2:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()

      ##print 'sig2 {}'.format(int(initial_sig1_entries))
      ##print 'sig2 {}'.format(int(proposed_sig1_entries))

      cutflow_line = '{qte} {log} {cut} & -{sig0_per}\% & -{sig1_per}\% & -{bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig0_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
          sig1_per = round((1 - proposed_sig1_entries / initial_sig1_entries)*100, 2), 
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
      )

    else:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()
      initial_sig2_entries = self.createHisto(self.files[3], 'sig', True, False).GetEntries()
      proposed_sig2_entries = self.createHisto(self.files[3], 'sig', True, True, self.proposed_cut).GetEntries()

      #print 'sig2 {}'.format(int(initial_sig1_entries))
      #print 'sig2 {}'.format(int(proposed_sig1_entries))
      #print 'sig3 {}'.format(int(initial_sig2_entries))
      #print 'sig3 {}'.format(int(proposed_sig2_entries))

      cutflow_line = '{qte} {log} {cut} {unit} & -{sig0_per}\% & -{sig1_per}\% & -{sig2_per}\% & -{bkg_per}\% \\\ '.format(
      #cutflow_line = '{qte} {log} {cut} {unit} & -{sig0_per}\% & -{sig1_per}\% & -{bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          unit = self.quantity.units,
          sig0_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2) if initial_sig_entries!=0. else '-', 
          sig1_per = round((1 - proposed_sig1_entries / initial_sig1_entries)*100, 2) if initial_sig1_entries!=0. else '-',
          sig2_per = round((1 - proposed_sig2_entries / initial_sig2_entries)*100, 2) if initial_sig2_entries!=0. else '-',
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
      )

    cutflow_line = cutflow_line.replace('#eta', '$\eta$')
    cutflow_line = cutflow_line.replace('#mu#mu#pi', '$\mu\mu\pi$')
    cutflow_line = cutflow_line.replace('#mu#pi', '$\mu\pi$')
    cutflow_line = cutflow_line.replace('#Delta', '$\Delta$')
    cutflow_line = cutflow_line.replace('#pi', '$\pi$')
    cutflow_line = cutflow_line.replace('#Phi', '$\Phi$')
    print cutflow_line


  def printEfficiencyLine(self):
    # get the entries without cuts at all. Used for the cumulative efficiency
    nocut_bkg_entries = self.createHisto(self.files[0], 'bkg', False, False).GetEntries()
    nocut_sig_entries = self.createHisto(self.files[1], 'sig', False, False).GetEntries()

    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()
    initial_sig_entries = self.createHisto(self.files[1], 'sig', True, False).GetEntries()
  
    proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()
    proposed_sig_entries = self.createHisto(self.files[1], 'sig', True, True, self.proposed_cut).GetEntries()

    if len(self.files)==1:
      cutflow_line = '{qte} {log} {cut} & {sig_per}\% & {bkg_per}\% \\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig_per = round((proposed_sig_entries / initial_sig_entries)*100, 1), 
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1),
      )

    elif len(self.files)==2:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()

      cutflow_line = '{qte} {log} {cut} & {sig0_per}\% & {sig1_per}\% & {bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig0_per = round((proposed_sig_entries / initial_sig_entries)*100, 1), 
          sig1_per = round((proposed_sig1_entries / initial_sig1_entries)*100, 1), 
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1),
      )

    else:
      nocut_sig1_entries = self.createHisto(self.files[2], 'sig', False, False).GetEntries()
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()
      nocut_sig2_entries = self.createHisto(self.files[3], 'sig', False, False).GetEntries()
      initial_sig2_entries = self.createHisto(self.files[3], 'sig', True, False).GetEntries()
      proposed_sig2_entries = self.createHisto(self.files[3], 'sig', True, True, self.proposed_cut).GetEntries()

      cutflow_line = '{qte} {log} {cut} {unit} & {sig0_per}\% & {sig0_cumul}\% & {sig1_per}\% & {sig1_cumul}\% & {sig2_per}\% & {sig2_cumul}\% & {bkg_per}\% & {bkg_cumul}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          unit = self.quantity.units,
          sig0_per = round((proposed_sig_entries / initial_sig_entries)*100, 1) if initial_sig_entries!=0. else '-', 
          sig1_per = round((proposed_sig1_entries / initial_sig1_entries)*100, 1) if initial_sig1_entries!=0. else '-',
          sig2_per = round((proposed_sig2_entries / initial_sig2_entries)*100, 1) if initial_sig2_entries!=0. else '-',
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1) if initial_bkg_entries!=0. else '-',
          sig0_cumul = round((proposed_sig_entries / nocut_sig_entries)*100, 1),
          sig1_cumul = round((proposed_sig1_entries / nocut_sig1_entries)*100, 1),
          sig2_cumul = round((proposed_sig2_entries / nocut_sig2_entries)*100, 1),
          bkg_cumul = round((proposed_bkg_entries / nocut_bkg_entries)*100, 1),
      )

    cutflow_line = cutflow_line.replace('#eta', '$\eta$')
    cutflow_line = cutflow_line.replace('#mu#mu#pi', '$\mu\mu\pi$')
    cutflow_line = cutflow_line.replace('#mu#pi', '$\mu\pi$')
    cutflow_line = cutflow_line.replace('#Delta', '$\Delta$')
    cutflow_line = cutflow_line.replace('#pi', '$\pi$')
    cutflow_line = cutflow_line.replace('#Phi', '$\Phi$')
    print cutflow_line



class FileCollection(object): 
  '''
  this class allows to associate to the studied files the corresponding 
  mass and coupling points. Useful when scaning through e.g masses
  '''
  def __init__(self, sample_name='', process='', signal_mass='', signal_ctau=''):
    self.sample_name = sample_name
    self.process = process
    self.signal_mass = signal_mass
    self.signal_ctau = signal_ctau
    # add label for reweighted?
    if self.process not in ['signal', 'background']:
      raise RuntimeError("Unknown process. Please choose among ['signal', 'background']")


class PreselectedQuantity(object):
  '''
  this class will allow to impose a preselection when studying 
  the impact of the cut on a given parameter
  '''
  def __init__(self, quantity, chosen_cut):
    self.quantity   = quantity
    self.chosen_cut = chosen_cut


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  file_background = FileCollection(
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH1_Run2018A/merged/bparknano_looseselection_dsaonly.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/bparknano_data_1file_looseselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/flat_bparknano_data_1file_looseselection.root',
      sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/merged/bparknano_loosepreselection_1file.root', # used for AN-v3
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22.root', # for baseline selection cutflow
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/bparknano_data_nocutatall_1000evts.root', # get the loose selection efficiency
      process = 'background',
      )

  file_m3 = FileCollection(
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V37/mass3.0_ctau100.0/nanoFiles/merged/bparknano_looseselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/bparknano_loosepreselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/bparknano_nj1.root',
      sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_5files.root', # used for AN-v3, preselection
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root', # for baseline selection cutflow
      process = 'signal',
      signal_mass = 3,
      signal_ctau = 100,
      )

  file_m4p5 = FileCollection(
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V37/mass4.5_ctau1.0/nanoFiles/merged/bparknano_looseselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau1.0/nanoFiles/merged/bparknano_loosepreselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n15/bparknano_nj1.root',
      sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_2files.root', # used for AN-v3, preselection
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root', # for baseline selection cutflow 
      process = 'signal',
      signal_mass = 4.5,
      signal_ctau = 1,
      )

  file_m1 = FileCollection(
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n38/merged/bparknano.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V37/mass1.0_ctau1000.0/nanoFiles/merged/bparknano_looseselection.root',
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n324/bparknano_nj1.root',
      sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30files.root', # used for AN-v3, preselection
      #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root', # for baseline selection cutflow
      process = 'signal',
      signal_mass = 1,
      signal_ctau = 1000,
      )

  file_V25 = FileCollection(
      sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_looseselection_muononly_alltrgmu.root',
      process = 'signal',
      signal_mass = 3,
      signal_ctau = 2000,
      )
  
  files = []
  files.append(file_background)
  files.append(file_m1)
  files.append(file_m3)
  files.append(file_m4p5)
  #files.append(file_V25)


  #trg_cond = Quantity('((BToMuMuPi_trg_mu_pt>7 && abs(BToMuMuPi_trg_mu_eta)<1.5) || (Muon_pt[BToMuMuPi_sel_mu_idx]>7 && abs( Muon_pt[BToMuMuPi_sel_mu_idx])<1.5))', '', 'trg_cond', 'Either muon triggers', '==', '', 0, 10)
  trg_cond = Quantity('((BToMuMuPi_trg_mu_pt>7 && abs(BToMuMuPi_trg_mu_eta)<1.5) || (Muon_pt[BToMuMuPi_sel_mu_idx]>7 && abs(Muon_eta[BToMuMuPi_sel_mu_idx])<1.5))', '', 'trg_cond', 'Either muon triggers', '==', '', 0, 10)
 
  b_mass = Quantity('BToMuMuPi_mass', 'b_mass', 'b_mass', '#mu#mu#pi invariant mass', '<', 'GeV', 5, 10)
  b_mass_m = Quantity('BToMuMuPi_mass', 'b_mass', 'b_mass', '#mu#mu#pi invariant mass', '>', 'GeV', 0, 3)
  b_pt = Quantity('BToMuMuPi_pt', 'b_pt', 'b_pt', '#mu#mu#pi pT', '>', 'GeV', 5, 15)
  b_eta = Quantity('abs(BToMuMuPi_eta)', 'abs(b_eta)', 'b_eta', '#mu#mu#pi |#eta|', '<', '', 0, 2.5) 

  sv_prob = Quantity('BToMuMuPi_sv_prob', 'sv_prob', 'sv_prob', 'SV probability', '>', '', 0, 0.05)
  sv_chi2 = Quantity('BToMuMuPi_sv_chi2', 'sv_chi2', 'sv_chi2', 'SV #chi^{2}', '<', '', 4, 14)  
  hnl_cos2D = Quantity('BToMuMuPi_hnl_cos2D', 'hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.93, 1) 
  hnl_charge = Quantity('', 'hnl_charge', 'hnl_charge', '\mu\pi charge', '==', '', 0, 1)
  sv_lxy_sig = Quantity('BToMuMuPi_sv_lxy_sig', 'sv_lxysig', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 100)
  sv_lxyz = Quantity('', 'sv_lxyz', 'sv_lxyz', '3D displacement of the SV', '<', '', 0, 100)
  sv_err = Quantity('', 'sv_lxy/sv_lxysig', 'sv_err', 'error on the SV 2D displacement', '>', '', 0, 100)

  trg_mu_pt = Quantity('BToMuMuPi_trg_mu_pt', 'trgmu_pt', 'trgmu_pt', 'primary muon pT', '>', 'GeV', 5, 13)
  trg_mu_eta = Quantity('fabs(BToMuMuPi_trg_mu_eta)', 'abs(trgmu_eta)', 'trgmu_eta', 'primary muon |#eta|', '<', '', 0, 2.5)
  #trg_mu_ip3d = Quantity('BToMuMuPi_trg_mu_ip3d', 'trgmu_ip3d', 'trigger muon IP (3D)', '<', 'cm', 5, 10)
  #trg_mu_sip3d = Quantity('BToMuMuPi_trg_mu_sip3d', 'trg_mu_sip3d', 'trigger muon IP significance (3D)', '>', '', 0, 3)
  trg_mu_dz = Quantity('abs(BToMuMuPi_trg_mu_dz)', 'abs(trgmu_dz)', 'abs(trg_mu_dz)', 'trigger muon |IP| on z', '>', 'cm', 0, 0.01)
  trg_mu_dxy = Quantity('abs(BToMuMuPi_trg_mu_dxy)', 'abs(trgmu_dxy)', 'abs(trg_mu_dxy)', 'trigger muon |IP| on xy', '>', 'cm', 0, 0.01)
  trg_mu_softid = Quantity('Muon_softId[BToMuMuPi_trg_mu_idx]', 'trgmu_softid', 'trg_mu_softid', 'trigger muon softId', '==', '', 0, 1)
  trg_mu_pfiso03rel = Quantity('Muon_pfiso03Rel_all[BToMuMuPi_trg_mu_idx]', 'trgmu_pfiso03rel', 'trg_mu_pfiso03rel', 'trigger muon relative PF iso03', '<', '', 0, 20)
  mu0_softid = Quantity('', 'mu0_softid', 'mu0_softid', 'primary muon soft ID', '==', '', 0, 1)
  
  mu_pt = Quantity('Muon_pt[BToMuMuPi_sel_mu_idx]', 'mu_pt', 'mu_pt', 'displaced muon pT', '>', 'GeV', 0, 10)
  mu_eta = Quantity('abs(Muon_eta[BToMuMuPi_sel_mu_idx])', 'abs(mu_eta)', 'mu_eta', 'displaced muon |#eta|', '<', '', 0, 2.5)
  mu_ip3d = Quantity('BToMuMuPi_sel_mu_ip3d', 'mu_ip3d', 'mu_ip3d', 'displaced muon IP (3D)', '>', 'cm', 0, 0.5)
  mu_sip3d = Quantity('BToMuMuPi_sel_mu_sip3d', 'mu_ip3dsig', 'mu_sip3d', 'displaced muon IP significance (3D)', '>', '', 0, 15)
  mu_dz = Quantity('abs(Muon_dz[BToMuMuPi_sel_mu_idx])', 'abs(mu_dz)', 'mu_dz', 'displaced muon |IP| on z', '>', 'cm', 0, 0.015)
  mu_dzS = Quantity('abs(Muon_dzS[BToMuMuPi_sel_mu_idx])', 'abs(mu_dzsig)', 'mu_dzS', 'displaced muon |IP| significance on z', '>', '', 0, 5)
  mu_dxy = Quantity('abs(Muon_dxy[BToMuMuPi_sel_mu_idx])', 'abs(mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 1)
  mu_dxy_BS = Quantity('abs(Muon_dxy_BS[BToMuMuPi_sel_mu_idx])', 'abs(mu_dxy)', 'mu_dxy_bs', 'displaced muon |IP| on xy', '>', 'cm', 0, 1)
  mu_dxyS = Quantity('abs(Muon_dxyS[BToMuMuPi_sel_mu_idx])', 'abs(mu_dxysig)', 'mu_dxyS', 'displaced muon |IP| significance on xy', '>', '', 0, 2)
  mu_dxyS_BS = Quantity('abs(Muon_dxyS_BS[BToMuMuPi_sel_mu_idx])', 'abs(mu_dxysig)', 'mu_dxyS', 'displaced muon |IP| significance on xy', '>', '', 0, 2)
  mu_pfiso03rel = Quantity('Muon_pfiso03Rel_all[BToMuMuPi_sel_mu_idx]', 'mu_pfiso03rel', 'mu_pfiso03rel', 'displaced muon relative PF iso03', '<', '', 0, 20)
  mu_looseid = Quantity('', 'mu_looseid', 'mu_looseid', 'displaced muon loose ID', '==', '', 0, 1)

  pi_pt = Quantity('BToMuMuPi_pi_pt', 'pi_pt', 'pi_pt', 'displaced pion pT', '>', 'GeV', 0, 6)
  pi_eta = Quantity('abs(BToMuMuPi_pi_eta)', 'abs(pi_eta)', 'pi_eta', 'displaced pion |#eta|', '<', '', 0, 2.5)
  pi_dz = Quantity('abs(BToMuMuPi_pi_dz)', 'abs(pi_dz)', 'pi_dz', 'displaced pion |IP| on z', '>', 'cm', 0, 0.015)
  pi_dxy = Quantity('abs(BToMuMuPi_pi_dxy)', 'abs(pi_dxy)', 'pi_dxy', 'displaced pion |IP| on xy', '>', 'cm', 0, 0.05)
  pi_dzS = Quantity('abs(BToMuMuPi_pi_dzS)', 'abs(pi_dzsig)', 'pi_dzS', 'displaced pion |IP| significance on z', '>', '', 0, 10)
  pi_dxyS = Quantity('abs(BToMuMuPi_pi_dxyS)', 'abs(pi_dxysig)', 'pi_dxyS', 'displaced pion |IP| significance on xy', '>', '', 0, 10)
  pi_DCASig = Quantity('BToMuMuPi_pi_DCASig', 'abs(pi_dcasig)', 'pi_dcasig', 'displaced pion DCA significance', '>', '', 0, 20)
  pi_packedcandhashighpurity = Quantity('', 'pi_packedcandhashighpurity', 'pi_packedcandhashighpurity', 'displaced pi has high purity track', '==', '', 0, 1)

  mu0mu_veto = Quantity('', '((mu0_charge!=mu_charge && abs(mu0_mu_mass-3.097)>0.15 && abs(mu0_mu_mass-3.686)>0.08 && abs(mu0_mu_mass-1.019)>0.01) || (mu0_charge==mu_charge))', 'mu0mu_veto', 'Resonances veto in \mu_0\mu mass spectrum', '==', '', 0, 1)
  mu0pi_veto = Quantity('', '((mu0_charge==mu_charge && abs(mu0_pi_mass-3.097)>0.05 && abs(mu0_pi_mass-1.76)>0.05) || (mu0_charge!=mu_charge))', 'mupi_veto', 'Resonances veto in \mu_0\pi mass spectrum', '==', '', 0, 1)

  #mu_Lxyz         = Quantity('BToMuMuPi_dimu_Lxyz'       , 'muons_Lxyz'     , '>', 'cm', 0, 0.02)
  #mu_trgmu_vzdiff = Quantity('BToMuMuPi_dimu_vzdiff'  , 'dimu_vzdiff', '#Delta vz(trgmu, mu)', '>', 'cm', 0, 0.02)
  #pi_trgmu_vzdiff = Quantity('BToMuMuPi_pi_mu_vzdiff'  , 'pi_trgmu_vzdiff', '#Delta vz(trgmu, pi)', '>', 'cm', 0, 0.02)

  #dr_mu_pi = Quantity('BToMuMuPi_dr_mu_pi', 'dr_mu_pi', '#DeltaR(mu, pi)', '<', '', 0, 2)
  #dr_trgmu_hnl = Quantity('BToMuMuPi_dr_trgmu_hnl', 'dr_trgmu_hnl', '#DeltaR(trgmu, hnl)', '<', '', 0, 2)

  hnl_mass_m = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass', 'hnl_mass', '#mu#pi invariant mass', '<', 'GeV', 4, 10) 
  hnl_mass_p = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass', 'hnl_mass', '#mu#pi invariant mass', '>', 'GeV', 0, 3) 
  hnl_pt = Quantity('BToMuMuPi_hnl_pt', 'hnl_pt', 'hnl_pt', '#mu#pi pT', '>', 'GeV', 2, 8) 
  hnl_eta = Quantity('abs(BToMuMuPi_hnl_eta)', 'abs(hnl_eta)', 'hnl_eta', '#mu#pi |#eta|', '<', '', 0, 2.5) 

  cos_theta_star_pion = Quantity('fabs(BToMuMuPi_cos_theta_star_pion)', 'fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.5, 1) 
  cos_theta_star_lepton = Quantity('fabs(BToMuMuPi_cos_theta_star_lepton)', 'fabs(cos_theta_star_lepton)', 'cos_theta_star_lepton', '|cos(#theta_{#mu}*)|', '<', '', 0.5, 1) 

  trgmu_pi_mass = Quantity('BToMuMuPi_trgmu_pi_mass', 'trgmu_pi_mass', 'trgmu_pi_mass', 'trg#mu#pi invariant mass', '<', 'GeV', 3, 8) 
  trgmu_mu_mass = Quantity('BToMuMuPi_trgmu_mu_mass', 'trgmu_mu_mass', 'trgmu_mu_mass', 'trg#mu#mu invariant mass', '<', 'GeV', 3, 8) 

  dpt_pi_fit_pi  = Quantity('abs(BToMuMuPi_dpt_pi_fit_pi)', 'abs(deltapt_pi_fit_pi)', 'dpt_pi_fit_pi', '|#Delta pT(#pi, fit#pi)|', '<', 'GeV', 0, 0.3) 
  deta_pi_fit_pi  = Quantity('abs(BToMuMuPi_deta_pi_fit_pi)', 'abs(deltaeta_pi_fit_pi)', 'deta_pi_fit_pi', '|#Delta #eta(#pi, fit#pi)|', '<', '', 0, 0.05) 
  dphi_pi_fit_pi  = Quantity('abs(BToMuMuPi_dphi_pi_fit_pi)', 'abs(deltaphi_pi_fit_pi)', 'dphi_pi_fit_pi', '|#Delta #Phi(#pi, fit#pi)|', '<', '', 0, 0.1) 
  
  dpt_mu_fit_mu  = Quantity('abs(BToMuMuPi_dpt_mu_fit_mu)', 'abs(deltapt_mu_fit_mu)', 'dpt_mu_fit_mu', '|#Delta pT(#mu, fit#mu)|', '<', 'GeV', 0, 0.8) 
  deta_mu_fit_mu  = Quantity('abs(BToMuMuPi_deta_mu_fit_mu)', 'abs(deltaeta_mu_fit_mu)', 'deta_mu_fit_mu', '|#Delta #eta(#mu, fit#mu)|', '<', '', 0, 0.1) 
  dphi_mu_fit_mu  = Quantity('abs(BToMuMuPi_dphi_mu_fit_mu)', 'abs(deltaphi_mu_fit_mu)', 'dphi_mu_fit_mu', '|#Delta #Phi(#mu, fit#mu)|', '<', '', 0, 0.1) 

  dr_trgmu_pi  = Quantity('abs(BToMuMuPi_dr_trgmu_pi)', 'abs(deltar_trgmu_pi)', 'dr_trgmu_pi', '|#DeltaR(trg#mu, #pi)|', '<', '', 0.3, 1.) 
  dr_trgmu_mu  = Quantity('abs(BToMuMuPi_dr_trgmu_mu)', 'abs(deltar_trgmu_mu)', 'dr_trgmu_mu', '|#DeltaR(trg#mu, #mu)|', '<', '', 0.3, 1.) 
  deta_mu_pi  = Quantity('abs(BToMuMuPi_deta_mu_pi)', 'abs(deltaeta_mu_pi)', 'deta_mu_pi', '|#Delta #eta(#mu, #pi)|', '<', '', 0, 1.8) 
  deta_trgmu_pi  = Quantity('abs(BToMuMuPi_deta_trgmu_pi)', 'abs(deltaeta_trgmu_pi)', 'deta_trgmu_pi', '|#Delta #eta(trg#mu, #pi)|', '<', '', 0, 1.8) 
  dphi_mu_pi  = Quantity('abs(BToMuMuPi_dphi_mu_pi)', 'abs(deltaphi_mu_pi)', 'dphi_mu_pi', '|#Delta #phi(#mu, #pi)|', '<', '', 0, 3) 
  dphi_trgmu_pi  = Quantity('abs(BToMuMuPi_dphi_trgmu_pi)', 'abs(deltaphi_trgmu_pi)', 'dphi_trgmu_pi', '|#Delta #phi(trg#mu, #pi)|', '<', '', 0, 3) 

  #'mu0_softid == 1',
  #'mu_looseid == 1',
  #'pi_packedcandhashighpurity == 1',
  #'((mu0_charge!=mu_charge && abs(mu0_mu_mass-3.097)>0.15 && abs(mu0_mu_mass-3.686)>0.08 && abs(mu0_mu_mass-1.019)>0.01) || (mu0_charge==mu_charge))',
  #'((mu0_charge==mu_charge && abs(mu0_pi_mass-3.097)>0.05 && abs(mu0_pi_mass-1.76)>0.05) || (mu0_charge!=mu_charge))',
  #'sv_lxyz < 100',
  #'sv_lxy/sv_lxysig>0',

  #soft_muon = Quantity('b_sel_mu_isSoft', 'mu_isSoft', '==', '', 0, 1)
 
  printCutflow = False
  printEfficiency = True
  printScan = False

  doPreselection = True
  doBaselineSelection = False
  doLooseSelection = False

  if doPreselection:

    preselection = [] 

    ##cut_trg_cond = 1
    ##if printScan: Selection(files, trg_cond, npoints=30, write_cut_analysis=False, proposed_cut=cut_trg_cond).getScanGraph()
    ###Selection(files, trg_cond, npoints=5, write_cut_analysis=False, proposed_cut=cut_trg_cond).getScanGraph()
    ###Selection(files, trg_cond, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, trg_cond, proposed_cut=cut_trg_cond).printCutflowLine()
    ###Selection(files, trg_cond, proposed_cut=cut_trg_cond).printCutflowLine()
    ###preselection.append(PreselectedQuantity(trg_cond, cut_trg_cond))

    #preselection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))

    #cut_trg_mu_pt = 1.5 #7
    #if printScan: Selection(files, trg_mu_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
    ##Selection(files, trg_mu_pt, npoints=5, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
    ##Selection(files, trg_mu_pt, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
    #if printEfficiency: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    #Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    #preselection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    #preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))



    cut_trg_mu_eta = 2 #1.5
    if printScan: Selection(files, trg_mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).getScanGraph()
    #Selection(files, trg_mu_eta, npoints=30).getROCGraph()
    if printCutflow: Selection(files, trg_mu_eta, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).printCutflowLine()
    if printEfficiency: Selection(files, trg_mu_eta, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).printEfficiencyLine()
    preselection.append(PreselectedQuantity(trg_mu_eta, cut_trg_mu_eta))

    ## for the moment we don't add it
    ##cut_trg_mu_softid = 1 
    ##if printScan: Selection(files, trg_mu_softid, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).getScanGraph()
    ##Selection(files, trg_mu_softid, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).getScanGraph()
    ##if printCutflow: Selection(files, trg_mu_softid, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).printCutflowLine()
    ##Selection(files, trg_mu_softid, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).printCutflowLine()
    ##preselection.append(PreselectedQuantity(trg_mu_softid, cut_trg_mu_softid))

    # depends on the mass of the signal, flat roc
    ##cut_trg_mu_pfiso03rel = 4 
    ##if printScan: Selection(files, trg_mu_pfiso03rel, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).getScanGraph()
    ###Selection(files, trg_mu_pfiso03rel, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).getScanGraph()
    ###Selection(files, trg_mu_pfiso03rel, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, trg_mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).printCutflowLine()
    ###Selection(files, trg_mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).printCutflowLine()
    ##preselection.append(PreselectedQuantity(trg_mu_pfiso03rel, cut_trg_mu_pfiso03rel))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    #preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    
    cut_pi_pt = 0.7 # this might be a bit tight
    if printScan: Selection(files, pi_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_pt).getScanGraph()
    #Selection(files, pi_pt, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_pt, preexisting_selection=preselection, proposed_cut=cut_pi_pt).printCutflowLine()
    #Selection(files, pi_pt, preexisting_selection=preselection, proposed_cut=cut_pi_pt).printCutflowLine()
    if printEfficiency: Selection(files, pi_pt, preexisting_selection=preselection, proposed_cut=cut_pi_pt).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_pt, cut_pi_pt))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    #preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    
    cut_pi_eta = 2
    if printScan: Selection(files, pi_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_eta).getScanGraph()
    #Selection(files, pi_eta, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_eta, preexisting_selection=preselection, proposed_cut=cut_pi_eta).printCutflowLine()
    #Selection(files, pi_eta, preexisting_selection=preselection, proposed_cut=cut_pi_eta).printCutflowLine()
    if printEfficiency: Selection(files, pi_eta, preexisting_selection=preselection, proposed_cut=cut_pi_eta).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_eta, cut_pi_eta))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    #preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    cut_pi_dz = 0.005
    if printScan: Selection(files, pi_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dz).getScanGraph()
    #Selection(files, pi_dz, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_dz, preexisting_selection=preselection, proposed_cut=cut_pi_dz).printCutflowLine()
    #Selection(files, pi_dz, preexisting_selection=preselection, proposed_cut=cut_pi_dz).printCutflowLine()
    if printEfficiency: Selection(files, pi_dz, preexisting_selection=preselection, proposed_cut=cut_pi_dz).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_dz, cut_pi_dz))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    #preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_pi_dxy = 0.005
    if printScan: Selection(files, pi_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).getScanGraph()
    #Selection(files, pi_dxy, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_dxy, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).printCutflowLine()
    #Selection(files, pi_dxy, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).printCutflowLine()
    if printEfficiency: Selection(files, pi_dxy, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_dxy, cut_pi_dxy))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    #preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_pi_dzS = 1.5
    if printScan: Selection(files, pi_dzS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).getScanGraph()
    #Selection(files, pi_dzS, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_dzS, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).printCutflowLine()
    #Selection(files, pi_dzS, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).printCutflowLine()
    if printEfficiency: Selection(files, pi_dzS, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_dzS, cut_pi_dzS))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    #preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_pi_dxyS = 3
    if printScan: Selection(files, pi_dxyS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).getScanGraph()
    #Selection(files, pi_dxyS, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_dxyS, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).printCutflowLine()
    #Selection(files, pi_dxyS, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).printCutflowLine()
    if printEfficiency: Selection(files, pi_dxyS, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_dxyS, cut_pi_dxyS))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    #preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    
    cut_pi_DCASig = 5
    if printScan: Selection(files, pi_DCASig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).getScanGraph()
    #Selection(files, pi_DCASig, npoints=30).getROCGraph()
    if printCutflow: Selection(files, pi_DCASig, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).printCutflowLine()
    #Selection(files, pi_DCASig, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).printCutflowLine()
    if printEfficiency: Selection(files, pi_DCASig, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).printEfficiencyLine()
    preselection.append(PreselectedQuantity(pi_DCASig, cut_pi_DCASig))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    #preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_mu_pt = 1.5
    if printScan: Selection(files, mu_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_pt).getScanGraph()
    #Selection(files, mu_pt, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_pt, preexisting_selection=preselection, proposed_cut=cut_mu_pt).printCutflowLine()
    #Selection(files, mu_pt, preexisting_selection=preselection, proposed_cut=cut_mu_pt).printCutflowLine()
    if printEfficiency: Selection(files, mu_pt, preexisting_selection=preselection, proposed_cut=cut_mu_pt).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_pt, cut_mu_pt))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    #preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_mu_eta = 2
    if printScan: Selection(files, mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_eta).getScanGraph()
    #Selection(files, mu_eta, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_eta, preexisting_selection=preselection, proposed_cut=cut_mu_eta).printCutflowLine()
    #Selection(files, mu_eta, preexisting_selection=preselection, proposed_cut=cut_mu_eta).printCutflowLine()
    if printEfficiency: Selection(files, mu_eta, preexisting_selection=preselection, proposed_cut=cut_mu_eta).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_eta, cut_mu_eta))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    #preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    cut_mu_dz = 0.0015
    if printScan: Selection(files, mu_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dz).getScanGraph()
    #Selection(files, mu_dz, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_dz, preexisting_selection=preselection, proposed_cut=cut_mu_dz).printCutflowLine()
    #Selection(files, mu_dz, preexisting_selection=preselection, proposed_cut=cut_mu_dz).printCutflowLine()
    if printEfficiency: Selection(files, mu_dz, preexisting_selection=preselection, proposed_cut=cut_mu_dz).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_dz, cut_mu_dz))

    # dxy computed wrt PV
    #cut_mu_dxy = 0.001
    #if printScan: Selection(files, mu_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).getScanGraph()
    ##Selection(files, mu_dxy, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, mu_dxy, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).printCutflowLine()
    ##Selection(files, mu_dxy, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_dxy, cut_mu_dxy))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    #preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    cut_mu_dxy_BS = 0.001
    if printScan: Selection(files, mu_dxy_BS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxy_BS).getScanGraph()
    #Selection(files, mu_dxy_BS, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_dxy_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxy_BS).printCutflowLine()
    #Selection(files, mu_dxy_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxy_BS).printCutflowLine()
    if printEfficiency: Selection(files, mu_dxy_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxy_BS).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_dxy_BS, cut_mu_dxy_BS))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    #preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    cut_mu_dzS = 1
    if printScan: Selection(files, mu_dzS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).getScanGraph()
    #Selection(files, mu_dzS, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_dzS, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).printCutflowLine()
    #Selection(files, mu_dzS, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).printCutflowLine()
    if printEfficiency: Selection(files, mu_dzS, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_dzS, cut_mu_dzS)) # not added for dsa muons

    # dxy computed wrt PV
    #cut_mu_dxyS = 1.5
    #if printScan: Selection(files, mu_dxyS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).getScanGraph()
    ##Selection(files, mu_dxyS, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, mu_dxyS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).printCutflowLine()
    ##Selection(files, mu_dxyS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_dxyS, cut_mu_dxyS)) # not added for dsa muons

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    #preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_mu_dxyS_BS = 1.5
    if printScan: Selection(files, mu_dxyS_BS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS_BS).getScanGraph()
    #Selection(files, mu_dxyS_BS, npoints=30).getROCGraph()
    if printCutflow: Selection(files, mu_dxyS_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS_BS).printCutflowLine()
    #Selection(files, mu_dxyS_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS_BS).printCutflowLine()
    if printEfficiency: Selection(files, mu_dxyS_BS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS_BS).printEfficiencyLine()
    preselection.append(PreselectedQuantity(mu_dxyS_BS, cut_mu_dxyS_BS)) # not added for dsa muons

    ## use the 2D significance instead
    ##cut_mu_sip3d = 7
    ##if printScan: Selection(files, mu_sip3d, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_sip3d).getScanGraph()
    ###Selection(files, mu_sip3d, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, mu_sip3d, preexisting_selection=preselection, proposed_cut=cut_mu_sip3d).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_sip3d, cut_mu_sip3d))

    ##cut_mu_pfiso03rel = 1.5
    ##if printScan: Selection(files, mu_pfiso03rel, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_pfiso03rel).getScanGraph()
    ###Selection(files, mu_pfiso03rel, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_pfiso03rel).getScanGraph()
    ###Selection(files, mu_pfiso03rel, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_mu_pfiso03rel).printCutflowLine()
    ##Selection(files, mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_mu_pfiso03rel).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_pfiso03rel, cut_mu_pfiso03rel))

    ## cuts applied post fit

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    #preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))

    cut_sv_prob = 0.001
    if printScan: Selection(files, sv_prob, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_prob).getScanGraph()
    #Selection(files, sv_prob, npoints=30).getROCGraph()
    if printCutflow: Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printCutflowLine()
    #Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printCutflowLine()
    if printEfficiency: Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printEfficiencyLine()
    preselection.append(PreselectedQuantity(sv_prob, cut_sv_prob))

    # not sufficiently discriminating
    ##cut_sv_chi2 = 9
    ##Selection(files, sv_chi2, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_chi2).getScanGraph()
    ##Selection(files, sv_chi2, npoints=30).getROCGraph()
    ##Selection(files, sv_chi2, preexisting_selection=preselection, proposed_cut=cut_sv_chi2).printCutflowLine()
    ##preselection.append(PreselectedQuantity(sv_chi2, cut_sv_chi2))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    #preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_hnl_cos2D = 0.995
    if printScan: Selection(files, hnl_cos2D, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).getScanGraph()
    #Selection(files, hnl_cos2D, npoints=30).getROCGraph()
    if printCutflow: Selection(files, hnl_cos2D, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).printCutflowLine()
    #Selection(files, hnl_cos2D, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).printCutflowLine()
    if printEfficiency: Selection(files, hnl_cos2D, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).printEfficiencyLine()
    preselection.append(PreselectedQuantity(hnl_cos2D, cut_hnl_cos2D))

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    #preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    preselection.append(PreselectedQuantity(b_mass, 8))


    cut_sv_lxy_sig = 15 
    if printScan: Selection(files, sv_lxy_sig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).getScanGraph()
    #Selection(files, b_mass, npoints=30).getROCGraph()
    if printCutflow: Selection(files, sv_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).printCutflowLine()
    #Selection(files, sv_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).printCutflowLine()
    if printEfficiency: Selection(files, sv_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).printEfficiencyLine()
    preselection.append(PreselectedQuantity(sv_lxy_sig, cut_sv_lxy_sig)) # not added for dsa muons

    preselection = []
    preselection.append(PreselectedQuantity(trg_mu_pt, 1.5))
    preselection.append(PreselectedQuantity(trg_mu_eta, 2))
    preselection.append(PreselectedQuantity(pi_pt, 0.7))
    preselection.append(PreselectedQuantity(pi_eta, 2))
    preselection.append(PreselectedQuantity(pi_dz, 0.005))
    preselection.append(PreselectedQuantity(pi_dxy, 0.005))
    preselection.append(PreselectedQuantity(pi_dzS, 1.5))
    preselection.append(PreselectedQuantity(pi_dxyS, 3))
    preselection.append(PreselectedQuantity(pi_DCASig, 5))
    preselection.append(PreselectedQuantity(mu_pt, 1.5))
    preselection.append(PreselectedQuantity(mu_eta, 2))
    preselection.append(PreselectedQuantity(mu_dz, 0.0015))
    preselection.append(PreselectedQuantity(mu_dxy_BS, 0.001))
    preselection.append(PreselectedQuantity(mu_dzS, 1)) # not added for dsa muons
    preselection.append(PreselectedQuantity(mu_dxyS_BS, 1.5)) # not added for dsa muons
    preselection.append(PreselectedQuantity(sv_prob, 0.01))
    preselection.append(PreselectedQuantity(hnl_cos2D, 0.995))
    preselection.append(PreselectedQuantity(sv_lxy_sig, 15)) # not added for dsa muons
    #preselection.append(PreselectedQuantity(b_mass, 8))

    cut_b_mass = 8 # move back to 8 if we want to have at hand this control region 
    if printScan: Selection(files, b_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_b_mass).getScanGraph()
    #Selection(files, b_mass, npoints=30).getROCGraph()
    if printCutflow: Selection(files, b_mass, preexisting_selection=preselection, proposed_cut=cut_b_mass).printCutflowLine()
    #Selection(files, b_mass, preexisting_selection=preselection, proposed_cut=cut_b_mass).printCutflowLine()
    if printEfficiency: Selection(files, b_mass, preexisting_selection=preselection, proposed_cut=cut_b_mass).printEfficiencyLine()
    preselection.append(PreselectedQuantity(b_mass, cut_b_mass))



    cut_trg_mu_pt = 1.5 #7
    if printScan: Selection(files, trg_mu_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
    #Selection(files, trg_mu_pt, npoints=5, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
    #Selection(files, trg_mu_pt, npoints=30).getROCGraph()
    if printCutflow: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
    if printEfficiency: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    Selection(files, trg_mu_pt, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    #preselection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))








    #cut_hnl_mass_m = 7 
    #if printScan: Selection(files, hnl_mass_m, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).getScanGraph()
    ##Selection(files, hnl_mass_m, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, hnl_mass_m, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).printCutflowLine()
    ##Selection(files, hnl_mass_m, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).printCutflowLine()
    #preselection.append(PreselectedQuantity(hnl_mass_m, cut_hnl_mass_m))

    #cut_cos_theta_star_pion = 0.9 
    #if printScan: Selection(files, cos_theta_star_pion, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_pion).getScanGraph()
    ##Selection(files, cos_theta_star_pion, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, cos_theta_star_pion, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_pion).printCutflowLine()
    ##Selection(files, cos_theta_star_pion, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_pion).printCutflowLine()
    #preselection.append(PreselectedQuantity(cos_theta_star_pion, cut_cos_theta_star_pion))

    ## to be added to the preselection?
    #cut_trgmu_pi_mass = 4.5 
    #if printScan: Selection(files, trgmu_pi_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trgmu_pi_mass).getScanGraph()
    ##Selection(files, trgmu_pi_mass, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, trgmu_pi_mass, preexisting_selection=preselection, proposed_cut=cut_trgmu_pi_mass).printCutflowLine()
    ##Selection(files, trgmu_pi_mass, preexisting_selection=preselection, proposed_cut=cut_trgmu_pi_mass).printCutflowLine()
    #preselection.append(PreselectedQuantity(trgmu_pi_mass, cut_trgmu_pi_mass))

    #cut_trgmu_mu_mass = 4.5 
    #if printScan: Selection(files, trgmu_mu_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trgmu_mu_mass).getScanGraph()
    ##Selection(files, trgmu_mu_mass, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, trgmu_mu_mass, preexisting_selection=preselection, proposed_cut=cut_trgmu_mu_mass).printCutflowLine()
    ##Selection(files, trgmu_mu_mass, preexisting_selection=preselection, proposed_cut=cut_trgmu_mu_mass).printCutflowLine()
    #preselection.append(PreselectedQuantity(trgmu_mu_mass, cut_trgmu_mu_mass))

    #cut_dr_trgmu_pi = 1.5 
    #if printScan: Selection(files, dr_trgmu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_pi).getScanGraph()
    ##Selection(files, dr_trgmu_pi, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, dr_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_pi).printCutflowLine()
    ##Selection(files, dr_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_pi).printCutflowLine()
    #preselection.append(PreselectedQuantity(dr_trgmu_pi, cut_dr_trgmu_pi))

    #cut_dr_trgmu_mu = 1. 
    #if printScan: Selection(files, dr_trgmu_mu, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_mu).getScanGraph()
    ##Selection(files, dr_trgmu_mu, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, dr_trgmu_mu, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_mu).printCutflowLine()
    ##Selection(files, dr_trgmu_mu, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_mu).printCutflowLine()
    #preselection.append(PreselectedQuantity(dr_trgmu_mu, cut_dr_trgmu_mu))

    ##cut_cos_theta_star_lepton = 0.9 
    ##if printScan: Selection(files, cos_theta_star_lepton, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_lepton).getScanGraph()
    ###Selection(files, cos_theta_star_lepton, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, cos_theta_star_lepton, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_lepton).printCutflowLine()
    ##Selection(files, cos_theta_star_lepton, preexisting_selection=preselection, proposed_cut=cut_cos_theta_star_lepton).printCutflowLine()
    ##preselection.append(PreselectedQuantity(cos_theta_star_lepton, cut_cos_theta_star_lepton))

    ## added for dsa
    ##cut_deta_mu_pi = 1 
    ##if printScan: Selection(files, deta_mu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_mu_pi).getScanGraph()
    ###Selection(files, deta_mu_pi, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, deta_mu_pi, preexisting_selection=preselection, proposed_cut=cut_deta_mu_pi).printCutflowLine()
    ###Selection(files, deta_mu_pi, preexisting_selection=preselection, proposed_cut=cut_deta_mu_pi).printCutflowLine()
    ##preselection.append(PreselectedQuantity(deta_mu_pi, cut_deta_mu_pi))

    ## added for dsa
    ##cut_deta_trgmu_pi = 0.8
    ##if printScan: Selection(files, deta_trgmu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_trgmu_pi).getScanGraph()
    ###Selection(files, deta_trgmu_pi, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, deta_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_deta_trgmu_pi).printCutflowLine()
    ###Selection(files, deta_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_deta_trgmu_pi).printCutflowLine()
    ##preselection.append(PreselectedQuantity(deta_trgmu_pi, cut_deta_trgmu_pi))

    ## added for dsa
    ##cut_dphi_mu_pi = 1 
    ##if printScan: Selection(files, dphi_mu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_pi).getScanGraph()
    ###Selection(files, dphi_mu_pi, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, dphi_mu_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_pi).printCutflowLine()
    ###Selection(files, dphi_mu_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_pi).printCutflowLine()
    ##preselection.append(PreselectedQuantity(dphi_mu_pi, cut_dphi_mu_pi))

    #cut_dphi_trgmu_pi = 0.7 
    #if printScan: Selection(files, dphi_trgmu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_trgmu_pi).getScanGraph()
    ##Selection(files, dphi_trgmu_pi, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, dphi_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_trgmu_pi).printCutflowLine()
    #Selection(files, dphi_trgmu_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_trgmu_pi).printCutflowLine()
    #preselection.append(PreselectedQuantity(dphi_trgmu_pi, cut_dphi_trgmu_pi))

    ##cut_hnl_mass_p = 0.3 
    ##if printScan: Selection(files, hnl_mass_p, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_p).getScanGraph()
    ###Selection(files, hnl_mass_p, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, hnl_mass_p, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_p).printCutflowLine()
    ##Selection(files, hnl_mass_p, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_p).printCutflowLine()
    ##preselection.append(PreselectedQuantity(hnl_mass_p, cut_hnl_mass_p))

    ##cut_dpt_pi_fit_pi = 0.1 
    ##if printScan: Selection(files, dpt_pi_fit_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dpt_pi_fit_pi).getScanGraph()
    ###Selection(files, dpt_pi_fit_pi, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dpt_pi_fit_pi).getScanGraph()
    ###Selection(files, dpt_pi_fit_pi, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, dpt_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dpt_pi_fit_pi).printCutflowLine()
    ##Selection(files, dpt_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dpt_pi_fit_pi).printCutflowLine()
    ##preselection.append(PreselectedQuantity(dpt_pi_fit_pi, cut_dpt_pi_fit_pi))

    #cut_deta_pi_fit_pi = 0.015 
    #if printScan: Selection(files, deta_pi_fit_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).getScanGraph()
    ##Selection(files, deta_pi_fit_pi, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).getScanGraph()
    ##Selection(files, deta_pi_fit_pi, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, deta_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).printCutflowLine()
    ##Selection(files, deta_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).printCutflowLine()
    #preselection.append(PreselectedQuantity(deta_pi_fit_pi, cut_deta_pi_fit_pi))

    #cut_dphi_pi_fit_pi = 0.03
    #if printScan: Selection(files, dphi_pi_fit_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).getScanGraph()
    ##Selection(files, dphi_pi_fit_pi, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).getScanGraph()
    ##Selection(files, dphi_pi_fit_pi, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, dphi_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).printCutflowLine()
    ##Selection(files, dphi_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).printCutflowLine()
    #preselection.append(PreselectedQuantity(dphi_pi_fit_pi, cut_dphi_pi_fit_pi))

    ##cut_dpt_mu_fit_mu = 0.1 
    ##if printScan: Selection(files, dpt_mu_fit_mu, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dpt_mu_fit_mu).getScanGraph()
    ###Selection(files, dpt_mu_fit_mu, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dpt_mu_fit_mu).getScanGraph()
    ###Selection(files, dpt_mu_fit_mu, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, dpt_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_dpt_mu_fit_mu).printCutflowLine()
    ##Selection(files, dpt_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_dpt_mu_fit_mu).printCutflowLine()
    ##preselection.append(PreselectedQuantity(dpt_mu_fit_mu, cut_dpt_mu_fit_mu))

    ##cut_deta_mu_fit_mu = 0.001
    ##if printScan: Selection(files, deta_mu_fit_mu, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_mu_fit_mu).getScanGraph()
    ###Selection(files, deta_mu_fit_mu, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_mu_fit_mu).getScanGraph()
    ###Selection(files, deta_mu_fit_mu, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, deta_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_deta_mu_fit_mu).printCutflowLine()
    ##Selection(files, deta_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_deta_mu_fit_mu).printCutflowLine()
    ##preselection.append(PreselectedQuantity(deta_mu_fit_mu, cut_deta_mu_fit_mu))

    ##cut_dphi_mu_fit_mu = 0.008 
    ##if printScan: Selection(files, dphi_mu_fit_mu, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_fit_mu).getScanGraph()
    ###Selection(files, dphi_mu_fit_mu, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_fit_mu).getScanGraph()
    ###Selection(files, dphi_mu_fit_mu, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, dphi_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_fit_mu).printCutflowLine()
    ##Selection(files, dphi_mu_fit_mu, preexisting_selection=preselection, proposed_cut=cut_dphi_mu_fit_mu).printCutflowLine()
    ##preselection.append(PreselectedQuantity(dphi_mu_fit_mu, cut_dphi_mu_fit_mu))

    # leave it to mva
    ##cut_dr_trgmu_hnl = 0.5
    #if printScan: Selection(files, dr_trgmu_hnl, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).getScanGraph()
    ##Selection(files, dr_trgmu_hnl, npoints=30).getROCGraph()
    #if printCutflow: Selection(files, dr_trgmu_hnl, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).printCutflowLine()
    ##Selection(files, dr_trgmu_hnl, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).printCutflowLine()
    #preselection.append(PreselectedQuantity(dr_trgmu_hnl, cut_dr_trgmu_hnl))
    
    ##cut_b_mass_m = 0.5
    ##if printScan: Selection(files, b_mass_m, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_b_mass_m).getScanGraph()
    ###Selection(files, b_mass, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, b_mass_m, preexisting_selection=preselection, proposed_cut=cut_b_mass_m).printCutflowLine()
    ##Selection(files, b_mass_m, preexisting_selection=preselection, proposed_cut=cut_b_mass_m).printCutflowLine()
    ##preselection.append(PreselectedQuantity(b_mass_m, cut_b_mass_m))


    # correlated
    ##cut_b_pt = 11
    ##if printScan: Selection(files, b_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_b_pt).getScanGraph()
    ###Selection(files, b_pt, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, b_pt, preexisting_selection=preselection, proposed_cut=cut_b_pt).printCutflowLine()
    ##preselection.append(PreselectedQuantity(b_pt, cut_b_pt))

    # correlated
    ##cut_b_eta = 1.7
    ##if printScan: Selection(files, b_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_b_eta).getScanGraph()
    ##Selection(files, b_eta, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, b_eta, preexisting_selection=preselection, proposed_cut=cut_b_eta).printCutflowLine()
    ##preselection.append(PreselectedQuantity(b_eta, cut_b_eta))

    ##cut_hnl_mass = 6.5 # signal dependent
    ##Selection(files, hnl_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_mass).getScanGraph()
    ##Selection(files, hnl_mass, npoints=30).getROCGraph()
    ##Selection(files, hnl_mass, preexisting_selection=preselection, proposed_cut=cut_hnl_mass).printCutflowLine()
    ##preselection.append(PreselectedQuantity(hnl_mass, cut_hnl_mass))

    # correlated
    ##cut_hnl_pt = 4
    ##if printScan: Selection(files, hnl_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_pt).getScanGraph()
    ###Selection(files, hnl_pt, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, hnl_pt, preexisting_selection=preselection, proposed_cut=cut_hnl_pt).printCutflowLine()
    ##preselection.append(PreselectedQuantity(hnl_pt, cut_hnl_pt))
    
    # correlated
    ##cut_hnl_eta = 1.8
    ##if printScan: Selection(files, hnl_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_eta).getScanGraph()
    ###Selection(files, hnl_eta, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, hnl_eta, preexisting_selection=preselection, proposed_cut=cut_hnl_eta).printCutflowLine()
    ##preselection.append(PreselectedQuantity(hnl_eta, cut_hnl_eta))

    # trigger muon not displaced
    ##cut_trg_mu_ip3d = 8
    ##Selection(files, trg_mu_ip3d, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_ip3d).getScanGraph()
    ##Selection(files, trg_mu_ip3d, npoints=30).getROCGraph()
    ##Selection(files, trg_mu_ip3d, preexisting_selection=preselection, proposed_cut=cut_trg_mu_ip3d).printCutflowLine()
    ##preselection.append(PreselectedQuantity(trg_mu_ip3d, cut_trg_mu_ip3d))

    ##cut_trg_mu_sip3d = 0.5
    ##Selection(files, trg_mu_sip3d, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_sip3d).getScanGraph()
    ##Selection(files, trg_mu_sip3d, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=3200).getROCGraph()
    ##preselection.append(PreselectedQuantity(trg_mu_sip3d, cut_trg_mu_sip3d))

    ##cut_trg_mu_dz = 0.0
    ##Selection(files, trg_mu_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_dz).getScanGraph()
    ##Selection(files, trg_mu_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=9).getROCGraph()
    ##preselection.append(PreselectedQuantity(trg_mu_dz, cut_trg_mu_dz))

    ##cut_trg_mu_dxy = 0.0
    ##Selection(files, trg_mu_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_dxy).getScanGraph()
    ##Selection(files, trg_mu_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0.15).getROCGraph()
    ##preselection.append(PreselectedQuantity(trg_mu_dxy, cut_trg_mu_dxy))

    ## no 3D if 2D
    ##cut_mu_ip3d = 0.01
    ##if printScan: Selection(files, mu_ip3d, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_ip3d).getScanGraph()
    ##Selection(files, mu_ip3d, npoints=30).getROCGraph()
    ##if printCutflow: Selection(files, mu_ip3d, preexisting_selection=preselection, proposed_cut=cut_mu_ip3d).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_ip3d, cut_mu_ip3d))
    ##cut_mu_Lxyz = 0.005
    ##Selection(files, mu_Lxyz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_Lxyz).getScanGraph()
    ##Selection(files, mu_Lxyz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_Lxyz).getROCGraph()
    ##preselection.append(PreselectedQuantity(mu_Lxyz, cut_mu_Lxyz))

    ##cut_pi_trgmu_vzdiff = 0.015
    ##Selection(files, pi_trgmu_vzdiff, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_trgmu_vzdiff).getScanGraph()
    ##Selection(files, pi_trgmu_vzdiff, npoints=30).getROCGraph()
    ##Selection(files, pi_trgmu_vzdiff, preexisting_selection=preselection, proposed_cut=cut_pi_trgmu_vzdiff).printCutflowLine()
    ##preselection.append(PreselectedQuantity(pi_trgmu_vzdiff, cut_pi_trgmu_vzdiff))

    ##cut_mu_trgmu_vzdiff = 0.005
    ##Selection(files, mu_trgmu_vzdiff, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_trgmu_vzdiff).getScanGraph()
    ##Selection(files, mu_trgmu_vzdiff, npoints=30).getROCGraph()
    ##Selection(files, mu_trgmu_vzdiff, preexisting_selection=preselection, proposed_cut=cut_mu_trgmu_vzdiff).printCutflowLine()
    ##preselection.append(PreselectedQuantity(mu_trgmu_vzdiff, cut_mu_trgmu_vzdiff))

    ##cut_dr_mu_pi = 1.8 # careful, signal dependent
    ##Selection(files, dr_mu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dr_mu_pi).getScanGraph()
    ##Selection(files, dr_mu_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0.8).getROCGraph()
    ##preselection.append(PreselectedQuantity(dr_mu_pi, cut_dr_mu_pi))

    #print Selection(files, b_mass, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0).getPreselectionString('bkg')



  if doBaselineSelection:

    baseline_selection = []

    cut_hnl_charge = 0
    if printCutflow: Selection(files, hnl_charge, proposed_cut=cut_hnl_charge).printCutflowLine()
    if printEfficiency: Selection(files, hnl_charge, proposed_cut=cut_hnl_charge).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(hnl_charge, cut_hnl_charge))

    cut_mu0_softid = 1
    if printCutflow: Selection(files, mu0_softid, preexisting_selection=baseline_selection, proposed_cut=cut_mu0_softid).printCutflowLine()
    if printEfficiency: Selection(files, mu0_softid, preexisting_selection=baseline_selection, proposed_cut=cut_mu0_softid).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(mu0_softid, cut_mu0_softid))

    cut_mu_looseid = 1
    if printCutflow: Selection(files, mu_looseid, preexisting_selection=baseline_selection, proposed_cut=cut_mu_looseid).printCutflowLine()
    if printEfficiency: Selection(files, mu_looseid, preexisting_selection=baseline_selection, proposed_cut=cut_mu_looseid).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(mu_looseid, cut_mu_looseid))

    cut_pi_packedcandhashighpurity = 1
    if printCutflow: Selection(files, pi_packedcandhashighpurity, preexisting_selection=baseline_selection, proposed_cut=cut_pi_packedcandhashighpurity).printCutflowLine()
    if printEfficiency: Selection(files, pi_packedcandhashighpurity, preexisting_selection=baseline_selection, proposed_cut=cut_pi_packedcandhashighpurity).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(pi_packedcandhashighpurity, cut_pi_packedcandhashighpurity))

    cut_mu0mu_veto = 1
    if printCutflow: Selection(files, mu0mu_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0mu_veto).printCutflowLine()
    if printEfficiency: Selection(files, mu0mu_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0mu_veto).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(mu0mu_veto, cut_mu0mu_veto))

    cut_mu0pi_veto = 1
    if printCutflow: Selection(files, mu0pi_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0pi_veto).printCutflowLine()
    if printEfficiency: Selection(files, mu0pi_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0pi_veto).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(mu0pi_veto, cut_mu0pi_veto))

    cut_sv_lxyz = 100
    if printCutflow: Selection(files, sv_lxyz, preexisting_selection=baseline_selection, proposed_cut=cut_sv_lxyz).printCutflowLine()
    if printEfficiency: Selection(files, sv_lxyz, preexisting_selection=baseline_selection, proposed_cut=cut_sv_lxyz).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(sv_lxyz, cut_sv_lxyz))

    cut_sv_err = 0
    if printCutflow: Selection(files, sv_err, preexisting_selection=baseline_selection, proposed_cut=cut_sv_err).printCutflowLine()
    if printEfficiency: Selection(files, sv_err, preexisting_selection=baseline_selection, proposed_cut=cut_sv_err).printEfficiencyLine()
    baseline_selection.append(PreselectedQuantity(sv_err, cut_sv_err))


  if doLooseSelection:

    loose_selection = []

    cut_trg_mu_pt = 0.3
    if printCutflow: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
    if printEfficiency: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
    loose_selection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))

    cut_trg_mu_eta = 2.8
    if printScan: Selection(files, trg_mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).getScanGraph()
    if printCutflow: Selection(files, trg_mu_eta, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).printCutflowLine()
    if printEfficiency: Selection(files, trg_mu_eta, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).printEfficiencyLine()
    loose_selection.append(PreselectedQuantity(trg_mu_eta, cut_trg_mu_eta))



