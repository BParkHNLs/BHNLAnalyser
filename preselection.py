import numpy as np
import ROOT
import os
import os.path
from os import path
from common import PlottingTools


def getLabel(str_): # move to utils
  new_str_ = str_
  if 'abs(' in str_: new_str_ = new_str_.replace('abs(', '')
  if ')'in str_: new_str_ = new_str_.replace(')', '')
  return new_str_


class Quantity(object): # make it inherit from utils Quantity and add logic
  def __init__(self, name_nano, label='', title='', logic='', units='', binMin=0., binMax=0.):
    self.name_nano = name_nano
    self.label = label
    self.title = title
    self.logic = logic
    self.units = units
    self.binMin = binMin
    self.binMax = binMax
                              

class Selection(PlottingTools):
  def __init__(self, files, quantity, preexisting_selection=None, npoints=50, write_cut_analysis=False, proposed_cut=None):
    # have separate files for signal and background ! for the moment keeping it like that
    self.files                 = files
    self.quantity              = quantity
    self.preexisting_selection = preexisting_selection
    self.npoints               = npoints
    self.write_cut_analysis    = write_cut_analysis
    self.proposed_cut          = proposed_cut


    # efficiency will be assessed based on the mass distributions
    # keep it like that?
    self.hnl_mass = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass','', binMin=0, binMax=10000) # with very loose selection some candidates have very large mupi invariant mass
    
    # baseline selection
    #self.baseline_selection = ' && '.join(['{}==0'.format('hnl_charge' if str_=='sig' else 'b_hnl_charge')])
    #                                  , 'b_mass<6.35']) 


  def createOutDir(self, outputdir): # this is to be moved to mother class
    if not path.exists(outputdir):
      os.system('mkdir {}'.format(outputdir))


  def getPreselectionString(self, str_):
    '''
    function to write the already-applied cuts in a string
    str_ is to make the distinction between nano and matched samples,
    where variables are named differently
    '''
    preselection_str = []
    for item, _ in enumerate(self.preexisting_selection):
      #name_variable = self.preexisting_selection[item].quantity.title if str_=='sig' else self.preexisting_selection[item].quantity.name_nano
      name_variable = self.preexisting_selection[item].quantity.name_nano
      preselection_str.append('{}{}{}'.format(name_variable,self.preexisting_selection[item].quantity.logic,self.preexisting_selection[item].chosen_cut))
    return ' && '.join(preselection_str)


  def computeLifetimeWeight(self, ctau0, ctau1, ct_reco): # to be moved to plotting tools
    return ctau0/ctau1 * math.exp(ct_reco*(1/ctau0 - 1/ctau1))


  def createHisto(self, file_, str_, with_extra_selection, cut=0): # move to common Plotter / PlottingTools class?
    '''
    str_ makes the difference between matched and nano sample
    '''

    filename = file_.sample_signal if str_=='sig' else file_.sample_background 
    tree_name = 'Events'
    #print filename

    #cut_variable = self.quantity.title if str_=='sig' else self.quantity.name_nano
    cut_variable = self.quantity.name_nano

    baseline_selection = 'BToMuMuPi_isMatched==1' if str_=='sig' else 'BToMuMuPi_isMatched>-99'
    
    c = ROOT.TCanvas()
    f = ROOT.TFile.Open(filename, 'READ')
    tree = f.Get(tree_name)
    
    # elements to fill the histogram
    distr_todraw = 'BToMuMuPi_hnl_mass' # self.hnl_mass.name_nano

    preselection = baseline_selection if self.preexisting_selection==None else baseline_selection + ' && ' + self.getPreselectionString(str_)
    if with_extra_selection:
      preselection += ' && {}{}{}'.format(cut_variable, self.quantity.logic, cut)

    #print preselection
    
    hist = ROOT.TH1D('hist', 'hist', 1500, self.hnl_mass.binMin, self.hnl_mass.binMax)
    #print 'begin draw'
    tree.Draw('{}>>hist'.format(distr_todraw), preselection)
    #print 'end draw'

    hist.SetDirectory(0)
    return hist


  #def getEff(self, file_, str_, cut):
  def getEff(self, entries_selected, entries_initial):
    #hist          = self.createHisto(file_, str_, False, cut)
    #hist_selected = self.createHisto(file_, str_, True, cut)
    eff = entries_selected / entries_initial if entries_initial!=0 else 0.
    return eff


  def getInfoBox(self, style, xmin, ymin, xmax, ymax, text, colour): # move to plottingtools as getTPaveText 
    box = ROOT.TPaveText(xmin, ymin, xmax, ymax, style)
    box.AddText(text)
    box.SetBorderSize(0)
    box.SetFillColor(ROOT.kWhite)
    box.SetTextColor(colour)
    box.SetTextSize(0.11)
    box.SetTextFont(42)
    box.SetTextAlign(11)
    return box

  def printCutFlowInformation(self):
    box2 = ROOT.TPaveText(0.05,0.7,0.4,0.98,"NDC")
    #box2.AddText('Cuts already applied: {}'.format('HNL charge==0' if self.preexisting_selection==None else 'HNL charge==0 && ' + self.getPreselectionString('sig')))
    #box2.AddText('Cuts already applied: all the previous ones')
    box2.AddText('Additional proposed cut: {}{}{}'.format(self.quantity.label, self.quantity.logic, self.proposed_cut))
    box2.GetListOfLines().Last().SetTextColor(ROOT.kRed)
    box2.SetBorderSize(0)
    box2.SetFillColor(ROOT.kWhite)
    box2.SetTextSize(0.11)
    box2.SetTextFont(42)
    box2.SetTextAlign(11)
    
    box21 = ROOT.TPaveText(0.02,0.3,0.08,0.4,"brNDC")
    box21.AddText('{}GeV'.format(self.files[0].signal_mass))
    #box21.AddText('{}'.format(self.files[0].v2))
    box21.SetBorderSize(0)
    box21.SetFillColor(ROOT.kWhite)
    #box21.SetTextColor(ROOT.kOrange+1)
    box21.SetTextSize(0.11)
    box21.SetTextFont(42)
    box21.SetTextAlign(11)

    box3 = ROOT.TPaveText(0.08,0.35,0.4,0.65,"brNDC")
    box3.AddText('N_matched ini: {}'.format(self.createHisto(self.files[0], 'sig', False).GetEntries()))
    box3.AddText('N_matched new: {}'.format(self.createHisto(self.files[0], 'sig', True, self.proposed_cut).GetEntries()))
    box3.SetBorderSize(0)
    box3.SetFillColor(ROOT.kWhite)
    box3.SetTextColor(ROOT.kOrange+1)
    box3.SetTextSize(0.11)
    box3.SetTextFont(42)
    box3.SetTextAlign(11)
    
    box4 = ROOT.TPaveText(0.35,0.35,0.5,0.65,"brNDC")
    box4.AddText('==> -{}%'.format(round((1 - self.createHisto(self.files[0], 'sig', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[0], 'sig', False).GetEntries())*100, 2)))
    print 'sig1',round((1 - self.createHisto(self.files[0], 'sig', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[0], 'sig', False).GetEntries())*100, 2)
    box4.SetBorderSize(0)
    box4.SetFillColor(ROOT.kWhite)
    box4.SetTextColor(ROOT.kOrange+1)
    box4.SetTextSize(0.11)
    box4.SetTextFont(42)
    box4.SetTextAlign(11)
    
    box5 = ROOT.TPaveText(0.08,0.0,0.4,0.3,"brNDC")
    box5.AddText('N_nano ini: {}'.format(self.createHisto(self.files[0], 'bkg', False).GetEntries()))
    box5.AddText('N_nano new: {}'.format(self.createHisto(self.files[0], 'bkg', True, self.proposed_cut).GetEntries()))
    box5.SetBorderSize(0)
    box5.SetFillColor(ROOT.kWhite)
    box5.SetTextColor(ROOT.kBlue+2)
    box5.SetTextSize(0.11)
    box5.SetTextFont(42)
    box5.SetTextAlign(11)
    
    box6 = ROOT.TPaveText(0.35,0.,0.5,0.3,"brNDC")
    box6.AddText('==> -{}%'.format(round((1 - self.createHisto(self.files[0], 'bkg', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[0], 'bkg', False).GetEntries())*100, 2)))
    print 'bkg1',round((1 - self.createHisto(self.files[0], 'bkg', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[0], 'bkg', False).GetEntries())*100, 2)
    box6.SetBorderSize(0)
    box6.SetFillColor(ROOT.kWhite)
    box6.SetTextColor(ROOT.kBlue+2)
    box6.SetTextSize(0.11)
    box6.SetTextFont(42)
    box6.SetTextAlign(11)
    
    box22 = ROOT.TPaveText(0.5,0.3,0.55,0.4,"brNDC")
    #box22.AddText('{}GeV'.format(self.files[1].signal_mass))
    #box22.AddText('{}'.format(self.files[1].v2))
    box22.SetBorderSize(0)
    box22.SetFillColor(ROOT.kWhite)
    #box21.SetTextColor(ROOT.kOrange+1)
    box22.SetTextSize(0.11)
    box22.SetTextFont(42)
    box22.SetTextAlign(11)

    box7 = ROOT.TPaveText(0.57,0.35,0.75,0.65,"brNDC")
    #box7.AddText('N_matched ini: {}'.format(self.createHisto(self.files[1], 'sig', False).GetEntries()))
    #box7.AddText('N_matched new: {}'.format(self.createHisto(self.files[1], 'sig', True, self.proposed_cut).GetEntries()))
    box7.SetBorderSize(0)
    box7.SetFillColor(ROOT.kWhite)
    box7.SetTextColor(ROOT.kOrange+1)
    box7.SetTextSize(0.11)
    box7.SetTextFont(42)
    box7.SetTextAlign(11)
    
    box8 = ROOT.TPaveText(0.85,0.35,0.95,0.65,"brNDC")
    #box8.AddText('==> -{}%'.format(round((1 - self.createHisto(self.files[1], 'sig', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[1], 'sig', False).GetEntries())*100, 2)))
    #print 'sig2',round((1 - self.createHisto(self.files[1], 'sig', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[1], 'sig', False).GetEntries())*100, 2)
    box8.SetBorderSize(0)
    box8.SetFillColor(ROOT.kWhite)
    box8.SetTextColor(ROOT.kOrange+1)
    box8.SetTextSize(0.11)
    box8.SetTextFont(42)
    box8.SetTextAlign(11)
  
    box9 = ROOT.TPaveText(0.57,0.0,0.75,0.3,"brNDC")
    #box9.AddText('N_nano ini: {}'.format(self.createHisto(self.files[1], 'bkg', False).GetEntries()))
    #box9.AddText('N_nano new: {}'.format(self.createHisto(self.files[1], 'bkg', True, self.proposed_cut).GetEntries()))
    box9.SetBorderSize(0)
    box9.SetFillColor(ROOT.kWhite)
    box9.SetTextColor(ROOT.kBlue+2)
    box9.SetTextSize(0.11)
    box9.SetTextFont(42)
    box9.SetTextAlign(11)
    
    box10 = ROOT.TPaveText(0.85,0.,0.95,0.3,"brNDC")
    #box10.AddText('==> -{}%'.format(round((1 - self.createHisto(self.files[1], 'bkg', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[1], 'bkg', False).GetEntries())*100, 2)))
    #print 'bkg2',round((1 - self.createHisto(self.files[1], 'bkg', True, self.proposed_cut).GetEntries() / self.createHisto(self.files[1], 'bkg', False).GetEntries())*100, 2)
    box10.SetBorderSize(0)
    box10.SetFillColor(ROOT.kWhite)
    box10.SetTextColor(ROOT.kBlue+2)
    box10.SetTextSize(0.11)
    box10.SetTextFont(42)
    box10.SetTextAlign(11)
  
    return [box2, box21, box3, box4, box5, box6, box22, box7, box8, box9, box10]


  def getScanGraph(self):

    '''
    plot the signal efficiency and background rejection as a function of the cut applied on the quantity of interest
    possibility to draw a line at the chosen cut
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints) 

    canv = ROOT.TCanvas('canv', 'canv', 900, 800)
    canv.SetGrid()
    
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
    pad_up.SetBottomMargin(0.1)
    pad_up.Draw()
    canv.cd()
    #pad_down = ROOT.TPad("pad_down","pad_down",0,0,0.5,0.25)
    pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
    pad_down.SetBottomMargin(0.15)
    pad_down.Draw()
    canv.cd()
    pad_leg = ROOT.TPad("pad_leg","pad_leg",0.5,0,1,0.25)
    pad_leg.SetBottomMargin(0.15)
    #pad_leg.Draw()

    if self.write_cut_analysis:
      #pad_up.cd()
      pad_down.cd()

    gs_sig = []
    gs_bkg = []

    # signal efficiency
    for ifile, file_ in enumerate(self.files): # do we want to keep that option?
      g_sig = ROOT.TGraph()

      initial_sig_entries = self.createHisto(file_, 'sig', False).GetEntries()

      for idx, cut in enumerate(points):
        selected_sig_entries = self.createHisto(file_, 'sig', True, cut).GetEntries()
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
    g_bkg = ROOT.TGraph()
    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', False).GetEntries()

    for idx, cut in enumerate(points):
      selected_bkg_entries = self.createHisto(self.files[0], 'bkg', True, cut).GetEntries()
      g_bkg.SetPoint(idx, cut, 1-self.getEff(selected_bkg_entries, initial_bkg_entries))
      g_bkg.SetLineWidth(0)
      g_bkg.SetMarkerColor(ROOT.kBlue+2)
      g_bkg.SetMarkerStyle(20)

    #gs_bkg.append(g_bkg)

    #for g_bkg in gs_bkg:
    g_bkg.Draw('P')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, 0, self.proposed_cut, 1)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      line.Draw('same')
        
    #legend = ROOT.TLegend(0.5, 0.37, 0.82, 0.61)
    #legend.SetTextSize(0.03);
    #legend.SetLineColor(0);
    #legend.SetFillColorAlpha(0, 0);
    #legend.SetBorderSize(0);
    legend = PlottingTools.getRootTLegend(self, xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    for ifile, file_ in enumerate(self.files):
        legend.AddEntry(gs_sig[ifile], 'sig efficiency ({}GeV, {}mm)'.format(file_.signal_mass,file_.signal_ctau))
    #legend.AddEntry(gs_bkg[ifile], 'bkg rejection')
    legend.AddEntry(g_bkg, 'bkg rejection')
    legend.Draw()

    # write cutflow information
    if self.write_cut_analysis:
      pad_down.cd()
      proposed_sig_entries = self.createHisto(self.files[0], 'sig', True, self.proposed_cut).GetEntries()
      proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, self.proposed_cut).GetEntries()

      box1 = self.getInfoBox("NDC", 0.05, 0.7, 0.4, 0.98, 'Additional proposed cut: {}{}{}'.format(self.quantity.label, self.quantity.logic, self.proposed_cut), ROOT.kRed)
      box2 = self.getInfoBox("brNDC", 0.02, 0.3, 0.08, 0.4, '{}GeV'.format(self.files[0].signal_mass), ROOT.kBlack)
      box3 = self.getInfoBox("brNDC", 0.08, 0.45, 0.4, 0.65, 'N_matched ini: {}'.format(initial_sig_entries), ROOT.kOrange+1)
      box4 = self.getInfoBox("brNDC", 0.08, 0.25, 0.4, 0.55, 'N_matched new: {}'.format(proposed_sig_entries), ROOT.kOrange+1)
      box5 = self.getInfoBox("brNDC", 0.35, 0.35, 0.5, 0.65, '==> -{}%'.format(round((1 - proposed_sig_entries / initial_sig_entries)*100, 2)), ROOT.kOrange+1)
      box6 = self.getInfoBox("brNDC", 0.08, 0.15, 0.4, 0.3, 'N_nano ini: {}'.format(initial_bkg_entries), ROOT.kBlue+2)
      box7 = self.getInfoBox("brNDC", 0.08, 0., 0.4, 0.15, 'N_nano new: {}'.format(proposed_bkg_entries), ROOT.kBlue+2)
      box8 = self.getInfoBox("brNDC", 0.35, 0., 0.5, 0.3, '==> -{}%'.format(round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2)), ROOT.kBlue+2)

      #[box1, box2, box3, box4, box5, box6, box7, box8, box9, box10, box11] = self.printCutFlowInformation()

      box1.Draw('same')
      box2.Draw('same')
      box3.Draw('same')
      box4.Draw('same')
      box5.Draw('same')
      box6.Draw('same')
      box7.Draw('same')
      box8.Draw('same')

    canv.cd()
    self.createOutDir('myPlots/selection')
    canv.SaveAs('myPlots/selection/scan_{}.png'.format(getLabel(self.quantity.label)))


  def getROCGraph(self):

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints)

    canv = ROOT.TCanvas('canv', 'canv', 1200, 1000)
    g = ROOT.TGraph2D()
    g.SetTitle('Cut on {} ({})'.format(self.quantity.title, self.quantity.logic))

    # do it for first file only
    file_ = self.files[0]
    
    initial_sig_entries = self.createHisto(file_, 'sig', False).GetEntries()
    initial_bkg_entries = self.createHisto(file_, 'bkg', False).GetEntries()

    for idx, cut in enumerate(points):
      selected_sig_entries = self.createHisto(file_, 'sig', True, cut).GetEntries()
      selected_bkg_entries = self.createHisto(file_, 'bkg', True, cut).GetEntries()
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
    

    self.createOutDir('myPlots/selection')
    canv.SaveAs('myPlots/selection/roc_{}.png'.format(getLabel(self.quantity.label)))


  def printCutflowLine(self):
    initial_sig_entries = self.createHisto(self.files[0], 'sig', False).GetEntries()
    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', False).GetEntries()
  
    proposed_sig_entries = self.createHisto(self.files[0], 'sig', True, self.proposed_cut).GetEntries()
    proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, self.proposed_cut).GetEntries()

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
      initial_sig1_entries = self.createHisto(self.files[1], 'sig', False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[1], 'sig', True, self.proposed_cut).GetEntries()

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
      initial_sig1_entries = self.createHisto(self.files[1], 'sig', False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[1], 'sig', True, self.proposed_cut).GetEntries()
      initial_sig2_entries = self.createHisto(self.files[2], 'sig', False).GetEntries()
      proposed_sig2_entries = self.createHisto(self.files[2], 'sig', True, self.proposed_cut).GetEntries()

      #print 'sig2 {}'.format(int(initial_sig1_entries))
      #print 'sig2 {}'.format(int(proposed_sig1_entries))
      #print 'sig3 {}'.format(int(initial_sig2_entries))
      #print 'sig3 {}'.format(int(proposed_sig2_entries))

      cutflow_line = '{qte} {log} {cut} {unit} & -{sig0_per}\% & -{sig1_per}\% & -{sig2_per}\% & -{bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          unit = self.quantity.units,
          sig0_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
          sig1_per = round((1 - proposed_sig1_entries / initial_sig1_entries)*100, 2), 
          sig2_per = round((1 - proposed_sig2_entries / initial_sig2_entries)*100, 2), 
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
      )

    cutflow_line = cutflow_line.replace('#eta', '$\eta$')
    cutflow_line = cutflow_line.replace('#mu#mu#pi', '$\mu\mu\pi$')
    cutflow_line = cutflow_line.replace('#mu#pi', '$\mu\pi$')
    cutflow_line = cutflow_line.replace('#Delta', '$\Delta$')
    cutflow_line = cutflow_line.replace('#pi', '$\pi$')
    cutflow_line = cutflow_line.replace('#Phi', '$\Phi$')
    print cutflow_line

    #print '{qte} {log} {cut} & -{sig_per}\% & -{bkg_per}\% \\ '.format(
    #    qte = self.quantity.title, 
    #    log = self.quantity.logic, 
    #    cut = self.proposed_cut, 
    #    #sig = int(proposed_sig_entries),
    #    sig_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
    #    #bkg = int(proposed_bkg_entries),
    #    bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
    #    )




class FileCollection(object): # move to utils? 
  '''
  this class allows to associate to the studied files the corresponding 
  mass and coupling points. Useful when scaning through e.g masses
  '''
  def __init__(self, sample_background='', sample_signal='', signal_mass='', signal_ctau=''):
    self.sample_background = sample_background
    self.sample_signal = sample_signal
    self.signal_mass = signal_mass
    #self.v2 = v2
    self.signal_ctau = signal_ctau
    # add label for reweighted?


class PreselectedQuantity(object): # keep here?
  '''
  this class will allow to impose a preselection when studying 
  the impact of the cut on a given parameter
  '''
  def __init__(self, quantity, chosen_cut):
    self.quantity   = quantity
    self.chosen_cut = chosen_cut


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  file_m3 = FileCollection(
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_withlooseselection.root',
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/Chunk0_n20/merged/bparknano_forselectionstudy.root',
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_forselectionstudy_full.root',
      sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH4_Run2018B/Chunk0_n10/bparknano_loosepreselection_v1_nj1.root',
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/Chunk0_n20/bparknano_forselectionstudy_nj1.root',
      #sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_withlooseselection_matched.root',
      #sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_alltrgmus_standardgenmatchingmatched.root',
      sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      signal_mass = 3,
      signal_ctau = 184,
      )

  file_m4p5 = FileCollection(
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/Chunk0_n20/bparknano_forselectionstudy_nj1.root',
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_forselectionstudy_full.root',
      #sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/bparknano_alltrgmus_standardgenmatchingmatched.root',
      sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      signal_mass = 4.5,
      signal_ctau = 1.2,
      )

  file_m1 = FileCollection(
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/Chunk0_n20/bparknano_forselectionstudy_nj1.root',
      #sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_forselectionstudy_full.root',
      #sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass1.0_ctau184.0/nanoFiles/merged/bparknano_alltrgmus_standardgenmatchingmatched.root',
      sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      signal_mass = 1,
      signal_ctau = 184,
      )

  file_V25 = FileCollection(
      sample_background = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/bparknano_forselectionstudy_full.root',
      sample_signal = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looseselection_muononly_alltrgmu.root',
      signal_mass = 3,
      signal_ctau = 2000,
      )
  
  files = []
  files.append(file_m3)
  files.append(file_m4p5)
  files.append(file_m1)
  #files.append(file_V25)

 
  b_mass = Quantity('BToMuMuPi_mass', 'b_mass', '#mu#mu#pi invariant mass', '<', 'GeV', 5, 10)
  b_mass_m = Quantity('BToMuMuPi_mass', 'b_mass', '#mu#mu#pi invariant mass', '>', 'GeV', 0, 3)
  b_pt = Quantity('BToMuMuPi_pt' , 'b_pt', '#mu#mu#pi pT', '>', 'GeV', 5, 15)
  b_eta = Quantity('abs(BToMuMuPi_eta)', 'b_eta', '#mu#mu#pi |#eta|', '<', '', 0, 2.5) 

  sv_prob = Quantity('BToMuMuPi_sv_prob', 'sv_prob', 'SV probability', '>', '', 0, 0.05)
  sv_chi2 = Quantity('BToMuMuPi_sv_chi2' , 'b_sv_chi2', 'SV #chi^{2}', '<', '', 4, 14)  
  hnl_cos2D = Quantity('BToMuMuPi_hnl_cos2D', 'hnl_cos2D', 'SV cos2D', '>', '', 0.93, 1) 
  sv_lxy_sig = Quantity('BToMuMuPi_sv_lxy_sig', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 100)

  trg_mu_pt = Quantity('BToMuMuPi_trg_mu_pt','trgmu_pt', 'trigger muon pT', '>', 'GeV', 5, 13)
  trg_mu_eta = Quantity('fabs(BToMuMuPi_trg_mu_eta)', 'trgmu_eta', 'trigger muon #eta', '<', '', 0, 2.5)
  #trg_mu_ip3d = Quantity('BToMuMuPi_trg_mu_ip3d', 'trgmu_ip3d', 'trigger muon IP (3D)', '>', 'cm', 0, 0.05)
  trg_mu_ip3d = Quantity('BToMuMuPi_trg_mu_ip3d', 'trgmu_ip3d', 'trigger muon IP (3D)', '<', 'cm', 5, 10)
  trg_mu_sip3d = Quantity('BToMuMuPi_trg_mu_sip3d', 'trg_mu_sip3d', 'trigger muon IP significance (3D)', '>', '', 0, 3)
  trg_mu_dz = Quantity('abs(BToMuMuPi_trg_mu_dz)', 'abs(trg_mu_dz)', 'trigger muon |IP| on z', '>', 'cm', 0, 0.01)
  trg_mu_dxy = Quantity('abs(BToMuMuPi_trg_mu_dxy)', 'abs(trg_mu_dxy)', 'trigger muon |IP| on xy', '>', 'cm', 0, 0.01)
  trg_mu_softid = Quantity('Muon_softId[BToMuMuPi_trg_mu_idx]', 'trg_mu_softid', 'trigger muon softId', '==', '', 0, 1)
  trg_mu_pfiso03rel = Quantity('Muon_pfiso03Rel_all[BToMuMuPi_trg_mu_idx]', 'trg_mu_pfiso03rel', 'trigger muon relative PF iso03', '<', '', 0, 20)
  
  #mu_pt = Quantity('BToMuMuPi_fit_mu_pt', 'mu_pt', 'displaced muon pT', '>', 'GeV', 0, 10)
  mu_pt = Quantity('Muon_pt[BToMuMuPi_sel_mu_idx]', 'mu_pt', 'displaced muon pT', '>', 'GeV', 0, 10)
  #mu_eta = Quantity('abs(BToMuMuPi_fit_mu_eta)', 'mu_eta', 'displaced muon |#eta|', '<', '', 0, 2.5)
  mu_eta = Quantity('abs(Muon_eta[BToMuMuPi_sel_mu_idx])', 'mu_eta', 'displaced muon |#eta|', '<', '', 0, 2.5)
  mu_ip3d = Quantity('BToMuMuPi_sel_mu_ip3d', 'mu_ip3d', 'displaced muon IP (3D)', '>', 'cm', 0, 0.5)
  mu_sip3d = Quantity('BToMuMuPi_sel_mu_sip3d', 'mu_sip3d', 'displaced muon IP significance (3D)', '>', '', 0, 15)
  mu_dz = Quantity('abs(BToMuMuPi_sel_mu_dz)', 'mu_dz', 'displaced muon |IP| on z', '>', 'cm', 0, 0.015)
  mu_dzS = Quantity('abs(Muon_dzS[BToMuMuPi_sel_mu_idx])', 'mu_dzS', 'displaced muon |IP| significance on z', '>', '', 0, 50)
  mu_dxy = Quantity('abs(BToMuMuPi_sel_mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 0.01)
  mu_dxyS = Quantity('abs(Muon_dxyS[BToMuMuPi_sel_mu_idx])', 'mu_dxyS', 'displaced muon |IP| significance on xy', '>', '', 0, 50)
  mu_pfiso03rel = Quantity('Muon_pfiso03Rel_all[BToMuMuPi_sel_mu_idx]', 'mu_pfiso03rel', 'displaced muon relative PF iso03', '<', '', 0, 20)

  #pi_pt = Quantity('BToMuMuPi_fit_pi_pt', 'pi_pt', 'displaced pion pT', '>', 0, 6)
  pi_pt = Quantity('ProbeTracks_pt[BToMuMuPi_pi_idx]', 'pi_pt', 'displaced pion pT', '>', 'GeV', 0, 6)
  #pi_eta = Quantity('abs(BToMuMuPi_fit_pi_eta)', 'pi_eta', 'displaced pion |#eta|', '<', '', 0, 2.5)
  pi_eta = Quantity('abs(ProbeTracks_eta[BToMuMuPi_pi_idx])', 'pi_eta', 'displaced pion |#eta|', '<', '', 0, 2.5)
  pi_dz = Quantity('abs(BToMuMuPi_pi_dz)', 'pi_dz', 'displaced pion |IP| on z', '>', 'cm', 0, 0.015)
  pi_dxy = Quantity('abs(BToMuMuPi_pi_dxy)', 'pi_dxy', 'displaced pion |IP| on xy', '>', 'cm', 0, 0.05)
  pi_dzS = Quantity('abs(BToMuMuPi_pi_dzS)', 'pi_dzS', 'displaced pion |IP| significance on z', '>', '', 0, 10)
  pi_dxyS = Quantity('abs(BToMuMuPi_pi_dxyS)', 'pi_dxyS', 'displaced pion |IP| significance on xy', '>', '', 0, 10)
  pi_DCASig = Quantity('BToMuMuPi_pi_DCASig', 'pi_dcasig', 'displaced pion DCA significance', '>', '', 0, 20)

  #mu_Lxyz        = Quantity('b_mu_Lxyz'       , 'muons_Lxyz'     , '<', 'cm', 0, 1.5)
  mu_Lxyz         = Quantity('BToMuMuPi_dimu_Lxyz'       , 'muons_Lxyz'     , '>', 'cm', 0, 0.02)
  #pi_trgmu_vzdiff = Quantity('b_pi_mu_vzdiff'  , 'pi_trgmu_vzdiff', '<', '', 0, 2)
  mu_trgmu_vzdiff = Quantity('BToMuMuPi_dimu_vzdiff'  , 'dimu_vzdiff', '#Delta vz(trgmu, mu)', '>', 'cm', 0, 0.02)
  pi_trgmu_vzdiff = Quantity('BToMuMuPi_pi_mu_vzdiff'  , 'pi_trgmu_vzdiff', '#Delta vz(trgmu, pi)', '>', 'cm', 0, 0.02)

  dr_mu_pi = Quantity('BToMuMuPi_dr_mu_pi', 'dr_mu_pi', '#DeltaR(mu, pi)', '<', '', 0, 2)
  dr_trgmu_hnl = Quantity('BToMuMuPi_dr_trgmu_hnl', 'dr_trgmu_hnl', '#DeltaR(trgmu, hnl)', '<', '', 0, 2)

  hnl_mass_m  = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass', '#mu#pi invariant mass', '<', 'GeV', 4, 10) 
  hnl_mass_p  = Quantity('BToMuMuPi_hnl_mass', 'hnl_mass', '#mu#pi invariant mass', '>', 'GeV', 0, 3) 
  hnl_pt    = Quantity('BToMuMuPi_hnl_pt', 'hnl_pt', '#mu#pi pT', '>', 'GeV', 2, 8) 
  hnl_eta   = Quantity('abs(BToMuMuPi_hnl_eta)', 'hnl_eta', '#mu#pi |#eta|', '<', '', 0, 2.5) 

  dpt_pi_fit_pi  = Quantity('abs(BToMuMuPi_dpt_pi_fit_pi)', 'dpt_pi_fit_pi', '|#Delta pT(#pi, fit#pi)|', '<', 'GeV', 0, 0.3) 
  deta_pi_fit_pi  = Quantity('abs(BToMuMuPi_deta_pi_fit_pi)', 'deta_pi_fit_pi', '|#Delta #eta(#pi, fit#pi)|', '<', '', 0, 0.05) 
  dphi_pi_fit_pi  = Quantity('abs(BToMuMuPi_dphi_pi_fit_pi)', 'dphi_pi_fit_pi', '|#Delta #Phi(#pi, fit#pi)|', '<', '', 0, 0.1) 
  
  dpt_mu_fit_mu  = Quantity('abs(BToMuMuPi_dpt_mu_fit_mu)', 'dpt_mu_fit_mu', '|#Delta pT(#mu, fit#mu)|', '<', 'GeV', 0, 0.8) 
  deta_mu_fit_mu  = Quantity('abs(BToMuMuPi_deta_mu_fit_mu)', 'deta_mu_fit_mu', '|#Delta #eta(#mu, fit#mu)|', '<', '', 0, 0.1) 
  dphi_mu_fit_mu  = Quantity('abs(BToMuMuPi_dphi_mu_fit_mu)', 'dphi_mu_fit_mu', '|#Delta #Phi(#mu, fit#mu)|', '<', '', 0, 0.1) 

  soft_muon = Quantity('b_sel_mu_isSoft', 'mu_isSoft', '==', '', 0, 1)
 
  print 'warning - some variables may have been renamed - check in case of crash'

  printCutflow = False
  printScan = True

  preselection = [] 

  cut_trg_mu_pt = 7
  if printScan: Selection(files, trg_mu_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
  #Selection(files, trg_mu_pt, npoints=5, write_cut_analysis=False, proposed_cut=cut_trg_mu_pt).getScanGraph()
  #Selection(files, trg_mu_pt, npoints=30).getROCGraph()
  if printCutflow: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
  #Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
  preselection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))

  cut_trg_mu_eta = 1.5 
  if printScan: Selection(files, trg_mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).getScanGraph()
  #Selection(files, trg_mu_eta, npoints=30).getROCGraph()
  if printCutflow: Selection(files, trg_mu_eta, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).printCutflowLine()
  #Selection(files, trg_mu_eta, preexisting_selection=preselection, proposed_cut=cut_trg_mu_eta).printCutflowLine()
  preselection.append(PreselectedQuantity(trg_mu_eta, cut_trg_mu_eta))

  # for the moment we don't add it
  cut_trg_mu_softid = 1 
  if printScan: Selection(files, trg_mu_softid, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).getScanGraph()
  #Selection(files, trg_mu_softid, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).getScanGraph()
  #if printCutflow: Selection(files, trg_mu_softid, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).printCutflowLine()
  #Selection(files, trg_mu_softid, preexisting_selection=preselection, proposed_cut=cut_trg_mu_softid).printCutflowLine()
  #preselection.append(PreselectedQuantity(trg_mu_softid, cut_trg_mu_softid))

  # depends on the mass of the signal, flat roc
  ##cut_trg_mu_pfiso03rel = 4 
  ##if printScan: Selection(files, trg_mu_pfiso03rel, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).getScanGraph()
  ###Selection(files, trg_mu_pfiso03rel, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).getScanGraph()
  ###Selection(files, trg_mu_pfiso03rel, npoints=30).getROCGraph()
  ##if printCutflow: Selection(files, trg_mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).printCutflowLine()
  ###Selection(files, trg_mu_pfiso03rel, preexisting_selection=preselection, proposed_cut=cut_trg_mu_pfiso03rel).printCutflowLine()
  ##preselection.append(PreselectedQuantity(trg_mu_pfiso03rel, cut_trg_mu_pfiso03rel))
  
  cut_pi_pt = 0.7 # this might be a bit tight
  if printScan: Selection(files, pi_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_pt).getScanGraph()
  #Selection(files, pi_pt, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_pt, preexisting_selection=preselection, proposed_cut=cut_pi_pt).printCutflowLine()
  #Selection(files, pi_pt, preexisting_selection=preselection, proposed_cut=cut_pi_pt).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_pt, cut_pi_pt))
  
  cut_pi_eta = 2
  if printScan: Selection(files, pi_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_eta).getScanGraph()
  #Selection(files, pi_eta, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_eta, preexisting_selection=preselection, proposed_cut=cut_pi_eta).printCutflowLine()
  #Selection(files, pi_eta, preexisting_selection=preselection, proposed_cut=cut_pi_eta).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_eta, cut_pi_eta))

  cut_pi_dz = 0.005
  if printScan: Selection(files, pi_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dz).getScanGraph()
  #Selection(files, pi_dz, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_dz, preexisting_selection=preselection, proposed_cut=cut_pi_dz).printCutflowLine()
  #Selection(files, pi_dz, preexisting_selection=preselection, proposed_cut=cut_pi_dz).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_dz, cut_pi_dz))

  cut_pi_dxy = 0.005
  if printScan: Selection(files, pi_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).getScanGraph()
  #Selection(files, pi_dxy, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_dxy, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).printCutflowLine()
  #Selection(files, pi_dxy, preexisting_selection=preselection, proposed_cut=cut_pi_dxy).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_dxy, cut_pi_dxy))

  cut_pi_dzS = 1.5
  if printScan: Selection(files, pi_dzS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).getScanGraph()
  #Selection(files, pi_dzS, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_dzS, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).printCutflowLine()
  #Selection(files, pi_dzS, preexisting_selection=preselection, proposed_cut=cut_pi_dzS).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_dzS, cut_pi_dzS))

  cut_pi_dxyS = 3
  if printScan: Selection(files, pi_dxyS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).getScanGraph()
  #Selection(files, pi_dxyS, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_dxyS, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).printCutflowLine()
  #Selection(files, pi_dxyS, preexisting_selection=preselection, proposed_cut=cut_pi_dxyS).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_dxyS, cut_pi_dxyS))
  
  cut_pi_DCASig = 5
  if printScan: Selection(files, pi_DCASig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).getScanGraph()
  #Selection(files, pi_DCASig, npoints=30).getROCGraph()
  if printCutflow: Selection(files, pi_DCASig, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).printCutflowLine()
  #Selection(files, pi_DCASig, preexisting_selection=preselection, proposed_cut=cut_pi_DCASig).printCutflowLine()
  preselection.append(PreselectedQuantity(pi_DCASig, cut_pi_DCASig))

  cut_mu_pt = 1.5
  if printScan: Selection(files, mu_pt, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_pt).getScanGraph()
  #Selection(files, mu_pt, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_pt, preexisting_selection=preselection, proposed_cut=cut_mu_pt).printCutflowLine()
  #Selection(files, mu_pt, preexisting_selection=preselection, proposed_cut=cut_mu_pt).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_pt, cut_mu_pt))

  cut_mu_eta = 2
  if printScan: Selection(files, mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_eta).getScanGraph()
  #Selection(files, mu_eta, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_eta, preexisting_selection=preselection, proposed_cut=cut_mu_eta).printCutflowLine()
  #Selection(files, mu_eta, preexisting_selection=preselection, proposed_cut=cut_mu_eta).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_eta, cut_mu_eta))

  cut_mu_dz = 0.0015
  if printScan: Selection(files, mu_dz, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dz).getScanGraph()
  #Selection(files, mu_dz, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_dz, preexisting_selection=preselection, proposed_cut=cut_mu_dz).printCutflowLine()
  #Selection(files, mu_dz, preexisting_selection=preselection, proposed_cut=cut_mu_dz).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_dz, cut_mu_dz))

  cut_mu_dxy = 0.001
  if printScan: Selection(files, mu_dxy, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).getScanGraph()
  #Selection(files, mu_dxy, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_dxy, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).printCutflowLine()
  #Selection(files, mu_dxy, preexisting_selection=preselection, proposed_cut=cut_mu_dxy).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_dxy, cut_mu_dxy))

  cut_mu_dzS = 1
  if printScan: Selection(files, mu_dzS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).getScanGraph()
  #Selection(files, mu_dzS, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_dzS, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).printCutflowLine()
  #Selection(files, mu_dzS, preexisting_selection=preselection, proposed_cut=cut_mu_dzS).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_dzS, cut_mu_dzS))

  cut_mu_dxyS = 1.5
  if printScan: Selection(files, mu_dxyS, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).getScanGraph()
  #Selection(files, mu_dxyS, npoints=30).getROCGraph()
  if printCutflow: Selection(files, mu_dxyS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).printCutflowLine()
  #Selection(files, mu_dxyS, preexisting_selection=preselection, proposed_cut=cut_mu_dxyS).printCutflowLine()
  preselection.append(PreselectedQuantity(mu_dxyS, cut_mu_dxyS))

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

  cut_sv_prob = 0.001
  if printScan: Selection(files, sv_prob, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_prob).getScanGraph()
  #Selection(files, sv_prob, npoints=30).getROCGraph()
  if printCutflow: Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printCutflowLine()
  #Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printCutflowLine()
  preselection.append(PreselectedQuantity(sv_prob, cut_sv_prob))

  # not sufficiently discriminating
  ##cut_sv_chi2 = 9
  ##Selection(files, sv_chi2, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_chi2).getScanGraph()
  ##Selection(files, sv_chi2, npoints=30).getROCGraph()
  ##Selection(files, sv_chi2, preexisting_selection=preselection, proposed_cut=cut_sv_chi2).printCutflowLine()
  ##preselection.append(PreselectedQuantity(sv_chi2, cut_sv_chi2))

  cut_hnl_cos2D = 0.99
  if printScan: Selection(files, hnl_cos2D, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).getScanGraph()
  #Selection(files, hnl_cos2D, npoints=30).getROCGraph()
  if printCutflow: Selection(files, hnl_cos2D, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).printCutflowLine()
  #Selection(files, hnl_cos2D, preexisting_selection=preselection, proposed_cut=cut_hnl_cos2D).printCutflowLine()
  preselection.append(PreselectedQuantity(hnl_cos2D, cut_hnl_cos2D))

  cut_sv_lxy_sig = 20 
  if printScan: Selection(files, sv_lxy_sig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).getScanGraph()
  #Selection(files, sv_lxy_sig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).getScanGraph()
  #Selection(files, b_mass, npoints=30).getROCGraph()
  if printCutflow: Selection(files, sv_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).printCutflowLine()
  #Selection(files, sv_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_sv_lxy_sig).printCutflowLine()
  preselection.append(PreselectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

  cut_b_mass = 5.4 # 8 # move back to 8 if we want to have at hand this control region 
  if printScan: Selection(files, b_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_b_mass).getScanGraph()
  #Selection(files, b_mass, npoints=30).getROCGraph()
  if printCutflow: Selection(files, b_mass, preexisting_selection=preselection, proposed_cut=cut_b_mass).printCutflowLine()
  #Selection(files, b_mass, preexisting_selection=preselection, proposed_cut=cut_b_mass).printCutflowLine()
  preselection.append(PreselectedQuantity(b_mass, cut_b_mass))

  ##cut_hnl_mass_m = 6.28 
  ##if printScan: Selection(files, hnl_mass_m, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).getScanGraph()
  ###Selection(files, hnl_mass_m, npoints=30).getROCGraph()
  ##if printCutflow: Selection(files, hnl_mass_m, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).printCutflowLine()
  ##Selection(files, hnl_mass_m, preexisting_selection=preselection, proposed_cut=cut_hnl_mass_m).printCutflowLine()
  ##preselection.append(PreselectedQuantity(hnl_mass_m, cut_hnl_mass_m))

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

  cut_deta_pi_fit_pi = 0.015 
  if printScan: Selection(files, deta_pi_fit_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).getScanGraph()
  #Selection(files, deta_pi_fit_pi, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).getScanGraph()
  #Selection(files, deta_pi_fit_pi, npoints=30).getROCGraph()
  if printCutflow: Selection(files, deta_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).printCutflowLine()
  #Selection(files, deta_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_deta_pi_fit_pi).printCutflowLine()
  preselection.append(PreselectedQuantity(deta_pi_fit_pi, cut_deta_pi_fit_pi))

  cut_dphi_pi_fit_pi = 0.03
  if printScan: Selection(files, dphi_pi_fit_pi, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).getScanGraph()
  #Selection(files, dphi_pi_fit_pi, npoints=10, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).getScanGraph()
  #Selection(files, dphi_pi_fit_pi, npoints=30).getROCGraph()
  if printCutflow: Selection(files, dphi_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).printCutflowLine()
  #Selection(files, dphi_pi_fit_pi, preexisting_selection=preselection, proposed_cut=cut_dphi_pi_fit_pi).printCutflowLine()
  preselection.append(PreselectedQuantity(dphi_pi_fit_pi, cut_dphi_pi_fit_pi))

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

  cut_dr_trgmu_hnl = 0.5
  if printScan: Selection(files, dr_trgmu_hnl, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).getScanGraph()
  #Selection(files, dr_trgmu_hnl, npoints=30).getROCGraph()
  if printCutflow: Selection(files, dr_trgmu_hnl, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).printCutflowLine()
  #Selection(files, dr_trgmu_hnl, preexisting_selection=preselection, proposed_cut=cut_dr_trgmu_hnl).printCutflowLine()
  preselection.append(PreselectedQuantity(dr_trgmu_hnl, cut_dr_trgmu_hnl))
  
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

  #print Selection(files, soft_muon, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0).getPreselectionString('sig')



