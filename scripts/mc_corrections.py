import os
import os.path
from os import path
import ROOT
from ROOT import gROOT, gStyle
from ROOT import RooFit
import math
import numpy as np
from array import array

from tools import Tools

import sys
sys.path.append('../objects')
from quantity import Quantity
from samples import data_samples, signal_samples
from baseline_selection import selection


class Plotter(Tools):
  def __init__(self, quantity='', data_files='', signal_files='', selection='', save_eos=False):
    self.tools = Tools()
    self.quantity = quantity
    self.data_files = data_files
    self.signal_files = signal_files
    self.selection = selection
    self.save_eos = save_eos


  def getQCDMCLabel(self, label1, label2, label3):
    '''
      based on MC label of the format ptXXtoYY 
    '''
    return label1[0:label1.find('to')] + label2[label2.find('to'):len(label2)] + label3[label3.find(' '):]


  def getDataLabel(self, label, version_label):
    idx = 0
    data_label = ''
    while idx < len(label):
      period = label[label.find('Parking', idx)+10:label.find('Parking', idx)+11]
      run = label[label.find('Run2018', idx)+7:label.find('Run2018', idx)+8]
      dataset = run+period
      if idx == 0: data_label = dataset # or idx>len(label)-20:
      else: data_label += '+{}'.format(dataset)
      idx = idx + 21
    data_label = data_label.replace('A1+A2+A3+A4+A5+A6', 'A*')
    data_label = data_label.replace('B1+B2+B3+B4+B5+B6', 'B*')
    data_label = data_label.replace('C1+C2+C3+C4+C5', 'C*')
    data_label = data_label.replace('D1+D2+D3+D4+D5', 'D*')
    if len(version_label) != 0:
      data_label += ' {}'.format(version_label)
    return data_label


  def getMaxRangeY(self, hist1, hist2, do_log=False, use_sig=False):
    margin = 0.15 if do_log==False else 0.5
    if use_sig: 
      max_hist1 = [hist.GetMaximum() for hist in hist1]
      the_max_hist1 = max(max_hist1)
      max_range = the_max_hist1+margin*the_max_hist1 if the_max_hist1>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    else:
      max_range = hist1.GetMaximum()+margin*hist1.GetMaximum() if hist1.GetMaximum()>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    return max_range


  def createHisto(self, tree, quantity, hist_name='hist', branchname='flat', weight=-99, selection='', scale=None, smear=None):
    ROOT.TH1.SetDefaultSumw2()
    hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)

    if scale != None and smear != None and 'signal' in hist_name:
      qte = '({})*({})*({})'.format(scale, smear, quantity.name_flat)
    elif scale != None and smear == None and 'signal' in hist_name:
      qte = '({})*({})'.format(scale, quantity.name_flat)
    else:
      qte = quantity.name_flat
      #if 'signal' not in hist_name:
      #  qte = quantity.name_flat
      #else:
      #  qte = 'double_gauss_corr_2*{}'.format(quantity.name_flat)
      if 'signal' not in hist_name and '_uncorrected' in qte: qte = qte.replace('_uncorrected', '')
      if 'signal' not in hist_name and '_corrected_weight' in qte: qte = qte.replace('_corrected_weight', '')
      if 'signal' not in hist_name and '_corrected_linearscale' in qte: qte = qte.replace('_corrected_linearscale', '')
    #print qte

    if quantity.name_flat != 'probe_dxy_bs_uncorrected':
      selection += ' && probe_dxy_bs > 0.015' #TODO apply correction once it's derived
    else:
      selection += ' && {} > 0.015'.format(qte)
    if selection == '' and weight == -99: selection_string = selection
    elif selection != '' and weight == -99: selection_string = '({})'.format(selection)
    elif selection == ''  and weight != -99: selection_string = '({})'.format(weight)
    else: selection_string = '({sel}) * ({wght})'.format(sel=selection, wght=weight)
    #print selection_string

    tree.Project(hist_name, qte, selection_string)

    hist.SetDirectory(0)
    return hist


  def computeWeight(self):
    '''
      Function that creates weights as the data/MC ratio
    '''
    # get the mc distribution
    signal_hists = [] # in principle not needed, since only one mc sample, but keep it generic
    for signal_file in self.signal_files:
      # get the tree
      f_sig = self.tools.getRootFile(signal_file.filename)
      tree_sig = self.tools.getTree(f_sig, 'tree')
    
      # define selection and weight
      matching_selection = 'ismatched==1'
      selection_signal = matching_selection if self.selection == '' else matching_selection + ' && ' + self.selection
      ## keep the mc non-weighted
      weight_sig = '(1)'

      # define the histogram
      hist_signal_name = 'hist_signal_{}'.format(self.quantity.label)
      hist_signal = self.createHisto(tree_sig, self.quantity, hist_name=hist_signal_name, branchname='flat', selection=selection_signal, weight=weight_sig)
      hist_signal.Sumw2()
    
      # overflow
      overflow_signal_tot = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
      error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
      hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
      hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
      hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
      hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)

      # normalise to unity
      int_signal = hist_signal.Integral()
      if int_signal != 0: hist_signal.Scale(1/int_signal)
      #hist_signal.SetFillColor(ROOT.kOrange+7)
      #hist_signal.SetFillStyle(3005)
      signal_hists.append(hist_signal)

    # get the data distribution (with sweights on)
    hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
    hist_data_tot.Sumw2()

    int_data_tot = 0.
    overflow_data_tot = 0
    error_overflow_data_tot = 0.

    for idata, data_file in enumerate(self.data_files): # again, not needed here, but kept for generality
      # get the tree
      f_data = self.tools.getRootFile(data_file.filename)
      tree_data = self.tools.getTree(f_data, 'tree')

      # get the histogram
      hist_data_name = 'hist_data_{}'.format(self.quantity.label)
      ## apply sweight on data
      hist_data = self.createHisto(tree_data, self.quantity, hist_name=hist_data_name, branchname='flat', selection=self.selection, weight='nbkg_sw')
      hist_data.Sumw2()

      # overflow
      overflow_data_tot = (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
      error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
      hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
      hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
      hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
      hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)
      
      # normalise to unity
      int_data_tot += hist_data.Integral()
      hist_data_tot.Add(hist_data)

      if hist_data.Integral() != 0: 
        hist_data.Scale(1./hist_data.Integral())

    if int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

    #hist_data_tot.SetMarkerStyle(20)

    # get the ratio
    hist_ratio = self.tools.getRatioHistogram(hist_data_tot, signal_hists[0])
    hist_ratio.Sumw2()

    # overflow (not needed)
    overflow_ratio = (hist_ratio.GetBinContent(hist_ratio.GetNbinsX()) + hist_ratio.GetBinContent(hist_ratio.GetNbinsX()+1))
    error_overflow_ratio = math.sqrt(math.pow(hist_ratio.GetBinError(hist_ratio.GetNbinsX()), 2) + math.pow(hist_ratio.GetBinError(hist_ratio.GetNbinsX()+1), 2)) 
    hist_ratio.SetBinContent(hist_ratio.GetNbinsX(), overflow_ratio)
    hist_ratio.SetBinError(hist_ratio.GetNbinsX(), error_overflow_ratio)
    hist_ratio.SetBinContent(hist_ratio.GetNbinsX()+1, 0)
    hist_ratio.SetBinError(hist_ratio.GetNbinsX()+1, 0)

    hist_ratio.SetLineWidth(2)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetTitle('')
    hist_ratio.GetXaxis().SetTitle(self.quantity.title)

    hist_ratio.GetYaxis().SetTitle('Data/MC')
    ROOT.gStyle.SetOptStat(0)

    # store the weight in a rootfile
    root_file = ROOT.TFile('mc_weight_{}.root'.format(self.quantity.label), 'RECREATE')
    #hist_data_tot.Write()
    #signal_hists[0].Write()
    hist_ratio.Write()
    root_file.Close()

    print '-> mc_weight_{}.root created'.format(self.quantity.label)


  def createScaleDistribution(self):
    '''
      Function that, in each bin, fetches the scale that makes the
      mc agree with data 

      NB: please do not mind the sketchy syntax...
    '''

    bins = np.arange(1., 61., 1)
    bin_min = 0.

    canv = self.tools.createTCanvas(name='canv', dimx=1200, dimy=1000)
    canv_graph = self.tools.createTCanvas(name='canv_graph', dimx=1200, dimy=1000)
    hist_scale = ROOT.TH1D('hist_scale', 'hist_scale', 30, 0.9, 1.5)

    hist_qte = ROOT.TH1D('hist_qte', 'hist_qte', 60, 0, 60)

    f_sig = self.tools.getRootFile(self.signal_files[0].filename)
    tree_sig = self.tools.getTree(f_sig, 'tree')

    f_data = self.tools.getRootFile(self.data_files[0].filename)
    tree_data = self.tools.getTree(f_data, 'tree')

    tree = ROOT.TTree('tree', 'tree')
    the_scale = array('d', [0])
    tree.Branch('scale', the_scale, 'scale/D')

    graph_scale = ROOT.TGraph()
    
    for bin_max in bins:
      print '\n{} {}'.format(bin_min, bin_max)
      quantity = Quantity(name_flat='probe_dxy_sig_bs_uncorrected', label='probe_dxy_sig_bs', title='probe #mu |d_{xy}| significance (BS)', nbins=1, bin_min=bin_min, bin_max=bin_max)
      bin_center = (bin_min + bin_max) / 2.
      bin_min = bin_max
      #if bin_min == 5.: continue
      scale = 0.99
      ratio = -1
      while (ratio < 0.985 or ratio > 1.015):
        # get the mc distribution
        signal_hists = [] # in principle not needed, since only one mc sample, but keep it generic
        for signal_file in self.signal_files:
          # get the tree
          #f_sig = self.tools.getRootFile(signal_file.filename)
          #tree_sig = self.tools.getTree(f_sig, 'tree')
        
          # define selection and weight
          matching_selection = 'ismatched==1'
          selection_signal = matching_selection if self.selection == '' else matching_selection + ' && ' + self.selection
          ## keep the mc non-weighted
          weight_sig = '(1)'

          # define the histogram
          hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(quantity.label, bin_min, bin_max, scale)
          hist_signal = self.createHisto(tree_sig, quantity, hist_name=hist_signal_name, branchname='flat', selection=selection_signal, weight=weight_sig, scale=scale)
          hist_signal.Sumw2()
        
          # overflow
          if bin_max == 60.:
            overflow_signal_tot = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
            error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
            hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
            hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
            hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
            hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)

          # normalise to unity
          int_signal = hist_signal.Integral()
          if int_signal != 0: hist_signal.Scale(1/int_signal)
          #hist_signal.SetFillColor(ROOT.kOrange+7)
          #hist_signal.SetFillStyle(3005)
          signal_hists.append(hist_signal)

        # get the data distribution (with sweights on)
        hist_data_tot_name = 'hist_data_tot_{}_{}_{}_{}'.format(quantity.label, bin_min, bin_max, scale)
        hist_data_tot = ROOT.TH1D(hist_data_tot_name, hist_data_tot_name, quantity.nbins, quantity.bin_min, quantity.bin_max)
        hist_data_tot.Sumw2()

        int_data_tot = 0.
        overflow_data_tot = 0
        error_overflow_data_tot = 0.

        for idata, data_file in enumerate(self.data_files): # again, not needed here, but kept for generality
          # get the tree
          #f_data = self.tools.getRootFile(data_file.filename)
          #tree_data = self.tools.getTree(f_data, 'tree')

          # get the histogram
          hist_data_name = 'hist_data_{}_{}_{}_{}'.format(quantity.label, bin_min, bin_max, scale)
          ## apply sweight on data
          hist_data = self.createHisto(tree_data, quantity, hist_name=hist_data_name, branchname='flat', selection=self.selection, weight='nbkg_sw')
          hist_data.Sumw2()

          # overflow
          if bin_max == 60.:
            overflow_data_tot = (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
            error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
            hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
            hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
            hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
            hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)
          
          # normalise to unity
          int_data_tot += hist_data.Integral()
          hist_data_tot.Add(hist_data)

          if hist_data.Integral() != 0: 
            hist_data.Scale(1./hist_data.Integral())

        if int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

        # get the ratio
        #hist_ratio = self.tools.getRatioHistogram(hist_data_tot, signal_hists[0])
        #hist_ratio.Sumw2()
        
        if signal_hists[0].GetBinContent(0) != 0:
          ratio = hist_data_tot.GetBinContent(0) / signal_hists[0].GetBinContent(0)
        else: ratio = 0
        print 'scale: {}'.format(scale)
        print 'ratio: {}'.format(ratio)
        scale = scale + 0.01
        if ratio == 0: break
        if scale > 1.25: break
        if ratio > 0.985 and ratio < 1.015: 
          print 'fill histo'
          hist_scale.Fill(scale)
          the_scale[0] = scale
          tree.Fill()
          point = graph_scale.GetN()
          graph_scale.SetPoint(point, bin_center, scale)
          graph_scale.GetXaxis().SetTitle('probe dxy significance')
          graph_scale.GetYaxis().SetTitle('scale')
          graph_scale.SetMarkerStyle(20)
          graph_scale.SetMarkerSize(2)
          hist_qte.SetBinContent(hist_qte.GetXaxis().FindBin(bin_center), scale)

    canv.cd()
    hist_scale.Draw()
    canv.SaveAs('myPlots/mc_corrections/scales.png')

    canv_graph.cd()
    graph_scale.Fit('pol1')
    ROOT.gStyle.SetOptFit() 
    graph_scale.Draw('AP')
    canv_graph.SaveAs('myPlots/mc_corrections/graph_scales.png')

    # store the weight in a rootfile
    root_file = ROOT.TFile('scales.root', 'RECREATE')
    #hist_scale.Write()
    tree.Write()
    root_file.Close()
    print '-> scales.root created'

    root_file_scale = ROOT.TFile('scale_graph.root', 'RECREATE')
    hist_qte.Write()
    root_file_scale.Close()
    print '-> scale_graph.root created'


  def performScaleFit(self):

    filename = 'scales_step1.root' 
    f = ROOT.TFile.Open(filename)
    tree = f.Get('tree')

    binMin = 0.9
    binMax = 1.3
    nbins = 30

    # define object to fit
    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    branch_name = 'scale'
    tree.Draw("{}>>hist".format(branch_name))

    scale = ROOT.RooRealVar("scale","scale", binMin, binMax)
    scale.setBins(nbins)

    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(scale), hist)

    # Define the PDF to fit: 
    mean_left = ROOT.RooRealVar("mean_left","mean_left", 1.04, 1., 1.1)
    sigma_left = ROOT.RooRealVar("sigma_left","sigma_left", 0.03, 0.005, 0.1)

    mean_right = ROOT.RooRealVar("mean_right","mean_right", 1.14, 1.1, 1.2)
    sigma_right = ROOT.RooRealVar("sigma_right","sigma_right", 0.03, 0.005, 0.1)

    gauss_left = ROOT.RooGaussian("gauss_left", "gauss_left", scale, mean_left, sigma_left)
    gauss_right = ROOT.RooGaussian("gauss_right", "gauss_right", scale, mean_right, sigma_right)

    # defines the relative importance of the two pdfs
    frac = ROOT.RooRealVar("frac","frac", 0.5, 0.0 ,1.0)

    # build the full model
    model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(gauss_left, gauss_right), ROOT.RooArgList(frac))

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    frame = scale.frame(ROOT.RooFit.Title(" "))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    result = model.fitTo(rdh)

    # plot the fit    
    model.plotOn(frame, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name('gauss_left'), ROOT.RooFit.Components('gauss_left'), ROOT.RooFit.LineStyle(ROOT.kDashed))
    model.plotOn(frame, ROOT.RooFit.LineColor(8), ROOT.RooFit.Name('gauss_right'), ROOT.RooFit.Components('gauss_right'), ROOT.RooFit.LineStyle(ROOT.kDashed))
    model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("model"), ROOT.RooFit.Components("model"))

    # and write the fit parameters
    model.paramOn(frame,   
         ROOT.RooFit.Layout(0.1, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # we compute the chisquare
    chisquare = frame.chiSquare("model","data")

    # and print it
    label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    label1.SetTextSize(0.03)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    qte = '#chi^{2}/ndof'
    label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    print "chisquare = {}".format(chisquare)

    canv.cd()

    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitle("correction scale")
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()
    label1.Draw()

    # save output
    canv.cd()
    outputdir = self.tools.getOutDir('./myPlots/mc_corrections', 'fits')
    canv.SaveAs("{}/fit.png".format(outputdir))
    canv.SaveAs("{}/fit.pdf".format(outputdir))

    # create workspace
    root_file = ROOT.TFile('workspace_mc_corrections.root', 'RECREATE')
    workspace = ROOT.RooWorkspace('workspace', 'workspace')
    getattr(workspace, 'import')(model)
    workspace.Write()
    workspace.Print()
    root_file.Close()


  def generateToys(self):
    f = ROOT.TFile.Open('workspace_mc_corrections.root')
    workspace = f.Get('workspace')
    workspace.Print()

    pdf = workspace.pdf('model')

    scale = ROOT.RooRealVar('scale', 'scale', 0.9, 1.5)
    data = pdf.generate(ROOT.RooArgSet(scale), 100)

    canv = self.tools.createTCanvas(name='canv', dimx=1200, dimy=1000)
    frame = scale.frame()
    data.plotOn(frame)
    frame.Draw()
    canv.SaveAs('myPlots/mc_corrections/toys.png')


  def createCompositeFunction(self):

    gauss_left = ROOT.TF1('gauss_left', 'gaus', 0.9, 1.3)
    gauss_left.SetParameters(1, 1.042, 0.034)
    gauss_left.SetLineColor(2)

    gauss_right = ROOT.TF1('gauss_right', 'gaus', 0.9, 1.3)
    gauss_right.SetParameters(1, 1.144, 0.026)
    gauss_right.SetLineColor(8)
    
    double_gauss = ROOT.TF1("double_gauss","[0]/sqrt(2.*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+(1-[0])/sqrt(2.*TMath::Pi())/[4]*exp(-(x-[3])*(x-[3])/2./[4]/[4])",0.9,1.3)
    double_gauss.SetParameter(0, 0.62)
    double_gauss.SetParameter(1, 1.042)
    double_gauss.SetParameter(2, 0.034)
    double_gauss.SetParameter(3, 1.144)
    double_gauss.SetParameter(4, 0.026)
    double_gauss.SetLineColor(4)

    canv = self.tools.createTCanvas(name='canv', dimx=1200, dimy=1000)
    #gauss_left.Draw()
    #gauss_right.Draw('same')
    double_gauss.Draw()
    canv.SaveAs('myPlots/mc_corrections/function.png')


  def plot(self, 
      title='', 
      branchname='flat', 
      treename='tree', 
      add_weight_hlt=False, 
      add_weight_pu=False, 
      add_weight_muid=False, 
      weight_hlt='', 
      weight_puqcd='',
      weight_pusig='',
      weight_mu0id='', 
      weight_muid='', 
      plot_data=True, 
      plot_sig=True, 
      plot_ratio=True, 
      do_shape=True, 
      do_luminorm=False, 
      do_stack=True, 
      do_log=False, 
      add_overflow=True, 
      add_CMSlabel=True, 
      CMS_tag='Preliminary', 
      do_tdrstyle=True,
      scale='',
      smear=''
      ):

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, scale, smear, do_log)
    canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
    canv.cd()

    # define the pads
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1) if plot_ratio else ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
    if plot_ratio: pad_up.SetBottomMargin(0.03)
    if do_log: pad_up.SetLogy()
    pad_up.Draw()
    canv.cd()
    if plot_ratio:
      pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
      pad_down.SetBottomMargin(0.25)
      pad_down.Draw()

    # prepare the legend
    if do_stack:
      legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.45, xmax=0.84, ymax=0.83, size=0.027)
    else:
      if not do_tdrstyle:
        legend = self.tools.getRootTLegend(xmin=0.47, ymin=0.65, xmax=0.84, ymax=0.83, size=0.027)
      else:
        legend = self.tools.getRootTLegend(xmin=0.42, ymin=0.57, xmax=0.79, ymax=0.79, size=0.038)

    pad_up.cd()

    # data
    if plot_data:
      hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
      hist_data_tot.Sumw2()

      int_data_tot = 0.
      overflow_data_tot = 0
      error_overflow_data_tot = 0.

      hist_data_stack = ROOT.THStack('hist_data_stack', '')
      data_label = ''

      for idata, data_file in enumerate(self.data_files):
        f_data = self.tools.getRootFile(data_file.filename)
        tree_data = self.tools.getTree(f_data, treename)
        hist_data_name = 'hist_data_{}_{}_{}'.format(self.quantity.label, scale, smear)
        # apply sweight on data
        hist_data = self.createHisto(tree_data, self.quantity, hist_name=hist_data_name, branchname=branchname, selection=self.selection, weight='nbkg_sw')
        hist_data.Sumw2()
        if add_overflow:
          overflow_data_tot = (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
          error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
          hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
          hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
          hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
          hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)
        if do_shape: 
          int_data_tot += hist_data.Integral()

        hist_data_tot.Add(hist_data)
        if idata == 0: data_label = data_file.label[:20]
        else: data_label += '_{}'.format(data_file.label[:20])

        version_label = data_file.label[data_file.label.find('('):data_file.label.find(')')+1]

        if do_shape and hist_data.Integral() != 0: 
          hist_data.Scale(1./hist_data.Integral())
          hist_data_stack.Add(hist_data)

      if do_shape and int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

      if not do_tdrstyle:
        legend.AddEntry(hist_data_tot, 'data - {}'.format(self.getDataLabel(data_label, version_label) if len(self.data_files)>1 else data_file.label))
      else:
        legend.AddEntry(hist_data_tot, 'data-driven background')

      ## set the style
      hist_data_tot.SetMarkerStyle(20)
      #hist_data_tot.SetFillColor(ROOT.kBlue-3)
      #hist_data_tot.SetFillStyle(3005)

    # signal
    if plot_sig:
      signal_hists = []
      for signal_file in self.signal_files:
        f_sig = self.tools.getRootFile(signal_file.filename)
        tree_sig = self.tools.getTree(f_sig, treename)

        hist_signal_name = 'hist_signal_{}_{}_{}'.format(self.quantity.label, scale, smear)
        matching_selection = 'ismatched==1' if branchname == 'flat' else 'BToMuMuPi_isMatched==1'
        selection_signal = matching_selection if self.selection == '' else matching_selection + ' && ' + self.selection
        weight_sig = '(1)'
        #weight_sig = '(weight_dxy_sig_bs)'
        if add_weight_hlt : weight_sig += ' * ({})'.format(weight_hlt)
        if add_weight_pu : weight_sig += ' * ({})'.format(weight_pusig)
        if add_weight_muid : weight_sig += ' * ({}) *({})'.format(weight_mu0id, weight_muid)

        hist_signal = self.createHisto(tree_sig, self.quantity, hist_name=hist_signal_name, branchname=branchname, selection=selection_signal, weight=weight_sig, scale=scale, smear=smear)
        hist_signal.Sumw2()
        if add_overflow:
          overflow_signal_tot = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
          error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
          hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
          hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
          hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
          hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)
        if do_shape: 
          int_signal = hist_signal.Integral()
          if int_signal != 0: hist_signal.Scale(1/int_signal)

        legend.AddEntry(hist_signal, 'signal - {}'.format(signal_file.label))

        #hist_signal.SetLineWidth(3)
        #hist_signal.SetLineColor(signal_file.colour)
        hist_signal.SetFillColor(ROOT.kOrange+7)
        hist_signal.SetFillStyle(3005)
        signal_hists.append(hist_signal)
      
    frame = hist_data_tot.Clone('frame')

    frame.SetTitle('')
    if not plot_ratio: 
      frame.GetXaxis().SetTitle(quantity.title)
      frame.GetXaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
      frame.GetXaxis().SetTitleSize(0.042)
      frame.GetXaxis().SetTitleOffset(1.1)
    if plot_ratio:
      frame.GetXaxis().SetLabelSize(0.0)
      frame.GetXaxis().SetTitleSize(0.0)
    frame.GetYaxis().SetTitle('Entries' if not do_shape else 'Normalised to unity')
    frame.GetYaxis().SetLabelSize(0.033 if not plot_ratio else 0.037)
    frame.GetYaxis().SetTitleSize(0.042)
    frame.GetYaxis().SetTitleOffset(1.3 if not plot_ratio else 1.1)
    frame.GetYaxis().SetRangeUser(1e-4, self.getMaxRangeY(signal_hists, hist_data_tot, do_log, use_sig=True))

    #ROOT.gStyle.SetPadLeftMargin(0.16) 
    ROOT.gStyle.SetOptStat(0)

    # draw the distributions
    frame.Draw()
    
    hist_data_tot.Draw('PE same')
    for hist_sig in signal_hists:
      hist_sig.Draw('histo same')
    hist_data_tot.Draw('PE same')

    hist_sig_err = signal_hists[0].Clone('hist_sig_err')
    hist_sig_err.SetLineWidth(0)
    hist_sig_err.SetFillStyle(3244)
    hist_sig_err.SetFillColor(ROOT.kGray+2)
    hist_sig_err.Draw('E2 same')

    # draw the legend
    legend.Draw('same')

    # add labels
    if not do_tdrstyle:
      self.tools.printLatexBox(0.65, 0.86, title, size=0.04 if plot_ratio else 0.036)
    else:
      self.tools.printLatexBox(0.60, 0.84, title, size=0.04 if plot_ratio else 0.038)
    if add_CMSlabel: self.tools.printCMSTag(pad_up, CMS_tag, size=0.55 if plot_ratio else 0.43)
    if do_luminorm:
      scale_text = ROOT.TPaveText(0.15, 0.83, 0.3, 0.88, "brNDC")
      scale_text.SetBorderSize(0)
      scale_text.SetFillColor(ROOT.kWhite)
      scale_text.SetTextSize(0.032)
      scale_text.SetTextAlign(31)
      scale_text.SetTextFont(42)
      scale_text.AddText('MC scaled by {}'.format(round(lumi_weight/self.tools.getDataLumi(self.data_files), 2)))
      scale_text.Draw()

    # plot the ratio
    if plot_ratio:
      pad_down.cd()

      hist_ratio = self.tools.getRatioHistogram(hist_data_tot, signal_hists[0])
      hist_ratio.Sumw2()

      hist_ratio.SetLineWidth(2)
      hist_ratio.SetMarkerStyle(20)
      hist_ratio.SetTitle('')
      hist_ratio.GetXaxis().SetTitle(self.quantity.title)

      hist_ratio.GetXaxis().SetLabelSize(0.1)
      hist_ratio.GetXaxis().SetTitleSize(0.13)
      hist_ratio.GetXaxis().SetTitleOffset(0.73)
      hist_ratio.GetYaxis().SetTitle('Data/MC')
      hist_ratio.GetYaxis().SetLabelSize(0.1)
      hist_ratio.GetYaxis().SetTitleSize(0.13)
      hist_ratio.GetYaxis().SetTitleOffset(0.345)
      val_min = hist_ratio.GetBinContent(hist_ratio.GetMinimumBin())
      val_max = hist_ratio.GetBinContent(hist_ratio.GetMaximumBin())
      hist_ratio.GetYaxis().SetRangeUser(val_min-0.15*val_min, val_max+0.15*val_max)

      hist_ratio.Draw('PE')

      ## draw line at ratio = 1
      line = ROOT.TLine(self.quantity.bin_min, 1, self.quantity.bin_max, 1)
      line.SetLineColor(4)
      line.SetLineWidth(2)
      line.Draw('same')

    if not self.save_eos:
      outputdir = self.tools.getOutDir('./myPlots/mc_corrections', outdirlabel=quantity.label, do_log=do_log)
    else:
      outputdir = self.tools.getOutDir('/eos/home-a/anlyon/www/BHNL/mc_corrections', outdirlabel=quantity.label, do_log=do_log)
    
    if scale != None and smear != None:
      name = '{}_scale{}_smear{}'.format(self.quantity.label, str(round(scale, 3)).replace('.', 'p'), smear[smear.rfind('_')+1:])
    else:
      name = '{}'.format(self.quantity.label)
    canv.SaveAs('{}/{}.png'.format(outputdir, name))
    canv.SaveAs('{}/{}.pdf'.format(outputdir, name))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  dxy = Quantity(name_flat='probe_dxy_bs_uncorrected', label='probe_dxy_bs', title='probe #mu |d_{xy}| (BS) [cm]', nbins=90, bin_min=0, bin_max=0.5)
  dxy_sig = Quantity(name_flat='probe_dxy_sig_bs_uncorrected', label='probe_dxy_sig_bs', title='probe #mu |d_{xy}| significance (BS)', nbins=90, bin_min=0, bin_max=60)
  lxy = Quantity(name_flat='lxy', label='lxy', title='#mu#mu l_{xy} [cm]', nbins=60, bin_min=0, bin_max=3)
  lxy_sig = Quantity(name_flat='lxy_sig', label='lxy_sig', title='#mu#mu l_{xy} significance', nbins=60, bin_min=0, bin_max=200)
  tag_pt = Quantity(name_flat='tag_pt', label='tag_pt', title='tag #mu pT [GeV]', nbins=60, bin_min=7, bin_max=20)
  probe_pt = Quantity(name_flat='probe_pt', label='probe_pt', title='probe #mu pT [GeV]', nbins=60, bin_min=7, bin_max=20)
  #dxy_sig = Quantity(name_flat='probe_dxy_sig_bs_corrected_linearscale', label='probe_dxy_sig_bs', title='probe #mu |d_{xy}| significance (BS)', nbins=90, bin_min=0, bin_max=60)
  #quantities = [dxy, dxy_sig]
  #quantities = [dxy, dxy_sig, lxy, lxy_sig]
  #quantities = [lxy_sig]
  quantities = [tag_pt, probe_pt]

  data_files = data_samples['tag_and_probe']
  signal_files = signal_samples['tag_and_probe_BToJPsiKstar']

  baseline_selection = selection['tag_and_probe'].flat

  save_eos = False

  scales = np.arange(1., 1.3, 0.01)
  smears = [
            'smeared_corr_0p01',
            #'smeared_corr_0p02',
            'smeared_corr_0p03',
            #'smeared_corr_0p04',
            'smeared_corr_0p05', 
            #'smeared_corr_0p06',
            #'smeared_corr_0p07',
            #'smeared_corr_0p08',
            'smeared_corr_0p1',
            #'smeared_corr_0p15',
            'smeared_corr_0p2',
            'smeared_corr_0p25',
            'smeared_corr_0p3',
            'smeared_corr_0p35',
            'smeared_corr_0p4',
            'smeared_corr_0p45',
            'smeared_corr_0p5',
           ]
  scales = [None]
  smears = [None]

  do_scanCorrections = True
  do_createWeightFile = False
  do_createScaleDistribution = False
  do_performScaleFit = False

  if do_scanCorrections:
    for quantity in quantities:
      plotter = Plotter(quantity=quantity, data_files=data_files, signal_files=signal_files, selection=baseline_selection, save_eos=save_eos)
      for scale in scales:
        for smear in smears:
          plotter.plot(scale=scale, smear=smear, do_log=True)
          plotter.plot(scale=scale, smear=smear, do_log=False)

  if do_createWeightFile:
    dxy = Quantity(name_flat='probe_dxy_bs_uncorrected', label='probe_dxy_bs', title='probe #mu |d_{xy}| (BS) [cm]', nbins=60, bin_min=0, bin_max=0.3)
    dxy_sig = Quantity(name_flat='probe_dxy_sig_bs_uncorrected', label='probe_dxy_sig_bs', title='probe #mu |d_{xy}| significance (BS)', nbins=60, bin_min=0, bin_max=60)

    plotter = Plotter(quantity=dxy, data_files=data_files, signal_files=signal_files, selection=baseline_selection)
    plotter.computeWeight()

    plotter = Plotter(quantity=dxy_sig, data_files=data_files, signal_files=signal_files, selection=baseline_selection)
    plotter.computeWeight()

  if do_createScaleDistribution:
    plotter = Plotter(data_files=data_files, signal_files=signal_files, selection=baseline_selection)
    plotter.createScaleDistribution()

  if do_performScaleFit:
    plotter = Plotter()
    plotter.performScaleFit()
    plotter.generateToys()
    plotter.createCompositeFunction()


