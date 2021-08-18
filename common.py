import os
import os.path
from os import path
import ROOT

from quantity import Quantity


class PlottingTools(object):
  def getRootFile(self, filename, with_ext=False):
    if not with_ext:
      f = ROOT.TFile.Open(filename, 'READ')
    else:
      f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+filename, 'READ')
    return f


  def getTree(self, rootfile, tree_name):
    tree = rootfile.Get(tree_name)
    if not tree:
      raise RuntimeError('Tree of name "{}" was not found in {}'.format(tree_name, rootfile.GetName()))
    return tree


  def createHisto(self, rootfile, tree_name, quantity, hist_name='hist', branchname='flat', weight=-99, selection='', cut=0):
    ROOT.TH1.SetDefaultSumw2()
    tree = self.getTree(rootfile, tree_name)
    hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)

    if selection == '' and weight == -99: selection_string = selection
    elif selection != '' and weight == -99: selection_string = '({})'.format(selection)
    elif selection == ''  and weight != -99: selection_string = '({})'.format(weight)
    else: selection_string = '({sel}) * ({wght})'.format(sel=selection, wght=weight)

    #print selection_string

    qte = quantity.name_flat if branchname != 'nano' else quantity.name_nano
    tree.Project(hist_name, qte, selection_string)

    hist.SetDirectory(0)
    return hist


  def createWeightedHistoQCDMC(self, qcd_files, white_list='', quantity='', selection=''):
    hist_mc_tot = ROOT.TH1D('hist_mc_tot', 'hist_mc_tot', quantity.nbins, quantity.bin_min, quantity.bin_max)
    hist_mc_tot.Sumw2()

    for ifile, qcd_file in enumerate(qcd_files):
      if qcd_file.label not in white_list: continue

      f_mc = self.getRootFile(qcd_file.filename)
      
      #weight_mc = self.computeQCDMCWeight(self.getTree(f_mc, 'signal_tree'), qcd_file.cross_section, qcd_file.filter_efficiency)
      weight_mc = self.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
      weight = '({}) * (weight_hlt)'.format(weight_mc)
      #weight = '({})'.format(weight_mc)
      hist_mc = self.createHisto(f_mc, 'signal_tree', quantity, branchname='flat', selection=selection, weight=weight) 
      hist_mc.Sumw2()
    
      hist_mc_tot.Add(hist_mc)

    return hist_mc_tot


  def getNminiAODEvts(self, rootfile):
    quantity = Quantity(name_nano='genEventCount', name_flat='geneventcount', nbins=100, bin_min=1, bin_max=1e9)
    hist = self.createHisto(rootfile=rootfile, tree_name='run_tree', quantity=quantity)
    n_reco_evts = hist.GetMean() * hist.GetEntries()
    return n_reco_evts


  #def computeQCDMCWeight(self, tree, cross_section, filter_efficiency):
  def computeQCDMCWeight(self, rootfile, cross_section, filter_efficiency):
    '''
      QCD MC weight defined as cross-section / nb of generated events
    '''
    #hist_name = 'hist_{}_{}'.format(cross_section, filter_efficiency)
    #quantity = Quantity(name_nano='genEventCount', name_flat='geneventcount', nbins=100, bin_min=1, bin_max=1e9)
    #hist = self.createHisto(rootfile=rootfile, tree_name='run_tree', quantity=quantity)
    #n_reco_evts = hist.GetMean() * hist.GetEntries()
    n_reco_evts = self.getNminiAODEvts(rootfile)
    ##n_reco_evts = self.getTree(rootfile, 'signal_tree').GetEntries()
    n_genevts = n_reco_evts / filter_efficiency
    weight = cross_section / n_genevts
    return weight


  def getRatioHistogram(self, hist1, hist2): 
    hist_ratio = hist1.Clone('hist_ratio')
    hist_ratio.Divide(hist2)
    return hist_ratio


  def createTCanvas(self, name, dimx=1200, dimy=1000):
    canv = ROOT.TCanvas(name, name, dimx, dimy)
    ROOT.SetOwnership(canv, False)
    return canv


  def getRootTLegend(self, xmin=0.65, ymin=0.7, xmax=0.85, ymax=0.9, size=0.03):
    legend = ROOT.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(size)
    legend.SetLineColor(0)
    legend.SetFillColorAlpha(0, 0)
    legend.SetBorderSize(0)
    return legend
    

  def getTextBox(self, style, xmin, ymin, xmax, ymax, text, colour):
    box = ROOT.TPaveText(xmin, ymin, xmax, ymax, style)
    box.AddText(text)
    box.SetBorderSize(0)
    box.SetFillColor(ROOT.kWhite)
    box.SetTextColor(colour)
    box.SetTextSize(0.11)
    box.SetTextFont(42)
    box.SetTextAlign(11)
    return box
  
  def getRootXAxis(self, hist, title='', label_size=0.037, title_size=0.042, offset=1.1, xmin=-99, xmax=-99): 
    ROOT.SetOwnership(hist, True)
    hist.GetXaxis().SetTitle(title)
    hist.GetXaxis().SetLabelSize(label_size)
    hist.GetXaxis().SetTitleSize(title_size)
    hist.GetXaxis().SetTitleOffset(offset)
    if xmin != -99 and xmax != -99:
      hist.GetXaxis().SetRangeUser(xmin, xmax)
    return hist


  def getRootYAxis(self, hist, title='', label_size=0.037, title_size=0.042, offset=1.1, ymin=-99, ymax=-99): 
    ROOT.SetOwnership(hist, False)
    hist.GetYaxis().SetTitle(title)
    hist.GetYaxis().SetLabelSize(label_size)
    hist.GetYaxis().SetTitleSize(title_size)
    hist.GetYaxis().SetTitleOffset(offset)
    if ymin != -99 and ymax != -99:
      hist.GetXaxis().SetRangeUser(ymin, ymax)
    return hist


  #def getRatioGraph(self, hist1, hist2):
  #  graph = ROOT.TGraphAsymmErrors()
  #  for ibin, bin_ in enumerate(bins):
  #    bin_min, bin_max = bin_
  #    x = (float(bin_min) + float(bin_max))/2.
  #    y = efficiency[ibin][0][imatch] if binning == 'displacement' else efficiency[0][ibin][imatch]
  #    err = error[ibin][0][imatch] if binning == 'displacement' else error[0][ibin][imatch]
  #    point = graph.GetN()
  #    graph.SetPoint(point, x, y)
  #    graph.SetPointError(point, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err, err)




  def getOutDir(self, maindir, outdirlabel, do_shape=False, do_luminorm=False, do_stack=False, do_log=False):
    if not path.exists(maindir):
      os.system('mkdir -p {}'.format(maindir))
    os.system('cp ./data/index.php {}'.format(maindir))
    dirlabel = outdirlabel

    outputdir = '{}/{}'.format(maindir, dirlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ./data/index.php {}'.format(outputdir))

    norm = None
    if do_shape: norm = 'shape'
    elif do_luminorm: norm = 'luminorm'

    #if not do_shape and not do_stack and not do_log: dirlabel += '/plain'
    if do_stack: 
      if norm==None and not do_log and not do_luminorm: dirlabel += '/stack'
      elif norm!=None and not do_log: dirlabel += '/stack_{}'.format(norm)
      else: dirlabel += '/stack_{}_log'.format(norm)
    else:
      if norm!=None and not do_log: dirlabel += '/{}'.format(norm)
      elif norm!=None and do_log: dirlabel += '/{}_log'.format(norm)
      elif norm==None and do_log: dirlabel += '/log'

    outputdir = '{}/{}'.format(maindir, dirlabel)
    
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ./data/index.php {}'.format(outputdir))

    return outputdir


  def getColour(self):
    return [
      ROOT.kBlue+2,
      ROOT.kBlue-11,
      ROOT.kBlue-9,
      ROOT.kBlue-10,
      ROOT.kRed-10,
      ROOT.kOrange-9,
      ROOT.kOrange+6,
      ROOT.kOrange-3,
      ROOT.kOrange+8
    ]
