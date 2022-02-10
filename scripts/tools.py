import os
import os.path
from os import path
import ROOT

import sys
sys.path.append('../objects')
from quantity import Quantity


class Tools(object):
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


  def createHisto(self, tree, quantity, hist_name='hist', branchname='flat', weight=-99, selection=''):
    ROOT.TH1.SetDefaultSumw2()
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


  def createWeightedHistoQCDMC(self, qcd_files, white_list='', quantity='', hist_name='hist_qcd_tot', selection='', add_weight_hlt=False, add_weight_pu=False, weight_hlt='', weight_puqcd=''):
    hist_mc_tot = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)
    hist_mc_tot.Sumw2()

    #lumi_tot = self.getQCDMCLumi(qcd_files, white_list)

    for ifile, qcd_file in enumerate(qcd_files):
      qcd_file_pthatrange = self.getPthatRange(qcd_file.label)
      if qcd_file_pthatrange not in white_list: continue

      f_mc = self.getRootFile(qcd_file.filename)
      tree_run = self.getTree(f_mc, 'run_tree')
      tree_mc = self.getTree(f_mc, 'signal_tree')

      weight_mc = self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
      #weight_mc = 1. / (self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency) * lumi_tot)
      weight_qcd = '({})'.format(weight_mc)
      if add_weight_hlt : weight_qcd += ' * ({})'.format(weight_hlt)
      if add_weight_pu : weight_qcd += ' * ({})'.format(weight_puqcd)
      hist_name = 'hist_qcd_{}'.format(qcd_file_pthatrange)
      hist_mc = self.createHisto(tree_mc, quantity, hist_name=hist_name, branchname='flat', selection=selection, weight=weight_qcd) 
      hist_mc.Sumw2()
    
      hist_mc_tot.Add(hist_mc)

    return hist_mc_tot


  def getPthatRange(self, label):
    '''
      Note that the qcd file label must contain XXtoYY
    '''
    return label[label.find('pt'):label.find(' ')]


  def getNminiAODEvts(self, tree_run):
    quantity = Quantity(name_nano='genEventCount', name_flat='geneventcount', nbins=100, bin_min=1, bin_max=1e9)
    #TODO apply weight here?
    hist = self.createHisto(tree=tree_run, quantity=quantity)
    n_reco_evts = hist.GetMean() * hist.GetEntries()
    return n_reco_evts


  def computeQCDMCWeight(self, tree, cross_section, filter_efficiency):
    '''
      QCD MC weight defined as cross-section / nb of generated events
    '''
    n_reco_evts = self.getNminiAODEvts(tree)
    n_genevts = n_reco_evts / filter_efficiency
    weight = cross_section * 1e03 / n_genevts # 1e03 factor for pb -> fb
    return weight


  def scaleToLumiWeight(self, data_files, qcd_files, white_list, selection, add_weight_hlt, add_weight_pu, weight_hlt, weight_puqcd):
    '''
      weight = lumi_data / lumi_mc = N_data * sigma_mc / (N_mc * sigma_data) estimated as N_data / N_mc
    '''

    quantity_forweight = Quantity(name_flat='hnl_mass', nbins=1, bin_min=0, bin_max=13000)
    n_obs_data = 0.
    n_err_data = 0.
    for data_file in data_files:
      f_data = self.getRootFile(data_file.filename)
      tree_data = self.getTree(f_data, 'signal_tree')
      hist_data = self.createHisto(tree_data, quantity_forweight, branchname='flat', selection=selection)
      hist_data.Sumw2()
      n_obs_data += hist_data.Integral()
      n_err_data += hist_data.GetBinError(1)

    hist_mc_tot = self.createWeightedHistoQCDMC(qcd_files=qcd_files, white_list=white_list, quantity=quantity_forweight, selection=selection, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, weight_hlt=weight_hlt, weight_puqcd=weight_puqcd)
    n_obs_mc = hist_mc_tot.Integral()
    n_err_mc = hist_mc_tot.GetBinError(1)

    weight = float(n_obs_data) / float(n_obs_mc) if float(n_obs_mc)!= 0. else 0.

    #if n_obs_data != 0 and n_obs_mc != 0:
    #  err = weight* (n_err_data / n_obs_data + n_err_mc / n_obs_mc)
    #else: 
    #  err = 0

    return weight


  def getSignalWeight(self, signal_file, sigma_B, lumi, isBc=False):
    '''
      weight = sigma_B * lumi * v_square * BR(B->muNX) * BR(N->mupi) * filter_eff / N_mini   
    '''
    # coupling square
    v_square = self.getVV(mass=signal_file.mass, ctau=signal_file.ctau, ismaj=True)

    # production branching ratio
    from decays import Decays 
    dec = Decays(mass=signal_file.mass, mixing_angle_square=1) # we factorise the mixing angle 
    if not isBc: BR_prod = dec.BR_tot_mu 
    else: BR_prod = dec.BR_Bc_mu 

    # decay branching ratio
    BR_NToMuPi = self.gamma_partial(mass=signal_file.mass, vv=v_square) / self.gamma_total(mass=signal_file.mass, vv=v_square)

    # number of generated events
    filter_efficiency = signal_file.filter_efficiency if not isBc else signal_file.filter_efficiency_Bc
    f = self.getRootFile(signal_file.filename)
    tree_run = self.getTree(f, 'run_tree')
    n_gen = self.getNminiAODEvts(tree_run)
    n_generated = 0.5 * n_gen / filter_efficiency # 0.5 factor since samples were produced 50% muon 50% electron  

    weight = sigma_B / 0.4 * lumi * v_square * BR_prod * BR_NToMuPi / n_generated 

    return weight


  def getCtauWeight(self, signal_file):
    filename = signal_file.filename # might need to be modified later on for Bc
    #original_ctau = filename[filename.find('ctau')+4:filename.find('/', filename.find('ctau')+1)]
    original_ctau = filename[filename.find('ctau')+4:filename.find('mm', filename.find('ctau')+1)]
    original_ctau = original_ctau.replace('p', '.')
    target_ctau = signal_file.ctau

    ctau_weight = 1. #-99
    if float(original_ctau) != float(target_ctau):
      ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

    return ctau_weight


  def getQCDMCLumi(self, qcd_files, white_list):
    lumi_tot = 0

    for qcd_file in qcd_files:
      qcd_file_pthatrange = self.getPthatRange(qcd_file.label)
      if qcd_file_pthatrange not in white_list: continue

      f_mc = self.getRootFile(qcd_file.filename)
      tree_run = self.getTree(f_mc, 'run_tree')
      
      lumi_inv = self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
      lumi = 1./lumi_inv
      lumi_tot += lumi

    return lumi_tot


  def getDataLumi(self, data_files):
    lumi_tot = 0.
    for data_file in data_files:
      lumi_tot += data_file.lumi

    return lumi_tot


  def getRatioHistogram(self, hist1, hist2): 
    hist_ratio = hist1.Clone('hist_ratio')
    hist_ratio.Divide(hist2)
    return hist_ratio


  def createTCanvas(self, name, dimx=1200, dimy=1000):
    canv = ROOT.TCanvas(name, name, dimx, dimy)
    ROOT.SetOwnership(canv, False)
    return canv


  def getRootTLegend(self, xmin=0.65, ymin=0.7, xmax=0.85, ymax=0.9, size=0.03, do_alpha=True):
    legend = ROOT.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(size)
    legend.SetLineColor(0)
    if do_alpha: legend.SetFillColorAlpha(0, 0)
    else: legend.SetFillColor(0)
    legend.SetBorderSize(0)
    return legend
    

  def getTextBox(self, style, xmin, ymin, xmax, ymax, text, size=0.11, colour=ROOT.kBlack):
    box = ROOT.TPaveText(xmin, ymin, xmax, ymax, style)
    box.AddText(text)
    box.SetBorderSize(0)
    box.SetFillColorAlpha(0, 0)
    box.SetTextColor(colour)
    box.SetTextSize(size)
    box.SetTextFont(42)
    box.SetTextAlign(11)
    return box


  def printLatexBox(self, x, y, text, size=0.04):
    box = ROOT.TLatex()
    box.SetNDC()
    box.SetTextFont(42)
    box.SetTextAlign(22) 
    box.SetTextSize(size)    
    box.DrawLatex(x, y, text)
  
  
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


  def printCMSTag(self, pad, cms_tag, size=0.55, offset=0.08):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin(), 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), 'CMS')
    # print CMS tag
    tag.SetTextFont(52)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    tag.DrawLatex(pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), cms_tag)      
    pad.Update()


  def printLumiTag(self, pad, lumi, size=0.43, offset=0.57):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    lumi_text = str(round(lumi, 2)) + ' fb^{-1} (13 TeV)'
    tag.SetTextFont(42)
    tag.SetTextAlign(11) 
    tag.SetTextSize(0.9*size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), lumi_text)
    pad.Update()


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




  def getOutDir(self, maindir, outdirlabel, do_shape=False, do_luminorm=False, do_stack=False, do_log=False, add_overflow=False):
    if not path.exists(maindir):
      os.system('mkdir -p {}'.format(maindir))
    os.system('cp ../data/index.php {}'.format(maindir))
    dirlabel = outdirlabel

    outputdir = '{}/{}'.format(maindir, dirlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ../data/index.php {}'.format(outputdir))
    os.system('cp ../data/index.php {}/..'.format(outputdir))

    norm = None
    if do_shape: norm = 'shape'
    elif do_luminorm: norm = 'luminorm'

    #if not do_shape and not do_stack and not do_log: dirlabel += '/plain'
    if do_stack: 
      if norm==None and not do_log and not do_luminorm: dirlabel += '/stack'
      elif norm==None and do_log:  dirlabel += '/stack/log'
      elif norm!=None and not do_log and not add_overflow: dirlabel += '/stack_{}'.format(norm)
      elif norm!=None and not do_log and add_overflow: dirlabel += '/stack_{}_overflow'.format(norm)
      elif norm!=None and do_log and not add_overflow: dirlabel += '/stack_{}_log'.format(norm)
      else: dirlabel += '/stack_{}_log_overflow'.format(norm)
    else:
      if norm==None and do_log: dirlabel += '/log'
      elif norm==None and not do_log: dirlabel += '/lin'
      elif norm!=None and not do_log and not add_overflow: dirlabel += '/{}'.format(norm)
      elif norm!=None and do_log and not add_overflow: dirlabel += '/{}_log'.format(norm)
      elif norm!=None and not do_log and add_overflow: dirlabel += '/{}_overflow'.format(norm)
      else: dirlabel += '/{}_log_overflow'.format(norm)

    outputdir = '{}/{}'.format(maindir, dirlabel)
    
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ../data/index.php {}'.format(outputdir))

    return outputdir


  def ctau_from_gamma(self, gamma):
    const_hbar = 6.582119569e-22 * 1e-03
    const_c = 299792458.                             
    tau_natural = 1. / gamma        # 1/GeV
    tau = tau_natural * const_hbar  # s
    ctau = tau * const_c * 1000     # mm
    return ctau


  def gamma_total(self, mass, vv):
    '''
    Total width for N (Dirac)
    '''
    from decays import HNLDecays
    gamma_total =  HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['tot']  # GeV
    return gamma_total


  def gamma_partial(self, mass, vv):
    '''
    Partial width for N->mupi (Dirac)
    '''
    from decays import HNLDecays
    gamma_partial = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['mupi'] # GeV
    return gamma_partial


  def getVV(self, mass=-99.,ctau=-99., ismaj=True):
    '''
    Helper function to go from ctau,m -> vv
    '''
    import numpy as np
    mult = 2. if ismaj else 1.
    ref_m = 1. # GeV
    ref_vv = 1. 
    ref_ctau = self.ctau_from_gamma(mult*self.gamma_total(mass=ref_m,vv=ref_vv))
    k = ref_ctau * np.power(ref_m,5) * ref_vv
    return k/(np.power(mass, 5) * ctau)


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)

