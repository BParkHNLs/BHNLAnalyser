import os
import os.path
from os import path
import ROOT

import sys
sys.path.append('../objects')
from quantity import Quantity
from samples import signal_samples


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


  def createWeightedHistoQCDMC(self, qcd_files, white_list='', quantity='', hist_name='hist_qcd_tot', selection='', add_weight_hlt=False, add_weight_pu=False, add_weight_muid=False, weight_hlt='', weight_puqcd='', weight_mu0id='', weight_muid='', treename='signal_tree'):
    hist_mc_tot = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)
    hist_mc_tot.Sumw2()

    #lumi_tot = self.getQCDMCLumi(qcd_files, white_list)

    for ifile, qcd_file in enumerate(qcd_files):
      qcd_file_pthatrange = self.getPthatRange(qcd_file.label)
      if qcd_file_pthatrange not in white_list: continue

      f_mc = self.getRootFile(qcd_file.filename)
      tree_run = self.getTree(f_mc, 'run_tree')
      tree_mc = self.getTree(f_mc, treename)

      weight_mc = self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency)
      #weight_mc = 1. / (self.computeQCDMCWeight(tree_run, qcd_file.cross_section, qcd_file.filter_efficiency) * lumi_tot)
      weight_qcd = '({})'.format(weight_mc)
      if add_weight_hlt : weight_qcd += ' * ({})'.format(weight_hlt)
      if add_weight_pu : weight_qcd += ' * ({})'.format(weight_puqcd)
      if add_weight_muid : weight_qcd += ' * ({}) * ({})'.format(weight_mu0id, weight_muid)
      hist_name = 'hist_qcd_{}'.format(qcd_file_pthatrange)
      hist_mc = self.createHisto(tree_mc, quantity, hist_name=hist_name, branchname='flat', selection=selection, weight=weight_qcd) 
      hist_mc.Sumw2()
    
      hist_mc_tot.Add(hist_mc)

    return hist_mc_tot


  def getSignalFileList(self, signal_label, mass, ctau, strategy='inclusive', is_bc=False):
    #TODO rename as getSignalSampleList?
    file_list = []
    generated_samples = signal_samples[signal_label] 

    # do not reweight in case their is a sample generated with that ctau
    is_generated_ctau = False
    for generated_sample in generated_samples:
      if generated_sample.mass != mass: continue
      filename = generated_sample.filename if not is_bc else generated_sample.filename_Bc
      if ctau == generated_sample.ctau and filename != None:
        file_list.append(generated_sample)
        is_generated_ctau = True

    # uncomment to not reweight points with generated ctau
    is_generated_ctau = False

    if not is_generated_ctau:
      file_list = []
      if strategy == 'inclusive':
        for signal_sample in generated_samples:
          if signal_sample.mass != mass: continue
          filename = signal_sample.filename if not is_bc else signal_sample.filename_Bc
          if filename == None: continue
          file_list.append(signal_sample)

      elif strategy == 'partial_inclusive':
        for signal_sample in generated_samples:
          if signal_sample.mass != mass: continue
          filename = signal_sample.filename if not is_bc else signal_sample.filename_Bc
          if filename == None: continue
          if signal_sample.ctau < ctau: 
          #  if len(file_list) == 0: # for mass 1, we reweight to even larger ctaus than the ones generated
          #    file_list.append(signal_sample)
            file_list.append(signal_sample)
            break
          file_list.append(signal_sample)

      elif strategy == 'exclusive_fromlargerctau':
        # only use the generated sample with the closest larger ctau
        for signal_sample in generated_samples:
          if signal_sample.mass != mass: continue
          filename = signal_sample.filename if not is_bc else signal_sample.filename_Bc
          if filename == None: continue
          if signal_sample.ctau >= ctau:
            file_list.append(signal_sample)
          elif signal_sample.ctau < ctau and len(file_list) == 0: 
            # for mass 1, we reweight to even larger ctaus than the ones generated
            file_list.append(signal_sample)

        if len(file_list) > 1:
          # only keep last element of the list
          file_list_tmp = []
          file_list_tmp.append(file_list.pop())
          file_list = file_list_tmp

      elif strategy == 'exclusive_fromsmallerctau':
        # only use the generated sample with the closest larger ctau
        for signal_sample in generated_samples:
          if signal_sample.mass != mass: continue
          filename = signal_sample.filename if not is_bc else signal_sample.filename_Bc
          if filename == None: continue
          if signal_sample.ctau <= ctau:
            file_list.append(signal_sample)
          elif signal_sample.ctau < ctau and len(file_list) == 0: 
            # for mass 1, we reweight to even larger ctaus than the ones generated
            file_list.append(signal_sample)

        if len(file_list) > 1:
          # only keep last element of the list
          file_list_tmp = []
          file_list_tmp.append(file_list[0])
          file_list = file_list_tmp

        if len(file_list) == 0:
          file_list.append(generated_samples.pop())

      elif strategy == 'unique':
        for signal_sample in generated_samples:
          if signal_sample.mass != mass: continue
          filename = signal_sample.filename if not is_bc else signal_sample.filename_Bc
          if filename == None: continue
          if float(mass) == 1. and 'ctau10p0mm' in signal_sample.filename: 
            file_list.append(signal_sample)
          if float(mass) == 3. and 'ctau10p0mm' in signal_sample.filename: 
            file_list.append(signal_sample)
          if float(mass) == 4.5 and 'ctau10p0mm' in signal_sample.filename: 
            file_list.append(signal_sample)

    return file_list


  def getPthatRange(self, label):
    '''
      Note that the qcd file label must contain XXtoYY
    '''
    if label.find(' ') != -1:
      return label[label.find('pt'):label.find(' ')]
    else:
      return label[label.find('pt'):]


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


  def scaleToLumiWeight(self, data_files, qcd_files, white_list, selection, add_weight_hlt, add_weight_pu, add_weight_muid, weight_hlt, weight_puqcd, weight_mu0id, weight_muid):
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

    hist_mc_tot = self.createWeightedHistoQCDMC(qcd_files=qcd_files, white_list=white_list, quantity=quantity_forweight, selection=selection, add_weight_hlt=add_weight_hlt, add_weight_pu=add_weight_pu, add_weight_muid=add_weight_muid, weight_hlt=weight_hlt, weight_puqcd=weight_puqcd, weight_mu0id=weight_mu0id, weight_muid=weight_muid)
    n_obs_mc = hist_mc_tot.Integral()
    n_err_mc = hist_mc_tot.GetBinError(1)

    weight = float(n_obs_data) / float(n_obs_mc) if float(n_obs_mc)!= 0. else 0.

    #if n_obs_data != 0 and n_obs_mc != 0:
    #  err = weight* (n_err_data / n_obs_data + n_err_mc / n_obs_mc)
    #else: 
    #  err = 0

    return weight


  def getSignalWeight(self, signal_files, mass, ctau, sigma_B, lumi, lhe_efficiency=0.08244, is_bu=False, is_bd=False, is_bs=False, is_bc=False):
    '''
      weight = sigma_B * lumi * v_square * BR(B->muNX) * BR(N->mupi) * filter_eff / N_mini 
    '''
    # coupling square
    v_square = self.getVV(mass=mass, ctau=ctau, ismaj=True)

    # production branching ratio
    from decays import Decays 
    dec = Decays(mass=mass, vv=1., fe=0., fu=1., ft=0.) # we factorise the mixing angle 
    if is_bu:
      BR_prod = dec.BR_B_mu
    elif is_bd:
      BR_prod = dec.BR_B0_mu
    elif is_bs:
      BR_prod = dec.BR_Bs_mu
    elif is_bc:
      BR_prod = dec.BR_Bc_mu_fq_weighted # fc included as is not accounted for in filter efficiency
    else: # previous normalisation strategy (faulty treatment of fq)
      BR_prod = dec.BR_tot_mu_fq_weighted

    # decay branching ratio
    BR_NToMuPi = self.gamma_partial(mass=mass, vv=1., fe=0., fu=1., ft=0.) / self.gamma_total(mass=mass, vv=1., fe=0., fu=1., ft=0.) # the coupling cancels in the ratio, valid for Majorana and Dirac

    # number of generated events (= n_gen / filter_efficiency = n_miniaod / filter_efficiency)
    n_gen_tot = 0
    if is_bu:
      for signal_file in signal_files:
        n_gen_tot += signal_file.n_miniaod_Bu
    elif is_bd:
      for signal_file in signal_files:
        n_gen_tot += signal_file.n_miniaod_Bd
    elif is_bs:
      for signal_file in signal_files:
        n_gen_tot += signal_file.n_miniaod_Bs
    elif is_bc:
      tree_run = ROOT.TChain('run_tree') 
      for signal_file in signal_files:
        tree_run.Add(signal_file.filename_Bc)
      n_gen_tot = self.getNminiAODEvts(tree_run)
    else: # previous normalisation strategy
      tree_run = ROOT.TChain('run_tree') 
      for signal_file in signal_files:
        signal_filename = signal_file.filename
        tree_run.Add(signal_filename)
      n_gen_tot = self.getNminiAODEvts(tree_run)

    # take the weighted average of the filter efficiencies of the generated samples
    filter_efficiency = 0.
    for signal_file in signal_files:
      if is_bu:
        the_filter_efficiency = signal_file.filter_efficiency_Bu
        n_gen = signal_file.n_miniaod_Bu
      elif is_bd:
        the_filter_efficiency = signal_file.filter_efficiency_Bd
        n_gen = signal_file.n_miniaod_Bd
      elif is_bs:
        the_filter_efficiency = signal_file.filter_efficiency_Bs
        n_gen = signal_file.n_miniaod_Bs
      elif is_bc:
        the_filter_efficiency = signal_file.filter_efficiency_Bc
        the_file = self.getRootFile(signal_file.filename_Bc)
        the_tree_run = self.getTree(the_file, 'run_tree')
        n_gen = self.getNminiAODEvts(the_tree_run)
      else: # previous normalisation strategy
        the_filter_efficiency = signal_file.filter_efficiency
        signal_filename = signal_file.filename
        the_file = self.getRootFile(signal_filename)
        the_tree_run = self.getTree(the_file, 'run_tree')
        n_gen = self.getNminiAODEvts(the_tree_run)

      filter_efficiency += n_gen * the_filter_efficiency
    filter_efficiency = filter_efficiency / n_gen_tot
    #print 'filter efficiency {}'.format(filter_efficiency)

    # in the case where the samples were produced with both the electron and muon channels, apply a correction
    # only consider events that were generated in the muon channel
    if not is_bc:
      corr = signal_file.muon_rate 
    else:
      corr = signal_file.muon_rate_Bc 

    efficiency = filter_efficiency if not is_bc else filter_efficiency * lhe_efficiency 
    n_generated = corr * n_gen_tot / efficiency
    print 'n_generated = corr * n_gen_tot / efficiency'
    print 'n_generated = {} * {} / {}'.format(corr, n_gen_tot, efficiency)

    f_u = 0.4 # B fragmentation fraction

    weight = sigma_B / f_u * lumi * v_square * BR_prod * BR_NToMuPi / n_generated 
    print 'sigma_B / f_u * lumi * v_square * BR_prod * BR_NToMuPi / n_generated'
    print 'weight = {} / 0.4 * {} * {} * {} * {} / {}'.format(sigma_B, lumi, v_square, BR_prod, BR_NToMuPi, n_generated)

    return weight


  def getSignalEfficiency(self, signal_files, selection=None, lhe_efficiency=0.08244, is_bc=False):
    #NOTE this method is not adapted to the exclusive normalisation strategy

    tree_signal = ROOT.TChain('signal_tree') 
    tree_run = ROOT.TChain('run_tree') 
    for signal_file in signal_files:
      if not is_bc:
        tree_run.Add(signal_file.filename)
        tree_signal.Add(signal_file.filename)
      else:
        tree_run.Add(signal_file.filename_Bc)
        tree_signal.Add(signal_file.filename_Bc)

    # number of selected events
    hist = ROOT.TH1D('hist', 'hist', 100, 0, 10)
    tree_signal.Project('hist', 'hnl_mass', selection)
    n_selected = hist.Integral()

    # number of generated events (= n_gen / filter_efficiency = n_miniaod / filter_efficiency)
    n_gen_tot = self.getNminiAODEvts(tree_run)
    # take the weighted average of the filter efficiencies of the generated samples
    filter_efficiency = 0.
    for signal_file in signal_files:
      if not is_bc:
        the_file = self.getRootFile(signal_file.filename)
        the_filter_efficiency = signal_file.filter_efficiency
      else:
        the_file = self.getRootFile(signal_file.filename_Bc)
        the_filter_efficiency = signal_file.filter_efficiency_Bc
      the_tree_run = self.getTree(the_file, 'run_tree')
      n_gen = self.getNminiAODEvts(the_tree_run)
      filter_efficiency += n_gen * the_filter_efficiency
    filter_efficiency = filter_efficiency / n_gen_tot

    # in the case where the samples were produced with both the electron and muon channels, apply a correction
    corr = signal_file.muon_rate # only consider events that were generated in the muon channel 
    #TODO fix the above if a different one is taken for Bc

    efficiency = filter_efficiency if not is_bc else filter_efficiency * lhe_efficiency 
    n_generated = corr * n_gen_tot / efficiency

    signal_efficiency = n_selected / n_generated 

    return signal_efficiency


  def getCtauWeight(self, signal_files, ctau, strategy='new', is_bu=False, is_bd=False, is_bs=False, is_bc=False):
    if strategy == 'new':
      # get the total number of gen (=miniaod) events
      n_miniaod_tot = 0
      if is_bu:
        for signal_file in signal_files:
          n_miniaod_tot += signal_file.n_miniaod_Bu
      elif is_bd:
        for signal_file in signal_files:
          n_miniaod_tot += signal_file.n_miniaod_Bd
      elif is_bs:
        for signal_file in signal_files:
          n_miniaod_tot += signal_file.n_miniaod_Bs
      elif is_bc:
        tree_run = ROOT.TChain('run_tree') 
        for signal_file in signal_files:
          tree_run.Add(signal_file.filename_Bc)
        n_miniaod_tot = self.getNminiAODEvts(tree_run)
      else: # previous normalisation strategy
        tree_run = ROOT.TChain('run_tree') 
        for signal_file in signal_files:
          tree_run.Add(signal_file.filename)
        n_miniaod_tot = self.getNminiAODEvts(tree_run)

      filter_efficiency_avg = 0.
      for signal_file in signal_files:
        if is_bu:
          the_filter_efficiency = signal_file.filter_efficiency_Bu
          n_gen = signal_file.n_miniaod_Bu
        elif is_bd:
          the_filter_efficiency = signal_file.filter_efficiency_Bd
          n_gen = signal_file.n_miniaod_Bd
        elif is_bs:
          the_filter_efficiency = signal_file.filter_efficiency_Bs
          n_gen = signal_file.n_miniaod_Bs
        elif is_bc:
          the_filter_efficiency = signal_file.filter_efficiency_Bc
          the_file = self.getRootFile(signal_file.filename_Bc)
          the_tree_run = self.getTree(the_file, 'run_tree')
          n_gen = self.getNminiAODEvts(the_tree_run)
        else: # previous normalisation strategy
          the_filter_efficiency = signal_file.filter_efficiency
          the_file = self.getRootFile(signal_file.filename)
          the_tree_run = self.getTree(the_file, 'run_tree')
          n_gen = self.getNminiAODEvts(the_tree_run)

        filter_efficiency_avg += n_gen * the_filter_efficiency
      filter_efficiency_avg = filter_efficiency_avg / n_miniaod_tot

      deno_weight = ''
      for ifile, signal_file in enumerate(signal_files):
        if is_bu:
          the_filter_efficiency = signal_file.filter_efficiency_Bu
          n_miniaod = signal_file.n_miniaod_Bu
        elif is_bd:
          the_filter_efficiency = signal_file.filter_efficiency_Bd
          n_miniaod = signal_file.n_miniaod_Bd
        elif is_bs:
          the_filter_efficiency = signal_file.filter_efficiency_Bs
          n_miniaod = signal_file.n_miniaod_Bs
        elif is_bc:
          the_filter_efficiency = signal_file.filter_efficiency_Bc
          the_file = self.getRootFile(signal_file.filename_Bc)
          the_tree_run = self.getTree(the_file, 'run_tree')
          n_miniaod = self.getNminiAODEvts(the_tree_run)
        else: # previous normalisation strategy
          the_filter_efficiency = signal_file.filter_efficiency
          the_file = self.getRootFile(signal_file.filename)
          the_tree_run = self.getTree(the_file, 'run_tree')
          n_miniaod = self.getNminiAODEvts(the_tree_run)

        if ifile == 0:
          deno_weight += ' {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod / the_filter_efficiency,
                ctau0 = signal_file.ctau,
                )
        else:
          deno_weight += ' + {n0} / {ctau0} * exp(-gen_hnl_ct / {ctau0})'.format(
                n0 = n_miniaod / the_filter_efficiency,
                ctau0 = signal_file.ctau,
                )
      #print 'ctau weight: ntot = {} / {} = {}'.format(n_miniaod_tot, filter_efficiency_avg, n_miniaod_tot / filter_efficiency_avg)
      weight_ctau = '({ntot} / {ctau1} * exp(-gen_hnl_ct / {ctau1}) * (1. / ({deno_weight})))'.format(
          ntot = n_miniaod_tot / filter_efficiency_avg,
          ctau1 = ctau,
          deno_weight = deno_weight,
          )

    elif strategy == 'old':
      signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
      if signal_file.is_private:
        original_ctau = signal_filename[signal_filename.find('ctau')+4:signal_filename.find('/', signal_filename.find('ctau')+1)]
      else:
        original_ctau = signal_filename[signal_filename.find('ctau')+4:signal_filename.find('mm', signal_filename.find('ctau')+1)]
      original_ctau = original_ctau.replace('p', '.')
      target_ctau = signal_file.ctau

      weight_ctau = 1. #-99
      if float(original_ctau) != float(target_ctau):
        ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

    return weight_ctau


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


  def printLatexBox(self, x, y, text, size=0.04, pos='center', font=42):
    box = ROOT.TLatex()
    box.SetNDC()
    box.SetTextFont(font)
    if pos == 'center': box.SetTextAlign(22) 
    else: box.SetTextAlign(12) 
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


  def printInnerCMSTag(self, pad, cms_tag, print_tag=False, x_pos=0.15, y_pos=0.83, size=0.55):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(x_pos, y_pos, 'CMS')
    ## print CMS tag
    tag.SetTextFont(52)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    x_pos_tag = x_pos
    y_pos_tag = y_pos - 0.06
    tag.DrawLatex(x_pos_tag, y_pos_tag, cms_tag)      
    pad.Update()


  def printCMSTagInFrame(self, pad, cms_tag, size=0.55, offset=0.11):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin()+0.2*pad.GetLeftMargin(), 1-pad.GetTopMargin()-0.8*pad.GetTopMargin(), 'CMS')
    # print CMS tag
    tag.SetTextFont(52)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    tag.DrawLatex(pad.GetLeftMargin()+0.2*pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()-0.8*pad.GetTopMargin(), cms_tag)      
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


  def gamma_total(self, mass, vv=1., fe=0., fu=1., ft=0.):
    '''
    Total width for N (Dirac)
    '''
    from decays import HNLDecays
    gamma_total = HNLDecays(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft).decay_rate['tot']  # GeV
    #gamma_total_e = 2*HNLDecays(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft).decay_rate['tot_el'] / (float(fe) * vv) if float(fe) != 0. else 0. # GeV
    #gamma_total_u = 2*HNLDecays(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft).decay_rate['tot_mu'] / (float(fu) * vv) if float(fu) != 0. else 0. # GeV
    #gamma_total_t = 2*HNLDecays(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft).decay_rate['tot_tau'] / (float(ft) * vv) if float(ft) != 0. else 0. # GeV
    #if fu == 0.5:
    #  print 'vv * (fe * gamma_total_e + fu * gamma_total_u + ft * gamma_total_t) = {} * ({} * {} + {} * {} + {} * {}) = {}'.format(vv, fe, gamma_total_e, fu, gamma_total_u, ft, gamma_total_t, 2*gamma_total)
    #  print 'gamm_tilde tot = {}'.format(2*gamma_total  / vv)
    #print gamma_total - vv * (fe * gamma_total_e + fu * gamma_total_u + ft * gamma_total_t)
    return gamma_total


  def gamma_partial(self, mass, vv=1., fe=0., fu=1., ft=0.):
    '''
    Partial width for N->mupi (Dirac)
    '''
    from decays import HNLDecays
    gamma_partial = HNLDecays(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft).decay_rate['mupi'] # GeV
    return gamma_partial


  def getVV(self, mass=-99.,ctau=-99., fe=0., fu=1., ft=0., ismaj=True):
    '''
    Helper function to go from ctau,m -> vv
    '''
    import numpy as np
    mult = 2. if ismaj else 1.
    ref_m = 1. # GeV
    ref_vv = 1. 
    #print '1/ct = maj * gamma_tot = {} * {}'.format(mult, self.gamma_total(mass=ref_m, vv=ref_vv, fe=fe, fu=fu, ft=ft))
    ref_ctau = self.ctau_from_gamma(mult*self.gamma_total(mass=ref_m, vv=ref_vv, fe=fe, fu=fu, ft=ft)) 
    k = ref_ctau * np.power(ref_m,5) * ref_vv
    return k/(np.power(mass, 5) * ctau)


  def getCtau(self, mass=-99, vv=-99, fe=0., fu=1., ft=0., ismaj=True):
    '''
    Helper function to go from vv,m(GeV) -> ctau (mm)
    '''
    mult = 2. if ismaj else 1.
    #if fu == 0.5:
    #  print '1/ct = mult * gamma_tot = {} * {}'.format(mult, self.gamma_total(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft))
    return self.ctau_from_gamma(mult*self.gamma_total(mass=mass, vv=vv, fe=fe, fu=fu, ft=ft))


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)


if __name__ == '__main__':
  tools = Tools()
  #print tools.getCtau(2.5, 1e-1)
  #print tools.getCtau(1.0, 6e-5, fe=0., fu=1., ft=0., ismaj=True)
  print tools.getVV(3., 1)
  #mass = 6.
  #print tools.gamma_partial(mass=mass, vv=1.) / tools.gamma_total(mass=mass, vv=1.)
  #print tools.gamma_partial(mass=mass, vv=1., fe=0., fu=1., ft=0.) / tools.gamma_total(mass=mass, vv=1., fe=0., fu=1., ft=0.) * 100.
  


