import os
import sys
from os import path
from time import time 

import math
from itertools import product

import ROOT
from root_pandas import read_root

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import array

from keras.models import load_model
from sklearn.metrics import roc_curve

sys.path.append('../')
sys.path.append('../../objects')
from tools import Tools
from mva_tools import MVATools

from samples import signal_samples, data_samples
from categories import categories
from baseline_selection import selection
from resolutions import resolutions
from quantity import Quantity as Qte


class TrainingInfo(object):
  '''
    Class that contains the training information:
      - model
      - scaler
      - features
  '''
  def __init__(self, dirname):
    self.dirname = dirname 
    self.indir = './outputs/{}'.format(self.dirname)
    self.model = self.loadModel()
    self.qt = self.loadScaler()
    self.features = self.loadFeatures()

  def loadModel(self):
    print '\n --> get the model'
    model_filename = '{}/net_model_weighted.h5'.format(self.indir)
    model = load_model(model_filename)
    return model


  def loadScaler(self):
    print '\n --> get the scaler'
    scaler_filename = '/'.join([self.indir, 'input_tranformation_weighted.pck'])
    qt = pickle.load(open(scaler_filename, 'rb'))
    return qt


  def loadFeatures(self):
    print '\n --> get the features'
    features_filename = '/'.join([self.indir, 'input_features.pck'])
    features = pickle.load(open(features_filename, 'rb'))
    return features



class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''
  def __init__(self, filename, selection, training_info, mass=None, ctau=None, colour=None, filter_efficiency=None, muon_rate=None):
    self.filename = filename
    self.selection = selection
    cutbased_selection_qte = ['hnl_pt', 'hnl_charge', 'sv_lxysig', 'hnl_cos2d'] 
    self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True, columns=training_info.features+cutbased_selection_qte)
    self.mass = mass
    self.ctau = ctau
    self.colour = colour
    self.filter_efficiency = filter_efficiency
    self.muon_rate = muon_rate

    print '\t selected events: {}'.format(len(self.df))



class MVAAnalyser(Tools, MVATools):
  def __init__(self, signal_files, data_files, dirname, baseline_selection, categories=None, plot_label=None, do_plotScore=False, do_createFiles=False, do_plotSigScan=False, do_plotROC=False, do_plotDistributions=False):
    self.tools = Tools()
    self.mva_tools = MVATools()
    self.signal_files = signal_files
    self.data_files = data_files
    self.dirname = dirname
    self.plot_label = plot_label
    self.baseline_selection = baseline_selection
    self.categories = categories
    self.do_plotScore = do_plotScore
    self.do_createFiles = do_createFiles
    self.do_plotSigScan = do_plotSigScan
    self.do_plotROC = do_plotROC
    self.do_plotDistributions = do_plotDistributions

    self.outdir = self.createOutDir()

    self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2'
    self.weight_pusig = 'weight_pu_sig_D'
    self.weight_mu0id = 'weight_mu0_softid'
    self.weight_muid = 'weight_mu_looseid'


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = './outputs/{}/analysis'.format(self.dirname)
    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))

    return outdir


  def getPandasQuery(self, selection):
    '''
      Converts selection to pandas query syntax
    '''
    query = selection.replace('&&', ' and ').replace('||', ' or ').replace('!=', ' not ').replace('!', ' not ')
    
    return query


  def getSamples(self, training_info, category):
    '''
      Function that fetches the samples into lists
    '''
    #TODO in case of more complex baseline selection, convert it to a format that is digestable by pandas query

    selection = self.baseline_selection + ' && ' + category.definition_flat
    #selection = 'hnl_charge==0'
    #query = self.getPandasQuery(selection)

    print('========> starting reading the trees')
    now = time()
    ## background 
    data_samples = []
    for data_file in self.data_files:
      #data_samples.append(Sample(filename=data_file.filename, selection=selection))
      data_samples.append(Sample(filename=data_file, selection=selection, training_info=training_info))

    ## signal
    mc_samples = []
    for signal_file in self.signal_files:
      mc_samples.append(Sample(filename=signal_file.filename, selection=selection, training_info=training_info, mass=signal_file.mass, ctau=signal_file.ctau, colour=signal_file.colour, filter_efficiency=signal_file.filter_efficiency, muon_rate=signal_file.muon_rate))
      
    print('========> it took %.2f seconds' %(time() - now))

    return data_samples, mc_samples


  def createDataframe(self, samples):
    '''
      Function that converts root samples to pandas dataframe
    '''
    df = pd.concat([idt.df for idt in samples], sort=False)

    return df


  def saveFig(self, plt, name):
    '''
      Save python figure
    '''
    plt.savefig('{}/{}.pdf'.format(self.outdir, name))    
    plt.savefig('{}/{}.png'.format(self.outdir, name))    
    print ' --> {}/{}.png created'.format(self.outdir, name)


  def predictScore(self, training_info, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=training_info.features)

    # apply the scaler
    xx = training_info.qt.transform(x[training_info.features])

    # predict
    score = training_info.model.predict(xx)

    return score


  def plotScore(self, training_info, mc_samples, data_samples, category, do_log):
    '''
      Plot the score distributions for signal and background
    '''
    canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
    if do_log:
      canv.SetLogy()
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)

    # background
    data_df = self.createDataframe(data_samples)
    bkg_score = self.predictScore(training_info, data_df)

    hist_bkg = ROOT.TH1F('bkg', 'bkg', 30, 0, 1)
    for score in bkg_score:
      hist_bkg.Fill(score)

    if hist_bkg.Integral()!=0: hist_bkg.Scale(1./hist_bkg.Integral())

    leg.AddEntry(hist_bkg, 'data-driven background')

    ROOT.gStyle.SetOptStat(0)
    hist_bkg.SetTitle(category.title)
    hist_bkg.SetFillColor(ROOT.kAzure-4)
    hist_bkg.SetLineColor(1)
    hist_bkg.GetXaxis().SetTitle('Score')
    hist_bkg.GetXaxis().SetLabelSize(0.033)
    hist_bkg.GetXaxis().SetTitleSize(0.042)
    hist_bkg.GetXaxis().SetTitleOffset(1.1)
    hist_bkg.GetYaxis().SetTitle('Normalised to unity')
    hist_bkg.GetYaxis().SetLabelSize(0.033) 
    hist_bkg.GetYaxis().SetTitleSize(0.042)
    hist_bkg.GetYaxis().SetTitleOffset(1.1)

    # signals
    hist_sigs = []
    for mc_sample in mc_samples:
      mc_df = self.createDataframe([mc_sample])
      sig_score = self.predictScore(training_info, mc_df)

      hist_name = 'sig_{}_{}'.format(mc_sample.mass, mc_sample.ctau)
      hist_sig = ROOT.TH1F(hist_name, hist_name, 30, 0, 1)
      hist_sig.SetDirectory(0)
      for score in sig_score:
        hist_sig.Fill(score)

      if hist_sig.Integral()!=0: hist_sig.Scale(1./hist_sig.Integral())

      leg.AddEntry(hist_sig, 'signal - mass {} GeV, ctau {} mm'.format(mc_sample.mass, mc_sample.ctau))

      hist_sig.SetLineWidth(3)
      hist_sig.SetLineColor(mc_sample.colour)
      hist_sigs.append(hist_sig)

    hist_bkg.Draw('histo')
    for hist_sig in hist_sigs:
      hist_sig.Draw('histo_same')
    leg.Draw('same')

    canv.cd()
    name = '{}/score_{}_{}'.format(self.outdir, category.label, self.plot_label) if self.plot_label != None else '{}/score_{}.png'.format(self.outdir, category.label)
    if do_log: name += '_log'
    canv.SaveAs(name + '.png')


  def plotROCCurve(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None
    # create dataframe
    data_df = self.createDataframe(data_samples)
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)


      mc_df = self.createDataframe([mc_sample])
      mc_df['is_signal'] = 1

      main_df = pd.concat([data_df, mc_df], sort=False)
      main_df.index = np.array(range(len(main_df)))
      main_df = main_df.sample(frac=1, replace=False, random_state=1986)

      Y = pd.DataFrame(main_df, columns=['is_signal'])
      score = self.predictScore(training_info=training_info, df=main_df)
      fpr, tpr, thresholds = roc_curve(Y, score) 

      plt.plot(fpr, tpr, linewidth=2, label='mva - ({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      plt.xlabel('False Positive Rate')
      plt.ylabel('True Positive Rate')

      # get cutbased performance
      if mc_sample.mass < 3.: cutbased_selection = category.cutbased_selection_lowmass
      elif mc_sample.mass >= 3. and mc_sample.mass < 4.5: cutbased_selection = category.cutbased_selection_mediummass
      elif mc_sample.mass >= 4.5: cutbased_selection = category.cutbased_selection_highmass
      if len(cutbased_selection) == 0:
        cutbased_selection = category.cutbased_selection

      main_df['cutbased_score'] = 0
      
      main_df_tmp1 = main_df.query(self.getPandasQuery(cutbased_selection))
      main_df_tmp1['cutbased_score'] = 1

      main_df_tmp2 = main_df.query(self.getPandasQuery('!({})'.format(cutbased_selection)))
      main_df_tmp2['cutbased_score'] = 0

      main_df_cutbased = pd.concat([main_df_tmp1, main_df_tmp2], sort=False)

      TP = len(main_df_cutbased.query('cutbased_score==1 & is_signal==1'))
      FN = len(main_df_cutbased.query('cutbased_score==0 & is_signal==1'))
      FP = len(main_df_cutbased.query('cutbased_score==1 & is_signal==0'))
      TN = len(main_df_cutbased.query('cutbased_score==0 & is_signal==0'))

      true_positive_rate = float(TP) / float(TP + FN)
      false_positive_rate = float(FP) / float(FP + TN)
      plt.plot(false_positive_rate, true_positive_rate, 'o', markersize=10, label='cutbased - ({}GeV, {}mm, {})'.format(mass, ctau, coupling))

      xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
      plt.plot(xy, xy, color='grey', linestyle='--')

      plt.title(r'{}'.format(category.title))
      plt.xlim(0, 1)
      plt.ylim(0, 1)

      if do_log:
        plt.xscale('log')
      plt.yscale('linear')

      if do_log:
        plt.legend(loc='upper left', framealpha=0.1)
      else:
        plt.legend(loc='lower right', framealpha=0.1)

    name = 'ROC_m{}_{}'.format(str(mc_samples[0].mass).replace('.', 'p'), category.label)
    if do_log: name += '_log'
    self.saveFig(plt, name)


  def plotScoreCurve(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None
    # create dataframe
    data_df = self.createDataframe(data_samples)
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      mc_df['is_signal'] = 1

      main_df = pd.concat([data_df, mc_df], sort=False)
      main_df.index = np.array(range(len(main_df)))
      main_df = main_df.sample(frac=1, replace=False, random_state=1986)

      Y = pd.DataFrame(main_df, columns=['is_signal'])
      score = self.predictScore(training_info=training_info, df=main_df)
      fpr, tpr, thresholds = roc_curve(Y, score) 

      plt.plot(thresholds, tpr-fpr, linewidth=2, label='mva - ({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      #plt.plot(thresholds, tpr / np.sqrt(fpr), linewidth=2, label='mva - ({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      plt.xlabel('Score')
      plt.ylabel('True Positive Rate - False Positive Rate')

      # get optimal score
      optimal_idx = np.argmax(tpr - fpr)
      optimal_score = thresholds[optimal_idx]
      optimal_performance = tpr[optimal_idx] - fpr[optimal_idx]
      plt.plot(optimal_score, optimal_performance, '*', markersize=10, label='Optimal - {}'.format(round(optimal_score, 2)))

      plt.title(r'{}'.format(category.title))
      plt.xlim(0, 1)
      plt.ylim(0, 1)

      if do_log:
        plt.xscale('log')
      plt.yscale('linear')

      if do_log:
        plt.legend(loc='upper left', framealpha=0.1)
      else:
        plt.legend(loc='lower right', framealpha=0.1)

    name = 'score_curve_m{}_{}'.format(str(mc_samples[0].mass).replace('.', 'p'), category.label)
    if do_log: name += '_log'
    self.saveFig(plt, name)


  def plotMVAPerformance(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None
    # create dataframe
    data_df = self.createDataframe(data_samples)
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      mc_df['is_signal'] = 1

      main_df = pd.concat([data_df, mc_df], sort=False)
      main_df.index = np.array(range(len(main_df)))
      main_df = main_df.sample(frac=1, replace=False, random_state=1986)

      # get mva performance
      Y = pd.DataFrame(main_df, columns=['is_signal'])
      score = self.predictScore(training_info=training_info, df=main_df)
      false_positive_rate_mva, true_positive_rate_mva, thresholds = roc_curve(Y, score) 

      # get cutbased performance
      if mc_sample.mass < 3.: cutbased_selection = category.cutbased_selection_lowmass
      elif mc_sample.mass >= 3. and mc_sample.mass < 4.5: cutbased_selection = category.cutbased_selection_mediummass
      elif mc_sample.mass >= 4.5: cutbased_selection = category.cutbased_selection_highmass
      if len(cutbased_selection) == 0:
        cutbased_selection = category.cutbased_selection

      main_df['cutbased_score'] = 0
      
      main_df_tmp1 = main_df.query(self.getPandasQuery(cutbased_selection))
      main_df_tmp1['cutbased_score'] = 1

      main_df_tmp2 = main_df.query(self.getPandasQuery('!({})'.format(cutbased_selection)))
      main_df_tmp2['cutbased_score'] = 0

      main_df_cutbased = pd.concat([main_df_tmp1, main_df_tmp2], sort=False)

      TP = len(main_df_cutbased.query('cutbased_score==1 & is_signal==1'))
      FN = len(main_df_cutbased.query('cutbased_score==0 & is_signal==1'))
      FP = len(main_df_cutbased.query('cutbased_score==1 & is_signal==0'))
      TN = len(main_df_cutbased.query('cutbased_score==0 & is_signal==0'))

      true_positive_rate_cutbased_val = float(TP) / float(TP + FN)
      false_positive_rate_cutbased_val = float(FP) / float(FP + TN)
      #TODO fill object of same dimension as those returned by roc

      true_positive_rate_cutbased = np.full(true_positive_rate_mva.shape, true_positive_rate_cutbased_val)
      false_positive_rate_cutbased = np.full(false_positive_rate_mva.shape, false_positive_rate_cutbased_val)

      relative_signal_efficiency = true_positive_rate_mva / true_positive_rate_cutbased
      relative_background_efficiency = false_positive_rate_mva / false_positive_rate_cutbased
      relative_significance = relative_signal_efficiency / np.sqrt(relative_background_efficiency)

      plt.plot(thresholds, relative_significance, linewidth=2, label='({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      plt.xlabel('Score')
      plt.ylabel('tpr(mva)/tpr(cut) / sqrt(fpr(mva)/fpr(cut))')

      ## get optimal score
      #optimal_idx = np.argmax(tpr - fpr)
      #optimal_score = score[optimal_idx]
      #optimal_performance = tpr[optimal_idx] - fpr[optimal_idx]
      #plt.plot(optimal_score, optimal_performance, '*', markersize=10, label='Optimal - {}'.format(round(optimal_score, 2)))

      plt.axhline(y = 1., color='black', linestyle = '--')
      plt.title(r'{}'.format(category.title))
      plt.xlim(0, 1)
      #plt.ylim(0, 1)

      if do_log:
        plt.xscale('log')
      plt.yscale('linear')

      if do_log:
        plt.legend(loc='upper left', framealpha=0.1)
      else:
        plt.legend(loc='lower right', framealpha=0.1)

    name = 'mva_performance_m{}_{}'.format(str(mc_samples[0].mass).replace('.', 'p'), category.label)
    if do_log: name += '_log'
    self.saveFig(plt, name)


  def computeSignalYields(self, tree_sig, mc_sample, selection, mass='', ctau='', lumi=0.774, sigma_B=472.8e9, add_weight_hlt=False, add_weight_pu=False, add_weight_muid=False, weight_hlt='weight_hlt_D', weight_pusig='weight_pu_sig_D', weight_mu0id='weight_mu0_softid', weight_muid='weight_mu_looseid', isBc=False):
    '''
      signal yields computed as sigma_HNL * lumi * efficiency
    '''

    # define ranges and binning (full window)
    fit_window_min = mass - 1.5
    fit_window_max = mass + 1.5 

    # define selection
    cond_sig = 'ismatched==1'
    selection_sig = cond_sig + ' && ' + selection

    # define signal weights
    weight_ctau = self.tools.getCtauWeight(signal_files=[mc_sample], ctau=ctau)
    lhe_efficiency = 0.08244 #FIXME
    weight_signal = self.tools.getSignalWeight(signal_files=[mc_sample], mass=mass, ctau=ctau, sigma_B=sigma_B, lumi=lumi, lhe_efficiency=lhe_efficiency)
    weight_sig = '({}) * ({})'.format(weight_signal, weight_ctau)
    if add_weight_hlt: weight_sig += ' * ({})'.format(weight_hlt)
    if add_weight_pu: weight_sig += ' * ({})'.format(weight_pusig)
    if add_weight_muid: weight_sig += ' * ({}) * ({})'.format(weight_mu0id, weight_muid)
    #print 'weight ',weight_sig

    # create histogram
    hist_name = 'hist_signal_{}'.format(isBc)
    hist = ROOT.TH1D(hist_name, hist_name, 100, fit_window_min, fit_window_max)
    branch_name = 'hnl_mass'
    tree_sig.Project(hist_name, branch_name , '({sel}) * ({wght})'.format(sel=selection_sig, wght=weight_sig))

    # get the number of yields
    n_sig = hist.Integral()

    return n_sig


  def plotSignificanceDiffScan(self, category, mc_samples, data_samples, n_points=30):
    '''
      Plot the gain in significance for different cuts on the score
      Indicates the performance of the cutbased method
    '''
    # compute signal yields
    scan_points = np.linspace(0, 1, n_points) 
      
    canv = self.tools.createTCanvas('canv', 900, 800)

    #ROOT.gStyle.SetPadLeftMargin(0.8)
    #ROOT.gStyle.SetPadRightMargin(0.8)
    
    gs_sig = []

    significance_diff_max = -1
    significance_diff_min = -1

    for mc_sample in mc_samples:
      # note: we don't use the lifetime reweigthing

      # get signal features
      signal_mass = mc_sample.mass
      signal_resolution = resolutions[signal_mass]
      signal_ctau = mc_sample.ctau
      #if signal_ctau != 1000.: continue

      # get signal tree, with score branch, without selection
      sig_label = 'sig_{}_{}_{}'.format(signal_mass, signal_ctau, category.label)
      tree_sig = self.mva_tools.getTreeWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, label=sig_label, treename='signal_tree', force_overwrite=True) 

      bkg_label = 'bkg_{}_{}'.format(signal_mass, category.label)
      tree_data = self.mva_tools.getTreeWithScore(files=data_samples, training_label='./outputs/'+self.dirname, label=bkg_label, treename='signal_tree', force_overwrite=True) 

      g_sig = ROOT.TGraph()
      quantity = Qte(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=signal_mass-2*signal_resolution, bin_max=signal_mass+2*signal_resolution)

      lumi_D1 = 5.302 # when computing the relative difference in significance, this value does not matter

      # compute the initial significance
      signal_selection_ini = 'ismatched == 1 && {} && {}'.format(self.baseline_selection, category.definition_flat) 
      signal_yields_ini = self.computeSignalYields(tree_sig=tree_sig, mc_sample=mc_sample, selection=signal_selection_ini, mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid) 

      background_selection_ini = '{} && {}'.format(self.baseline_selection, category.definition_flat)
      hist_data_name = 'hist_data_ini_{}_{}'.format(signal_mass, category.label)
      hist_data = self.tools.createHisto(tree_data, quantity, hist_name=hist_data_name, branchname='flat', selection=background_selection_ini)
      background_yields_ini = hist_data.Integral()

      if background_yields_ini != 0:
        significance_ini = signal_yields_ini / math.sqrt(background_yields_ini)
      else:
        significance_ini = -99.

      for idx, cut in enumerate(scan_points):
        print 'score > {}'.format(cut)
        # compute signal yields
        signal_selection_sel = 'ismatched == 1 && {} && {} && score>{}'.format(self.baseline_selection, category.definition_flat, cut) 
        signal_yields_sel = self.computeSignalYields(tree_sig=tree_sig, mc_sample=mc_sample, selection=signal_selection_sel, mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid) 

        # get background yields from unblinded data
        background_selection_sel = '{} && {} && score>{}'.format(self.baseline_selection, category.definition_flat, cut)
        hist_data_name = 'hist_data_sel_{}_{}'.format(signal_mass, category.label)
        hist_data = self.tools.createHisto(tree_data, quantity, hist_name=hist_data_name, branchname='flat', selection=background_selection_sel)
        background_yields_sel = hist_data.Integral()

        # compute significance
        if background_yields_sel != 0:
          significance_sel = signal_yields_sel / math.sqrt(background_yields_sel)
        else:
          significance_sel = -99.

        # compute relative difference in significance
        if significance_ini != -99. and significance_sel != -99.:
          significance_diff = ((significance_sel / significance_ini) - 1) * 100
        else:
          significance_diff = -1

        if significance_diff_max == -1 or significance_diff > significance_diff_max:
          significance_diff_max = significance_diff

        if significance_diff_min == -1 or significance_diff < significance_diff_min:
          significance_diff_min = significance_diff

        if significance_diff_min < -50: significance_diff_min = -50

        g_sig.SetPoint(idx, cut, significance_diff)
        g_sig.SetLineWidth(2)
        g_sig.SetLineStyle(9)
        g_sig.SetLineColor(mc_sample.colour)
        g_sig.SetMarkerColor(mc_sample.colour)
        g_sig.SetMarkerStyle(20)

      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on MVA score')
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetTitle('Relative difference in significance [%]')
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetTitleSize(0.045)
      gs_sig[0].GetYaxis().SetTitleOffset(1.15)
      gs_sig[0].GetYaxis().SetRangeUser(significance_diff_min-0.15*significance_diff_min, significance_diff_max+0.15*significance_diff_max)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('APL')
        else:
          g_sig.Draw('PL')

      # get cutbased performance
      if signal_mass < 3: cutbased_selection = category.cutbased_selection_lowmass
      elif signal_mass >= 3 and signal_mass < 4.5: cutbased_selection = category.cutbased_selection_mediummass
      elif signal_mass >= 4.5: cutbased_selection = category.cutbased_selection_highmass
      if len(cutbased_selection) == 0:
        cutbased_selection = category.cutbased_selection

      signal_selection_cutbased = 'ismatched == 1 && {} && {} && {}'.format(self.baseline_selection, category.definition_flat, cutbased_selection) 
      signal_yields_cutbased = self.computeSignalYields(tree_sig=tree_sig, mc_sample=mc_sample, selection=signal_selection_cutbased, mass=signal_mass, ctau=signal_ctau, lumi=lumi_D1, sigma_B=472.8e9, isBc=False, add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt=self.weight_hlt, weight_pusig=self.weight_pusig, weight_mu0id=self.weight_mu0id, weight_muid=self.weight_muid) 

      background_selection_cutbased = '{} && {} && {}'.format(self.baseline_selection, category.definition_flat, cutbased_selection)
      hist_data_name = 'hist_data_cutbased_{}_{}'.format(signal_mass, category.label)
      hist_data = self.tools.createHisto(tree_data, quantity, hist_name=hist_data_name, branchname='flat', selection=background_selection_cutbased)
      background_yields_cutbased = hist_data.Integral()

      if background_yields_cutbased != 0:
        significance_cutbased = signal_yields_cutbased / math.sqrt(background_yields_cutbased) 
        significance_cutbased_diff = ((significance_cutbased / significance_ini) - 1) * 100

      if background_yields_cutbased != 0:
        line = ROOT.TLine(0, significance_cutbased_diff, 1, significance_cutbased_diff)
        line.SetLineColor(mc_sample.colour)
        line.SetLineStyle(9)
        line.SetLineWidth(3)
        line.Draw('same')

    self.tools.printLatexBox(0.67, 0.25, category.title, size=0.041)
          
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.25, xmax=0.82, ymax=0.4, size=0.04)
    for ifile, mc_sample in enumerate(mc_samples):
      legend.AddEntry(gs_sig[ifile], '{}GeV, {}mm'.format(mc_sample.mass, mc_sample.ctau))
    legend.Draw()

    canv.cd()
    canv.SaveAs('{}/significance_diff_scan_{}_{}.png'.format(self.outdir, category.label, self.plot_label))


  def plotDistributions(self, quantity, sample, category, cut_score, do_shape=False):
    '''
      Plot the distribution shapes with and without cut on score
    '''

    canv = self.tools.createTCanvas('canv', 800, 700) 
    label = 'tree'

    selection_noscore = self.baseline_selection + ' && ' + category.definition_flat
    #selection_score = selection_noscore + ' && score > {}'.format(cut_score)

    tree = self.mva_tools.getTreeWithScore(files=[sample], training_label='./outputs/'+self.dirname, selection=selection_noscore, label=label, treename='signal_tree', force_overwrite=True) 
    #branches = tree.GetListOfBranches()
    #for i in range(0, len(branches)):
    #  print branches.At(i).GetName()

    hist1 = self.tools.createHisto(tree, quantity, hist_name='hist1', branchname='flat', selection='score > -99.')
    hist2 = self.tools.createHisto(tree, quantity, hist_name='hist2', branchname='flat', selection='score > {}'.format(cut_score))

    if do_shape:
      if hist1.Integral() != 0: hist1.Scale(1./hist1.Integral())
      if hist2.Integral() != 0: hist2.Scale(1./hist2.Integral())

    hist1.SetLineColor(ROOT.kBlue)
    hist1.SetLineWidth(3)
    hist2.SetLineColor(ROOT.kMagenta)
    hist2.SetLineWidth(3)
    
    hist1.SetTitle(category.title)
    hist1.GetXaxis().SetTitle(quantity.title)
    hist1.GetXaxis().SetLabelSize(0.033)
    hist1.GetXaxis().SetTitleSize(0.042)
    hist1.GetXaxis().SetTitleOffset(1.1)
    hist1.GetYaxis().SetTitle('Normalised to unity')
    hist1.GetYaxis().SetLabelSize(0.033)
    hist1.GetYaxis().SetTitleSize(0.042)
    hist1.GetYaxis().SetTitleOffset(1.3)

    hist1.Draw('histo')
    hist2.Draw('histo same')
  
    leg = self.tools.getRootTLegend(xmin=0.65, ymin=0.7, xmax=0.85, ymax=0.9, size=0.03)
    leg.AddEntry(hist1, 'No cut on score')
    leg.AddEntry(hist2, 'Score > {}'.format(cut_score))
    leg.Draw()

    ROOT.gStyle.SetOptStat(0)

    if sample.mass != None:
      name = '{}_m{}_ctau{}_{}'.format(quantity.label, sample.mass, sample.ctau, category.label)
    else:
      name = '{}_{}'.format(quantity.label, category.label)

    if not path.exists('{}/plots'.format(self.outdir)):
      os.system('mkdir -p {}/plots'.format(self.outdir))
    canv.SaveAs('{}/plots/{}.png'.format(self.outdir, name))

    os.system('rm tree.root')



  def createRootFile(self, training_info, samples, name):
    '''
      Add the score to the dataframe and save the root file with all the branches
    '''
    filenames = [] 
    for sample in samples:
      # get the score
      df = self.createDataframe([sample])
      score = self.predictScore(training_info, df) 

      # add the score
      df.insert(len(df.columns), 'score', score)

      # create outfile
      if sample.mass != None and sample.ctau != None:
        root_filename = '{}/{}_m{}_ctau{}.root'.format(self.outdir, name, str(sample.mass).replace('.', 'p'), str(sample.ctau).replace('.', 'p')) #NOTE for the analysis code, save it in scratch (i.e ./)
      else:
        root_filename = '{}/{}.root'.format(self.outdir, name) #NOTE for the analysis code, save it in scratch (i.e ./)
      df.to_root(root_filename, key='signal_tree', store_index=False)

      print ' --> {} created'.format(root_filename)

      filenames.append(root_filename)
    
    return filenames


  def createFriendTree(self, training_info, samples, name):
    '''
      Create tree that contains the score and store it in a root file
    '''
    for sample in samples:
      # get the score
      df = self.createDataframe([sample])
      score = self.predictScore(training_info, df) 

      # create file
      if sample.mass != None and sample.ctau != None:
        root_filename = '{}/{}_m{}_ctau{}.root'.format(self.outdir, name, str(sample.mass).replace('.', 'p'), str(sample.ctau).replace('.', 'p')) #NOTE for the analysis code, save it in scratch (i.e ./)
      else:
        root_filename = '{}/{}.root'.format(self.outdir, name) #NOTE for the analysis code, save it in scratch (i.e ./)
      friend_file = ROOT.TFile(root_filename, 'RECREATE')

      # create tree
      tree = ROOT.TTree('score_tree', 'score_tree')
      the_score = array('d',[0])
      tree.Branch('score', the_score, 'score/D')
      for entry in range(len(score)):
        the_score[0] = score[entry][0]
        #print score[entry][0]
        tree.Fill()
      #tree.Fill()
      tree.Write()
      friend_file.Close()

      print ' --> {} created'.format(root_filename)
    

  def compareTrees(self):
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root'
    friend_filename = './outputs/test_20Aug2022_13h40m03s/analysis/score_mc_m4p5_ctau100p0.root' 
    score_filename = './outputs/test_20Aug2022_13h40m03s/analysis/sample_mc_m4p5_ctau100p0.root'

    f = self.tools.getRootFile(filename)
    f_friend = self.tools.getRootFile(friend_filename)
    f_score = self.tools.getRootFile(score_filename)

    tree = self.tools.getTree(f, 'signal_tree')
    friend_tree = self.tools.getTree(f_friend, 'score_tree')
    score_tree = self.tools.getTree(f_score, 'signal_tree')

    tree.AddFriend(friend_tree)

    canv = self.tools.createTCanvas('canv', 800, 700) 
    hist1 = ROOT.TH1D('hist1', 'hist1', 30, 0, 1)
    hist1.SetLineColor(ROOT.kRed)
    #tree.Draw('score >> hist1')
    tree.Project('hist1', 'score', 'score>0.9')
    hist1.Draw()

    hist2 = ROOT.TH1D('hist2', 'hist2', 30, 0, 1)
    hist2.SetLineColor(ROOT.kBlue)
    score_tree.Draw('score >> hist2')

    hist2.Draw()
    hist1.Draw('same')

    name = '{}/score_comparison.png'.format(self.outdir)
    canv.SaveAs(name)


  def process(self):
    print '---- MVA Analyser ----'

    print '\n -> get the training information'
    training_info = TrainingInfo(self.dirname)

    for category in self.categories:
      print '\n -> get the samples'
      data_samples, mc_samples = self.getSamples(training_info=training_info, category=category)

      if self.do_plotScore:
        print '\n -> plot the score'
        self.plotScore(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)
        self.plotScore(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=True)

      if self.do_plotROC:
        if category.label == 'incl': continue
        self.plotROCCurve(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)
        self.plotROCCurve(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=True)
        self.plotScoreCurve(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)
        self.plotMVAPerformance(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)

      if self.do_createFiles:
        print '\n -> create rootfiles'
        self.createRootFile(training_info, mc_samples, 'sample_mc')
        self.createRootFile(training_info, data_samples, 'sample_data')
       
        self.createFriendTree(training_info, mc_samples, 'score_mc')
        self.createFriendTree(training_info, data_samples, 'score_data')

        #self.compareTrees()
  
      if self.do_plotSigScan:
        self.plotSignificanceDiffScan(category=category, mc_samples=mc_samples, data_samples=data_samples, n_points=11)

      hnl_mass = Qte(name_nano='BToMuMuPi_hnl_mass', name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=80, bin_min=0., bin_max=5.4)

      if self.do_plotDistributions:
        self.plotDistributions(quantity=hnl_mass, sample=data_samples[0], category=category, cut_score=0.9, do_shape=True) 
        #for mc_sample in mc_samples:
        #  self.plotDistributions(quantity=hnl_mass, sample=mc_sample, category=category, cut_score=0.9, do_shape=True) 



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  #dirname = 'test_20Aug2022_13h40m03s' # large ntuples
  #dirname = 'test_22Aug2022_16h08m33s'  # re-indexed
  #dirname = 'test_22Aug2022_16h17m44s' # not re-indexed
  #dirname = 'test_29Aug2022_11h52m47s' # trained on m1 ctau100
  #dirname = 'test_08Sep2022_21h04m04s' # more features, more stat
  #dirname = 'test_08Sep2022_21h56m30s' # larger learning rate
  #dirname = 'test_08Sep2022_23h11m59s' # trained on ctau100 and ctau1
  #dirname = 'test_08Sep2022_23h37m00s' # mass 1.5
  dirname = 'test_09Sep2022_00h12m12s' # mass 4.5

  #baseline_selection = 'hnl_charge==0'
  baseline_selection = selection['baseline_08Aug22'].flat + ' && hnl_charge==0'

  do_plotScore = True
  do_createFiles = False
  do_plotSigScan = False
  do_plotROC = True
  do_plotDistributions = False

  signal_labels = ['V12_08Aug22_m4p5']
  #signal_labels = ['V12_08Aug22_m1', 'V12_08Aug22_m1p5', 'V12_08Aug22_m2', 'V12_08Aug22_m3', 'V12_08Aug22_m4p5']
  plot_label = 'm4p5'

  signal_files = []
  for signal_label in signal_labels:
    for signal_file in signal_samples[signal_label]:
      signal_files.append(signal_file)

  data_files = []
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj1.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_sr.root')
  data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk1_n500/flat/flat_bparknano_08Aug22_sr.root')

  categories = categories['V12_08Aug22_permass']
  #categories = categories['inclusive']

  analyser = MVAAnalyser(
      signal_files = signal_files,
      data_files = data_files,
      dirname = dirname,
      plot_label = plot_label,
      baseline_selection = baseline_selection,
      categories = categories,
      do_plotScore = do_plotScore,
      do_createFiles = do_createFiles,
      do_plotSigScan = do_plotSigScan,
      do_plotROC = do_plotROC,
      do_plotDistributions = do_plotDistributions,
      )

  analyser.process()





