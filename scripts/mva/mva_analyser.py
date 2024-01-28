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
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from array import array

from keras.models import load_model
from sklearn.metrics import roc_curve
from sklearn.metrics import auc as sklearn_auc

sys.path.append('../')
sys.path.append('../../objects')
from tools import Tools
from mva_tools import MVATools

from samples import signal_samples, data_samples
from categories import categories
from baseline_selection import selection
from quantity import Quantity as Qte
from quantity import quantities
quantities_preselection = quantities['preselection']
quantities_trackid = quantities['trackId']
quantities_muonid = quantities['muonId']
quantities_small = quantities['small']
quantities_pNN = quantities['pNN_features']


class TrainingInfo(object):
  '''
    Class that contains the training information:
      - model
      - scaler
      - features
  '''
  def __init__(self, dirname, label):
    self.dirname = dirname 
    self.indir = './outputs/{}'.format(self.dirname)
    self.label = label
    self.model = self.loadModel()
    self.qt = self.loadScaler()
    self.features = self.loadFeatures()

  def loadModel(self):
    print '\n --> get the model'
    model_filename = '{}/net_model_weighted_{}.h5'.format(self.indir, self.label)
    model = load_model(model_filename)
    return model


  def loadScaler(self):
    print '\n --> get the scaler'
    scaler_filename = '/'.join([self.indir, 'input_tranformation_weighted_{}.pck'.format(self.label)])
    qt = pickle.load(open(scaler_filename, 'rb'))
    return qt


  def loadFeatures(self):
    print '\n --> get the features'
    features_filename = '/'.join([self.indir, 'input_features_{}.pck'.format(self.label)])
    features = pickle.load(open(features_filename, 'rb'))
    if 'mass_key' in features: features.remove('mass_key')
    return features



class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''
  def __init__(self, filename, selection, training_info, do_parametric=False, mass=None, ctau=None, colour=None, filter_efficiency=None, muon_rate=None):
    self.filename = filename
    self.selection = selection
    #self.selection = 'hnl_charge==0'
    cutbased_selection_qte = ['hnl_pt', 'hnl_charge', 'sv_lxysig', 'hnl_cos2d', 'hnl_mass'] 
    self.extra_branches = [
      'ismatched',
      'mu0_softid',
      'mu_looseid', 
      'pi_packedcandhashighpurity',
      'mu0_charge',
      'mu_charge',
      'b_mass',
      'mu0_mu_mass',
      'mu0_pi_mass',
      'sv_lxyz',
      'hnl_charge',
      'sv_lxy',
      'sv_chi2',
      #'weight_hlt_D1',
      'mu0_pt', 
      'mu_pt',
      'pi_pt',
      'b_pt',
      'mu0_mu_mass',
      'mu0_pi_mass',
      'deltar_mu0_mu',
      'deltar_mu0_pi',
      'hnl_cos2d',
      'sv_prob',
      'sv_lxysig',
      'pi_dcasig',
      'pi_dcasig_corr',
      'mu0_pfiso03_rel',
      'mu_pfiso03_rel',
      'pi_numberoftrackerlayers',
      'pi_numberofpixellayers',
      'mu_numberoftrackerlayers',
      'mu_numberofpixellayers',
      'mu0_numberoftrackerlayers',
      'mu0_numberofpixellayers',
      'mu0_eta',
      'mu_eta',
      'pi_eta',
      ]
    if mass != None and ctau != None: 
      self.extra_branches.append('gen_hnl_ct')
      #self.extra_branches.append('weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable')
      #self.extra_branches.append('weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2')
      self.extra_branches.append('weight_pu_sig_D')
      self.extra_branches.append('weight_mu0_softid')
      self.extra_branches.append('weight_mu_looseid')
    #print self.selection
    try:
      #self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True, columns=training_info.features+cutbased_selection_qte)
      self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True, columns=training_info.features+cutbased_selection_qte+self.extra_branches)
      #self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True)
    except:
      self.df = pd.DataFrame()
    self.mass = mass
    self.ctau = ctau
    self.colour = colour
    self.filter_efficiency = filter_efficiency
    self.muon_rate = muon_rate

    print '\t selected events: {}'.format(len(self.df))



class MVAAnalyser(Tools, MVATools):
  def __init__(self, signal_files, data_files, dirname, baseline_selection, categories=None, do_parametric=False, do_addCutbased=False, do_plotScore=False, do_createFiles=False, do_plotSigScan=False, do_plotROC=False, do_plotScoreCurve=False, do_plotMVAPerformance=False, do_plotWPScan=False, do_plotAUC=False, do_plotAUCvsLifetime=False, do_plotMass=False, do_plotPNNComparison=False, do_plotPreselection=False, do_plotSignalBackgroundComparison=False, do_plotDistributionComparison=False, do_compareROC=False, do_studyDisplacedTracks=False, do_plotScoreNorm=False):
    self.tools = Tools()
    self.mva_tools = MVATools()
    self.signal_files = signal_files
    self.data_files = data_files
    self.dirname = dirname
    self.baseline_selection = baseline_selection
    self.categories = categories
    self.do_parametric = do_parametric
    self.do_addCutbased = do_addCutbased
    self.do_plotScore = do_plotScore
    self.do_createFiles = do_createFiles
    self.do_plotSigScan = do_plotSigScan
    self.do_plotROC = do_plotROC
    self.do_plotScoreCurve = do_plotScoreCurve
    self.do_plotMVAPerformance = do_plotMVAPerformance
    self.do_plotWPScan = do_plotWPScan
    self.do_plotAUC = do_plotAUC
    self.do_plotMass = do_plotMass
    self.do_plotAUCvsLifetime = do_plotAUCvsLifetime
    self.do_plotPNNComparison = do_plotPNNComparison 
    self.do_plotPreselection = do_plotPreselection
    self.do_plotSignalBackgroundComparison = do_plotSignalBackgroundComparison
    self.do_plotDistributionComparison = do_plotDistributionComparison
    self.do_compareROC = do_compareROC
    self.do_studyDisplacedTracks = do_studyDisplacedTracks
    self.do_plotScoreNorm = do_plotScoreNorm

    self.outdir = self.createOutDir()

    self.resolution_p0 = 6.98338e-04
    self.resolution_p1 = 7.78382e-03

    self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2'
    #self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_efftable'
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
    print selection
    #selection = 'hnl_charge==0'
    #query = self.getPandasQuery(selection)

    print('========> starting reading the trees')
    now = time()
    ## background 
    data_samples = []
    for data_file in self.data_files:
      #data_samples.append(Sample(filename=data_file.filename, selection=selection))
      data_sample = Sample(filename=data_file, selection=selection, training_info=training_info)
      if len(data_sample.df) != 0:
        data_samples.append(data_sample)

    ## signal
    mc_samples = []
    for signal_file in self.signal_files:
      signal_filename = signal_file.filename if not category.is_bc else signal_file.filename_Bc
      filter_efficiency = signal_file.filter_efficiency if not category.is_bc else signal_file.filter_efficiency_Bc
      mc_sample = Sample(filename=signal_filename, selection=selection, training_info=training_info, mass=signal_file.mass, ctau=signal_file.ctau, colour=signal_file.colour, filter_efficiency=filter_efficiency, muon_rate=signal_file.muon_rate)
      if len(mc_sample.df) != 0:
        mc_samples.append(mc_sample)
      
    print('========> it took %.2f seconds' %(time() - now))

    return data_samples, mc_samples


  def getMCSamples(self, signal_files, training_info, category):
    '''
      Function that fetches the samples into lists
    '''

    selection = self.baseline_selection + ' && ' + category.definition_flat

    mc_samples = []
    for signal_file in signal_files:
      signal_filename = signal_file.filename if not category.is_bc else signal_file.filename_Bc
      filter_efficiency = signal_file.filter_efficiency if not category.is_bc else signal_file.filter_efficiency_Bc
      mc_sample = Sample(filename=signal_filename, selection=selection, training_info=training_info, mass=signal_file.mass, ctau=signal_file.ctau, colour=signal_file.colour, filter_efficiency=filter_efficiency, muon_rate=signal_file.muon_rate)
      if len(mc_sample.df) != 0:
        mc_samples.append(mc_sample)

    return mc_samples


  def createDataframe(self, samples):
    '''
      Function that converts root samples to pandas dataframe
    '''
    df = pd.concat([idt.df for idt in samples], sort=False)

    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

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
    if self.do_parametric:
      x = pd.DataFrame(df, columns=training_info.features + ['mass_key'])#, 'weight_hlt_D1'])
    else:
      x = pd.DataFrame(df, columns=training_info.features)

    # apply the scaler
    if self.do_parametric:
      xx = training_info.qt.transform(x[training_info.features + ['mass_key']])#, 'weight_hlt_D1']])
    else:
      xx = training_info.qt.transform(x[training_info.features])

    # predict
    score = training_info.model.predict(xx)

    return score


  def plotScore(self, training_info, mc_samples, data_samples, category, do_log):
    '''
      Plot the score distributions for signal and background
    '''
    pd.options.mode.chained_assignment = None
    canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
    if do_log:
      canv.SetLogy()
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    print masses
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    # background
    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
    if self.do_parametric:
      data_df['mass_key'] = mass
    bkg_score = self.predictScore(training_info, data_df)

    hist_bkg = ROOT.TH1F('bkg', 'bkg', 30, 0, 1)
    for score in bkg_score:
      hist_bkg.Fill(score)

    if hist_bkg.Integral()!=0: hist_bkg.Scale(1./hist_bkg.Integral())

    leg.AddEntry(hist_bkg, 'data-driven background')

    ROOT.gStyle.SetOptStat(0)
    hist_bkg.SetTitle(' ')
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
      if self.do_parametric:
        mc_df['mass_key'] = mass
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
    name = '{}/score_m{}_{}'.format(self.outdir, str(mass).replace('.', 'p'), category.label)
    if do_log: name += '_log'
    canv.SaveAs(name + '.png')


  def plotScoreNorm(self, training_info, mc_samples, data_samples, category, do_log):
    '''
      Plot the score distributions for signal and background
    '''
    pd.options.mode.chained_assignment = None

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    #print masses
    #if len(masses) != 1:
    #  raise RuntimeError('Please provide signal samples of the same mass')

    for mass in masses:
      canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
      if do_log:
        canv.SetLogy()
      leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)

      resolution = self.resolution_p0 + self.resolution_p1 * mass

      # background
      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
      if self.do_parametric:
        data_df['mass_key'] = mass
      bkg_score = self.predictScore(training_info, data_df)

      hist_bkg = ROOT.TH1F('bkg', 'bkg', 30, 0, 1)
      for score in bkg_score:
        hist_bkg.Fill(score)

      leg.AddEntry(hist_bkg, 'data (5.3 /fb)')

      ROOT.gStyle.SetOptStat(0)
      hist_bkg.SetTitle(category.title)
      hist_bkg.SetFillColor(ROOT.kAzure-4)
      hist_bkg.SetLineColor(1)
      hist_bkg.GetXaxis().SetTitle('Score')
      hist_bkg.GetXaxis().SetLabelSize(0.033)
      hist_bkg.GetXaxis().SetTitleSize(0.042)
      hist_bkg.GetXaxis().SetTitleOffset(1.1)
      hist_bkg.GetYaxis().SetTitle('Entries')
      hist_bkg.GetYaxis().SetLabelSize(0.033) 
      hist_bkg.GetYaxis().SetTitleSize(0.042)
      hist_bkg.GetYaxis().SetTitleOffset(1.1)
      hist_bkg.GetYaxis().SetRangeUser(1e-1, hist_bkg.GetMaximum() + 0.15*hist_bkg.GetMaximum())

      # signals
      hist_sigs = []
      for mc_sample in mc_samples:
        if mc_sample.mass != mass: continue
        mc_df = self.createDataframe([mc_sample])
        if self.do_parametric:
          mc_df['mass_key'] = mass
        sig_score = self.predictScore(training_info, mc_df)

        hist_name = 'sig_{}_{}'.format(mc_sample.mass, mc_sample.ctau)
        hist_sig = ROOT.TH1F(hist_name, hist_name, 30, 0, 1)
        hist_sig.SetDirectory(0)
        for score in sig_score:
          hist_sig.Fill(score)

        # previous normalisation strategy
        weight_signal = self.tools.getSignalWeight(signal_files=[mc_sample], mass=mass, ctau=mc_sample.ctau, sigma_B=572.0e9, lumi=5.302, lhe_efficiency=0.08244, is_bu=False, is_bd=False, is_bs=False, is_bc=False)
        hist_sig.Scale(weight_signal)
        #if hist_sig.Integral()!=0: hist_sig.Scale(1./hist_sig.Integral())

        leg.AddEntry(hist_sig, 'signal - mass {} GeV, ctau {} mm'.format(mc_sample.mass, mc_sample.ctau))

        hist_sig.SetLineWidth(3)
        hist_sig.SetLineColor(mc_sample.colour)
        hist_sigs.append(hist_sig)

      hist_bkg.Draw('histo')
      for hist_sig in hist_sigs:
        hist_sig.Draw('histo_same')
      leg.Draw('same')

      canv.cd()
      name = '{}/score_norm_m{}_{}'.format(self.outdir, str(mass).replace('.', 'p'), category.label)
      if do_log: name += '_log'
      canv.SaveAs(name + '.png')



  def plotROCCurve(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    # create dataframe
    data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
    if self.do_parametric:
      data_df['mass_key'] = mass
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      if self.do_parametric:
        mc_df['mass_key'] = mass
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
      if self.do_addCutbased:
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

      plt.title(r'{}'.format(category.label))
      if do_log:
        plt.xlim(1e-5, 1)
      else:
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


  def plotAUCGraph(self, training_info, mc_samples, data_samples, category):
    pd.options.mode.chained_assignment = None

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)

    graph = ROOT.TGraph()

    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      resolution = self.resolution_p0 + self.resolution_p1 * mass

      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      # create dataframe
      data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
      if self.do_parametric:
        data_df['mass_key'] = mass
      data_df['is_signal'] = 0

      mc_df = self.createDataframe([mc_sample])
      if self.do_parametric:
        mc_df['mass_key'] = mass
      mc_df['is_signal'] = 1

      main_df = pd.concat([data_df, mc_df], sort=False)
      main_df.index = np.array(range(len(main_df)))
      main_df = main_df.sample(frac=1, replace=False, random_state=1986)

      Y = pd.DataFrame(main_df, columns=['is_signal'])
      score = self.predictScore(training_info=training_info, df=main_df)
      fpr, tpr, thresholds = roc_curve(Y, score) 
      auc = sklearn_auc(fpr, tpr)

      point = graph.GetN() 
      graph.SetPoint(point, mass, auc)
      graph.SetMarkerColor(ROOT.kBlue)
      graph.SetMarkerStyle(20)
      graph.GetXaxis().SetTitle('Signal mass [GeV]')
      graph.GetXaxis().SetLabelSize(0.037)
      graph.GetXaxis().SetTitleSize(0.042)
      graph.GetXaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetTitle('AUC')
      graph.GetYaxis().SetLabelSize(0.037)
      graph.GetYaxis().SetTitleSize(0.042)
      graph.GetYaxis().SetTitleOffset(1.1)
      graph.GetYaxis().SetRangeUser(0.7, 1)

    canv = self.tools.createTCanvas('canv', 800, 700) 
    graph.Draw('AP')

    name = 'AUC_{}'.format(category.label)
    canv.SaveAs('{}/{}.png'.format(self.outdir, name))


  def plotAUCvsLifetimeGraph(self, training_info, mc_samples, data_samples, category):
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.13)

    #masses = []
    #for mc_sample in mc_samples:
    #  if mc_sample.mass not in masses: masses.append(mc_sample.mass)

    #graph = ROOT.TGraph()

    ctaus = [1, 10, 100, 1000]
    colours = [ROOT.kRed-9, ROOT.kRed-7, ROOT.kRed+2, ROOT.kRed+4]

    graphs = []

    leg = self.tools.getRootTLegend(xmin=0.5, ymin=0.15, xmax=0.82, ymax=0.35, size=0.05)

    used_ctaus = []
    for ictau, ctau in enumerate(ctaus):
      graph = ROOT.TGraph()

      for mc_sample in mc_samples:
        if mc_sample.ctau != ctau: continue 

        mass = mc_sample.mass
        resolution = self.resolution_p0 + self.resolution_p1 * mass

        # consider the 10 sigma window around the signal mass
        window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
        # create dataframe
        data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
        if self.do_parametric:
          data_df['mass_key'] = mass
        data_df['is_signal'] = 0

        mc_df = self.createDataframe([mc_sample])
        if self.do_parametric:
          mc_df['mass_key'] = mass
        mc_df['is_signal'] = 1

        main_df = pd.concat([data_df, mc_df], sort=False)
        main_df.index = np.array(range(len(main_df)))
        main_df = main_df.sample(frac=1, replace=False, random_state=1986)

        Y = pd.DataFrame(main_df, columns=['is_signal'])
        score = self.predictScore(training_info=training_info, df=main_df)
        fpr, tpr, thresholds = roc_curve(Y, score) 
        auc = sklearn_auc(fpr, tpr)

        #graph = ROOT.TGraph()
        point = graph.GetN() 
        graph.SetLineColor(colours[ictau])
        graph.SetLineWidth(2)
        graph.SetLineStyle(9)
        graph.SetPoint(point, mass, auc)
        graph.SetMarkerColor(colours[ictau])
        graph.SetMarkerStyle(20)
        #graph.SetMarkerSize(2)

        if ctau not in used_ctaus:
          leg.AddEntry(graph, '{c} = {ctau} mm'.format(c='#it{c}#tau_{N}', ctau=ctau))
          used_ctaus.append(ctau)

        #graph.GetXaxis().SetTitle('Signal mass [GeV]')
        #graph.GetXaxis().SetLabelSize(0.037)
        #graph.GetXaxis().SetTitleSize(0.042)
        #graph.GetXaxis().SetTitleOffset(1.1)
        #graph.GetXaxis().SetRangeUser(0.5, 5)
        #graph.GetYaxis().SetTitle('AUC')
        #graph.GetYaxis().SetLabelSize(0.037)
        #graph.GetYaxis().SetTitleSize(0.042)
        #graph.GetYaxis().SetTitleOffset(1.1)
        #graph.GetYaxis().SetRangeUser(0.7, 1)
        graphs.append(graph)

    frame = ROOT.TGraph()
    range_min = 1.
    #range_max = 4.5
    range_max = 3.
    frame.SetPoint(0, range_min, 0.94)
    frame.SetPoint(1, range_max, 1.012)
    frame.GetXaxis().SetTitle('#it{m}_{N} (GeV)')
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.047)
    frame.GetXaxis().SetTitleOffset(1.)
    frame.GetYaxis().SetTitle('AUC')
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.047)
    frame.GetYaxis().SetTitleOffset(1.4)
      
    canv = self.tools.createTCanvas('canv', 800, 700) 
    pad = ROOT.TPad("pad","pad",0,0,1,1)
    pad.Draw()
    pad.cd()

    frame.Draw('AP')
    for igraph, graph in enumerate(graphs):
      graph.Draw('PL same')

    line = ROOT.TLine(frame.GetXaxis().GetXmin(), 1, frame.GetXaxis().GetXmax(), 1)
    line.SetLineColor(1)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw('same')

    #self.tools.printCMSTagInFrame(pad, 'Preliminary', size=0.55)
    self.tools.printLatexBox(0.22, 0.82, 'CMS', size=0.06, font=61)
    #self.tools.printLatexBox(0.4, 0.817, 'Preliminary', size=0.055, font=52)
    self.tools.printLumiTag(pad, 5.2, size=0.5, offset=0.52)
    self.tools.printLatexBox(0.66, 0.4, category.title, size=0.05, font=22)

    leg.Draw()

    name = 'AUC_vs_lifetime_{}'.format(category.label)
    canv.cd()
    canv.SaveAs('{}/{}.png'.format(self.outdir, name))
    canv.SaveAs('{}/{}.pdf'.format(self.outdir, name))
    canv.SaveAs('{}/{}.C'.format(self.outdir, name))


  def plotPNNComparisonGraph(self, mc_samples, data_samples, category):
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.13)

    #dirname_1 = 'test_2022Nov29_09h26m28s' # adding muon isolation
    #dirname_2 = 'test_2022Oct12_15h12m37s' # mass 3 only (without muon isolation)
    dirname_1 = 'training_Aug23'
    dirname_2 = 'V13_06Feb23_m3_2023Nov28_20h32m09s'
      
    training_info_1 = TrainingInfo(dirname_1, category.label)
    training_info_2 = TrainingInfo(dirname_2, category.label)

    #signal_labels_trained_1 = ['V42_08Aug22_m1p0', 'V42_08Aug22_m1p5', 'V42_08Aug22_m2p0', 'V42_08Aug22_m3p0', 'V42_08Aug22_m4p5']
    #signal_labels_nottrained_1 = ['V42_08Aug22_m1p26', 'V42_08Aug22_m1p77', 'V42_08Aug22_m2p5', 'V42_08Aug22_m3p4', 'V42_08Aug22_m4p1']
    #signal_labels_trained_2 = ['V42_08Aug22_m3p0']
    #signal_labels_nottrained_2 = ['V42_08Aug22_m1p0', 'V42_08Aug22_m1p5', 'V42_08Aug22_m2p0', 'V42_08Aug22_m4p5', 'V42_08Aug22_m1p26', 'V42_08Aug22_m1p77', 'V42_08Aug22_m2p5', 'V42_08Aug22_m3p4', 'V42_08Aug22_m4p1']
    #signal_labels_tot = ['V42_08Aug22_m1p0', 'V42_08Aug22_m1p26', 'V42_08Aug22_m1p5', 'V42_08Aug22_m1p77', 'V42_08Aug22_m2p0', 'V42_08Aug22_m2p5', 'V42_08Aug22_m3p0', 'V42_08Aug22_m3p4', 'V42_08Aug22_m4p1', 'V42_08Aug22_m4p5']

    #signal_labels_trained_1 = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p5', 'V42_06Feb23_m2p0', 'V42_06Feb23_m3p0', 'V42_06Feb23_m4p5']
    #signal_labels_nottrained_1 = ['V42_06Feb23_m1p26', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p5', 'V42_06Feb23_m3p4', 'V42_06Feb23_m4p1']
    #signal_labels_trained_2 = ['V42_06Feb23_m3p0']
    #signal_labels_nottrained_2 = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p5', 'V42_06Feb23_m2p0', 'V42_06Feb23_m4p5', 'V42_06Feb23_m1p26', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p5', 'V42_06Feb23_m3p4', 'V42_06Feb23_m4p1']
    #signal_labels_tot = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p26', 'V42_06Feb23_m1p5', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p0', 'V42_06Feb23_m2p5', 'V42_06Feb23_m3p0', 'V42_06Feb23_m3p4', 'V42_06Feb23_m4p1', 'V42_06Feb23_m4p5']
    signal_labels_trained_1 = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p5', 'V42_06Feb23_m2p0', 'V42_06Feb23_m3p0']
    signal_labels_nottrained_1 = ['V42_06Feb23_m1p26', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p5']
    signal_labels_trained_2 = ['V42_06Feb23_m3p0']
    signal_labels_nottrained_2 = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p5', 'V42_06Feb23_m2p0', 'V42_06Feb23_m4p5', 'V42_06Feb23_m1p26', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p5']
    signal_labels_tot = ['V42_06Feb23_m1p0', 'V42_06Feb23_m1p26', 'V42_06Feb23_m1p5', 'V42_06Feb23_m1p77', 'V42_06Feb23_m2p0', 'V42_06Feb23_m2p5', 'V42_06Feb23_m3p0']

    signal_files_trained_1 = []
    for signal_label in signal_labels_trained_1:
      for signal_file in signal_samples[signal_label]:
        signal_files_trained_1.append(signal_file)

    #for sig in signal_files_trained_1:
    #  print sig.filename

    signal_files_trained_2 = []
    for signal_label in signal_labels_trained_2:
      for signal_file in signal_samples[signal_label]:
        signal_files_trained_2.append(signal_file)

    signal_files_nottrained_1 = []
    for signal_label in signal_labels_nottrained_1:
      for signal_file in signal_samples[signal_label]:
        signal_files_nottrained_1.append(signal_file)

    signal_files_nottrained_2 = []
    for signal_label in signal_labels_nottrained_2:
      for signal_file in signal_samples[signal_label]:
        signal_files_nottrained_2.append(signal_file)

    signal_files_tot = []
    for signal_label in signal_labels_tot:
      for signal_file in signal_samples[signal_label]:
        signal_files_tot.append(signal_file)
      
    mc_samples_trained_1 = self.getMCSamples(signal_files=signal_files_trained_1, training_info=training_info_1, category=category)
    mc_samples_trained_2 = self.getMCSamples(signal_files=signal_files_trained_2, training_info=training_info_2, category=category)
    mc_samples_nottrained_1 = self.getMCSamples(signal_files=signal_files_nottrained_1, training_info=training_info_1, category=category)
    mc_samples_nottrained_2 = self.getMCSamples(signal_files=signal_files_nottrained_2, training_info=training_info_2, category=category)
    mc_samples_tot = self.getMCSamples(signal_files=signal_files_tot, training_info=training_info_1, category=category)

    mc_samples_1 = [mc_samples_trained_1, mc_samples_nottrained_1]
    mc_samples_2 = [mc_samples_trained_2, mc_samples_nottrained_2]

    colours = [ROOT.kBlue, ROOT.kRed]
    markerstyle = [20, 24]

    graphs = []
    legends = []

    leg = self.tools.getRootTLegend(xmin=0.3, ymin=0.15, xmax=0.8, ymax=0.35, size=0.05)

    for itraining, training_info in enumerate([training_info_1, training_info_2]):
      if itraining == 0:
        mc_samples_list = mc_samples_1
      else:
        mc_samples_list = mc_samples_2

      for ilist, mc_sample_list in enumerate(mc_samples_list):
        graph = ROOT.TGraph()
        for imc, mc_sample in enumerate(mc_sample_list):
          mass = mc_sample.mass
          resolution = self.resolution_p0 + self.resolution_p1 * mass

          # consider the 10 sigma window around the signal mass
          window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
          # create dataframe
          data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
          if self.do_parametric:
            data_df['mass_key'] = mass
          data_df['is_signal'] = 0

          mc_df = self.createDataframe([mc_sample])
          if self.do_parametric:
            mc_df['mass_key'] = mass
          mc_df['is_signal'] = 1

          main_df = pd.concat([data_df, mc_df], sort=False)
          main_df.index = np.array(range(len(main_df)))
          main_df = main_df.sample(frac=1, replace=False, random_state=1986)

          Y = pd.DataFrame(main_df, columns=['is_signal'])
          score = self.predictScore(training_info=training_info, df=main_df)
          fpr, tpr, thresholds = roc_curve(Y, score) 
          auc = sklearn_auc(fpr, tpr)

          point = graph.GetN() 
          graph.SetPoint(point, mass, auc)
          graph.SetLineColor(colours[itraining])
          graph.SetLineWidth(2)
          graph.SetLineStyle(9)
          graph.SetMarkerColor(colours[itraining])
          graph.SetMarkerStyle(markerstyle[ilist])
          graph.SetMarkerSize(2)
          if itraining == 0: 
            training_label = 'pNN'
          else:
            #training_label = 'NN trained on 3 GeV'
            training_label = 'NN'
          if ilist == 0:
            point_label = '(point trained on)'
          else:
            point_label = '(point not trained on)'

          legend = '{} {}'.format(training_label, point_label)
          if legend not in legends:
            leg.AddEntry(graph, legend)
            legends.append(legend)

          graphs.append(graph)

    mc_samples_merged_1 = mc_samples_trained_1 + mc_samples_nottrained_1
    mc_samples_merged_2 = mc_samples_trained_2 + mc_samples_nottrained_2

    colours = [ROOT.kBlue, ROOT.kRed]
    markerstyle = [20, 24]

    graphs_merged = []
    #masses = [1, 1.26, 1.5, 1.77, 2., 2.5, 3., 3.4, 4.1, 4.5]

    for itraining, training_info in enumerate([training_info_1, training_info_2]):
      graph = ROOT.TGraph()
      for ilist, mc_sample in enumerate(mc_samples_tot):
        mass = mc_sample.mass
        resolution = self.resolution_p0 + self.resolution_p1 * mass

        # consider the 10 sigma window around the signal mass
        window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
        # create dataframe
        data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
        if self.do_parametric:
          data_df['mass_key'] = mass
        data_df['is_signal'] = 0

        mc_df = self.createDataframe([mc_sample])
        if self.do_parametric:
          mc_df['mass_key'] = mass
        mc_df['is_signal'] = 1

        main_df = pd.concat([data_df, mc_df], sort=False)
        main_df.index = np.array(range(len(main_df)))
        main_df = main_df.sample(frac=1, replace=False, random_state=1986)

        Y = pd.DataFrame(main_df, columns=['is_signal'])
        score = self.predictScore(training_info=training_info, df=main_df)
        fpr, tpr, thresholds = roc_curve(Y, score) 
        auc = sklearn_auc(fpr, tpr)

        point = graph.GetN() 
        graph.SetPoint(point, mass, auc)
        graph.SetLineColor(colours[itraining])
        graph.SetLineWidth(2)
        graph.SetLineStyle(9)
        #graph.SetMarkerColor(colours[itraining])
        #graph.SetMarkerStyle(markerstyle[ilist])

        graphs_merged.append(graph)

    frame = ROOT.TGraph()
    range_min = 1.
    #range_max = 4.5
    range_max = 3.
    frame.SetPoint(0, range_min, 0.77)
    frame.SetPoint(1, range_max, 1.05)
    frame.GetXaxis().SetTitle('#it{m}_{N} (GeV)')
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.047)
    frame.GetXaxis().SetTitleOffset(1.)
    frame.GetYaxis().SetTitle('AUC')
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.047)
    frame.GetYaxis().SetTitleOffset(1.4)
      
    canv = self.tools.createTCanvas('canv', 800, 700) 
    pad = ROOT.TPad("pad","pad",0,0,1,1)
    pad.Draw()
    pad.cd()

    frame.Draw('AP')
    for graph in graphs:
      graph.Draw('P same')

    for graph in graphs_merged:
      graph.Draw('L same')

    #self.tools.printCMSTagInFrame(pad, 'Preliminary', size=0.55)
    self.tools.printLatexBox(0.22, 0.82, 'CMS', size=0.06, font=61)
    #self.tools.printLatexBox(0.4, 0.817, 'Preliminary', size=0.055, font=52)
    self.tools.printLumiTag(pad, 5.2, size=0.5, offset=0.52)
    self.tools.printLatexBox(0.66, 0.4, category.title, size=0.05, font=22)

    leg.Draw()

    name = 'pNN_comparison_{}'.format(category.label)
    canv.SaveAs('{}/{}.png'.format(self.outdir, name))
    canv.SaveAs('{}/{}.pdf'.format(self.outdir, name))
    canv.SaveAs('{}/{}.C'.format(self.outdir, name))


  def plotScoreCurve(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    # create dataframe
    data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
    if self.do_parametric:
      data_df['mass_key'] = mass
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      if self.do_parametric:
        mc_df['mass_key'] = mass
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

      plt.title(r'{}'.format(category.label))
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

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    # create dataframe
    data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
    if self.do_parametric:
      data_df['mass_key'] = mass
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      if self.do_parametric:
        mc_df['mass_key'] = mass
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
      plt.title(r'{}'.format(category.label))
      #plt.xlim(0, 1)
      plt.xlim(0.8, 1)
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


  def plotWPScan(self, training_info, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None

    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    # create dataframe
    data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
    if self.do_parametric:
      data_df['mass_key'] = mass
    data_df['is_signal'] = 0

    plt.clf()
    for mc_sample in mc_samples:
      mass = mc_sample.mass
      ctau = mc_sample.ctau
      v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
      coupling = self.tools.getCouplingLabel(v2)

      mc_df = self.createDataframe([mc_sample])
      if self.do_parametric:
        mc_df['mass_key'] = mass
      mc_df['is_signal'] = 1

      main_df = pd.concat([data_df, mc_df], sort=False)
      main_df.index = np.array(range(len(main_df)))
      main_df = main_df.sample(frac=1, replace=False, random_state=1986)

      # get mva performance
      Y = pd.DataFrame(main_df, columns=['is_signal'])
      score = self.predictScore(training_info=training_info, df=main_df)
      false_positive_rate_mva, true_positive_rate_mva, thresholds = roc_curve(Y, score) 

      # get performance without cut on score
      true_positive_rate_ini = 1.#true_positive_rate_mva[len(true_positive_rate_mva)]
      false_positive_rate_ini = 1.#false_positive_rate_mva[len(false_positive_rate_ini)]
      #print thresholds
      #print true_positive_rate_ini 
      #print '\n\n\n'
      #print false_positive_rate_ini

      relative_signal_efficiency = true_positive_rate_mva / true_positive_rate_ini
      relative_background_efficiency = false_positive_rate_mva / false_positive_rate_ini
      relative_significance = relative_signal_efficiency / np.sqrt(relative_background_efficiency)

      plt.plot(thresholds, relative_significance, linewidth=2, label='({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      #plt.plot(thresholds, false_positive_rate_mva, linewidth=2, label='({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      #plt.plot(thresholds, true_positive_rate_mva, linewidth=2, label='({}GeV, {}mm, {})'.format(mass, ctau, coupling))
      plt.xlabel('Score')
      plt.ylabel('tpr(mva)/tpr(0) / sqrt(fpr(mva)/fpr(0))')

      plt.axhline(y = 1., color='black', linestyle = '--')
      plt.title(r'{}'.format(category.label))
      #plt.xlim(0, 1)
      plt.xlim(0.8, 1)
      #plt.ylim(0, 1)

      if do_log:
        plt.xscale('log')
      plt.yscale('linear')

      if do_log:
        plt.legend(loc='upper left', framealpha=0.1)
      else:
        plt.legend(loc='lower right', framealpha=0.1)

    name = 'working_point_m{}_{}'.format(str(mc_samples[0].mass).replace('.', 'p'), category.label)
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
    scan_points = np.linspace(0.9, 1, n_points) 
      
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
      signal_resolution = self.resolution_p0 + self.resolution_p1 * signal_mass
      signal_ctau = mc_sample.ctau
      #if signal_ctau != 1000.: continue

      # get signal tree, with score branch, without selection
      extra_branches = [
        'ismatched',
        'mu0_softid',
        'mu_looseid', 
        'pi_packedcandhashighpurity',
        'mu0_charge',
        'mu_charge',
        'mu0_mu_mass',
        'mu0_pi_mass',
        'sv_lxyz',
        'hnl_charge',
        'sv_lxy',
        ]

      sig_label = 'sig_{}_{}_{}'.format(signal_mass, signal_ctau, category.label)
      filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=signal_mass, selection=self.baseline_selection, weights=extra_branches+['gen_hnl_ct', self.weight_hlt, self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=False) 
      file_sig = self.tools.getRootFile(filename_sig)
      tree_sig = self.tools.getTree(file_sig, 'signal_tree')

      bkg_label = 'bkg_{}_{}'.format(signal_mass, category.label)
      filename_data = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=signal_mass, selection=self.baseline_selection, weights=extra_branches, label=bkg_label, treename='signal_tree', force_overwrite=False) 
      file_data = self.tools.getRootFile(filename_data)
      tree_data = self.tools.getTree(file_data, 'signal_tree')

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

      if self.do_addCutbased:
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

    #os.system('rm *m{}*{}.root'.format(signal_mass, category.label))

    canv.cd()
    canv.SaveAs('{}/significance_diff_scan_m{}_{}.png'.format(self.outdir, str(signal_mass).replace('.', 'p'), category.label))


  def plotMass(self, training_info, mc_samples, data_samples, category, object_mass):
    '''
      Plot the background mass distribution for given cut on scores
    '''
    if object_mass not in ['hnl_mass', 'b_mass']:
      raise RuntimeError('Unknown object_mass. Choose among ["hnl_mass", "b_mass"]')

    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetOptStat(0)
    canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)
    
    masses = []
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    if len(masses) != 1:
      raise RuntimeError('Please provide signal samples of the same mass')

    #scores = [0.95, 0.]
    scores = [0.99, 0.]
    colours = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]

    mass = mc_samples[0].mass
    resolution = self.resolution_p0 + self.resolution_p1 * mass

    if object_mass == 'hnl_mass':
      quantity = Qte(name_nano='BToMuMuPi_hnl_mass', name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=80, bin_min=mass-10*resolution, bin_max=mass+10*resolution)
    elif object_mass == 'b_mass':
      quantity = Qte(name_nano='BToMuMuPi_mass', name_flat='b_mass', label='b_mass', title='#mu#mu#pi invariant mass [GeV]', nbins=80, bin_min=0, bin_max=8)

    # consider the 10 sigma window around the signal mass
    window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
    bkg_label = 'bkg_{}_{}'.format(mass, category.label)
    extra_branch = ['b_mass']
    selection_bkg = self.baseline_selection + ' && ' + window + ' && ' + category.definition_flat
    #TODO understand why the category selection should be added while the samples are selected with that selection applied
    print selection_bkg
    #filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=self.baseline_selection + ' && ' + window, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
    filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=selection_bkg, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
    file_bkg = self.tools.getRootFile(filename_bkg)
    tree_bkg = self.tools.getTree(file_bkg, 'signal_tree')

    hists = []
    for iscore, score in enumerate(scores):
      hist_bkg = self.tools.createHisto(tree_bkg, quantity, hist_name='hist_bkg'+str(score), branchname='flat', selection='score > {}'.format(score))
      if hist_bkg.Integral() != 0: hist_bkg.Scale(1./hist_bkg.Integral())

      if score == 0.:
        hist_bkg.SetFillColor(ROOT.kBlue-3)
        hist_bkg.SetFillStyle(3005)
      else:
        hist_bkg.SetMarkerStyle(20)
        hist_bkg.SetMarkerColor(colours[iscore])
      
      hist_bkg.SetTitle(category.title)
      if object_mass == 'hnl_mass':
        x_label = '#mu#pi invariant mass [GeV]'
      elif object_mass == 'b_mass':
        x_label = '#mu#mu#pi invariant mass [GeV]'
      hist_bkg.GetXaxis().SetTitle(x_label)
      hist_bkg.GetXaxis().SetLabelSize(0.033)
      hist_bkg.GetXaxis().SetTitleSize(0.042)
      hist_bkg.GetXaxis().SetTitleOffset(1.1)
      hist_bkg.GetYaxis().SetTitle('Normalised to unity')
      hist_bkg.GetYaxis().SetLabelSize(0.033)
      hist_bkg.GetYaxis().SetTitleSize(0.042)
      hist_bkg.GetYaxis().SetTitleOffset(1.1)
      hists.append(hist_bkg)
      if score == 0:
        leg_entry = 'no cut on the pNN score'
      else:
        leg_entry = 'pNN score > {}'.format(score)
      leg.AddEntry(hist_bkg, leg_entry)

    for ihist, hist in enumerate(hists):
      if ihist == 0: hist.Draw('')
      elif ihist != 0 and ihist != len(scores)-1: hist.Draw('same')
      else: hist.Draw('hist same')
    leg.Draw('same')

    canv.cd()
    name = '{}/{}_m{}_{}'.format(self.outdir, object_mass, str(mass).replace('.', 'p'), category.label)
    canv.SaveAs(name + '.png')


  def plotPreselectionQuantity(self, training_info, mc_samples, quantity, category):
    '''
      Plot the background mass distribution for given cut on scores
    '''
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetOptStat(0)
    canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
    leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.45, xmax=0.65, ymax=0.83, size=0.04)
    
    #masses = []
    #for mc_sample in mc_samples:
    #  if mc_sample.mass not in masses: masses.append(mc_sample.mass)
    #if len(masses) != 1:
    #  raise RuntimeError('Please provide signal samples of the same mass')

    #scores = [0.99, 0.]
    scores = [0.95]
    colours = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]

    #mass = mc_samples[0].mass
    #resolution = self.resolution_p0 + self.resolution_p1 * mass

    #if object_mass == 'hnl_mass':
    #  quantity = Qte(name_nano='BToMuMuPi_hnl_mass', name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=80, bin_min=mass-10*resolution, bin_max=mass+10*resolution)
    #elif object_mass == 'b_mass':
    #  quantity = Qte(name_nano='BToMuMuPi_mass', name_flat='b_mass', label='b_mass', title='#mu#mu#pi invariant mass [GeV]', nbins=80, bin_min=0, bin_max=8)

    # consider the 10 sigma window around the signal mass
    #window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)

    hists = []
    for mc_sample in mc_samples:
      sig_label = 'sig_{}_{}_{}'.format(mc_sample.mass, mc_sample.ctau, category.label)
      extra_branch = [quantity.name_flat]
      #filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=self.baseline_selection + ' && ' + window, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
      #filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=self.baseline_selection, weights=extra_branch+['gen_hnl_ct', self.weight_hlt, self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=True) 
      filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=self.baseline_selection, weights=extra_branch+['gen_hnl_ct', self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=True) 
      file_sig = self.tools.getRootFile(filename_sig)
      tree_sig = self.tools.getTree(file_sig, 'signal_tree')

      #hists = []
      for iscore, score in enumerate(scores):
        #score = 0.99
        hist_sig = self.tools.createHisto(tree_sig, quantity, hist_name='hist_sig'+str(score), branchname='flat', selection='score > {}'.format(score))
        if hist_sig.Integral() != 0: hist_sig.Scale(1./hist_sig.Integral())

        if score == 0.:
          hist_sig.SetFillColor(mc_sample.colour)#ROOT.kBlue-3)
          hist_sig.SetFillStyle(3005)
        else:
          hist_sig.SetMarkerStyle(20)
          hist_sig.SetMarkerColor(mc_sample.colour)#colours[iscore])
        
        hist_sig.SetTitle(category.title)
        #if object_mass == 'hnl_mass':
        #  x_label = '#mu#pi invariant mass [GeV]'
        #elif object_mass == 'b_mass':
        #  x_label = '#mu#mu#pi invariant mass [GeV]'
        hist_sig.GetXaxis().SetTitle(quantity.title)#x_label)
        hist_sig.GetXaxis().SetLabelSize(0.033)
        hist_sig.GetXaxis().SetTitleSize(0.042)
        hist_sig.GetXaxis().SetTitleOffset(1.1)
        hist_sig.GetYaxis().SetTitle('Normalised to unity')
        hist_sig.GetYaxis().SetLabelSize(0.033)
        hist_sig.GetYaxis().SetTitleSize(0.042)
        hist_sig.GetYaxis().SetTitleOffset(1.1)
        hists.append(hist_sig)
        leg_entry = '{} GeV {} mm score>{}'.format(mc_sample.mass, mc_sample.ctau, score)
        #if score == 0:
        #  leg_entry = '{} GeV {} mm score>0'.format(mc_sample.mass, mc_sample.ctau)
        #else:
        #  leg_entry = 'pNN score > {}'.format(score)
        leg.AddEntry(hist_sig, leg_entry)

    for ihist, hist in enumerate(hists):
      if ihist == 0: hist.Draw('')
      else: hist.Draw('same')
      #elif ihist != 0 and ihist != len(scores)-1: hist.Draw('same')
      #else: hist.Draw('hist same')
    leg.Draw('same')

    canv.cd()
    os.system('mkdir -p {}/preselection'.format(self.outdir))
    name = '{}/preselection/{}_{}'.format(self.outdir, quantity.label, category.label)
    canv.SaveAs(name + '.png')


  def plotSignalBackgroundComparison(self, training_info, mc_samples, data_samples, category, quantity, score):
    '''
      Plot the quantity distribution for signal and background
    '''
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetOptStat(0)

    masses = [] # [1.0, 1.5, 2.0, 3.0, 4.5]
    #masses = [3.0]
    for mc_sample in mc_samples:
      if mc_sample.mass not in masses:
        masses.append(mc_sample.mass)

    for mass in masses:
      canv = self.tools.createTCanvas('canv'+category.label+str(mass), 800, 700) 
      leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)
      
      #score = 0.99
      colours = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]

      resolution = self.resolution_p0 + self.resolution_p1 * mass

      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      bkg_label = 'bkg_{}_{}'.format(mass, category.label)
      #extra_branch = [quantity.name_flat]
      extra_branch = [
        'hnl_mass',
        'b_mass',
        'mu0_pt', 
        'mu_pt',
        'pi_pt',
        'mu0_mu_mass',
        'mu0_pi_mass',
        'deltar_mu0_mu',
        'deltar_mu0_pi',
        'hnl_cos2d',
        'sv_prob',
        'sv_lxysig',
        'pi_dcasig',
        'mu0_pfiso03_rel',
        'mu_pfiso03_rel',
        'pi_numberoftrackerlayers',
        'pi_numberofpixellayers',
        'mu_numberoftrackerlayers',
        'mu_numberofpixellayers',
        'mu0_numberoftrackerlayers',
        'mu0_numberofpixellayers',
        'sv_lxy',
        'mu0_eta',
        'mu_eta',
        'pi_eta',
      ]
      selection_bkg = self.baseline_selection + ' && ' + window + ' && ' + category.definition_flat
      print selection_bkg
      #filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=self.baseline_selection + ' && ' + window, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
      #filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=selection_bkg, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
      filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=selection_bkg, label=bkg_label, treename='signal_tree', force_overwrite=True, weights=extra_branch) 
      file_bkg = self.tools.getRootFile(filename_bkg)
      tree_bkg = self.tools.getTree(file_bkg, 'signal_tree')

      hist_bkg = self.tools.createHisto(tree_bkg, quantity, hist_name='hist_bkg'+str(score), branchname='flat', selection='score > {}'.format(score))
      if hist_bkg.Integral() != 0: hist_bkg.Scale(1./hist_bkg.Integral())

      hist_bkg.SetFillColor(ROOT.kBlue-3)
      hist_bkg.SetFillStyle(3005)
        
      hist_bkg.SetTitle(category.title)
      x_label = quantity.title
      hist_bkg.GetXaxis().SetTitle(x_label)
      hist_bkg.GetXaxis().SetLabelSize(0.033)
      hist_bkg.GetXaxis().SetTitleSize(0.042)
      hist_bkg.GetXaxis().SetTitleOffset(1.1)
      hist_bkg.GetYaxis().SetTitle('Normalised to unity')
      hist_bkg.GetYaxis().SetLabelSize(0.033)
      hist_bkg.GetYaxis().SetTitleSize(0.042)
      hist_bkg.GetYaxis().SetTitleOffset(1.1)
      leg.AddEntry(hist_bkg, 'background')

      # signal
      hists_sig = []
      for mc_sample in mc_samples:
        if mc_sample.mass != mass: continue
        sig_label = 'sig_{}_{}_{}'.format(mc_sample.mass, mc_sample.ctau, category.label)
        selection_sig = self.baseline_selection + ' && ' + category.definition_flat
        #extra_branch = [quantity.name_flat]
        extra_branch = [
          'mu0_pt', 
          'mu_pt',
          'pi_pt',
          'mu0_mu_mass',
          'mu0_pi_mass',
          'deltar_mu0_mu',
          'deltar_mu0_pi',
          'hnl_cos2d',
          'sv_prob',
          'sv_lxysig',
          'pi_dcasig',
          'mu0_pfiso03_rel',
          'mu_pfiso03_rel',
          'pi_numberoftrackerlayers',
          'pi_numberofpixellayers',
          'mu_numberoftrackerlayers',
          'mu_numberofpixellayers',
          'mu0_numberoftrackerlayers',
          'mu0_numberofpixellayers',
          'sv_lxy',
          'mu0_eta',
          'mu_eta',
          'pi_eta',
        ]
        #filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=selection_sig, weights=extra_branch+['gen_hnl_ct', self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=True) 
        filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=selection_sig, weights=extra_branch+['gen_hnl_ct', self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=True) 
        file_sig = self.tools.getRootFile(filename_sig)
        tree_sig = self.tools.getTree(file_sig, 'signal_tree')

        hist_sig = self.tools.createHisto(tree_sig, quantity, hist_name='hist_sig'+str(score), branchname='flat', selection='score > {}'.format(score))
        if hist_sig.Integral() != 0: hist_sig.Scale(1./hist_sig.Integral())

        hist_sig.SetMarkerStyle(20)
        hist_sig.SetMarkerColor(mc_sample.colour)#colours[iscore])
        hists_sig.append(hist_sig)
        leg.AddEntry(hist_sig, 'signal {} GeV {} mm'.format(mc_sample.mass, mc_sample.ctau))

      hist_bkg.Draw('hist')
      for hist in hists_sig:
        hist.Draw('same')
      leg.Draw('same')

      canv.cd()
      outputdir = '{}/distributions'.format(self.outdir)
      if not path.exists(outputdir):
        os.system('mkdir -p {}'.format(outputdir))
      name = '{}/{}_m{}_{}_score{}'.format(outputdir, quantity.label, str(mass).replace('.', 'p'), category.label, str(score).replace('.', 'p'))
      canv.SaveAs(name + '.png')


  def plotDistributionComparison(self, training_info, mc_samples, data_samples, category, quantity, score):
    '''
      Plot the quantity distribution for signal and background
    '''
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.13)

    #masses = [1.0, 2.0, 4.5]
    masses = [2.0]
    #colours = [ROOT.kOrange+0, ROOT.kRed+1, ROOT.kRed+4]
    colours = [ROOT.kBlue+1]
    #for mc_sample in mc_samples:
    #  if mc_sample.mass not in masses:
    #    masses.append(mc_sample.mass)

    canv = self.tools.createTCanvas('canv'+category.label, 800, 700) 
    pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
    pad.Draw()
    pad.cd()
    leg = self.tools.getRootTLegend(xmin=0.4, ymin=0.55, xmax=0.65, ymax=0.67, size=0.04)

    hists_sig = []
    hists_bkg = []

    max_val = -99.

    for imass, mass in enumerate(masses):
      #canv = self.tools.createTCanvas('canv'+category.label+str(mass), 800, 700) 
      #leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)

      resolution = self.resolution_p0 + self.resolution_p1 * mass

      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      bkg_label = 'bkg_{}_{}'.format(mass, category.label)
      extra_branch = [
        'hnl_mass',
        'b_mass',
        'mu0_pt', 
        'mu_pt',
        'pi_pt',
        'mu0_mu_mass',
        'mu0_pi_mass',
        'deltar_mu0_mu',
        'deltar_mu0_pi',
        'hnl_cos2d',
        'sv_prob',
        'sv_lxysig',
        'pi_dcasig',
        'mu0_pfiso03_rel',
        'mu_pfiso03_rel',
        'pi_numberoftrackerlayers',
        'pi_numberofpixellayers',
        'mu_numberoftrackerlayers',
        'mu_numberofpixellayers',
        'mu0_numberoftrackerlayers',
        'mu0_numberofpixellayers',
        'sv_lxy',
        'mu0_eta',
        'mu_eta',
        'pi_eta',
      ]
      selection_bkg = self.baseline_selection + ' && ' + window + ' && ' + category.definition_flat
      filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=selection_bkg, label=bkg_label, treename='signal_tree', weights=extra_branch, force_overwrite=False) 
      file_bkg = self.tools.getRootFile(filename_bkg)
      tree_bkg = self.tools.getTree(file_bkg, 'signal_tree')

      hist_bkg = self.tools.createHisto(tree_bkg, quantity, hist_name='hist_bkg'+str(score), branchname='flat', selection='score > {}'.format(score))
      if hist_bkg.Integral() != 0: hist_bkg.Scale(1./hist_bkg.Integral())

      hist_bkg.SetFillColor(ROOT.kBlue-3)#colours[imass-1])
      hist_bkg.SetFillStyle(3005)

      if hist_bkg.GetMaximum() > max_val: max_val = hist_bkg.GetMaximum()
        
      hist_bkg.SetTitle(' ')
      x_label = quantity.title
      hist_bkg.GetXaxis().SetTitle(x_label)
      hist_bkg.GetXaxis().SetLabelSize(0.04)
      hist_bkg.GetXaxis().SetTitleSize(0.045)
      hist_bkg.GetXaxis().SetLabelOffset(0.015)
      hist_bkg.GetXaxis().SetTitleOffset(1.2)
      hist_bkg.GetYaxis().SetTitle('Normalised to unity')
      hist_bkg.GetYaxis().SetLabelSize(0.04)
      hist_bkg.GetYaxis().SetTitleSize(0.045)
      hist_bkg.GetYaxis().SetLabelOffset(0.01)
      hist_bkg.GetYaxis().SetTitleOffset(1.4)
      leg.AddEntry(hist_bkg, 'data')

      hists_bkg.append(hist_bkg)

      # signal
      #hists_sig = []
      for mc_sample in mc_samples:
        if mc_sample.mass != mass: continue
        sig_label = 'sig_{}_{}_{}'.format(mc_sample.mass, mc_sample.ctau, category.label)
        selection_sig = self.baseline_selection + ' && ' + category.definition_flat
        extra_branch = [
          'mu0_pt', 
          'mu_pt',
          'pi_pt',
          'b_mass',
          'mu0_mu_mass',
          'mu0_pi_mass',
          'deltar_mu0_mu',
          'deltar_mu0_pi',
          'hnl_cos2d',
          'sv_prob',
          'sv_lxysig',
          'pi_dcasig',
          'mu0_pfiso03_rel',
          'mu_pfiso03_rel',
          'pi_numberoftrackerlayers',
          'pi_numberofpixellayers',
          'mu_numberoftrackerlayers',
          'mu_numberofpixellayers',
          'mu0_numberoftrackerlayers',
          'mu0_numberofpixellayers',
          'sv_lxy',
          'mu0_eta',
          'mu_eta',
          'pi_eta',
        ]
        filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=selection_sig, weights=extra_branch+['gen_hnl_ct', self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=False) 
        file_sig = self.tools.getRootFile(filename_sig)
        tree_sig = self.tools.getTree(file_sig, 'signal_tree')

        hist_sig = self.tools.createHisto(tree_sig, quantity, hist_name='hist_sig'+str(score), branchname='flat', selection='score > {}'.format(score))
        if hist_sig.Integral() != 0: hist_sig.Scale(1./hist_sig.Integral())

        if hist_sig.GetMaximum() > max_val: max_val = hist_sig.GetMaximum()

        hist_sig.SetMarkerStyle(20)
        hist_sig.SetMarkerColor(colours[imass])
        hists_sig.append(hist_sig)
        leg.AddEntry(hist_sig, 'signal - {} GeV, {} mm'.format(round(mc_sample.mass, 1), round(mc_sample.ctau, 1)))

    do_log = False
    if quantity.label == 'hnl_cos2d': do_log = True
    if do_log:
      pad.SetLogy()
      #pad.SetLogx()

    for i, hist in enumerate(hists_bkg):
      if i == 0:
        if quantity.label == 'b_mass': 
          hist.GetYaxis().SetRangeUser(0, max_val + 1.*max_val)
        else:
          hist.GetYaxis().SetRangeUser(0, max_val + 0.3*max_val) if not do_log else hist.GetYaxis().SetRangeUser(1e-4, max_val + 30*max_val)
        if quantity.label == 'hnl_cos2d': 
          hist.GetXaxis().SetNdivisions(3)
        hist.Draw('hist')
        hist.Draw('PE same')
      else:
        hist.Draw('hist same')
        hist.Draw('PE same')
    for hist in hists_sig:
      hist.Draw('same')
    leg.Draw('same')

    print_tag = True
    #CMS_tag = 'Preliminary'
    CMS_tag = ''
    self.tools.printInnerCMSTag(pad, CMS_tag, print_tag, x_pos=0.17, y_pos=0.83, size=0.55)
    self.tools.printLatexBox(0.41, 0.84, category.title, size=0.04, pos='left', font=42)
    if 'Bc' in category.label:
      b_mass_label = '#mu{P}#mu^{#pm}#pi^{#mp} mass > 5.7 GeV'
    else:
      b_mass_label = '#mu_{P}#mu^{#pm}#pi^{#mp} mass #leq 5.7 GeV'
    self.tools.printLatexBox(0.41, 0.78, b_mass_label, size=0.04, pos='left', font=42)
    self.tools.printLatexBox(0.41, 0.73, 'dimuon channel', size=0.04, pos='left', font=42)
    self.tools.printLumiTag(pad, 5.2, size=0.5, offset=0.52)

    canv.cd()
    outputdir = '{}/pNN_features'.format(self.outdir)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    name = '{}/{}_{}_score{}'.format(outputdir, quantity.label, category.label, str(score).replace('.', 'p'))
    canv.SaveAs(name + '.png')
    canv.SaveAs(name + '.pdf')
    canv.SaveAs(name + '.C')


  def compareROCCurve(self, mc_samples, data_samples, category, do_log=False):
    pd.options.mode.chained_assignment = None

    masses = [1.0, 2.0, 3.0, 4.5, 1.5, 5.5]

    #training_info1 = TrainingInfo('V13_06Feb23_2023Apr06_14h13m31s', category.label) # used for unblinding
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jun11_22h15m20s', category.label) # track id
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jun15_23h04m36s', category.label) # displaced mu id
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jun28_12h19m48s', category.label) # pi dcasig, no sv chi2
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jun29_10h53m14s', category.label) # mu0 id
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul04_19h20m03s', category.label) # no b pt/mass, no deltaR
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul04_19h21m47s', category.label) # no b pt/mass, no invmass
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul04_21h09m09s', category.label) # no b pt/mass, no layers
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul04_22h10m04s', category.label) # modified patience
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul05_16h43m39s', category.label) # with weight
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul05_17h29m50s', category.label) # layers instead of hits, 32 nodes
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul05_21h30m59s', category.label) # scaling weight
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul05_21h30m59s', category.label) # scaling weight
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul12_11h17m28s', category.label) # weight 1
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul12_11h17m59s', category.label) # trigger weight
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul14_16h33m54s', category.label) # on track info on mu, early stopping
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul14_16h33m54s', category.label) # early stopping, patience 5, period 100
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul17_21h35m41s', category.label) # early stopping, patience 10, period 100
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul17_23h01m06s', category.label) # no early stopping, period 100
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul17_23h02m14s', category.label) # no early stopping, period 60
    #training_info1 = TrainingInfo('V13_06Feb23_2023Jul18_15h36m17s', category.label) # adding back mu layers (used for LLPJul23 talk)
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul24_21h21m38s', category.label) # more Bc stat
    training_info1 = TrainingInfo('training_Jul23', category.label) # mix of the two above
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul27_22h29m08s', category.label) # removing mu iso
    #training_info2 = TrainingInfo('V13_06Feb23_2023Jul31_11h57m54s', category.label) # removing mu iso v2
    #training_info2 = TrainingInfo('V13_06Feb23_2023Aug04_20h56m51s', category.label) # no iso, pi_dcasig_corr
    #training_info2 = TrainingInfo('V13_06Feb23_2023Aug05_00h06m25s', category.label) # pi_dcasig_corr
    #training_info2 = TrainingInfo('V13_06Feb23_2023Aug05_00h06m44s', category.label) # pi_dcasig_corr
    training_info2 = TrainingInfo('V13_06Feb23_2023Aug05_00h07m05s', category.label) # pi_dcasig_corr

    #leg1 = 'training w/o track ID'
    #leg2 = 'training w/ track ID'

    #leg1 = 'training w/ track ID'
    #leg2 = 'training w/ track + muon ID'

    #leg1 = 'training w/o ID'
    #leg2 = 'training w/ track + muon ID'

    #leg1 = 'training w/o pi dcasig'
    #leg2 = 'training w/ pi dcasig'

    #leg1 = 'training w/o mu0 ID'
    #leg2 = 'training w/ mu0 ID'

    #leg1 = 'Previous training'
    #leg2 = 'Updated training'

    leg1 = 'Reference training'
    leg2 = 'Updated training'

    #leg1 = 'Early stopping patience 5'
    #leg2 = 'Early stopping patience 10'

    #leg1 = 'Early stopping'
    #leg2 = 'No early stopping'

    #leg1 = 'No early stopping, period 100'
    #leg2 = 'No early stopping, period 60'

    #leg1 = 'w/o weight'
    #leg2 = 'w/ weight'

    for mass in masses:

      resolution = self.resolution_p0 + self.resolution_p1 * mass

      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      # create dataframe
      data_df = self.createDataframe(data_samples).query(self.getPandasQuery(window))
      if self.do_parametric:
        data_df['mass_key'] = mass
      data_df['is_signal'] = 0

      plt.clf()
      for mc_sample in mc_samples:
        if mc_sample.mass != mass: continue
        ctau = mc_sample.ctau
        v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
        coupling = self.tools.getCouplingLabel(v2)

        mc_df = self.createDataframe([mc_sample])
        if self.do_parametric:
          mc_df['mass_key'] = mass
        mc_df['is_signal'] = 1

        main_df = pd.concat([data_df, mc_df], sort=False)
        main_df.index = np.array(range(len(main_df)))
        main_df = main_df.sample(frac=1, replace=False, random_state=1986)

        Y = pd.DataFrame(main_df, columns=['is_signal'])

        score1 = self.predictScore(training_info=training_info1, df=main_df)
        fpr1, tpr1, thresholds1 = roc_curve(Y, score1) 
        idx1 = 0
        for idx, thrs in enumerate(thresholds1):
          if thrs >= 0.99: idx1 = idx
          if thrs < 0.99: break

        score2 = self.predictScore(training_info=training_info2, df=main_df)
        fpr2, tpr2, thresholds2 = roc_curve(Y, score2) 
        idx2 = 0
        for idx, thrs in enumerate(thresholds2):
          if thrs >= 0.99: idx2 = idx
          if thrs < 0.99: break

        plt.plot(fpr1, tpr1, linewidth=2, label=leg1)
        plt.plot(fpr2, tpr2, linewidth=2, label=leg2)
        plt.plot(fpr1[idx1], tpr1[idx1], '*', markersize=10)
        plt.plot(fpr2[idx2], tpr2[idx2], '*', markersize=10)

        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

        xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
        plt.plot(xy, xy, color='grey', linestyle='--')

        plt.title(r'{}, {} GeV, {} mm'.format(category.label, mass, ctau))
        if do_log:
          plt.xlim(1e-5, 1)
        else:
          plt.xlim(0, 1)
        plt.ylim(0, 1)

        if do_log:
          plt.xscale('log')
        plt.yscale('linear')

        if do_log:
          plt.legend(loc='upper left', framealpha=0.1)
        else:
          plt.legend(loc='lower right', framealpha=0.1)

      name = 'comparison_ROC_m{}_{}'.format(str(mass).replace('.', 'p'), category.label)
      if do_log: name += '_log'
      self.saveFig(plt, name)


  def studyDisplacedTracks(self, training_info, mc_samples, data_samples, category):
    '''
      Plot the quantity distribution for signal and background
    '''
    pd.options.mode.chained_assignment = None
    ROOT.gStyle.SetOptStat(0)

    masses = [1.0, 1.5, 2.0, 3.0, 4.5]

    for mass in masses:
      print 'mass {}'.format(mass)

      canv = self.tools.createTCanvas('canv'+category.label+str(mass), 800, 700) 
      leg = self.tools.getRootTLegend(xmin=0.2, ymin=0.65, xmax=0.65, ymax=0.83, size=0.04)
      
      score = 0.99
      colours = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+3]

      resolution = self.resolution_p0 + self.resolution_p1 * mass
      quantity = Qte(name_nano='BToMuMuPi_hnl_mass', name_flat='hnl_mass', label='hnl_mass', title='#mu#pi invariant mass [GeV]', nbins=80, bin_min=mass-10*resolution, bin_max=mass+10*resolution)

      # consider the 10 sigma window around the signal mass
      window = 'hnl_mass > {} && hnl_mass < {}'.format(mass-10*resolution, mass+10*resolution)
      bkg_label = 'bkg_{}_{}'.format(mass, category.label)
      extra_branch = [
        'mu0_pt', 
        'mu_pt',
        'pi_pt',
        'mu0_mu_mass',
        'mu0_pi_mass',
        'deltar_mu0_mu',
        'deltar_mu0_pi',
        'hnl_cos2d',
        'sv_prob',
        'sv_lxysig',
        'pi_dcasig',
        'mu0_pfiso03_rel',
        'mu_pfiso03_rel',
        'pi_numberoftrackerlayers',
        'pi_numberofpixellayers',
        'mu_numberoftrackerlayers',
        'mu_numberofpixellayers',
        'mu0_numberoftrackerlayers',
        'mu0_numberofpixellayers',
        'sv_lxy',
      ]
      selection_bkg = self.baseline_selection + ' && ' + window + ' && ' + category.definition_flat
      filename_bkg = self.mva_tools.getFileWithScore(files=data_samples, training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mass, selection=selection_bkg, label=bkg_label, treename='signal_tree', force_overwrite=False, weights=extra_branch) 
      file_bkg = self.tools.getRootFile(filename_bkg)
      tree_bkg = self.tools.getTree(file_bkg, 'signal_tree')

      selection_bkg_pNN_1 = 'score > {}'.format(score)
      selection_bkg_pNN_2 = '((sv_lxy>18.0 && pi_numberofpixellayers<=0) || sv_lxy<=18) && ((sv_lxy>13.0 && pi_numberofpixellayers<=1) || sv_lxy<=13) && ((sv_lxy>8.0 && pi_numberofpixellayers<=2) || sv_lxy<=8) && ((sv_lxy>5.0 && pi_numberofpixellayers<=3) || sv_lxy<=5) && score > {}'.format(score)

      hist_bkg_1 = self.tools.createHisto(tree_bkg, quantity, hist_name='hist_bkg'+str(score), branchname='flat', selection=selection_bkg_pNN_1)
      hist_bkg_2 = self.tools.createHisto(tree_bkg, quantity, hist_name='hist_bkg'+str(score), branchname='flat', selection=selection_bkg_pNN_2)

      print 'bkg yields with track sel / without track sel = {} / {} = {}'.format(hist_bkg_2.Integral(), hist_bkg_1.Integral(), hist_bkg_2.Integral()/hist_bkg_1.Integral())

      # signal
      hists_sig = []
      for mc_sample in mc_samples:
        if mc_sample.mass != mass: continue
        sig_label = 'sig_{}_{}_{}'.format(mc_sample.mass, mc_sample.ctau, category.label)
        selection_sig = self.baseline_selection + ' && ' + category.definition_flat
        extra_branch = [
          'mu0_pt', 
          'mu_pt',
          'pi_pt',
          'mu0_mu_mass',
          'mu0_pi_mass',
          'deltar_mu0_mu',
          'deltar_mu0_pi',
          'hnl_cos2d',
          'sv_prob',
          'sv_lxysig',
          'pi_dcasig',
          'mu0_pfiso03_rel',
          'mu_pfiso03_rel',
          'pi_numberoftrackerlayers',
          'pi_numberofpixellayers',
          'mu_numberoftrackerlayers',
          'mu_numberofpixellayers',
          'mu0_numberoftrackerlayers',
          'mu0_numberofpixellayers',
          'sv_lxy',
        ]
        filename_sig = self.mva_tools.getFileWithScore(files=[mc_sample], training_label='./outputs/'+self.dirname, category_label=category.label, do_parametric=self.do_parametric, mass=mc_sample.mass, selection=selection_sig, weights=extra_branch+['gen_hnl_ct', self.weight_pusig, self.weight_mu0id, self.weight_muid], label=sig_label, treename='signal_tree', force_overwrite=True) 
        file_sig = self.tools.getRootFile(filename_sig)
        tree_sig = self.tools.getTree(file_sig, 'signal_tree')

        selection_sig_pNN_1 = 'score > {}'.format(score)
        selection_sig_pNN_2 = '((sv_lxy>18.0 && pi_numberofpixellayers<=0) || sv_lxy<=18) && ((sv_lxy>13.0 && pi_numberofpixellayers<=1) || sv_lxy<=13) && ((sv_lxy>8.0 && pi_numberofpixellayers<=2) || sv_lxy<=8) && ((sv_lxy>5.0 && pi_numberofpixellayers<=3) || sv_lxy<=5) && score > {}'.format(score)

        hist_sig_1 = self.tools.createHisto(tree_sig, quantity, hist_name='hist_sig'+str(score), branchname='flat', selection=selection_sig_pNN_1)
        hist_sig_2 = self.tools.createHisto(tree_sig, quantity, hist_name='hist_sig'+str(score), branchname='flat', selection=selection_sig_pNN_2)

        print 'sig yields with track sel / without track sel = {} / {} = {}'.format(hist_sig_2.Integral(), hist_sig_1.Integral(), hist_sig_2.Integral()/hist_sig_1.Integral())


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

    for category in self.categories:
      if category.label == 'incl': continue
      #if category.label != 'lxysiggt150_SS' and category.label != 'lxysig50to150_SS': continue
      #if category.label != 'lxysig0to50_OS' and category.label != 'lxysig50to150_SS': continue
      if category.label != 'lxysiggt150_OS': continue# and category.label != 'lxysig0to50_SS': continue
      #if 'Bc' in category.label: continue

      print '\n -> get the training information'
      training_info = TrainingInfo(self.dirname, category.label)

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

      if self.do_plotScoreCurve:
        self.plotScoreCurve(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)

      if self.do_plotMVAPerformance:
        self.plotMVAPerformance(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)

      if self.do_plotWPScan:
        self.plotWPScan(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)

      if self.do_plotAUC:
        self.plotAUCGraph(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category)

      if self.do_plotAUCvsLifetime:
        self.plotAUCvsLifetimeGraph(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category)

      if self.do_plotPNNComparison:
        self.plotPNNComparisonGraph(mc_samples=mc_samples, data_samples=data_samples, category=category)

      if self.do_createFiles:
        print '\n -> create rootfiles'
        self.createRootFile(training_info, mc_samples, 'sample_mc')
        self.createRootFile(training_info, data_samples, 'sample_data')
       
        self.createFriendTree(training_info, mc_samples, 'score_mc')
        self.createFriendTree(training_info, data_samples, 'score_data')

        #self.compareTrees()
  
      if self.do_plotSigScan:
        self.plotSignificanceDiffScan(category=category, mc_samples=mc_samples, data_samples=data_samples, n_points=30)

      if self.do_plotSignalBackgroundComparison:
        for quantity in quantities_pNN:
          #self.plotSignalBackgroundComparison(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, quantity=quantity, category=category, score=0.0)
          self.plotSignalBackgroundComparison(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, quantity=quantity, category=category, score=0.99)

      if self.do_plotDistributionComparison:
        for quantity in quantities_pNN:
          self.plotDistributionComparison(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, quantity=quantity, category=category, score=0.0)

      if self.do_plotMass:
        self.plotMass(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, object_mass='hnl_mass')
        self.plotMass(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, object_mass='b_mass')

      if self.do_plotPreselection:
        #quantitites = quantities['preselection']
        #for quantity in quantities_preselection:
        for quantity in quantities_trackid:
          print quantity.name_flat
          self.plotPreselectionQuantity(training_info=training_info, mc_samples=mc_samples, quantity=quantity, category=category)

      if self.do_compareROC:
        if category.label == 'incl': continue
        self.compareROCCurve(mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)
        self.compareROCCurve(mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=True)

      if self.do_studyDisplacedTracks:
        self.studyDisplacedTracks(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category)

      if self.do_plotScoreNorm:
        print '\n -> plot the score'
        self.plotScoreNorm(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=False)
        self.plotScoreNorm(training_info=training_info, mc_samples=mc_samples, data_samples=data_samples, category=category, do_log=True)





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
  #dirname = 'test_09Sep2022_00h12m12s' # mass 4.5
  #dirname = 'test_12Sep2022_14h47m55s' # mass 4.5, whnl nn
  #dirname = 'test_12Sep2022_22h23m42s' # mass 3, categories
  #dirname = 'test_13Sep2022_10h23m24s' # mass 3, categories, mass window
  #dirname = 'test_13Sep2022_21h38m47s'
  #dirname = 'test_14Sep2022_10h06m37s' # mass 3, non parametric
  #dirname = 'test_14Sep2022_10h08m03s' # mass 3, parametric extended to full spectrum
  #dirname = 'test_14Sep2022_11h24m19s' # trained on masses 1.5, 3
  #dirname = 'test_14Sep2022_11h40m11s' # trained on masses 1, 1.5, 2, 3, 4.5
  #dirname = 'test_14Sep2022_14h26m36s' # as above, scaled key
  #dirname = 'test_14Sep2022_14h56m50s' # 30 nsigma
  #dirname = 'test_14Sep2022_20h33m22s' # 10 nsigma, different features with ctau and mass dependency, more balanced training
  #dirname = 'test_15Sep2022_00h37m43s' # balanced
  #dirname = 'test_15Sep2022_01h01m18s' # mixed ctau signal
  #dirname = 'test_19Sep2022_16h33m33s' # mixed ctau signal, all samples
  #dirname = 'test_19Sep2022_17h05m13s' # all masses apart from mass 3
  #dirname = 'test_21Sep2022_13h36m30s' # parametric as above trained on mass 3 only
  #dirname = 'test_2022Sep21_15h58m10s' # all categories, but lxygt5_SS
  #dirname = 'test_2022Sep21_17h43m05s' # category lxygt5_SS
  #dirname = 'test_2022Sep21_17h54m44s' # all categories
  #dirname = 'test_2022Sep27_10h58m03s' # 0_150_300
  #dirname = 'test_2022Sep27_16h35m02s' # 0 50 150
  #dirname = 'test_2022Sep29_09h11m21s' # without pi dcasig
  #dirname = 'test_2022Oct06_17h34m15s' # new feature set
  #dirname = 'test_2022Oct10_16h40m58s' # invmass and deltaR added
  #dirname = 'test_2022Oct12_15h12m37s' # mass 3 only
  #dirname = 'test_2022Nov29_09h26m28s' # adding muon isolation
  #dirname = 'test_2023Jan25_12h38m31s' # test statistics
  #dirname = 'test_Bc_2023Jan25_22h42m38s' # Bc
  #dirname = 'test_2023Jan30_17h10m34s' # test statistics ML review
  #dirname = 'V13_06Feb23_2023Apr06_14h13m31s' # used for unblinding
  #dirname = 'V13_06Feb23_2023Jun11_22h15m20s' # track id
  #dirname = 'V13_06Feb23_2023Jun15_23h04m36s' # displaced mu id
  #dirname = 'V13_06Feb23_2023Jun28_12h19m48s' #pi dcasig, no sv_chi2
  #dirname = 'V13_06Feb23_2023Jun29_10h53m14s' # with mu0 id
  #dirname = 'V13_06Feb23_2023Jul05_16h43m39s' # with weight
  #dirname = 'V13_06Feb23_2023Jul05_22h51m32s' # layers instead of hits, no weight
  #dirname = 'V13_06Feb23_2023Jul11_08h46m17s' # weight 1
  #dirname = 'V13_06Feb23_2023Jul12_09h53m40s' # test
  #dirname = 'V13_06Feb23_2023Jul12_11h17m28s' # weight 1
  #dirname = 'V13_06Feb23_2023Jul12_11h17m59s' # trigger weight
  #dirname = 'V13_06Feb23_2023Jul14_16h33m54s' # early stopping, patience 5
  #dirname = 'V13_06Feb23_2023Jul17_21h35m41s' # early stopping, patience 10
  #dirname = 'V13_06Feb23_2023Jul17_23h01m06s' # no early stopping, period 100
  #dirname = 'V13_06Feb23_2023Jul17_23h02m14s' # no early stopping, period 60
  #dirname = 'V13_06Feb23_2023Jul18_15h36m17s' # adding back mu layers
  #dirname = 'V13_06Feb23_2023Jul24_21h21m38s' # Bc stat
  #dirname = 'training_Jul23' # mix of the two previous ones (Bc + gt150 taken from previous one)
  #dirname = 'V13_06Feb23_2023Aug04_20h56m51s' # no iso, pi_dcasig_corr
  #dirname = 'V13_06Feb23_2023Aug05_00h06m25s' # pi_dcasig_corr
  #dirname = 'V13_06Feb23_2023Aug05_00h06m44s' # pi_dcasig_corr
  #dirname = 'V13_06Feb23_2023Aug05_00h07m05s' # pi_dcasig_corr
  dirname = 'training_Aug23' # mix of three previous ones

  baseline_selection = selection['baseline_06Feb23'].flat + ' && hnl_charge==0'
  #categories = categories['V12_08Aug22_permass']
  categories = categories['categories_0_50_150_Bc']

  do_parametric = True

  data_files = []
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_partial.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_norm.root')
  data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_31Jul23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root')

  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_06Feb23_31Jul23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk1_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk2_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk3_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk4_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk5_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk6_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk7_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk8_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk9_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk10_n500/flat/flat_bparknano_06Feb23_15Jun23.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk11_n500/flat/flat_bparknano_06Feb23_15Jun23.root')

  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk1_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk2_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk3_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk4_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk5_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk6_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk7_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk8_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk9_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk10_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk11_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk12_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk13_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk14_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk15_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk16_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk17_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk18_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk19_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk20_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk21_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk22_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk23_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk24_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk25_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk26_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk27_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk28_n500/flat/flat_bparknano_08Aug22_sr.root')
  #data_files.append('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk29_n500/flat/flat_bparknano_08Aug22_sr.root')

  do_analyseMVA = False    # assess performance of mva
  do_compareMVA = False   # compare mva performance to that of the cutbased method

  if do_analyseMVA:
    #signal_labels = ['V12_08Aug22_m1', 'V12_08Aug22_m1p5', 'V12_08Aug22_m2', 'V12_08Aug22_m3', 'V12_08Aug22_m4p5']
    #signal_labels = ['V12_08Aug22_m1', 'V12_08Aug22_m1p5', 'V12_08Aug22_m2', 'V12_08Aug22_m3', 'V12_08Aug22_m4p5']
    signal_labels = ['V13_06Feb23_m1', 'V13_06Feb23_m1p5', 'V13_06Feb23_m2', 'V13_06Feb23_m3', 'V13_06Feb23_m4p5']
    #masses = ['m1', 'm1p5', 'm2', 'm3', 'm4p5']
    masses = ['m2']

    for mass in masses:
      signal_files = []
      for signal_label in signal_labels:
        if mass not in signal_label: continue
        if mass == 'm1' and 'm1p5' in signal_label: continue
        for signal_file in signal_samples[signal_label]:
          #print signal_file.filename
          signal_files.append(signal_file)

      analyser = MVAAnalyser(
          signal_files = signal_files,
          data_files = data_files,
          dirname = dirname,
          baseline_selection = baseline_selection,
          categories = categories,
          do_parametric = do_parametric,
          do_addCutbased = False,
          do_plotScore = True,
          do_createFiles = False,
          do_plotSigScan = False,
          do_plotROC = True,
          do_plotScoreCurve = False,
          do_plotMVAPerformance = False,
          do_plotWPScan = False,
          do_plotAUC = False,
          do_plotMass = True,
          do_plotPreselection = False,
          do_plotSignalBackgroundComparison = False,
          do_compareROC = False,
          )

      analyser.process()

  elif do_compareMVA:
    signal_labels = ['V12_08Aug22_m3']

    signal_files = []
    for signal_label in signal_labels:
      for signal_file in signal_samples[signal_label]:
        signal_files.append(signal_file)


    analyser = MVAAnalyser(
        signal_files = signal_files,
        data_files = data_files,
        dirname = dirname,
        baseline_selection = baseline_selection,
        categories = categories,
        do_parametric = do_parametric,
        do_addCutbased = True,
        do_plotScore = False,
        do_createFiles = False,
        do_plotSigScan = False,
        do_plotROC = True,
        do_plotScoreCurve = False,
        do_plotMVAPerformance = True,
        do_plotWPScan = False,
        do_plotAUC = False,
        do_plotMass = False,
        do_plotPreselection = False,
        do_plotSignalBackgroundComparison = False,
        do_compareROC = False,
        )

    analyser.process()

  else:

    do_plotScore = False
    do_createFiles = False
    do_plotSigScan = False
    do_plotScoreCurve = False
    do_plotMVAPerformance = False
    do_plotWPScan = False
    do_plotROC = False
    do_plotAUC = False
    do_plotMass = False
    do_plotAUCvsLifetime = False
    do_plotPNNComparison = True
    do_plotPreselection = False
    do_plotSignalBackgroundComparison = False
    do_plotDistributionComparison = False
    do_compareROC = False
    do_studyDisplacedTracks = False
    do_plotScoreNorm = False

    #signal_labels = ['V12_08Aug22_m1', 'V12_08Aug22_m1p5', 'V12_08Aug22_m2', 'V12_08Aug22_m3', 'V12_08Aug22_m4p5']
    #signal_labels = ['V13_06Feb23_m1', 'V13_06Feb23_m1p5', 'V13_06Feb23_m2', 'V13_06Feb23_m3', 'V13_06Feb23_m4p5']
    signal_labels = ['V13_06Feb23_m1', 'V13_06Feb23_m1p5', 'V13_06Feb23_m2', 'V13_06Feb23_m3']
    #signal_labels = ['V12_08Aug22_sensitivity']
    #signal_labels = ['V13_06Feb23_trackid']

    #signal_labels = [
    #  #'V42_08Aug22_m0p5',
    #  'V42_08Aug22_m0p6',
    #  'V42_08Aug22_m0p7',
    #  'V42_08Aug22_m0p8',
    #  'V42_08Aug22_m0p9',
    #  'V42_08Aug22_m1p02',
    #  'V42_08Aug22_m1p04',
    #  'V42_08Aug22_m1p06',
    #  'V42_08Aug22_m1p08',
    #  'V42_08Aug22_m1p1',
    #  'V42_08Aug22_m1p12',
    #  #'V42_08Aug22_m1p14',
    #  'V42_08Aug22_m1p16',
    #  'V42_08Aug22_m1p18',
    #  'V42_08Aug22_m1p2',
    #  'V42_08Aug22_m1p22',
    #  'V42_08Aug22_m1p24',
    #  'V42_08Aug22_m1p26',
    #  #'V42_08Aug22_m1p28',
    #  'V42_08Aug22_m1p3',
    #  'V42_08Aug22_m1p32',
    #  'V42_08Aug22_m1p34',
    #  'V42_08Aug22_m1p36',
    #  'V42_08Aug22_m1p38',
    #  'V42_08Aug22_m1p4',
    #  'V42_08Aug22_m1p42',
    #  'V42_08Aug22_m1p44',
    #  'V42_08Aug22_m1p46',
    #  'V42_08Aug22_m1p48',
    #  'V42_08Aug22_m1p5',
    #  'V42_08Aug22_m1p53',
    #  'V42_08Aug22_m1p56',
    #  'V42_08Aug22_m1p59',
    #  #'V42_08Aug22_m1p62',
    #  'V42_08Aug22_m1p65',
    #  #'V42_08Aug22_m1p68',
    #  'V42_08Aug22_m1p71',
    #  'V42_08Aug22_m1p74',
    #  'V42_08Aug22_m1p77',
    #  'V42_08Aug22_m1p8',
    #  'V42_08Aug22_m1p83',
    #  #'V42_08Aug22_m1p86',
    #  #'V42_08Aug22_m1p89',
    #  'V42_08Aug22_m1p92',
    #  'V42_08Aug22_m1p95',
    #  #'V42_08Aug22_m1p98',
    #  #'V42_08Aug22_m2p05',
    #  'V42_08Aug22_m2p0',
    #  'V42_08Aug22_m2p1',
    #  'V42_08Aug22_m2p15',
    #  #'V42_08Aug22_m2p2',
    #  'V42_08Aug22_m2p25',
    #  'V42_08Aug22_m2p3',
    #  'V42_08Aug22_m2p35',
    #  #'V42_08Aug22_m2p4',
    #  'V42_08Aug22_m2p45',
    #  'V42_08Aug22_m2p5',
    #  'V42_08Aug22_m2p55',
    #  'V42_08Aug22_m2p6',
    #  'V42_08Aug22_m2p65',
    #  'V42_08Aug22_m2p7',
    #  #'V42_08Aug22_m2p75',
    #  'V42_08Aug22_m2p8',
    #  'V42_08Aug22_m2p85',
    #  'V42_08Aug22_m2p9',
    #  'V42_08Aug22_m2p95',
    #  'V42_08Aug22_m3p0',
    #  'V42_08Aug22_m3p05',
    #  'V42_08Aug22_m3p1',
    #  'V42_08Aug22_m3p15',
    #  'V42_08Aug22_m3p2',
    #  'V42_08Aug22_m3p25',
    #  'V42_08Aug22_m3p3',
    #  'V42_08Aug22_m3p35',
    #  'V42_08Aug22_m3p4',
    #  'V42_08Aug22_m3p45',
    #  'V42_08Aug22_m3p5',
    #  'V42_08Aug22_m3p55',
    #  'V42_08Aug22_m3p6',
    #  'V42_08Aug22_m3p65',
    #  'V42_08Aug22_m3p7',
    #  'V42_08Aug22_m3p75',
    #  'V42_08Aug22_m3p8',
    #  'V42_08Aug22_m3p85',
    #  'V42_08Aug22_m3p9',
    #  'V42_08Aug22_m3p95',
    #  'V42_08Aug22_m4p0',
    #  'V42_08Aug22_m4p1',
    #  'V42_08Aug22_m4p2',
    #  #'V42_08Aug22_m4p3',
    #  #'V42_08Aug22_m4p4',
    #  'V42_08Aug22_m4p5',
    #  'V42_08Aug22_m4p6',
    #  #'V42_08Aug22_m4p7',
    #  #'V42_08Aug22_m4p8',
    #]

    signal_files = []
    for signal_label in signal_labels:
      for signal_file in signal_samples[signal_label]:
        signal_files.append(signal_file)

    analyser = MVAAnalyser(
        signal_files = signal_files,
        data_files = data_files,
        dirname = dirname,
        baseline_selection = baseline_selection,
        categories = categories,
        do_parametric = do_parametric,
        do_plotScore = do_plotScore,
        do_createFiles = do_createFiles,
        do_plotSigScan = do_plotSigScan,
        do_plotROC = do_plotROC,
        do_plotScoreCurve = do_plotScoreCurve,
        do_plotMVAPerformance = do_plotMVAPerformance,
        do_plotWPScan = do_plotWPScan,
        do_plotAUC = do_plotAUC,
        do_plotMass = do_plotMass,
        do_plotAUCvsLifetime = do_plotAUCvsLifetime,
        do_plotPNNComparison = do_plotPNNComparison,
        do_plotPreselection = do_plotPreselection,
        do_plotSignalBackgroundComparison = do_plotSignalBackgroundComparison,
        do_plotDistributionComparison = do_plotDistributionComparison,
        do_compareROC = do_compareROC,
        do_studyDisplacedTracks = do_studyDisplacedTracks,
        do_plotScoreNorm = do_plotScoreNorm,
        )

    analyser.process()





