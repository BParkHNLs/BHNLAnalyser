import os
import sys
from os import path

import ROOT

from root_pandas import read_root

import pickle
import numpy as np
import pandas as pd
from array import array

from keras.models import load_model


class TrainingInfo(object):
  '''
    Class that contains the training information:
      - model
      - scaler
      - features
  '''
  def __init__(self, training_label):
    self.training_label = training_label
    #self.indir = './mva/outputs/{}'.format(self.training_label)
    self.indir = '/t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNLAnalyser/scripts/mva/outputs/{}'.format(self.training_label) #FIXME
    self.model = self.loadModel()
    self.qt = self.loadScaler()
    self.features = self.loadFeatures()

    if not path.exists(self.indir):
      raise RuntimeError('Training set with label "{}" was not to be found in "{}". Please check training label'.format(self.training_label, self.indir))


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



class MVATools(object):
  '''
    This class contains tools to allow the mva selection
  '''

  def convertRootToDF(self, sample, training_info, signal_treename, selection):
    '''
      Function that converts root samples to pandas dataframe
    '''
    #TODO make sure that the selection is readable by pandas query
    #df = read_root(sample, signal_treename, where=selection, warn_missing_tree=True)
    df = read_root(sample, signal_treename, where=selection, warn_missing_tree=True, columns=training_info.features)

    return df


  def createDataframe(self, training_info, samples_filename, selection, signal_treename):
    '''
      Function that returns a dataframe out of a list of samples
    '''
    df = pd.concat([self.convertRootToDF(idt, training_info, signal_treename, selection) for idt in samples_filename], sort=False)

    return df


  def getTrainingInfo(self, training_label):
    '''
      Get training information
    '''
    training_info = TrainingInfo(training_label=training_label)

    return training_info


  #def predictScore(self, training_label, df):
  def predictScore(self, training_info, df):
    '''
      Return score with scaled input features
    '''

    #training_info = self.getTrainingInfo(training_label)

    x = pd.DataFrame(df, columns=training_info.features) #TODO think of selection

    # apply the scaler
    xx = training_info.qt.transform(x[training_info.features])

    # predict
    score = training_info.model.predict(xx)

    return score


  def saveTreeInMemory(self, tree):
    '''
      Only way that was found so far to deal with PyROOT memory issues...
      tree.SetDirectory(0) alone would not work
    '''
    c = ROOT.TCanvas()
    tree.Draw('score')
    tree.SetDirectory(0)


  def createScoreTree(self, training_label, samples_filename, selection, label, signal_treename, score_treename):
    '''
      Create tree that contains the score and store it in a root file
    '''
    # get the training information
    training_info = self.getTrainingInfo(training_label)

    # get the score
    df = self.createDataframe(training_info=training_info, samples_filename=samples_filename, selection=selection, signal_treename=signal_treename)
    #score = self.predictScore(training_label=training_label, df=df) 
    score = self.predictScore(training_info=training_info, df=df) 

    # create file
    root_filename = './{}.root'.format(label.replace('.', 'p'))
    friend_file = ROOT.TFile(root_filename, 'RECREATE')

    # create tree
    tree = ROOT.TTree(score_treename, score_treename)
    the_score = array('d',[0])
    tree.Branch('score', the_score, 'score/D')
    for entry in range(len(score)):
      the_score[0] = score[entry][0]
      tree.Fill()
    tree.Write()
    friend_file.Close()

    print ' --> {} created'.format(root_filename)

    
  def getScoreTree(self, training_label, samples_filename, selection, label, signal_treename, score_treename):
    '''
      Function that returns the score tree
    '''
    root_filename = './{}.root'.format(label.replace('.', 'p')) #TODO is this correct? also for signal?
    if not path.exists(root_filename):
      self.createScoreTree(training_label=training_label, samples_filename=samples_filename, selection=selection, label=label, signal_treename=signal_treename, score_treename=score_treename)
    f = ROOT.TFile.Open(root_filename)
    tree = f.Get(score_treename)
    self.saveTreeInMemory(tree)

    return tree


  def getTreeWithScore(self, files=None, training_label='', selection='hnl_charge>-99', label='', signal_treename='signal_tree', score_treename='score_tree'): 
    '''
      This function returns the signal tree with the score added via a friend tree
    '''
    #FIXME how to deal with selection
    #TODO implement different treatment between data and signal
    #NOTE should we first predict the score and then apply the selection or vice versa? make a check. This may be particularly important if we have 
    # per category training

    # get signal tree
    samples_filename = []
    signal_tree = ROOT.TChain(signal_treename)
    for file_ in files:
      signal_tree.Add(file_.filename) 
      samples_filename.append(file_.filename)

    # get the score tree 
    score_tree = self.getScoreTree(samples_filename=samples_filename, training_label=training_label, selection=selection, label=label, signal_treename=signal_treename, score_treename=score_treename)
    self.saveTreeInMemory(score_tree)

    #root_filename = './{}.root'.format(label.replace('.', 'p')) #TODO is this correct? also for signal?
    #if not path.exists(root_filename):
    #  self.createScoreTree(training_label=training_label, samples_filename=samples_filename, selection=selection, label=label, signal_treename=signal_treename, score_treename=score_treename)
    #f = ROOT.TFile.Open(root_filename)
    #score_tree = f.Get(score_treename)
    #self.saveTreeInMemory(score_tree)

    if(signal_tree.GetEntries() != score_tree.GetEntries()):
      raise RuntimeError('[MVA Tools] [getTreeWithScore] The signal and score trees have different number of entries. Please check') 

    # add score tree to signal tree as a friend
    signal_tree.AddFriend(score_tree)

    #signal_tree.SetDirectory(0)
    self.saveTreeInMemory(signal_tree)

    return signal_tree


