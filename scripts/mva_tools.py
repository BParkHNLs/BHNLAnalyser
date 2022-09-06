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
    self.indir = './{}'.format(self.training_label)
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

  def getSelection(self, selection):
    '''
      Function that removes the cut on the score if any from the selection string
    '''
    if 'score' not in selection:
      selection_string = selection
    else:
      idx = selection.find('score')
      idx_in = selection.rfind('&&')
      if idx_in < idx:
        selection_string = selection[:idx_in-1]
      else:
        selection_string = selection[:idx-1] + selection[idx_in+2:]
    
    return selection_string


  def getScoreSelection(self, selection):
    '''
      Function that only keeps the the cut on the score
    '''
    if 'score' not in selection:
      raise RuntimeError('Please specify which cut on the score to apply when running with the mva method')
    else:
      idx = selection.find('score')
      idx_in = selection.rfind('&&')
      if idx_in < idx:
        selection_string = selection[idx:]
      else:
        selection_string = selection[idx:idx_in-1]
    
    return selection_string


  def convertRootToDF(self, sample, training_info, signal_treename, selection, weights):
    '''
      Function that converts root samples to pandas dataframe
      The only kept branches are the hnl mass, score and signal weights
    '''
    if weights == None:
      extra_columns = ['hnl_mass']
    else:
      extra_columns = ['hnl_mass'] + weights

    df = read_root(sample, signal_treename, where=self.getSelection(selection), warn_missing_tree=True, columns=training_info.features+extra_columns)

    return df


  def createDataframe(self, training_info, samples_filename, selection, weights, signal_treename):
    '''
      Function that returns a dataframe out of a list of samples
    '''
    df = pd.concat([self.convertRootToDF(idt, training_info, signal_treename, selection, weights) for idt in samples_filename], sort=False)

    return df


  def getTrainingInfo(self, training_label):
    '''
      Get training information
    '''
    training_info = TrainingInfo(training_label=training_label)

    return training_info


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


  def createFileWithAnalysisTree(self, training_label, samples_filename, selection, weights, label, treename):
    '''
      Create tree that contains the hnl mass and the score and store it in a root file
      This tree is going to be used for the rest of the analysis (see fitter)
      Note that both the baseline selection and category definition have to be applied
    '''
    # get the training information
    training_info = self.getTrainingInfo(training_label)

    # create dataframe
    print '\n --> create dataframe'
    df = self.createDataframe(training_info=training_info, samples_filename=samples_filename, selection=selection, weights=weights, signal_treename=treename)

    # get the score
    print '\n --> predict the score'
    score = self.predictScore(training_info=training_info, df=df) 

    # get other quantities to fill the branches with
    hnl_mass = df['hnl_mass']
    weight_val = {}
    if weights != None:
      for weight in weights:
        weight_val[weight] = df[weight]

    # create file
    root_filename = './{}.root'.format(label.replace('.', 'p'))
    out_file = ROOT.TFile(root_filename, 'RECREATE')

    # create tree
    print ' --> create analysis tree'
    tree = ROOT.TTree(treename, treename)

    # initialise branches
    the_score = array('d', [0])
    the_hnl_mass = array('d', [0])
    the_weight = {}
    if weights != None:
      for weight in weights:
        the_weight[weight] = array('d', [0])

    tree.Branch('score', the_score, 'score/D')
    tree.Branch('hnl_mass', the_hnl_mass, 'hnl_mass/D')
    if weights != None:
      for weight in weights:
        tree.Branch(weight, the_weight[weight], '{}/D'.format(weight))

    # fill the tree
    for entry in range(len(score)):
      the_score[0] = score[entry][0]
      the_hnl_mass[0] = hnl_mass[entry]
      if weights != None:
        for weight in weights:
          the_weight[weight][0] = weight_val[weight][entry]

      tree.Fill()

    tree.Write()
    out_file.Close()

    print ' --> {} created'.format(root_filename)

    
  def getFileWithScore(self, files=None, training_label='', selection='hnl_charge>-99', weights=None, label='', treename='signal_tree', force_overwrite=False): 
    '''
      This function returns the file with the analysis tree that contains the hnl mass, the score and other quantities used for signal reweighting
      The argument 'weights' is the list of branches that will need to be added to the tree for the reweighting at analysis level
      Note that the baseline selection and the category definition are to be applied, as only the cut on the score is kept at analysis level (see fitter)
      The argument 'selection' can still contain the cut on the score, as it is going to be removed automatically when creating the dataframes
    '''
    samples_filename = []
    for file_ in files:
      samples_filename.append(file_.filename)
      
    root_filename = './{}.root'.format(label.replace('.', 'p'))
    if force_overwrite or not path.exists(root_filename):
      self.createFileWithAnalysisTree(training_label=training_label, samples_filename=samples_filename, selection=selection, weights=weights, label=label, treename=treename)

    return root_filename

