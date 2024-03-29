import os 
import sys
from os import path

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import product
from time import time 
from datetime import datetime
#import seaborn as sns
import glob

import ROOT
import root_pandas
from root_numpy import root2array
from root_pandas import read_root

from keras.models import Sequential, Model
from keras.layers import Dense, Input, Dropout, BatchNormalization
from keras.utils import plot_model
from keras.callbacks import EarlyStopping, Callback, ReduceLROnPlateau, ModelCheckpoint
from keras.constraints import unit_norm
from keras.optimizers import SGD, Adam

import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.preprocessing import RobustScaler, StandardScaler

sys.path.append('../')
sys.path.append('../../objects')
from tools import Tools
from samples import signal_samples
from baseline_selection import selection
from categories import categories

#TODO add gen-matching
#TODO make sure that only gen matched events are saved in bc

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to perform the NN training', add_help=True)
  parser.add_argument('--category_batch', type=str, dest='category_batch', help='category label'                               , default=None)
  parser.add_argument('--outdir'        , type=str, dest='outdir'        , help='output directory'                             , default=None)
  parser.add_argument('--process'       , dest='process'                 , help='run the process function', action='store_true', default=False)
  return parser.parse_args()


class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''
  def __init__(self, filename, selection):
    self.filename = filename
    self.selection = selection
    #print self.selection
    self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True)
    try:
      self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True)
    except:
      print 'No entry was found with the requested selection'
      self.df = pd.DataFrame()
    print '\t selected events: {}'.format(len(self.df))




class Trainer(object):
  def __init__(self, features, epochs, batch_size, learning_rate, scaler_type, do_early_stopping, do_reduce_lr, do_parametric, signal_label, data_pl, data_tagnano, data_tagflat, nsigma, dirname, baseline_selection, categories, category_batch, outdir):
    self.features = features
    self.epochs = epochs
    self.batch_size = batch_size
    self.learning_rate = learning_rate
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%Y%b%d_%Hh%Mm%Ss')

    self.baseline_selection = baseline_selection
    self.categories = categories
    self.category_batch = category_batch
    self.signal_label = signal_label
    self.signal_files = signal_samples[self.signal_label]

    self.outdir = outdir

    self.username = 'anlyon'
    self.data_pl = data_pl
    self.data_tagnano = data_tagnano
    self.data_tagflat = data_tagflat

    self.do_parametric = do_parametric
    self.nsigma = nsigma

    self.target_branch = 'is_signal'
    self.do_scale_key = True

    self.resolution_p0 = 6.98338e-04
    self.resolution_p1 = 7.78382e-03


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = './outputs/{}'.format(self.dirname)
    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))

    return outdir


  def writeSubmitter(self, category):
    '''
      Write bash submitter
    '''
    content = '\n'.join([
        '#!/bin/bash',
        'homedir="$PWD"',
        'workdir=/scratch/{}/training_{}'.format(self.username, category.label),
        'mkdir -p $workdir',
        'cp -r ../*py $workdir', #TODO adapt if changing directory
        'cp -r ./*py $workdir',
        'cp -r ../../objects/*py $workdir',
        'cd $workdir',
        'DATE_START=`date +%s`',
        'python trainer.py --category_batch {cat} --outdir {out} --process'.format(cat=category.label, out=self.outdir),
        'DATE_END=`date +%s`',
        'runtime=$((DATE_END-DATE_START))',
        'echo " --> Wallclock running time: $runtime s"',
        'cp -r *h5 $homedir/{out}'.format(out=self.outdir),
        'cp -r *pck $homedir/{out}'.format(out=self.outdir),
        'cp -r *png $homedir/{out}'.format(out=self.outdir),
        'cp -r *pdf $homedir/{out}'.format(out=self.outdir),
        'cp -r *txt $homedir/{out}'.format(out=self.outdir),
        'cp trainer.py $homedir/{out}'.format(out=self.outdir),
        'cd $homedir',
        'rm -r $workdir',
        ])
    submitter_name = '{}/submitter_{}.sh'.format(self.outdir, category.label)
    submitter = open(submitter_name, 'w+')
    submitter.write(content)
    submitter.close()
    '\n -> {} created'.format(submitter_name)


  def submit(self, category):
    '''
      Submit bash script on slurm
    '''
    command = 'sbatch -p standard --account t3 -o {out}/log_{cat}.txt -e {out}/log_{cat}.txt --job-name=trainer_{cat} --mem {mem} {out}/submitter_{cat}.sh'.format(
        out = self.outdir,
        cat = category.label,
        mem = 20000,
        )
    print '\n --> submitting category {}'.format(category.label)
    os.system(command)


  def saveFig(self, plt, name):
    '''
      Save python figure
    '''
    plt.savefig('{}/{}.pdf'.format(self.outdir, name))    
    plt.savefig('{}/{}.png'.format(self.outdir, name))    
    print ' --> {}/{}.png created'.format(self.outdir, name)


  def getDataSamples(self, extra_selection=None, max_files=-1):
    '''
      Function that fetches the samples into lists
    '''
    print('========> starting reading the trees')
    now = time()
    ## data 
    path = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/data'.format(self.username)
    if self.data_tagflat != None:
      filename = '{path}/{pl}/ParkingBPH1_Run2018D/Chunk*/flat/flat_bparknano_{tagnano}_{tagflat}.root'.format(path=path, pl=self.data_pl, tagnano=self.data_tagnano, tagflat=self.data_tagflat)
    else:
      filename = '{path}/{pl}/ParkingBPH1_Run2018D/Chunk*/flat/flat_bparknano_{tagnano}.root'.format(path=path, pl=self.data_pl, tagnano=self.data_tagnano)
    data_filenames =  [f for f in glob.glob(filename)] 
    data_filenames = sorted(data_filenames, key=lambda file_: float(file_[file_.find('Chunk')+5:file_.find('_', file_.find('Chunk')+5)]))

    data_samples = []
    for ifile, data_filename in enumerate(data_filenames):
      #if ifile<10: continue # those datasets will be used for the validation
      #if ifile<30: continue # those datasets will be used for the validation
      if max_files != -1 and ifile > max_files+9: continue
      print data_filename
      data_samples.append(Sample(filename=data_filename, selection=self.baseline_selection + ' && ' + extra_selection))

    print('========> it took %.2f seconds' %(time() - now))

    return data_samples


  def getMCSamples(self, extra_selection=None, max_files=-1):
    '''
      Function that fetches the samples into lists
    '''
    print('========> starting reading the trees')
    now = time()
    filename_mc_1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root'
    filename_mc_2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root'
    filename_mc_3 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root'

    mc_samples = [
      Sample(filename=filename_mc_1, selection=self.baseline_selection + ' && ' + extra_selection),
      Sample(filename=filename_mc_2, selection=self.baseline_selection + ' && ' + extra_selection),
      Sample(filename=filename_mc_3, selection=self.baseline_selection + ' && ' + extra_selection),
    ]
      
    print('========> it took %.2f seconds' %(time() - now))

    return mc_samples


  def getPandasQuery(self, selection):
    '''
      Converts selection to pandas query syntax
    '''
    query = selection.replace('&&', ' and ').replace('||', ' or ').replace('!=', ' not ').replace('!', ' not ')
    
    return query


  def getMassList(self, is_bc=False):
    '''
      Get the list of signal masses used in the training
    '''
    masses = []
    for signal_file in self.signal_files:
      mass = signal_file.mass
      if mass not in masses:
        if not is_bc:
          masses.append(mass)
        elif is_bc and mass >= 3.:
          masses.append(mass)

    masses.sort()

    return masses


  def removeInfs(self, df):
    '''
      Remove rows with infs/nan in at least one of the features
    '''
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=self.features)

    return df


  def createDataframe(self, data_samples, mc_samples):
    '''
      Function that converts root samples to pandas dataframe
    '''
    data_df = pd.concat([idt.df for idt in data_samples], sort=False)
    mc_df   = pd.concat([imc.df for imc in mc_samples]  , sort=False)

    data_df = self.removeInfs(data_df)
    mc_df = self.removeInfs(mc_df)

    return data_df, mc_df


  def createParametrisedDataframe(self, extra_selection, data_type, is_bc=False, statistics=None, max_files=-1):
    '''
      Function to create the dataframe with the training features and parameter
      The statistics is balanced between the signal and background, among the signal
      hypotheses and the background windows
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data type "{}". Please check'.format(data_type))

    if data_type == 'data' and statistics == None:
      raise RuntimeError('Please first create mc dataframe in order to get the available statistics')

    pd.options.mode.chained_assignment = None

    masses = self.getMassList(is_bc=is_bc)
    print 'masses ',masses
    
    if data_type == 'mc':
      stats = dict.fromkeys(masses, 0)
      n_ctaus = dict.fromkeys(masses, 0)
      dfs = dict.fromkeys(masses, []) 

      for signal_file in self.signal_files: 
        # get sample
        signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
        if is_bc and signal_file.filename_Bc == None: continue # skip points for which there is no Bc sample
        print signal_filename
        sample = Sample(filename=signal_filename, selection=self.baseline_selection + ' && ' + extra_selection)

        # get dataframe
        the_df = sample.df

        # get the statistics per mass point
        stats[signal_file.mass] += len(the_df)

        # get number of ctau samples per mass point
        n_ctaus[signal_file.mass] += 1

        # set the mass parameter to the exact value
        the_df['mass_key'] = signal_file.mass

        # compute the signal weight (unused for the time being)
        # lumi set to one as the normalisation does not matter
        #the_df['weight'] = Tools().getSignalWeight(signal_files=[signal_file], mass=signal_file.mass, ctau=signal_file.ctau, sigma_B=472.8e9, lumi=1, lhe_efficiency=0.08244) #TODO what about isBc?

        dfs[signal_file.mass] = dfs[signal_file.mass] + [the_df]

      # get the minimum statistics 
      #statistics = min(stats[min(stats, key=stats.get)], 5000)
      #aimed_statistics = stats[min(stats, key=stats.get)]
      statistics = stats[min(stats, key=stats.get)]
      
      df = pd.DataFrame()
      for mass in masses:
        df_tmp = pd.concat([(idt.sample(statistics/n_ctaus[mass]) if statistics/n_ctaus[mass]<len(idt) else idt) for idt in dfs[mass]], sort=False)
        df = pd.concat([df, df_tmp], sort=False) 
      
        
    elif data_type == 'data':
      # get the sample
      samples = self.getDataSamples(extra_selection=extra_selection, max_files=max_files) 

      df_full = [] # will be used for plotting purposes later on
      dfs = dict.fromkeys(masses, []) 

      for sample in samples:
        # get the dataframe
        the_df = sample.df
        df_full.append(the_df)

        # first set an invalid parameter for the full spectrum. Ranges with this parameter will be removed from the training
        the_df['mass_key'] = -1 

        # and define the windows of nsigma around the mass hypothesis
        # set the mass parameter for a nsigma window range
        window_check = -1
        for mass in masses:
          resolution = self.resolution_p0 + self.resolution_p1 * mass
          window_min = mass - self.nsigma * resolution
          window_max = mass + self.nsigma * resolution
          # prevent windows to overlap
          #if window_check != -1 and window_min <= window_check:
          #  raise RuntimeError('There is an overlap of window around mass {} GeV with the previous window. Please adapt the window size.'.format(mass))

          window = 'hnl_mass > {} && hnl_mass < {}'.format(window_min, window_max)
          df_window = the_df.query(self.getPandasQuery(window))
          df_window['mass_key'] = mass

          # set the weight
          df_window['weight'] = 1

          dfs[mass] = dfs[mass] + [df_window]

          window_check = window_max

      df_full = pd.concat([idt for idt in df_full], sort=False)

      df = pd.DataFrame()
      for mass in masses:
        requested_stat = statistics
        stat_permass = 0
        for the_df in dfs[mass]:
          stat_permass += len(the_df)
        if stat_permass < statistics:
          print 'Requesting {} events from a sample of {} events. Please provide a larger input data sample'.format(statistics, stat_permass)
          requested_stat = stat_permass
        df_tmp = pd.concat([idt for idt in dfs[mass]], sort=False)
        df_tmp = df_tmp.sample(requested_stat)
        df = pd.concat([df, df_tmp], sort=False) 
    
    df = self.removeInfs(df)

    if data_type == 'mc':
      return df, statistics
    else:
      return df, df_full


  def assignTarget(self, df, branch, target):
    '''
      Add the target branch to the dataframes
    '''
    pd.options.mode.chained_assignment = None
    df[branch] = target

    return df


  def doScaling(self, X):
      '''
        Normalise the input features with a keras scaler 
      '''
      if self.scaler_type == 'robust':
        qt = RobustScaler()
      elif self.scaler_type == 'standard':
        qt = StandardScaler()
      else:
        raise RuntimeError('Unknown scaler "{}" - Aborting...'.format(self.scaler_type))

      if not self.do_parametric:
        qt.fit(X[self.features])
        xx = qt.transform(X[self.features])
      else:
        if self.do_scale_key:
          qt.fit(X[self.features + ['mass_key']])
          xx = qt.transform(X[self.features + ['mass_key']])
        else:
          qt.fit(X[self.features])
          xx = qt.transform(X[self.features])
          xx = np.insert(xx, (len(self.features)), X['mass_key'], axis=1)

      return xx, qt


  def preprocessing(self, data_df, mc_df, label):
    '''
      Preprocessing of data before training/testing the NN
      This includes:
        - building the main_df
        - building the scaler
        - get the scaled features xx
        - get the target Y
    '''

    # concatenate the events and shuffle
    main_df = pd.concat([data_df, mc_df], sort=False)

    # re-index
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=1986) # of course, keep R's seed ;)

    # X and Y
    if not self.do_parametric:
      X = pd.DataFrame(main_df, columns=list(set(self.features)))
    else:
      X = pd.DataFrame(main_df, columns=list(set(self.features)) + ['mass_key'])
    Y = pd.DataFrame(main_df, columns=[self.target_branch])

    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)

    # and save the scaler, which will have to be used throughout the full process, even at evaluation time
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted_{}.pck'.format(label)])
    pickle.dump(qt,open(scaler_filename, 'wb'))
    print ' --> {} created'.format(scaler_filename)

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features_{}.pck'.format(label)])
    if not self.do_parametric:
      pickle.dump(self.features, open(features_filename, 'wb' ))
    else:
      pickle.dump(self.features + ['mass_key'], open(features_filename, 'wb' ))
    print ' --> {} created'.format(features_filename)

    return main_df, qt, xx, Y


  def defineModel(self):
    '''
      Define the NN
    '''
    #NOTE for the moment, everything is hardcoded

    activation = 'relu'
    
    # define the net
    n_input = len(self.features) if not self.do_parametric else len(self.features) + 1 
    input  = Input((n_input,))
    layer  = Dense(64, activation=activation, name='dense1', kernel_constraint=unit_norm())(input)
    output = Dense(1 , activation='sigmoid' , name='output', )(layer)

    # Define outputs of your model
    model = Model(input, output)
    optimiser = Adam(lr=self.learning_rate)
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['mae', 'acc'])
    
    print(model.summary())

    return model


  def defineCallbacks(self, label):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=10)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.00001, cooldown=10, verbose=True)
    
    # save the model every now and then
    # kept only during excecution time and removed afterwards
    filepath = '/'.join([self.outdir, 'saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    
    callbacks = [save_model]
    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    return callbacks


  def prepareInputs(self, xx, Y):
    '''
      Separate the main dataframe into training and validation sets
      Note: the input xx should arlready be scaled
    '''

    x_train, x_val, y_train, y_val = train_test_split(xx, Y, test_size=0.2, shuffle=True)
    
    # the x should only contain the features and not all the branches of main_df
    if not self.do_parametric:
      x_train = pd.DataFrame(x_train, columns=list(set(self.features))) # alternative to X_train[self.features[:]]
      x_val = pd.DataFrame(x_val, columns=list(set(self.features)))
    else:
      x_train = pd.DataFrame(x_train, columns=list(set(self.features)) + ['mass_key'])
      x_val = pd.DataFrame(x_val, columns=list(set(self.features)) + ['mass_key'])

    x_train = x_train.reset_index(drop=True)
    x_val = x_val.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)
    y_val = y_val.reset_index(drop=True)

    return x_train, x_val, y_train, y_val


  def prepareScaledInputs(self, main_df, Y, qt):
    '''
      Separate the main dataframe into training and validation sets
      The features are normalised according to the scaler qt

      NB: function not adapted to the parametric case
    '''

    x_train, x_val, y_train, y_val = train_test_split(main_df, Y, test_size=0.2, shuffle=True)

    # scale the features
    x_train = qt.transform(x_train[self.features])
    x_val = qt.transform(x_val[self.features])

    return x_train, x_val, y_train, y_val


  def train(self, model, x_train, y_train, x_val, y_val, callbacks):
    '''
      Perform the training
    '''
    history = model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True)
    #TODO apply weight if any

    return history
    

  def predictScore(self, model, df, label):
    '''
      Return score with scaled input features
    '''
    if not self.do_parametric or not self.do_scale_key:
      x = pd.DataFrame(df, columns=self.features)
    else:
      x = pd.DataFrame(df, columns=self.features + ['mass_key'])

    # apply the scaler
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted_{}.pck'.format(label)])
    qt = pickle.load(open(scaler_filename, 'rb'))

    # if not scaling key
    if not self.do_parametric or not self.do_scale_key:
      xx = qt.transform(x[self.features])
    else:
      xx = qt.transform(x[self.features + ['mass_key']])

    # predict
    score = model.predict(xx)

    return score


  def saveModel(self, model, label):
    '''
      Save the model
    '''
    model_filename = '{}/net_model_weighted_{}.h5'.format(self.outdir, label)
    model.save(model_filename)
    print ' --> {} created'.format(model_filename)


  def plotParametrisedMass(self, df, df_full=None, data_type='', label=''):
    '''
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data type "{}". Please check'.format(data_type))

    masses = self.getMassList()
    bin_min = 0
    bin_max = 6.5
    nbins = 80.

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #if not df_full.empty:
    if isinstance(df_full, pd.DataFrame):
      ax.hist(df_full['hnl_mass'], bins=np.arange(bin_min, bin_max, (bin_max-bin_min)/nbins), color='whitesmoke', label='unused')

    for mass in masses:
      resolution = self.resolution_p0 + self.resolution_p1 * mass
      window_min = mass - self.nsigma * resolution
      window_max = mass + self.nsigma * resolution
      window = 'hnl_mass > {} && hnl_mass < {}'.format(window_min, window_max)
      df_window = df.query(self.getPandasQuery(window))
      ax.hist(df_window['hnl_mass'], bins=np.arange(bin_min, bin_max, (bin_max-bin_min)/nbins), label='mass key = {}'.format(mass))

    ax.legend(loc='upper right',prop={'size': 12})
    ax.set_title(label, fontsize=20)
    ax.set_xlabel('Mass spectrum',fontsize=18)
    self.saveFig(fig, 'mass_{}_{}'.format(data_type, label))
    plt.clf()


  def plotLoss(self, history, label):
    '''
      Plot the loss for training and validation sets
    '''
    loss_train = history.history['loss']
    loss_val = history.history['val_loss']
    epochs_range = range(1, self.epochs+1)
    epochs = epochs_range
    plt.plot(epochs, loss_train, 'g', label='Training loss')
    plt.plot(epochs, loss_val, 'b', label='Validation loss')
    plt.title('Training and Validation Loss - {}'.format(label))
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'loss_{}'.format(label))
    plt.clf()


  def plotAccuracy(self, history, label):
    '''
      Plot the accuracy for training and validation sets
    '''
    acc_train = history.history['acc']
    acc_val = history.history['val_acc']
    epochs_range = range(1, self.epochs+1)
    epochs = epochs_range
    plt.plot(epochs, acc_train, 'g', label='Training accuracy')
    plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
    plt.title('Training and Validation Accuracy - {}'.format(label))
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend(loc='lower left')
    self.saveFig(plt, 'accuracy_{}'.format(label))
    plt.clf()


  def plotScore(self, model, mc_test_df, data_test_df, label):
    '''
      Plot the score distributions for signal and background
    '''
    # get the score for the test dataframe
    sig_score = self.predictScore(model, mc_test_df, label)
    bkg_score = self.predictScore(model, data_test_df, label)
    # add the score to the dataframes
    mc_test_df['score'] = sig_score
    data_test_df['score'] = bkg_score

    # plot the score distributions
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bkg_score = [data_test_df['score']]
    bkg_name=['data-driven background']
    ax.hist(bkg_score, bins=np.arange(0,1.025,0.025), color='blue', alpha=0.8, label=bkg_name)
    ax.hist(mc_test_df['score'], bins=np.arange(0,1.025,0.025), color='darkorange', alpha=1, label='signal', histtype='step', linewidth=2)
    ax.legend(loc='upper left',prop={'size': 12})
    ax.set_title("Score distribution for testing set - {}".format(label), fontsize=20)
    ax.set_xlabel('Score',fontsize=18)
    self.saveFig(fig, 'score_{}'.format(label))
    plt.clf()


  def plotROC(self, model, x_train, y_train, x_val, y_val, label):
    '''
      Plot the ROC curve
      Note that the x inputs are already scaled
    '''
    score_train = model.predict(x_train)
    fpr, tpr, wps = roc_curve(y_train, score_train) 

    plt.plot(fpr, tpr, label='train ROC')
    #print("AUC train",sk.metrics.auc(fpr,tpr))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    score_val = model.predict(x_val)
    fpr, tpr, wps = roc_curve(y_val, score_val) 
    plt.plot(fpr, tpr, label='test ROC')
    plt.title('ROC - {}'.format(label))
    #print("AUC test",sk.metrics.auc(fpr,tpr))

    xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.yscale('linear')
    plt.legend()
    self.saveFig(plt, 'ROC_{}'.format(label))
    plt.clf()


  #def plotCorrelations(self, model, df, data_type, label):
  #  '''
  #    Plot the correlation matrix based on the training set
  #  '''
  #  if data_type not in ['data', 'mc']:
  #    raise RuntimeError('Unknown data_type "{}". Aborting'.format(data_type))

  #  # get the score for the test dataframe
  #  score = self.predictScore(model, df, label)

  #  # add the score to the dataframes
  #  df['score'] = score

  #  corr = df[self.features + ['score']].corr()
  #  
  #  # Set up the matplotlib figure
  #  f, ax = plt.subplots(figsize=(11, 9))
  #  
  #  # Generate a custom diverging colormap
  #  cmap = sns.diverging_palette(220, 10, as_cmap=True)

  #  # Draw the heatmap with the mask and correct aspect ratio
  #  g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
  #                  square=True, linewidths=.8, cbar_kws={"shrink": .8})

  #  # rotate axis labels
  #  g.set_xticklabels(self.features+['score'], rotation='vertical')
  #  g.set_yticklabels(self.features+['score'], rotation='horizontal')

  #  plt.title('Linear Correlation Matrix - {}'.format(data_type))
  #  plt.tight_layout()
  #  self.saveFig(plt, 'correlations_{}_{}'.format(data_type, label))
  #  plt.clf()


  def plotKSTest(self, model, x_train, x_val, y_train, y_val, data_type, label):
    '''
      Plot the outcome of the Kolmogorov test
      Statistical test of compatibility in shape between two histograms
      Used to test the overfitting
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data_type "{}". Aborting'.format(data_type))


    # only keep the data or mc components of the features
    # does it mess up with the normalisation?
    if data_type == 'data':
      x_train_part = x_train.drop(y_train.query('is_signal==1').index)
      x_val_part = x_val.drop(y_val.query('is_signal==1').index)
    else:
      x_train_part = x_train.drop(y_train.query('is_signal==0').index)
      x_val_part = x_val.drop(y_val.query('is_signal==0').index)

    score_train = model.predict(x_train_part)
    score_val = model.predict(x_val_part)

    h1 = ROOT.TH1F('train', 'train', 30, -1, 1)
    h2 = ROOT.TH1F('val', 'val', 30, -1, 1)
    for t, b in zip(score_train, score_val):
      h1.Fill(t)
      h2.Fill(b)

    c1=ROOT.TCanvas()
    if h1.Integral()!=0: h1.Scale(1./h1.Integral())
    if h2.Integral()!=0: h2.Scale(1./h2.Integral())
    c1.Draw()
    h1.SetTitle(label)
    h1.Draw("hist")
    h2.SetLineColor(ROOT.kRed)
    h2.SetFillColor(ROOT.kWhite)
    h1.SetFillColor(ROOT.kWhite)
    h2.Draw("hist SAME")
    c1.BuildLegend()
    ks_score = h1.KolmogorovTest(h2)
    ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
    ks_value.AddText('KS score {} = {}'.format(data_type, round(ks_score, 3)))
    ks_value.SetFillColor(0)
    ks_value.Draw('EP same')

    c1.SaveAs('{}/KS_test_{}_{}.png'.format(self.outdir, data_type, label))


  def process(self):
    print '---- MVA Trainer ----'
    
    # create output directory
    if self.outdir == None:
      print '\n -> create output directory'
      self.outdir = self.createOutDir()
    else:
      self.outdir = '.'

    # copy script
    os.system('cp trainer.py {}'.format(self.outdir))
        
    for category in self.categories:
      if category.label == 'incl': continue
      #if category.label != 'lxysiggt150_OS_Bc': continue
      print '\n-.-.-'
      print 'category: {}'.format(category.label)
      print '-.-.-'

      if self.category_batch != None and category.label != category_batch: continue 

      # get the samples
      if not self.do_parametric:
        print '\n -> get the samples'
        data_samples = self.getDataSamples(extra_selection=category.definition_flat)
        mc_samples = self.getMCSamples(extra_selection=category.definition_flat)

        # create dataframes
        print '\n -> create the dataframes'
        data_df, mc_df = self.createDataframe(data_samples, mc_samples)

      else:
        print '\n -> create the dataframes'
        # do not load too much files in case of large statistics
        #if category.label in ['lxysig0to50_OS', 'lxysig0to50_SS', 'lxysig0to50_OS_Bc', 'lxysig0to50_SS_Bc']: max_files = 10 
        #elif category.label in ['lxysig50to150_OS', 'lxysig50to150_SS', 'lxysig50to150_OS_Bc', 'lxysig50to150_SS_Bc']: max_files = 30 
        #else: max_files = -1
        max_files = -1

        mc_df, statistics = self.createParametrisedDataframe(extra_selection=category.definition_flat, data_type='mc', is_bc=category.is_bc)
        data_df, data_df_full = self.createParametrisedDataframe(extra_selection=category.definition_flat, data_type='data', statistics=statistics, max_files=max_files, is_bc=category.is_bc)

      # assign the signal tag
      data_df = self.assignTarget(data_df, self.target_branch, 0) 
      mc_df = self.assignTarget(mc_df, self.target_branch, 1) 

      # preprocessing the dataframes
      print '\n -> preprocessing the dataframes' 
      main_df, qt, xx, Y = self.preprocessing(data_df, mc_df, category.label)
      #TODO add the scatter plots?

      if self.do_parametric:
        self.plotParametrisedMass(df=data_df, df_full=data_df_full, data_type='data', label=category.label)
        self.plotParametrisedMass(df=mc_df, data_type='mc', label=category.label)

      # define the NN
      print '\n -> defining the model' 
      model = self.defineModel()

      # define the callbacks
      print '\n -> defining the callbacks' 
      callbacks = self.defineCallbacks(category.label)

      # out of the main_df, define which data chunks to 
      # train and test on. Make sure that the features
      # are scaled
      print '\n -> prepare the inputs' 
      #x_train, x_val, y_train, y_val = self.prepareScaledInputs(main_df, Y, qt)
      x_train, x_val, y_train, y_val = self.prepareInputs(xx, Y)

      # create statistics file
      stat_filename = '{}/statistics_{}.txt'.format(self.outdir, category.label)
      stat_file = open(stat_filename, 'w+')
      stat_file.write('\nAimed statistics per mass parameter: {}'.format(statistics))
      stat_file.write('\nFull statistics for data: {}'.format(len(data_df)))
      stat_file.write('\nFull statistics for mc: {}'.format(len(mc_df)))
      stat_file.write('\nFull statistics for data+mc: {}'.format(len(main_df)))
      stat_file.write('\nData+mc statistics used for training: {}'.format(len(x_train)))
      stat_file.write('\nData+mc statistics used for validation: {}'.format(len(x_val)))
      stat_file.close()
      print ' --> {} created'.format(stat_filename) 

      # do the training
      print '\n -> training...' 
      history = self.train(model, x_train, y_train, x_val, y_val, callbacks)

      # save the model
      print '\n -> save the model' 
      self.saveModel(model, category.label)

      # plotting
      print '\n -> plotting...' 
      if self.do_parametric:
        self.plotParametrisedMass(df=data_df, df_full=data_df_full, data_type='data', label=category.label)
        self.plotParametrisedMass(df=mc_df, data_type='mc', label=category.label)
      self.plotLoss(history, category.label)
      self.plotAccuracy(history, category.label)
      #self.plotScore(model, mc_test_df, data_test_df, category.label)
      #self.plotScore(model, mc_df, data_df, category.label)
      self.plotROC(model, x_train, y_train, x_val, y_val, category.label)
      #self.plotCorrelations(model, data_df, 'data', category.label)
      #self.plotCorrelations(model, mc_df, 'mc', category.label)
      self.plotKSTest(model, x_train, x_val, y_train, y_val, 'data', category.label)
      self.plotKSTest(model, x_train, x_val, y_train, y_val, 'mc', category.label)

      # cleaning
      print '\n -> cleaning'
      os.system('rm -r {}/saved-model*h5'.format(self.outdir))

      print '\n --- Done ---'


  def process_batch(self):
    print '---- MVA Trainer ----'

    # create output directory
    print '\n -> create output directory'
    self.outdir = self.createOutDir()

    for category in self.categories:
      if category.label == 'incl': continue
      print '\n-.-.-'
      print 'category: {}'.format(category.label)
      print '-.-.-'
      self.writeSubmitter(category=category)
      self.submit(category=category)
      
    print '\n --- Done ---'



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  #features = ['pi_pt', 'b_mass', 'pi_dcasig', 'hnl_cos2d', 'sv_prob']
  #features = ['pi_pt', 'pi_eta', 'mu_pt', 'mu_eta', 'mu0_pt', 'mu0_eta', 'b_pt', 'b_eta', 'b_mass', 'mu0_mu_mass', 'mu0_mu_pt', 'mu0_pi_mass', 'mu0_pi_pt', 'deltar_mu_pi', 'deltar_mu0_hnl', 'deltar_mu0_mu', 'deltar_mu0_pi', 'sv_chi2', 'sv_prob']
  #features = ['pi_pt','mu_pt', 'mu0_pt','b_mass', 'mu0_mu_mass', 'mu0_pi_mass', 'deltar_mu_pi', 'deltar_mu0_mu', 'deltar_mu0_pi', 'sv_prob']
  #features = ['pi_pt','mu_pt', 'mu0_pt','b_mass', 'mu0_mu_mass', 'mu0_pi_mass', 'sv_prob'] # remove the deltaRs but keep masses as they can learn about possible sm resonances?
  #features = ['pi_pt','mu_pt', 'mu0_pt','b_mass', 'hnl_cos2d', 'pi_dcasig', 'sv_lxysig', 'sv_prob']
  features = ['pi_pt','mu_pt', 'mu0_pt','b_mass', 'hnl_cos2d', 'sv_lxysig', 'sv_prob', 'sv_chi2', 'b_pt', 'mu0_mu_mass', 'mu0_pi_mass', 'deltar_mu0_mu', 'deltar_mu0_pi', 'mu0_pfiso03_rel', 'mu_pfiso03_rel']
  #features = ['pi_pt', 'pi_dcasig']
  epochs = 60
  batch_size = 32
  learning_rate = 0.01
  scaler_type = 'robust'
  do_early_stopping = False
  do_reduce_lr = True
  dirname = 'V13_06Feb23'
  baseline_selection = 'hnl_charge==0 && ' + selection['baseline_08Aug22'].flat 
  categories = categories['categories_0_50_150_Bc']
  category_batch = getOptions().category_batch
  outdir = getOptions().outdir
  #NOTE add optimiser, learning rate etc? 

  submit_batch = True

  do_parametric = True
  nsigma = 10
  signal_label = 'V13_06Feb23_training_large'
  data_pl = 'V13_06Feb23'
  data_tagnano = '06Feb23'
  data_tagflat = 'partial'
  #signal_label = 'V12_08Aug22_training_large'
  #data_pl = 'V12_08Aug22'
  #data_tagnano = '08Aug22'
  #data_tagflat = 'sr'

  trainer = Trainer(
      features = features, 
      epochs = epochs,
      batch_size = batch_size,
      learning_rate = learning_rate,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      do_parametric = do_parametric,
      signal_label = signal_label, 
      data_pl = data_pl,
      data_tagnano = data_tagnano,
      data_tagflat = data_tagflat,
      nsigma = nsigma,
      dirname = dirname,
      baseline_selection = baseline_selection,
      categories = categories,
      category_batch = category_batch,
      outdir = outdir,
      )

  if submit_batch and not getOptions().process:
    trainer.process_batch()
  else:
    trainer.process()




