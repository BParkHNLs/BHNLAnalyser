import os 
import sys
from os import path

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from time import time 
from datetime import datetime
import seaborn as sns

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

class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''
  def __init__(self, filename, selection):
    self.filename = filename
    self.selection = selection
    self.df = read_root(self.filename, 'signal_tree', where=self.selection, warn_missing_tree=True)
    print '\t selected events: {}'.format(len(self.df))




class Trainer(object):
  def __init__(self, features, epochs, batch_size, scaler_type, do_early_stopping, do_reduce_lr, dirname, baseline_selection):
    self.features = features
    self.epochs = epochs
    self.batch_size = batch_size
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')

    self.baseline_selection = baseline_selection

    self.target_branch = 'is_signal'


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = './outputs/{}'.format(self.dirname)
    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))

    return outdir


  def saveFig(self, plt, name):
    '''
      Save python figure
    '''
    plt.savefig('{}/{}.pdf'.format(self.outdir, name))    
    plt.savefig('{}/{}.png'.format(self.outdir, name))    
    print ' --> {}/{}.png created'.format(self.outdir, name)


  def getSamples(self):
    '''
      Function that fetches the samples into lists
    '''
    #FIXME this should not be hardcoded!
    #TODO in case of more complex baseline selection, convert it to a format that is digestable by pandas query

    print('========> starting reading the trees')
    now = time()
    ## data 
    # used for training
    filename_data_1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj1.root'
    filename_data_2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj2.root'
    filename_data_3 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj3.root'
    data_samples = [
      Sample(filename=filename_data_1, selection=self.baseline_selection),
      Sample(filename=filename_data_2, selection=self.baseline_selection),
      Sample(filename=filename_data_3, selection=self.baseline_selection),
    ]

    # used for testing
    filename_data_4 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj4.root'
    filename_data_5 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj5.root'
    filename_data_6 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj6.root'
    data_test_samples = [
      Sample(filename=filename_data_4, selection=self.baseline_selection),
      Sample(filename=filename_data_5, selection=self.baseline_selection),
      Sample(filename=filename_data_6, selection=self.baseline_selection),
    ]

    ## signal
    # used for training
    filename_mc_1 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj1.root'
    filename_mc_2 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj2.root'
    filename_mc_3 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj3.root'
    filename_mc_4 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj4.root'
    filename_mc_5 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj5.root'
    filename_mc_6 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj6.root'
    mc_samples = [
      Sample(filename=filename_mc_1, selection=self.baseline_selection),
      Sample(filename=filename_mc_2, selection=self.baseline_selection),
      Sample(filename=filename_mc_3, selection=self.baseline_selection),
      Sample(filename=filename_mc_4, selection=self.baseline_selection),
      Sample(filename=filename_mc_5, selection=self.baseline_selection),
      Sample(filename=filename_mc_6, selection=self.baseline_selection),
    ]

    # used for testing
    filename_mc_4 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj4.root'
    filename_mc_5 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj5.root'
    filename_mc_6 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj6.root'
    #filename_mc_4 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n15/flat/flat_bparknano_08Aug22_nj4.root'
    #filename_mc_5 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n15/flat/flat_bparknano_08Aug22_nj5.root'
    #filename_mc_6 = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n15/flat/flat_bparknano_08Aug22_nj6.root'
    mc_test_samples = [
      Sample(filename=filename_mc_4, selection=self.baseline_selection),
      Sample(filename=filename_mc_5, selection=self.baseline_selection),
      Sample(filename=filename_mc_6, selection=self.baseline_selection),
    ]
      
    print('========> it took %.2f seconds' %(time() - now))

    return data_samples, data_test_samples, mc_samples, mc_test_samples


  def createDataframe(self, data_samples, mc_samples):
    '''
      Function that converts root samples to pandas dataframe
    '''
    data_df = pd.concat([idt.df for idt in data_samples], sort=False)
    mc_df   = pd.concat([imc.df for imc in mc_samples]  , sort=False)

    return data_df, mc_df


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

      qt.fit(X[features])
      xx = qt.transform(X[features])

      return xx, qt


  def preprocessing(self, data_df, mc_df):
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
    X = pd.DataFrame(main_df, columns=list(set(self.features)))
    Y = pd.DataFrame(main_df, columns=[self.target_branch])

    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)

    # and save the scaler, which will have to be used throughout the full process, even at evaluation time
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    pickle.dump(qt,open(scaler_filename, 'wb'))
    print ' --> {} created'.format(scaler_filename)

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print ' --> {} created'.format(features_filename)

    return main_df, qt, xx, Y


  def defineModel(self):
    '''
      Define the NN
    '''
    #NOTE for the moment, everything is hardcoded

    activation = 'relu'
    
    # define the net
    input  = Input((len(features),))
    layer  = Dense(64, activation=activation, name='dense1', kernel_constraint=unit_norm())(input)
    output = Dense(1 , activation='sigmoid' , name='output', )(layer)

    # Define outputs of your model
    model = Model(input, output)
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['mae', 'acc'])
    
    print(model.summary())

    return model


  def defineCallbacks(self):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=50)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.00001, cooldown=10, verbose=True)
    
    # save the model every now and then
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
    x_train = pd.DataFrame(x_train, columns=list(set(self.features))) # alternative to X_train[self.features[:]]
    x_val = pd.DataFrame(x_val, columns=list(set(self.features)))

    return x_train, x_val, y_train, y_val


  def prepareScaledInputs(self, main_df, Y, qt):
    '''
      Separate the main dataframe into training and validation sets
      The features are normalised according to the scaler qt
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

    return history
    

  def plotLoss(self, history):
    '''
      Plot the loss for training and validation sets
    '''
    loss_train = history.history['loss']
    loss_val = history.history['val_loss']
    epochs_range = range(1, self.epochs+1)
    epochs = epochs_range
    plt.plot(epochs, loss_train, 'g', label='Training loss')
    plt.plot(epochs, loss_val, 'b', label='Validation loss')
    plt.title('Training and Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'loss')
    plt.clf()


  def plotAccuracy(self, history):
    '''
      Plot the accuracy for training and validation sets
    '''
    acc_train = history.history['acc']
    acc_val = history.history['val_acc']
    epochs_range = range(1, self.epochs+1)
    epochs = epochs_range
    plt.plot(epochs, acc_train, 'g', label='Training accuracy')
    plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
    plt.title('Training and Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'accuracy')
    plt.clf()


  def predictScore(self, model, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=self.features)

    # apply the scaler
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    qt = pickle.load(open(scaler_filename, 'rb'))
    xx = qt.transform(x[self.features])

    # predict
    score = model.predict(xx)

    return score


  def plotScore(self, model, mc_test_df, data_test_df):
    '''
      Plot the score distributions for signal and background
    '''
    # get the score for the test dataframe
    sig_score = self.predictScore(model, mc_test_df)
    bkg_score = self.predictScore(model, data_test_df)

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
    ax.set_title("Score distribution for testing set", fontsize=20)
    ax.set_xlabel('Score',fontsize=18)
    fig.savefig('outputs/score.pdf')
    fig.savefig('outputs/score.png')
    self.saveFig(fig, 'score')
    plt.clf()


  def plotROC(self, model, x_train, y_train, x_val, y_val):
    '''
      Plot the ROC curve
      Note that the x inputs are already scaled
    '''
    score_train = model.predict(x_train)
    fpr, tpr, wps = roc_curve(y_train, score_train) 

    plt.clf()
    plt.plot(fpr, tpr, label='train ROC')
    #print("AUC train",sk.metrics.auc(fpr,tpr))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    score_val = model.predict(x_val)
    fpr, tpr, wps = roc_curve(y_val, score_val) 
    plt.plot(fpr, tpr, label='test ROC')
    #print("AUC test",sk.metrics.auc(fpr,tpr))

    xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.yscale('linear')
    plt.legend()
    self.saveFig(plt, 'ROC')
    plt.clf()


  def plotCorrelations(self, model, df, label):
    '''
      Plot the correlation matrix based on the training set
    '''
    if label not in ['data', 'mc']:
      raise RuntimeError('Unknown label "{}". Aborting'.format(label))

    # get the score for the test dataframe
    score = self.predictScore(model, df)

    # add the score to the dataframes
    df['score'] = score

    corr = df[self.features + ['score']].corr()
    
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))
    
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
                    square=True, linewidths=.8, cbar_kws={"shrink": .8})

    # rotate axis labels
    g.set_xticklabels(self.features+['score'], rotation='vertical')
    g.set_yticklabels(self.features+['score'], rotation='horizontal')

    plt.title('Linear Correlation Matrix - {}'.format(label))
    plt.tight_layout()
    self.saveFig(plt, 'correlations_{}'.format(label))
    plt.clf()


  def plotKSTest(self, model, x_train, x_val, y_train, y_val, label):
    '''
      Plot the outcome of the Kolmogorov test
      Used to test the overfitting
    '''
    if label not in ['data', 'mc']:
      raise RuntimeError('Unknown label "{}". Aborting'.format(label))


    # only keep the data or mc components of the features
    # does it mess up with the normalisation?
    if label == 'data':
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
    h1.Draw("hist")
    h2.SetLineColor(ROOT.kRed)
    h2.SetFillColor(ROOT.kWhite)
    h1.SetFillColor(ROOT.kWhite)
    h2.Draw("hist SAME")
    c1.BuildLegend()
    ks_score = h1.KolmogorovTest(h2)
    ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
    ks_value.AddText('KS score {} = {}'.format(label, round(ks_score, 3)))
    ks_value.SetFillColor(0)
    ks_value.Draw('EP same')

    c1.SaveAs(self.outdir + 'KS_test_{}.png'.format(label))
    #print("KS score: ",ks_score, len(train_pred),len(test_pred))



  def process(self):
    print '---- MVA Trainer ----'
    
    # create output directory
    print '\n -> create output directory'
    self.outdir = self.createOutDir()

    # get the samples
    print '\n -> get the samples'
    data_samples, data_test_samples, mc_samples, mc_test_samples = self.getSamples()

    # create dataframes
    print '\n -> create the dataframes'
    data_df, mc_df = self.createDataframe(data_samples, mc_samples)
    data_test_df, mc_test_df = self.createDataframe(data_test_samples, mc_test_samples)

    # assign the signal tag
    data_df = self.assignTarget(data_df, self.target_branch, 0) 
    mc_df = self.assignTarget(mc_df, self.target_branch, 1) 

    # here in the rjpsi, there would be some index engineering to isolate the testing
    # set from the df. Instead, we want to test on different samples, which do not
    # necessarily have the same signal signature, so we keep the full df and will
    # use for testing the above defined *_test_df instead

    # preprocessing the dataframes
    print '\n -> preprocessing the dataframes' 
    main_df, qt, xx, Y = self.preprocessing(data_df, mc_df)
    #TODO add the scatter plots?

    # define the NN
    print '\n -> defining the model' 
    model = self.defineModel()

    # define the callbacks
    print '\n -> defining the callbacks' 
    callbacks = self.defineCallbacks()

    # out of the main_df, define which data chunks to 
    # train and test on. Make sure that the features
    # are scaled
    print '\n -> prepare the inputs' 
    #x_train, x_val, y_train, y_val = self.prepareScaledInputs(main_df, Y, qt)
    x_train, x_val, y_train, y_val = self.prepareInputs(xx, Y)

    # do the training
    print '\n -> training...' 
    history = self.train(model, x_train, y_train, x_val, y_val, callbacks)

    # plotting
    print '\n -> plotting...' 
    self.plotLoss(history)
    self.plotAccuracy(history)
    self.plotScore(model, mc_test_df, data_test_df)
    self.plotROC(model, x_train, y_train, x_val, y_val)
    self.plotCorrelations(model, data_df, 'data')
    self.plotCorrelations(model, mc_df, 'mc')
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 'data')
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 'mc')




if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  features = ['pi_pt', 'b_mass', 'pi_dcasig', 'hnl_cos2d']
  #features = ['pi_pt', 'pi_dcasig']
  epochs = 50
  batch_size = 32
  scaler_type = 'robust'
  do_early_stopping = False
  do_reduce_lr = True
  dirname = 'test'
  baseline_selection = 'hnl_charge==0'
  #NOTE add optimiser, learning rate etc? 

  trainer = Trainer(
      features = features, 
      epochs = epochs,
      batch_size = batch_size,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      dirname = dirname,
      baseline_selection = baseline_selection,
      )

  trainer.process()




