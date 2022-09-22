'''
  The config class gathers all the entries of the analysis cfgs/*_cfg.py
  The checkConfig() function is also defined
'''

class Config(object):
  def __init__(self, 
               # general
               data_label=None, 
               qcd_label=None, 
               signal_labels=None, 
               points_label=None,
               sample_type=None,
               tree_name=None, 
               reweighting_strategy=None,
               categories_label=None, 
               baseline_selection_label=None, 
               do_cutbased=None,
               do_mva=None,
               training_label=None,
               do_parametric=None,
               cut_score=None,
               quantities_label=None, 
               weight_label=None,
               qcd_white_list=None,
               add_weight_hlt=None,
               add_weight_pu=None,
               add_weight_muid=None,
               branch_weight_hlt=None,
               branch_weight_puqcd=None,
               branch_weight_pusig=None,
               branch_weight_mu0id=None,
               branch_weight_muid=None,

               # for plotter
               plot_CR=None, 
               plot_SR=None, 
               plot_dataSig=None,
               plot_ratio=None,
               do_shape=None,
               do_luminorm=None,
               do_stack=None,
               do_log=None, 
               add_overflow=None,
               add_CMSlabel=None,
               add_lumilabel=None,
               CMStag=None,
               do_tdrstyle=None,

               # for datacards
               ABCD_label=None,
               do_ABCD=None,
               do_ABCDHybrid=None,
               do_TF=None,
               do_realData=None,
               do_counting=None,
               do_shape_analysis=None,
               do_shape_TH1=None,
               use_discrete_profiling=None,
               signal_model_label=None, 
               background_model_label=None,
               do_binned_fit=None, 
               do_blind=None, 
               mass_window_size=None,
               fit_window_size=None,
               nbins=None, 
               plot_pulls=None,
               do_categories=None,
               add_Bc=None,
               plot_prefit=None,
               lumi_target=None,
               sigma_B=None,
               lhe_efficiency=None,
               sigma_mult_window=None,

               # for limits
               run_blind=None,
  ):

    self.data_label = data_label
    self.qcd_label = qcd_label
    self.signal_labels = signal_labels
    self.points_label = points_label
    self.sample_type = sample_type
    self.tree_name = tree_name
    self.categories_label = categories_label
    self.reweighting_strategy = reweighting_strategy
    self.baseline_selection_label = baseline_selection_label
    self.do_cutbased = do_cutbased
    self.do_mva = do_mva
    self.training_label = training_label
    self.do_parametric = do_parametric
    self.cut_score = cut_score
    self.quantities_label = quantities_label
    self.weight_label = weight_label
    self.qcd_white_list = qcd_white_list
    self.add_weight_hlt = add_weight_hlt
    self.add_weight_pu = add_weight_pu
    self.add_weight_muid = add_weight_muid
    self.branch_weight_hlt = branch_weight_hlt
    self.branch_weight_puqcd = branch_weight_puqcd
    self.branch_weight_pusig = branch_weight_pusig
    self.branch_weight_mu0id = branch_weight_mu0id
    self.branch_weight_muid = branch_weight_muid
    self.plot_CR = plot_CR
    self.plot_SR = plot_SR
    self.plot_dataSig = plot_dataSig
    self.plot_ratio = plot_ratio
    self.do_shape = do_shape
    self.do_luminorm = do_luminorm
    self.do_stack = do_stack
    self.do_log = do_log
    self.add_overflow = add_overflow
    self.add_CMSlabel = add_CMSlabel
    self.add_lumilabel = add_lumilabel
    self.CMStag = CMStag
    self.do_tdrstyle = do_tdrstyle
    self.ABCD_label = ABCD_label
    self.do_ABCD = do_ABCD
    self.do_ABCDHybrid = do_ABCDHybrid
    self.do_TF = do_TF
    self.do_realData = do_realData
    self.do_counting = do_counting
    self.do_shape_analysis = do_shape_analysis
    self.do_shape_TH1 = do_shape_TH1
    self.use_discrete_profiling = use_discrete_profiling
    self.signal_model_label = signal_model_label
    self.background_model_label = background_model_label
    self.do_binned_fit = do_binned_fit
    self.do_blind = do_blind
    self.mass_window_size = mass_window_size
    self.fit_window_size = fit_window_size
    self.nbins = nbins
    self.plot_pulls = plot_pulls
    self.do_categories = do_categories
    self.add_Bc = add_Bc
    self.plot_prefit = plot_prefit
    self.lumi_target = lumi_target
    self.sigma_B = sigma_B
    self.lhe_efficiency = lhe_efficiency
    self.sigma_mult_window = sigma_mult_window
    self.run_blind = run_blind


  def checkConfig(self):
    import ROOT 

    from samples import data_samples, qcd_samples, signal_samples
    if self.data_label != None and self.data_label not in data_samples.keys():
      raise RuntimeError('The data sample label (data_label) is not valid. Please choose amongst {}.'.format(data_samples.keys()))
    else: 
      print '       ---> Data sample label OK'
      for sample in data_samples[self.data_label]:
        if not ROOT.TFile.Open(sample.filename, 'READ'):
          raise RuntimeError('The sample "{}" does not exist. Please check the filename.'.format(sample.filename))
        else:
          if not ROOT.TFile.Open(sample.filename, 'READ').Get(self.tree_name):
            raise RuntimeError('The tree "{}" does not exist in "{}". Please check the tree name.'.format(self.tree_name, sample.filename))
    print '       ---> Data samples OK'

    if self.qcd_label != None and self.qcd_label not in qcd_samples.keys():
      raise RuntimeError('The qcd sample label (qcd_label) is not valid. Please choose amongst {}.'.format(qcd_samples.keys()))
    else: 
      print '       ---> QCD sample label OK'
      if self.qcd_label != None:
        for sample in qcd_samples[self.qcd_label]:
          if 'pt' not in sample.label or 'to' not in sample.label:
            raise RuntimeError('Please make sure that the qcd sample label contains "ptXXtoYY" for the white list to work')
          if not ROOT.TFile.Open(sample.filename, 'READ'):
            raise RuntimeError('The sample "{}" does not exist. Please check the filename.'.format(sample.filename))
          else:
            if not ROOT.TFile.Open(sample.filename, 'READ').Get(self.tree_name):
              raise RuntimeError('The tree "{}" does not exist in "{}". Please check the tree name.'.format(self.tree_name, sample.filename))
      print '       ---> QCD samples OK'

    for signal_label in self.signal_labels:
      if signal_label != None and signal_label not in signal_samples.keys():
        raise RuntimeError('The signal sample label "{}" is not valid. Please choose amongst {}.'.format(signal_label, signal_samples.keys()))
    print '       ---> Signal sample label OK'
    for signal_label in self.signal_labels:
        for sample in signal_samples[signal_label]:
          if not ROOT.TFile.Open(sample.filename, 'READ'):
            raise RuntimeError('The sample "{}" does not exist. Please check the filename.'.format(sample.filename))
          else:
            if not ROOT.TFile.Open(sample.filename, 'READ').Get(self.tree_name):
              raise RuntimeError('The tree "{}" does not exist in "{}". Please check the tree name.'.format(self.tree_name, sample.filename))
    print '       ---> Signal samples OK'

    from points import points
    if self.points_label != None and self.points_label not in points.keys():
      raise RuntimeError('The ctau points label (points_label) is not valid. Please choose amongst {}.'.format(points.keys()))
    else: 
      print '       ---> Lifetime points label OK'

    from qcd_white_list import white_list
    if self.qcd_white_list != None and self.qcd_white_list not in white_list.keys():
      raise RuntimeError('The qcd white list label (qcd_white_list) is not valid. Please choose amongst {}.'.format(white_list.keys()))
    else: 
      print '       ---> QCD pthat white list OK'

    reweighting_strategies = ['inclusive', 'partial_inclusive', 'exclusive_fromlargerctau', 'exclusive_fromsmallerctau']
    if self.reweighting_strategy not in reweighting_strategies:
      raise RuntimeError('Unknown reweighting strategy "{}". Please choose amongst {}'.format(self.reweighting_strategy, reweighting_strategies))
    else:
      print '       ---> Reweighting strategy OK'

    from categories import categories
    if self.categories_label != None and self.categories_label not in categories.keys():
      raise RuntimeError('The categories label (categories_label) is not valid. Please choose amongst {}.'.format(categories.keys()))
    else: 
      print '       ---> Categories label OK'

    from baseline_selection import selection
    if self.baseline_selection_label != None and self.baseline_selection_label not in selection.keys():
      raise RuntimeError('The selection label (baseline_selection_label) is not valid. Please choose amongst {}.'.format(selection.keys()))
    else: 
      print '       ---> Baseline selection label OK'

    if self.do_cutbased + self.do_mva != 1:
      raise RuntimeError('Choose either the cutbased or mva selection methods (do_cutbased or do_mva)')
    else:
      print '       ---> Selection method OK'

    from os import path
    if self.do_mva:
      if not path.exists('./scripts/mva/outputs/{}'.format(self.training_label)):
        raise RuntimeError('The training "{}" was not found. Please check'.format(self.training_label))
      else:
        print '       ---> Training label OK'

    from quantity import quantities
    for quantity_label in self.quantities_label:
      if quantity_label != None and quantity_label not in quantities.keys():
        raise RuntimeError('The quantities label (quantities_label) is not valid. Please choose amongst {}.'.format(quantities.keys()))
    print '       ---> Quantities label OK'

    if self.add_weight_hlt and self.branch_weight_hlt == None:
      raise RuntimeError('Please indicate the name of the branch for the hlt weight (branch_weight_hlt)')
    #elif self.add_weight_hlt and self.branch_weight_hlt not in ['weight_hlt_A1', 'weight_hlt_A1_6', 'weight_hlt_A1_6_B1', 'weight_hlt_HLT_Mu9_IP6_A1_6']:
    #  raise RuntimeError('Unrecognised branch "{}" for the hlt weight. Please check.'.format(self.branch_weight_hlt))
    else:
      print '       ---> HLT weight OK'

    #if self.add_weight_pu and self.branch_weight_pu == None:
    #  raise RuntimeError('Please indicate the name of the branch for the pile-up weight (branch_weight_pu)')
    #elif self.add_weight_pu and self.branch_weight_pu not in ['weight_pu_qcd_ntrueint', 'weight_pu_qcd_ntrueint_weighted']:
    #  raise RuntimeError('Unrecognised branch "{}" for the pile-up weight. Please check.'.format(self.branch_weight_pu))
    #else:
    #  print '       ---> Pile-up weight OK'

    if self.add_weight_muid and (self.branch_weight_muid == None or self.branch_weight_mu0id == None):
      raise RuntimeError('Please indicate the name of the branch for the muid and mu0id weights (branch_weight_muid and branch_weight_mu0id)')
    else:
      print '       ---> Muon ID weight OK'

    if self.do_shape == self.do_luminorm:
      raise RuntimeError('Invalid arguments for --do_shape ({}) and --do_luminorm ({}) \nPlease only set exactly one option to True at a time'.format(self.do_shape, self.do_luminorm))
    else: 
      print '       ---> Plotter request OK'

    if self.do_ABCD + self.do_ABCDHybrid + self.do_TF + self.do_realData != 1:
      raise RuntimeError('Please choose only one of the following background estimation method "[do_ABCD, do_ABCDHybrid, do_TF, do_realData]"')

    if self.sigma_B not in [472.8e9, 327e9]:
      raise RuntimeError('The value of sigma_B not in [472.8e9, 327e9]. Please check')

    print '       ---> Datacard request OK'

    if self.do_counting + self.do_shape_analysis + self.do_shape_TH1 != 1:
      raise RuntimeError("Please specify which analysis strategy to use among ['do_counting', 'do_shape_analysis', 'do_shape_TH1']")

    signal_model_list = ['voigtian']
    if self.do_shape_analysis and self.signal_model_label not in signal_model_list:
      raise RuntimeError('Please make sure to enter a valid signal model. Choose among {}'.format(signal_model_list))

    background_model_list = ['chebychev']
    if self.do_shape_analysis  and not self.use_discrete_profiling and self.background_model_label not in background_model_list:
      raise RuntimeError('Please make sure to enter a valid background model. Choose among {}'.format(background_model_list))
      

  # check that either lumi or shape, or can process both in parallel?


