
class Config(object):
  def __init__(self, 
               data_label=None, 
               qcd_label=None, 
               signal_label=None, 
               sample_type=None,
               tree_name=None, 
               categories_label=None, 
               selection_label=None, 
               quantities_label=None, 
               weight_label=None,
               qcd_white_list=None,
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
               CMStag=None,
               # for datacards
               ABCD_label=None,
               do_categories=None,
               add_Bc=None,
  ):

    self.data_label = data_label
    self.qcd_label = qcd_label
    self.signal_label = signal_label
    self.sample_type = sample_type
    self.tree_name = tree_name
    self.categories_label = categories_label
    self.selection_label = selection_label
    self.quantities_label = quantities_label
    self.weight_label = weight_label
    self.qcd_white_list = qcd_white_list
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
    self.CMStag = CMStag
    self.ABCD_label = ABCD_label
    self.do_categories = do_categories
    self.add_Bc = add_Bc
    # add lumi, sigma_B? 


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
      for sample in qcd_samples[self.qcd_label]:
        if 'pt' not in sample.label or 'to' not in sample.label:
          raise RuntimeError('Please make sure that the qcd sample label contains "ptXXtoYY" for the white list to work')
        if not ROOT.TFile.Open(sample.filename, 'READ'):
          raise RuntimeError('The sample "{}" does not exist. Please check the filename.'.format(sample.filename))
        else:
          if not ROOT.TFile.Open(sample.filename, 'READ').Get(self.tree_name):
            raise RuntimeError('The tree "{}" does not exist in "{}". Please check the tree name.'.format(self.tree_name, sample.filename))
    print '       ---> QCD samples OK'

    if self.signal_label != None and self.signal_label not in signal_samples.keys():
      raise RuntimeError('The signal sample label (signal_label) is not valid. Please choose amongst {}.'.format(signal_samples.keys()))
    else: 
      print '       ---> Signal sample label OK'
      for sample in signal_samples[self.signal_label]:
        if not ROOT.TFile.Open(sample.filename, 'READ'):
          raise RuntimeError('The sample "{}" does not exist. Please check the filename.'.format(sample.filename))
        else:
          if not ROOT.TFile.Open(sample.filename, 'READ').Get(self.tree_name):
            raise RuntimeError('The tree "{}" does not exist in "{}". Please check the tree name.'.format(self.tree_name, sample.filename))
    print '       ---> Signal samples OK'

    from qcd_white_list import white_list
    if self.qcd_white_list != None and self.qcd_white_list not in white_list.keys():
      raise RuntimeError('The qcd white list label (qcd_white_list) is not valid. Please choose amongst {}.'.format(white_list.keys()))
    else: 
      print '       ---> QCD pthat white list OK'

    from categories import categories
    if self.categories_label != None and self.categories_label not in categories.keys():
      raise RuntimeError('The categories label (categories_label) is not valid. Please choose amongst {}.'.format(categories.keys()))
    else: 
      print '       ---> Categories label OK'

    from selection import selection
    if self.selection_label != None and self.selection_label not in selection.keys():
      raise RuntimeError('The selection label (selection_label) is not valid. Please choose amongst {}.'.format(selection.keys()))
    else: 
      print '       ---> Selection label OK'

    from quantity import quantities
    if self.quantities_label != None and self.quantities_label not in quantities.keys():
      raise RuntimeError('The quantities label (quantities_label) is not valid. Please choose amongst {}.'.format(quantities.keys()))
    else: 
      print '       ---> Quantities label OK'

  # check that either lumi or shape, or can process both in parallel?


