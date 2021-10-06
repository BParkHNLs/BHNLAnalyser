
class Config(object):
  def __init__(self, 
               data_label='', 
               qcd_label='', 
               signal_label='', 
               sample_type='',
               tree_name='', 
               categories_label='', 
               selection_label='', 
               quantities_label='', 
               weight_label='',
               # for plotter
               plot_CR='', 
               plot_SR='', 
               plot_dataSig='',
               plot_ratio='',
               do_shape='',
               do_luminorm='',
               do_stack='',
               do_log='', 
               add_overflow='',
               add_CMSlabel='',
               CMStag='',
               # for datacards
               ABCD_label='',
               do_categories='',
               add_Bc='',
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

  # add checkConfig function?


