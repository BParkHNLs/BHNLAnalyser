
class TernaryStyle(object):
  '''
    Class that sets the style of the ternary plots

    Parameters:
      - fontsize: sets generic fontsize
      - log: plot z axis (colorbar) in log scale
      - fixed_range: True fixes the range between [range_min, range_max], 
                          False sets the range to the min/max values for that given mass
      - range_min, range_max: values of the range if fixed_range is True
      - scientific: True sets scientific notation in z axis (colorbar)
      - exponent: if scientific is True, sets the exponent to which the values in the z-axis are reported
                  format: 1e-X
  '''
  def __init__(self, fontsize=10., log=False, scientific=False, fixed_range=False, range_min=1e-5, range_max=1e-3, exponent=1e-4):
    self.fontsize = fontsize
    self.log = log
    self.fixed_range = fixed_range
    self.range_min = range_min
    self.range_max = range_max
    self.scientific = scientific
    self.exponent = exponent

    if self.log and not self.scientific: 
      print 'WARNING: options log(True) and scientific(False) is not well supported'

