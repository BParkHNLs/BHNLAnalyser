'''
Class to define the vetoes in the displaced lepton-pion invariant mass
'''

class Veto(object):
  def __init__(self, mass=None, range_min=None, range_max=None):
    self.mass = mass
    self.range_min = range_min
    self.range_max = range_max



vetoes = [
  Veto(
      mass = 1.76,
      range_min = 1.74,
      range_max = 1.80,
      ),
  Veto(
      mass = 3.686,
      range_min = 3.65,
      range_max = 3.75,
      ),
  Veto(
      mass = 3.097,
      range_min = 3.05,
      range_max = 3.15,
      ),
  ]
