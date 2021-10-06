'''
  This class defines extra selection to be applied at analysis level
  Both the 'flat' and 'nano' syntaxes can be reported
'''

class Selection(object):
  def __init__(self, flat='', nano=''):
    self.flat = flat
    self.nano = nano


selection = {}
selection['standard'] = Selection(
    flat = ' && '.join([
      'mu_isdsa !=1',
      'trgmu_softid == 1',
      'mu_looseid == 1',
      'mu_intimemuon == 1',
      'mu_trackerhighpurityflag == 1'
    ])
)
