'''
  Class that defines the coupling flavour fractions (fe, fu, ft)
'''
class CouplingScenario(object):
  def __init__(self, fe, fu, ft):
    self.fe = float(fe)
    self.fu = float(fu)
    self.ft = float(ft)

    #if float(self.fe + self.fu + self.ft) != 1.0:
    #  print self.fe + self.fu + self.ft
    #  raise RuntimeError('Invalid coupling scenario (fe={}, fu={}, ft={}). Make sure that the sum of the couplings is equal to 1.'.format(self.fe, self.fu, self.ft))



coupling_scenarios = {}

coupling_scenarios['muon_only'] = (
    CouplingScenario(fe=0., fu=1., ft=0.),
    )

coupling_scenarios['benchmark'] = (
    CouplingScenario(fe=0., fu=1., ft=0.),
    CouplingScenario(fe=0., fu=0.5, ft=0.5),
    CouplingScenario(fe=0.5, fu=0.5, ft=0.),
    CouplingScenario(fe=1./3., fu=1./3., ft=1./3.),
    )

coupling_scenarios['triangle'] = (
    CouplingScenario(fe=0., fu=1., ft=0.),
    CouplingScenario(fe=0., fu=0.9, ft=0.1),
    CouplingScenario(fe=0., fu=0.8, ft=0.2),
    CouplingScenario(fe=0., fu=0.7, ft=0.3),
    CouplingScenario(fe=0., fu=0.6, ft=0.4),
    CouplingScenario(fe=0., fu=0.5, ft=0.5),
    CouplingScenario(fe=0., fu=0.4, ft=0.6),
    CouplingScenario(fe=0., fu=0.3, ft=0.7),
    CouplingScenario(fe=0., fu=0.2, ft=0.8),
    CouplingScenario(fe=0., fu=0.1, ft=0.9),
    CouplingScenario(fe=0., fu=0., ft=1.),
    CouplingScenario(fe=0.1, fu=0.9, ft=0.),
    CouplingScenario(fe=0.1, fu=0.8, ft=0.1),
    CouplingScenario(fe=0.1, fu=0.7, ft=0.2),
    CouplingScenario(fe=0.1, fu=0.6, ft=0.3),
    CouplingScenario(fe=0.1, fu=0.5, ft=0.4),
    CouplingScenario(fe=0.1, fu=0.4, ft=0.5),
    CouplingScenario(fe=0.1, fu=0.3, ft=0.6),
    CouplingScenario(fe=0.1, fu=0.2, ft=0.7),
    CouplingScenario(fe=0.1, fu=0.1, ft=0.8),
    CouplingScenario(fe=0.1, fu=0., ft=0.9),
    CouplingScenario(fe=0.2, fu=0.8, ft=0.),
    CouplingScenario(fe=0.2, fu=0.7, ft=0.1),
    CouplingScenario(fe=0.2, fu=0.6, ft=0.2),
    CouplingScenario(fe=0.2, fu=0.5, ft=0.3),
    CouplingScenario(fe=0.2, fu=0.4, ft=0.4),
    CouplingScenario(fe=0.2, fu=0.3, ft=0.5),
    CouplingScenario(fe=0.2, fu=0.2, ft=0.6),
    CouplingScenario(fe=0.2, fu=0.1, ft=0.7),
    CouplingScenario(fe=0.2, fu=0., ft=0.8),
    CouplingScenario(fe=0.3, fu=0.7, ft=0.),
    CouplingScenario(fe=0.3, fu=0.6, ft=0.1),
    CouplingScenario(fe=0.3, fu=0.5, ft=0.2),
    CouplingScenario(fe=0.3, fu=0.4, ft=0.3),
    CouplingScenario(fe=0.3, fu=0.3, ft=0.4),
    CouplingScenario(fe=0.3, fu=0.2, ft=0.5),
    CouplingScenario(fe=0.3, fu=0.1, ft=0.6),
    CouplingScenario(fe=0.3, fu=0., ft=0.7),
    CouplingScenario(fe=0.4, fu=0.6, ft=0.),
    CouplingScenario(fe=0.4, fu=0.5, ft=0.1),
    CouplingScenario(fe=0.4, fu=0.4, ft=0.2),
    CouplingScenario(fe=0.4, fu=0.3, ft=0.3),
    CouplingScenario(fe=0.4, fu=0.2, ft=0.4),
    CouplingScenario(fe=0.4, fu=0.1, ft=0.5),
    CouplingScenario(fe=0.4, fu=0., ft=0.6),
    CouplingScenario(fe=0.5, fu=0.5, ft=0.),
    CouplingScenario(fe=0.5, fu=0.4, ft=0.1),
    CouplingScenario(fe=0.5, fu=0.3, ft=0.2),
    CouplingScenario(fe=0.5, fu=0.2, ft=0.3),
    CouplingScenario(fe=0.5, fu=0.1, ft=0.4),
    CouplingScenario(fe=0.5, fu=0., ft=0.5),
    CouplingScenario(fe=0.6, fu=0.4, ft=0.),
    CouplingScenario(fe=0.6, fu=0.3, ft=0.1),
    CouplingScenario(fe=0.6, fu=0.2, ft=0.2),
    CouplingScenario(fe=0.6, fu=0.1, ft=0.3),
    CouplingScenario(fe=0.6, fu=0., ft=0.4),
    CouplingScenario(fe=0.7, fu=0.3, ft=0.),
    CouplingScenario(fe=0.7, fu=0.2, ft=0.1),
    CouplingScenario(fe=0.7, fu=0.1, ft=0.2),
    CouplingScenario(fe=0.7, fu=0., ft=0.3),
    CouplingScenario(fe=0.8, fu=0.2, ft=0.),
    CouplingScenario(fe=0.8, fu=0.1, ft=0.1),
    CouplingScenario(fe=0.8, fu=0., ft=0.2),
    CouplingScenario(fe=0.9, fu=0.1, ft=0.),
    CouplingScenario(fe=0.9, fu=0., ft=0.1),
    CouplingScenario(fe=1., fu=0., ft=0.),
    )
