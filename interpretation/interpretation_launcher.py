import os
import sys
from os import path
sys.path.append('../objects')
sys.path.append('../scripts')
sys.path.append('../limits')
from coupling_scenarios import coupling_scenarios
from reweight_datacard import DatacardReweighter
from combine_flavour_channels import FlavourChannelsCombinator
from produce_limits import LimitProducer
from tools import Tools


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Run the flavour combined limits for the given coupling scenario', add_help=True)
  parser.add_argument('--homedir'       , type=str , dest='homedir'               , help='name of the homedir'                               , default=None)
  parser.add_argument('--outdirlabel'   , type=str , dest='outdirlabel'           , help='name of the outdir'                                , default=None)
  parser.add_argument('--subdirlabel'   , type=str , dest='subdirlabel'           , help='name of the subdir'                                , default=None)
  parser.add_argument('--mass'          , type=str , dest='mass'                  , help='mass'                                              , default='1.0')
  parser.add_argument('--ctau'          , type=str , dest='ctau'                  , help='ctau'                                              , default=None)
  parser.add_argument('--fe'            , type=str , dest='fe'                    , help='electron coupling fraction'                        , default='1.0')
  parser.add_argument('--fu'            , type=str , dest='fu'                    , help='muon coupling fraction'                            , default='1.0')
  parser.add_argument('--ft'            , type=str , dest='ft'                    , help='tau coupling fraction'                             , default='1.0')
  parser.add_argument('--muon_label'    , type=str , dest='muon_label'            , help='muon_label'                                        , default='')
  parser.add_argument('--electron_label', type=str , dest='electron_label'        , help='electron_label'                                    , default='')
  parser.add_argument('--use_discrete_profiling'   , dest='use_discrete_profiling', help='use discrete profiling method', action='store_true', default=False)
  parser.add_argument('--do_blind'                 , dest='do_blind'              , help='run blinded?'                 , action='store_true', default=False)
  return parser.parse_args()


class InterpretationLauncher(object):
  '''
    This class monitors the signal normalisation reweighting to a given coupling scenario,
    the combination of the datacards between the two flavour channels and the limit production
  '''
  def __init__(self, mass, ctau, fe, fu, ft, muon_label, electron_label, homedir, outdirlabel, subdirlabel, do_blind, use_discrete_profiling):
    self.tools = Tools()
    self.mass = float(mass)
    self.ctau = float(ctau)
    self.fe = float(fe)
    self.fu = float(fu)
    self.ft = float(ft)
    self.muon_label = muon_label
    self.electron_label = electron_label
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.do_blind = do_blind
    self.use_discrete_profiling = use_discrete_profiling

    self.templatename_muon = 'datacard_combined_Majorana_m_{mass}_ctau_{ctau}_v2_{v2}.txt' 
    #if self.mass not in [3., 4.5]:
    #  self.templatename_electron = 'HNL_m_{mass}_ctau{ctau}_PF_combined.txt'
    #else:
    #  self.templatename_electron = 'HNL_m_{mass}_ctau{ctau}_PF_combined_B_Bc.txt'
    self.templatename_electron = 'HNL_m_{mass}_ctau{ctau}_PF_combined.txt'

    #self.path_motherdir = '/eos/home-a/anlyon/BHNLDatacards'
    self.path_motherdir = '/work/anlyon/BHNLDatacards/BHNLDatacards/'


  def getSignalCoupling(self, mass, ctau):
    signal_v2 = self.tools.getVV(mass=mass, ctau=ctau, ismaj=True)
    signal_coupling = self.tools.getCouplingLabel(signal_v2)

    return signal_coupling


  def checkDatacard(self, datacard_name, flavour_channel):
    if flavour_channel not in ['muon', 'electron']:
      raise RuntimeError('Unknown flavour channel. Please choose amongst ["muon", "electron"]')

    if flavour_channel == 'muon':
      datacard = '{}/muon/{}/{}'.format(self.path_motherdir, self.muon_label, datacard_name)
    elif flavour_channel == 'electron':
      datacard = '{}/electron/{}/{}'.format(self.path_motherdir, self.electron_label, datacard_name)

    if not path.exists(datacard):
      print ' -> Datacard "{}" not found'.format(datacard)
      return False
    else:
      return True


  def process(self):
    print '\n --------------------------'
    print ' Interpretation Launcher  '
    print ' --------------------------'

    #for coupling_scenario in self.coupling_scenarios:
    print '\n### '
    print 'coupling scenario: (fe={}, fu={}, ft={})'.format(self.fe, self.fu, self.ft)
    print '### '

    fe = str(round(self.fe, 1)).replace('.', 'p')
    fu = str(round(self.fu, 1)).replace('.', 'p')
    ft = str(round(self.ft, 1)).replace('.', 'p')

    self.v2 = self.getSignalCoupling(mass=self.mass, ctau=self.ctau)

    datacard_name_muon = self.templatename_muon.format(mass=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.v2).replace('.', 'p').replace('-', 'm'))
    if not self.checkDatacard(datacard_name=datacard_name_muon, flavour_channel='muon'):
      raise RuntimeError('Muon datacard "{}" not found'.format(datacard_name_muon))

    datacard_name_electron = self.templatename_electron.format(mass='{:.2f}'.format(self.mass).replace('.','p'), ctau='{:.3f}'.format(self.ctau).replace('.','p'))

    if not self.checkDatacard(datacard_name=datacard_name_electron, flavour_channel='electron'):
      raise RuntimeError('Electron datacard "{}" not found'.format(datacard_name_electron))

    print '\n -> reweight signal rate in the muon datacard'
    reweighter_muon = DatacardReweighter(
        datacard_name = datacard_name_muon,
        mass = self.mass,
        ctau = self.ctau,
        fe = self.fe,
        fu = self.fu,
        ft = self.ft,
        path_motherdir = self.path_motherdir,
        indirlabel = self.muon_label, 
        outdirlabel = self.outdirlabel,
        subdirlabel = self.subdirlabel,
        flavour_channel = 'muon',
        )
    reweighter_muon.process()

    print '\n -> reweight signal rate in the electron datacard'
    reweighter_electron = DatacardReweighter(
        datacard_name = datacard_name_electron,
        mass = self.mass,
        ctau = self.ctau,
        fe = self.fe,
        fu = self.fu,
        ft = self.ft,
        path_motherdir = self.path_motherdir,
        indirlabel = self.electron_label, 
        outdirlabel = self.outdirlabel,
        subdirlabel = self.subdirlabel,
        flavour_channel = 'electron',
        )
    reweighter_electron.process()

    print '\n -> combine the datacards between the flavour channels'

    combinator = FlavourChannelsCombinator(
        datacard_name_muon = datacard_name_muon,
        datacard_name_electron = datacard_name_electron,
        fe = self.fe,
        fu = self.fu,
        ft = self.ft,
        outdirlabel = self.outdirlabel,
        subdirlabel = self.subdirlabel,
        )
    combinator.process()

    print '\n -> produce the limit'

    limit_producer = LimitProducer(
        mass = self.mass,
        ctau = self.ctau,
        fe = self.fe,
        fu = self.fu,
        ft = self.ft,
        homedir = self.homedir,
        indirlabel = './outputs/{}/{}/weighted_datacards/coupling_{}_{}_{}/flavour_combined'.format(self.outdirlabel, self.subdirlabel, fe, fu, ft),
        outdirlabel = self.outdirlabel,
        subdirlabel = self.subdirlabel,
        do_blind = self.do_blind,
        use_discrete_profiling = self.use_discrete_profiling,
        scenario = 'Majorana',
        )
    limit_producer.process()

    print '\n -- Done --'



if __name__ == '__main__':

  opt = getOptions()
  launcher = InterpretationLauncher(
      mass = opt.mass,
      ctau = opt.ctau,
      fe = opt.fe,
      fu = opt.fu,
      ft = opt.ft,
      muon_label = opt.muon_label,
      electron_label = opt.electron_label,
      homedir = opt.homedir,
      outdirlabel = opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      use_discrete_profiling = True if opt.use_discrete_profiling else False,
      do_blind = True if opt.do_blind else False,
      )
  launcher.process()

