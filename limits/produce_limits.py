import sys 
import os
from os import path
import re
import subprocess
from decimal import Decimal
from tools import Tools


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Run limit on a single mass/coupling point', add_help=True)
  parser.add_argument('--homedir'    , type=str, dest='homedir'                , help='name of the homedir'                               , default=None)
  parser.add_argument('--indirlabel' , type=str , dest='indirlabel'            , help='name of the indir'                                 , default=None)
  parser.add_argument('--outdirlabel', type=str , dest='outdirlabel'           , help='name of the outdir'                                , default=None)
  parser.add_argument('--subdirlabel', type=str , dest='subdirlabel'           , help='name of the subdir'                                , default=None)
  parser.add_argument('--mass'       , type=str , dest='mass'                  , help='mass'                                              , default='1.0')
  parser.add_argument('--ctau'       , type=str , dest='ctau'                  , help='ctau'                                              , default=None)
  parser.add_argument('--use_discrete_profiling', dest='use_discrete_profiling', help='use discrete profiling method', action='store_true', default=False)
  parser.add_argument('--run_blind'             , dest='run_blind'             , help='run blinded?'                 , action='store_true', default=False)
  return parser.parse_args()

  # add FitDiagnostics?
  # combine -M FitDiagnostics datacard_bhnl_m_2p0_v2_1p7em03_cat_lxysiggt150_SS.txt --plots --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams


class LimitProducer(object):
  #def __init__(self, opt):
  def __init__(self, mass, ctau, homedir, indirlabel, outdirlabel, subdirlabel, run_blind, use_discrete_profiling, fe=None, fu=None, ft=None):
    self.tools = Tools()
    self.mass = mass
    self.ctau = ctau
    self.homedir = homedir
    self.indirlabel = indirlabel
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.run_blind = run_blind
    self.use_discrete_profiling = use_discrete_profiling
    self.fe = fe
    self.fu = fu
    self.ft = ft


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)


  def process(self):
    v2 = self.tools.getVV(mass=float(self.mass), ctau=float(self.ctau), ismaj=True)
    coupling = self.getCouplingLabel(v2)

    command = 'combine -M AsymptoticLimits {i}/datacard_combined_m_{m}_v2_{c}.txt'.format(i=self.indirlabel, m=str(self.mass).replace('.', 'p'), c=str(coupling).replace('.', 'p').replace('-', 'm'))
    if self.run_blind:
      command += ' --run blind'
    if self.use_discrete_profiling:
      command += ' --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams' 
    
    print '\t\t',command
     
    results = subprocess.check_output(command.split())
    
    if self.fe == None and self.fu == None and self.ft == None:
      outputdir = '{}/outputs/{}/limits/{}/results'.format(self.homedir, self.outdirlabel, self.subdirlabel)
    else:
      fe = str(round(self.fe, 1)).replace('.', 'p')
      fu = str(round(self.fu, 1)).replace('.', 'p')
      ft = str(round(self.ft, 1)).replace('.', 'p')
      outputdir = '{}/outputs/{}/limits/{}/results_{}_{}_{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))    

    result_file_name = '{}/result_m_{}_v2_{}.txt'.format(outputdir, str(self.mass), str(coupling)) 
    with open(result_file_name, 'w') as ff:
        print >> ff, results


if __name__ == "__main__":

  opt = getOptions()
  producer = LimitProducer(
      mass = opt.mass,
      ctau = opt.ctau,
      homedir = opt.homedir,
      indirlabel = opt.indirlabel,
      outdirlabel = opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      run_blind = opt.run_blind,
      use_discrete_profiling = opt.use_discrete_profiling,
      )
  producer.process()


