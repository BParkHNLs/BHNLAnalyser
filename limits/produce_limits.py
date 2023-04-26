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
  parser.add_argument('--scenario'   , type=str , dest='scenario'              , help='signal under consideration'                        , default='Majorana', choices=['Majorana', 'Dirac'])
  parser.add_argument('--use_discrete_profiling', dest='use_discrete_profiling', help='use discrete profiling method', action='store_true', default=False)
  parser.add_argument('--do_blind'               , dest='do_blind'             , help='run blinded?'                 , action='store_true', default=False)
  return parser.parse_args()


class LimitProducer(object):
  def __init__(self, mass, ctau, scenario, homedir, indirlabel, outdirlabel, subdirlabel, do_blind, use_discrete_profiling, fe=None, fu=None, ft=None):
    self.tools = Tools()
    self.mass = mass
    self.ctau = ctau
    self.scenario = scenario
    self.homedir = homedir
    self.indirlabel = indirlabel
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.do_blind = do_blind
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
    self.coupling = self.getCouplingLabel(v2)

    # produce limits
    command = 'combine -M AsymptoticLimits {i}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt'.format(i=self.indirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario)
    if self.do_blind:
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

    result_file_name = '{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(outputdir, self.scenario, str(self.mass), str(self.ctau), str(self.coupling)) 
    with open(result_file_name, 'w') as ff:
        print >> ff, results

    # produce prefit plots
    os.system('mkdir -p  {}/outputs/{}/datacards/{}/prefit_postfit_plots_{}_m_{}_ctau_{}_v2_{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.scenario, str(self.mass).replace('.', 'p'), str(self.ctau).replace('.', 'p'), str(self.coupling).replace('.', 'p').replace('-', 'm')))
    os.system('cp ../data/index.php {}/outputs/{}/datacards/{}/prefit_postfit_plots_{}_m_{}_ctau_{}_v2_{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.scenario, str(self.mass).replace('.', 'p'), str(self.ctau).replace('.', 'p'), str(self.coupling).replace('.', 'p').replace('-', 'm')))
    command_prefit_plot = 'combine -M FitDiagnostics {i}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt --plots --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams --out {hm}/outputs/{out}/datacards/{sub}/prefit_postfit_plots_{sc}_m_{m}_ctau_{ctau}_v2_{v2}'.format(i=self.indirlabel, hm=self.homedir, out=self.outdirlabel, sub=self.subdirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario)
    os.system(command_prefit_plot)

    # compute p-value
    command_pvalue = 'combine -M Significance {i}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt --pval'.format(i=self.indirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario)
    if self.use_discrete_profiling:
      command_pvalue += ' --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams' 
    
    print '\t\t',command_pvalue
     
    results_pvalue = subprocess.check_output(command_pvalue.split())
    
    result_pvalue_file_name = '{}/pvalue_{}_m_{}_ctau_{}_v2_{}.txt'.format(outputdir, self.scenario, str(self.mass), str(self.ctau), str(self.coupling)) 
    with open(result_pvalue_file_name, 'w') as ff:
        print >> ff, results_pvalue


if __name__ == "__main__":

  opt = getOptions()
  producer = LimitProducer(
      mass = opt.mass,
      ctau = opt.ctau,
      scenario = opt.scenario,
      homedir = opt.homedir,
      indirlabel = opt.indirlabel,
      outdirlabel = opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      do_blind = opt.do_blind,
      use_discrete_profiling = opt.use_discrete_profiling,
      )
  producer.process()


