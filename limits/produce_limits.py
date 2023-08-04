import sys 
import os
from os import path
import re
import subprocess
from decimal import Decimal
from tools import Tools

import sys
sys.path.append('../objects')
from categories import categories


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Run limit on a single mass/coupling point', add_help=True)
  parser.add_argument('--homedir'         , type=str, dest='homedir'               , help='name of the homedir'                               , default=None)
  parser.add_argument('--indirlabel'      , type=str, dest='indirlabel'            , help='name of the indir'                                 , default=None)
  parser.add_argument('--outdirlabel'     , type=str, dest='outdirlabel'           , help='name of the outdir'                                , default=None)
  parser.add_argument('--subdirlabel'     , type=str, dest='subdirlabel'           , help='name of the subdir'                                , default=None)
  parser.add_argument('--mass'            , type=str, dest='mass'                  , help='mass'                                              , default='1.0')
  parser.add_argument('--ctau'            , type=str, dest='ctau'                  , help='ctau'                                              , default=None)
  parser.add_argument('--scenario'        , type=str, dest='scenario'              , help='signal under consideration'                        , default='Majorana', choices=['Majorana', 'Dirac'])
  parser.add_argument('--categories_label', type=str, dest='categories_label'      , help='label of the list of categories'                   , default='standard')
  parser.add_argument('--use_discrete_profiling'    , dest='use_discrete_profiling', help='use discrete profiling method', action='store_true', default=False)
  parser.add_argument('--do_blind'                  , dest='do_blind'              , help='run blinded?'                 , action='store_true', default=False)
  return parser.parse_args()


class LimitProducer(object):
  def __init__(self, mass, ctau, scenario, homedir, indirlabel, outdirlabel, subdirlabel, do_blind, use_discrete_profiling, categories_label=None, fe=None, fu=None, ft=None, do_fitdiagnostics=True, do_pvalue=True):
    self.tools = Tools()
    self.mass = mass
    self.ctau = ctau
    self.scenario = scenario
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.do_blind = do_blind
    self.use_discrete_profiling = use_discrete_profiling
    self.fe = fe
    self.fu = fu
    self.ft = ft
    self.do_limits_percategory = False
    if categories_label != None:
      self.categories = categories[categories_label]
    self.indirlabel = indirlabel if not self.do_limits_percategory else './'
    self.do_fitdiagnostics = do_fitdiagnostics
    self.do_pvalue = do_pvalue


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)


  def process(self):
    v2 = self.tools.getVV(mass=float(self.mass), ctau=float(self.ctau), ismaj=True)
    self.coupling = self.getCouplingLabel(v2)

    if not self.do_limits_percategory:
      # create outputdir
      if self.fe == None and self.fu == None and self.ft == None:
        outputdir = '{}/outputs/{}/limits/{}/results'.format(self.homedir, self.outdirlabel, self.subdirlabel)
      else:
        fe = str(round(self.fe, 1)).replace('.', 'p')
        fu = str(round(self.fu, 1)).replace('.', 'p')
        ft = str(round(self.ft, 1)).replace('.', 'p')
        outputdir = '{}/outputs/{}/limits/{}/results_{}_{}_{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft)
      if not path.exists(outputdir):
        os.system('mkdir -p {}'.format(outputdir))    

      # produce limits
      command = 'combine -M AsymptoticLimits {i}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt'.format(i=self.indirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario)
      if self.do_blind:
        command += ' --run blind'
      if self.use_discrete_profiling:
        command += ' --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams' 
      
      print '\t\t',command
       
      results = subprocess.check_output(command.split())
      
      result_file_name = '{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(outputdir, self.scenario, str(self.mass), str(self.ctau), str(self.coupling)) 
      with open(result_file_name, 'w') as ff:
          print >> ff, results

    else: # do limits per category
      for category in self.categories:
        if 'incl' in category.label: continue
        if self.fe == None and self.fu == None and self.ft == None:
          outputdir = '{}/outputs/{}/limits/{}/results/{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, category.label)
        else:
          fe = str(round(self.fe, 1)).replace('.', 'p')
          fu = str(round(self.fu, 1)).replace('.', 'p')
          ft = str(round(self.ft, 1)).replace('.', 'p')
          outputdir = '{}/outputs/{}/limits/{}/results_{}_{}_{}/{}'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft, category.label)
        if not path.exists(outputdir):
          os.system('mkdir -p {}'.format(outputdir))    

        # produce limits
        command = 'combine -M AsymptoticLimits {i}/datacard_bhnl_m_{m}_ctau_{ctau}_v2_{v2}_cat_{cat}.txt'.format(i=self.indirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario, cat=category.label)
        if self.do_blind:
          command += ' --run blind'
        if self.use_discrete_profiling:
          command += ' --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams' 
        
        print '\t\t',command
         
        try:
          results = subprocess.check_output(command.split())
          
          result_file_name = '{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(outputdir, self.scenario, str(self.mass), str(self.ctau), str(self.coupling)) 
          with open(result_file_name, 'w') as ff:
              print >> ff, results
        except:
          continue

    # produce fit diagnostics
    if self.do_fitdiagnostics:
      command_fits = 'combine -M FitDiagnostics {i}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt --plots --cminDefaultMinimizerStrategy=0 --X-rtd MINIMIZER_freezeDisassociatedParams'.format(i=self.indirlabel, hm=self.homedir, out=self.outdirlabel, sub=self.subdirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario)
      os.system(command_fits)
      os.system('mv fitDiagnosticsTest.root {hm}/outputs/{out}/limits/{sub}/results/fitDiagnostics_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.root'.format(hm=self.homedir, out=self.outdirlabel, sub=self.subdirlabel, m=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'), v2=str(self.coupling).replace('.', 'p').replace('-', 'm'), sc=self.scenario))

    # compute p-value
    if self.do_pvalue:
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
      categories_label = opt.categories_label,
      do_blind = opt.do_blind,
      use_discrete_profiling = opt.use_discrete_profiling,
      )
  producer.process()


