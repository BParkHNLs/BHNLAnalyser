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
  parser.add_argument('--outdirlabel', type=str, dest='outdirlabel' , help='name of the outdir'               , default=None)
  parser.add_argument('--subdirlabel', type=str, dest='subdirlabel' , help='name of the subdir'               , default=None)
  parser.add_argument('--mass'       , type=str, dest='mass'        , help='mass'                             , default='1.0')
  parser.add_argument('--ctau'       , type=str, dest='ctau'        , help='ctau'                             , default=None)
  parser.add_argument('--run_blind'            , dest='run_blind'   , help='run blinded?', action='store_true', default=False)
  return parser.parse_args()


def getCouplingLabel(v2):
  coupling = "{:e}".format(v2)
  part1 = coupling[:coupling.find('e')]
  part1 = str(round(float(part1), 1))
  part2 = coupling[coupling.find('e'):]
  return (part1+part2)


if __name__ == "__main__":

        opt = getOptions()
        mass = opt.mass
        v2 = Tools().getVV(mass=float(opt.mass), ctau=float(opt.ctau), ismaj=True)
        coupling = getCouplingLabel(v2)

        command = 'combine -M AsymptoticLimits {s}/datacard_combined_m_{m}_v2_{c}.txt'.format(s=opt.subdirlabel, m=mass, c=coupling)
        if opt.run_blind:
          command += ' --run blind'
        
        print '\t\t',command
         
        results = subprocess.check_output(command.split())
        
        outputdir = './outputs/{}/limits/{}/results'.format(opt.outdirlabel, opt.subdirlabel)
        if not path.exists(outputdir):
          os.system('mkdir -p {}'.format(outputdir))    

        result_file_name = '{}/result_m_{}_v2_{}.txt'.format(outputdir, str(mass), str(coupling)) 
        with open(result_file_name, 'w') as ff:
            print >> ff, results



