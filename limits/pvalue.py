import sys
import os
from os import path
import glob
import re
import pickle
import math
import numpy as np
import otherExp_limits as db 
from collections import OrderedDict
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from intersection import intersection
from utils import getMassList


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Run limit on a single mass/coupling point', add_help=True)
  parser.add_argument('--scenario'          , type=str, dest='scenario'          , help='signal under consideration'                   , default='Majorana', choices=['Majorana', 'Dirac'])
  parser.add_argument('--homedir'           , type=str, dest='homedir'           , help='name of the homedir'                          , default=None)
  parser.add_argument('--outdirlabel'       , type=str, dest='outdirlabel'       , help='name of the outdir'                           , default=None)
  parser.add_argument('--subdirlabel'       , type=str, dest='subdirlabel'       , help='name of the subdir'                           , default=None)
  parser.add_argument('--mass_whitelist'    , type=str, dest='mass_whitelist'    , help='allowed values for masses'                    , default='None')
  parser.add_argument('--mass_blacklist'    , type=str, dest='mass_blacklist'    , help='values for masses to skip'                    , default='None')
  parser.add_argument('--coupling_whitelist', type=str, dest='coupling_whitelist', help='allowed values for couplings'                 , default='None')
  parser.add_argument('--coupling_blacklist', type=str, dest='coupling_blacklist', help='values for couplings to skip'                 , default='None')
  parser.add_argument('--fe'                , type=str, dest='fe'                , help='electron coupling fraction'                   , default='1.0')
  parser.add_argument('--fu'                , type=str, dest='fu'                , help='muon coupling fraction'                       , default='1.0')
  parser.add_argument('--ft'                , type=str, dest='ft'                , help='tau coupling fraction'                        , default='1.0')
  parser.add_argument('--do_blind'                    , dest='do_blind'          , help='run blinded or unblinded', action='store_true', default=False)
  return parser.parse_args()


class LimitPlotter(object):
  def __init__(self, scenario, homedir, outdirlabel, subdirlabel, mass_whitelist, mass_blacklist, coupling_whitelist, coupling_blacklist, do_blind, fe, fu, ft):
    self.scenario = scenario
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.mass_whitelist = mass_whitelist
    self.mass_blacklist = mass_blacklist
    self.coupling_whitelist = coupling_whitelist
    self.coupling_blacklist = coupling_blacklist
    self.do_blind = do_blind
    self.fe = fe
    self.fu = fu
    self.ft = ft
    if self.fe != 'None' and self.fu != 'None' and self.ft != 'None':
      self.do_coupling_scenario = True
    else:
      self.do_coupling_scenario = False
    if self.do_coupling_scenario:
      self.fe = str(round(float(fe), 1)).replace('.', 'p')
      self.fu = str(round(float(fu), 1)).replace('.', 'p')
      self.ft = str(round(float(ft), 1)).replace('.', 'p')


  def sortList(self, input):
    return float(input)


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2)


  def get_coupling(self, couplings, values, crossing=1):
    '''
      Function that returns the coupling at which the limit intersects with 1 (crossing) in the log-log plane
    '''
    # first, search for couplings whose associated limit is the closest to up and down 1
    coupling_up = 0
    coupling_down = 1e9
    value_up = 0
    value_down = 0
    for icoupling, coupling in enumerate(couplings):
      if icoupling+1 >= len(couplings): continue
      #print '{} {} {} {}'.format(values[icoupling], values[icoupling+1], crossing, values[icoupling] > crossing and values[icoupling+1] < crossing)
      if values[icoupling] >= crossing and values[icoupling+1] <= crossing:
        coupling_up = couplings[icoupling+1]
        value_up = values[icoupling+1]
        coupling_down = couplings[icoupling]
        value_down = values[icoupling]
        break

    if abs(value_down -1) < abs(value_up -1):
      coupling = coupling_down
    else:
      coupling = coupling_up

    return coupling


  def process(self):
    #lumi =  '5.3 fb'+r'$^{-1}$'+' projected to 41.6 fb'+r'$^{-1}$'
    lumi =  '40.0 fb'+r'$^{-1}$'

    # get the files 
    if not self.do_coupling_scenario:
      pathToResults = '{}/outputs/{}/limits/{}/results/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    else:
      pathToResults = '{}/outputs/{}/limits/{}/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.fe, self.fu, self.ft) 

    fileName = 'pvalue*{}*.txt'.format(self.scenario)

    print pathToResults+fileName

    files = [f for f in glob.glob(pathToResults+fileName)]
    for f in files:
      print f
   
    # get the list of the masses from the fileNames
    masses = getMassList(files)
    masses.sort(key=self.sortList)
   
    # needed for the 2D limit plot
    limits2D = OrderedDict()

    # create directory
    plotDir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    os.system('mkdir -p {d}'.format(d=plotDir)) 

    the_masses = []
    the_pvalues = []

    for mass in masses:
      # get white/black listed mass
      if self.mass_whitelist != 'None':
        if mass not in self.mass_whitelist.split(','): continue
      
      if self.mass_blacklist != 'None':
        if mass in self.mass_blacklist.split(','): continue

      #if float(mass) < 1.2 or float(mass)>2.: continue

      print '\nmass {}'.format(mass)

      v2s       = []
      obs       = []
      minus_two = []
      minus_one = []
      central   = []
      plus_one  = []
      plus_two  = []

      for limitFile in files:
        if 'm_{}_'.format(mass) not in limitFile: continue
     
        # for each mass, get the list of the couplings and ctaus from the file name
        ctau = limitFile[limitFile.find('ctau_')+5:limitFile.find('_', limitFile.find('ctau_')+5)]
        coupling = limitFile[limitFile.find('v2_')+3:limitFile.find('.txt')]
        val_coupling = float(coupling)
      
        # get white/black listed coupling
        if self.coupling_whitelist!='None':
          if str(val_coupling) not in self.coupling_whitelist.split(','): continue 
        
        if self.coupling_blacklist!='None':
          if str(val_coupling) in self.coupling_blacklist.split(','): continue 
        
        try:
          thefile = open('{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, self.scenario, mass, ctau, coupling), 'r')
          
          # get the necessary information from the result files
          val_obs       = None
          val_minus_two = None
          val_minus_one = None
          val_central   = None
          val_plus_one  = None
          val_plus_two  = None

          content = thefile.readlines()
          for line in content:
            if 'Observed' in line:
              values = re.findall(r'\d+', line)
              val_obs = values[0] + '.' + values[1]
          obs.append(float(val_obs))
          v2s.append(val_coupling)

        except:
          print 'Cannot open {}result_{}_m_{}_v2_{}.txt'.format(pathToResults, self.scenario, mass, coupling)


      graph = zip(v2s, obs)
      graph.sort(key = lambda x : float(x[0])) # sort by coupling
    
      v2s = [jj[0] for jj in graph]
      obs = [jj[1] for jj in graph]
      
      # find the coupling the closest to the exclusion
      coupling_obs = self.get_coupling(v2s, obs)
      coupling_label = self.getCouplingLabel(coupling_obs)

      # get the corresponding p-value
      p_value_filename =  glob.glob('{}/outputs/{}/limits/{}/results/pvalue_{}_m_{}_ctau_*_v2_{}.txt'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.scenario, mass, coupling_label))[0]
      p_value_file = open(p_value_filename)
      lines = p_value_file.readlines()
      for line in lines:
        if 'p-value of background' not in line: continue
        p_value = line[line.find('p-value of background')+23:len(line)-1]
      print p_value
      the_masses.append(round(float(mass), 2))
      the_pvalues.append(float(p_value))


    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    coupling_scenario = r'(f$_{e}$={fe}, f$_{mu}$={fu}, f$_{tau}$={ft})'.format(e='e', fe=self.fe.replace('p', '.'), mu=r'\mu', fu=self.fu.replace('p', '.'), tau=r'\tau', ft=self.ft.replace('p', '.'))
    plt.plot(the_masses, the_pvalues, color='black', label='p-value', linewidth=2)
    plt.axhline(y=1.5e-1, color='red', linestyle='-')
    plt.axhline(y=2.1e-2, color='red', linestyle='-')
    plt.axhline(y=1.1e-3, color='red', linestyle='-')

    ax.text(0.97, 0.84, r'$1\sigma$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, fontweight='bold', color='red')
    ax.text(0.97, 0.56, r'$2\sigma$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, fontweight='bold', color='red')
    ax.text(0.97, 0.15, r'$3\sigma$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, fontweight='bold', color='red')

    plt.ylim(5e-4, 0.6)

    plt.yticks(fontsize=17)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yscale('log')
    plt.xticks(fontsize=17)
    plt.xlabel(r'$m_N$ (GeV)', fontsize=23)
    plt.ylabel('Local p-value', fontsize=23)

    plt.title(lumi + ' (13 TeV)', loc='right', fontsize=23)
    #plt.legend()
    if not self.do_coupling_scenario:
      name_log = 'p_value_{}'.format(self.scenario) 
    else:
      name_log = 'p_value_{}_scenario_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    print '--> {}/{}.png created'.format(plotDir, name_log)
    plt.savefig('{}/{}.pdf'.format(plotDir, name_log))
    plt.savefig('{}/{}.png'.format(plotDir, name_log))
    


if __name__ == "__main__":

  opt=getOptions()

  plotter = LimitPlotter(
      scenario = opt.scenario,
      homedir = opt.homedir,
      outdirlabel = opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      mass_whitelist = opt.mass_whitelist,
      mass_blacklist = opt.mass_blacklist,
      coupling_whitelist = opt.coupling_whitelist,
      coupling_blacklist = opt.coupling_blacklist,
      do_blind = True if opt.do_blind else False,
      fe = opt.fe,
      fu = opt.fu,
      ft = opt.ft,
      )

  plotter.process()


