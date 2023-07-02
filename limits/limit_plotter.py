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
  #parser.add_argument('--signal_type'       , type=str, dest='signal_type'       , help='signal under consideration'                   , default='majorana', choices=['majorana', 'dirac'])
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


  def get_intersection(self, couplings, values, crossing=1):
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

    #print 'coupling up: {} {}'.format(coupling_up, value_up)
    #print 'coupling down: {} {}'.format(coupling_down, value_down)

    # do not proceed if crossing is not found
    if coupling_up == 0 or coupling_down == 1e9:
      intersection = -99
    else:
      # then, search for the intersection with 1 in the log-log plane
      # powerlaw y = kx^m behaves as a linear law in the log-log plane: log(y) = mlog(x) + log(k)
      # get the slope m
      m = (math.log(value_up) - math.log(value_down)) / (math.log(coupling_up) - math.log(coupling_down))

      # and the coordinate k
      k = value_up / math.pow(coupling_up, m)

      # to finally compute the coupling where there is the intersection
      intersection = math.pow(float(crossing)/float(k), 1./float(m)) 

    return intersection


  def smooth(self, y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


  def process(self):
    #signal_type = self.signal_type 
    #lumi =  '41.6 fb'+r'$^{-1}$'
    #lumi =  '5.3 fb'+r'$^{-1}$'+' projected to 41.6 fb'+r'$^{-1}$'
    lumi =  '40.0 fb'+r'$^{-1}$'
    #lumi =  '41.5 fb'+r'$^{-1}$'

    # get the files 
    if not self.do_coupling_scenario:
      pathToResults = '{}/outputs/{}/limits/{}/results/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    else:
      pathToResults = '{}/outputs/{}/limits/{}/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.fe, self.fu, self.ft) 

    fileName = 'result*{}*.txt'.format(self.scenario)

    files = [f for f in glob.glob(pathToResults+fileName)]
   
    # get the list of the masses from the fileNames
    masses = getMassList(files)
    masses.sort(key=self.sortList)
   
    # needed for the 2D limit plot
    limits2D = OrderedDict()

    # create directory
    plotDir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    os.system('mkdir -p {d}'.format(d=plotDir)) 

    for mass in masses:
      # get white/black listed mass
      if self.mass_whitelist != 'None':
        if mass not in self.mass_whitelist.split(','): continue
      
      if self.mass_blacklist != 'None':
        if mass in self.mass_blacklist.split(','): continue

      #if float(mass) < 3. or float(mass) > 4.: continue

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
        coupling = limitFile[limitFile.rfind('v2_')+3:limitFile.find('.txt')]
        val_coupling = float(coupling)
      
        # get white/black listed coupling
        if self.coupling_whitelist!='None':
          if str(val_coupling) not in self.coupling_whitelist.split(','): continue 
        
        if self.coupling_blacklist!='None':
          if str(val_coupling) in self.coupling_blacklist.split(','): continue 
        
        try:
          thefile = open('{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, self.scenario, mass, ctau, coupling), 'r')
          #print '{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, self.scenario, mass, ctau, coupling)
          #thefile = open('{}/result_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, mass, ctau, coupling), 'r')
          #thefile = open('{}/result_m_{}_v2_{}.txt'.format(pathToResults, mass, coupling), 'r')
          
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
            elif 'Expected  2.5' in line: 
              values = re.findall(r'\d+', line)
              val_minus_two = values[2] + '.' + values[3]
            elif 'Expected 16' in line: 
              values = re.findall(r'\d+', line)
              val_minus_one = values[2] + '.' + values[3]
            elif 'Expected 50' in line: 
              values = re.findall(r'\d+', line)
              val_central = values[2] + '.' + values[3]
            elif 'Expected 84' in line: 
              values = re.findall(r'\d+', line)
              val_plus_one = values[2] + '.' + values[3]
            elif 'Expected 97.5' in line: 
              values = re.findall(r'\d+', line)
              val_plus_two = values[2] + '.' + values[3]
          
          if all([jj is None for jj in [val_minus_two, val_minus_one, val_central, val_plus_one, val_plus_two] ]): continue
          if not self.do_blind:
            if val_obs is None: 
              print 'WARNING: cannot plot unblinded if limits were produced blinded'
              print '--> Aborting'
              continue
         
          v2s.append(val_coupling)
          minus_two.append(float(val_minus_two))
          minus_one.append(float(val_minus_one))
          central.append(float(val_central))
          plus_one.append(float(val_plus_one))
          plus_two.append(float(val_plus_two))
          if not self.do_blind: obs.append(float(val_obs))

        except:
          print 'Cannot open {}result_{}_m_{}_v2_{}.txt'.format(pathToResults, self.scenario, mass, coupling)


      print '-> will plot 1D limit for mass {}'.format(mass)
      
      if not self.do_blind:
          graph = zip(v2s, minus_two, minus_one, central, plus_one, plus_two, obs)
      else:
        graph = zip(v2s, minus_two, minus_one, central, plus_one, plus_two)    
      graph.sort(key = lambda x : float(x[0])) # sort by coupling
    

      v2s       = [jj[0] for jj in graph]
      minus_two = [jj[1] for jj in graph]
      minus_one = [jj[2] for jj in graph]
      central   = [jj[3] for jj in graph]
      plus_one  = [jj[4] for jj in graph]
      plus_two  = [jj[5] for jj in graph]
      if not self.do_blind:
          obs = [jj[6] for jj in graph]
      
      print '\n-> Plots will be saved in {}'.format(plotDir)
          
      plt.clf()
      ax = plt.axes()
      print '   couplings: {}'.format(v2s)
      coupling_scenario = r'(f$_{e}$={fe}, f$_{mu}$={fu}, f$_{tau}$={ft})'.format(e='e', fe=self.fe.replace('p', '.'), mu=r'\mu', fu=self.fu.replace('p', '.'), tau=r'\tau', ft=self.ft.replace('p', '.'))
      plt.fill_between(v2s, minus_two, plus_two, color='gold', label=r'$\pm 2 \sigma$')
      plt.fill_between(v2s, minus_one, plus_one, color='forestgreen' , label=r'$\pm 1 \sigma$')
      plt.plot(v2s, central, color='red', label='central expected', linewidth=2)
      if not self.do_blind:
          plt.plot(v2s, obs, color='black', label='observed')    
      
      plt.axhline(y=1, color='black', linestyle='-')
      plt.xlabel(r'$|V|^2$')
      plt.ylabel('exclusion limit 95% CL')

      #plt.title('HNL m = %s GeV %s' %(mass, signal_type))
      if self.do_coupling_scenario:
        plt.title(r'Majorana HN, m = %s GeV, %s' %(mass, coupling_scenario))
      else:
        if str(mass) == '1': 
          mass = '1.0'
        if str(mass) == '2': 
          mass = '2.0'
        if str(mass) == '3': 
          mass = '3.0'
        #plt.title(r'Majorana HN, m = %s GeV' %(mass,))
        plt.title('Majorana HN, m = {} GeV'.format(mass))
      plt.legend()
      plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
      #ax.text(0.5, 0.5, coupling_scenario, style='normal', bbox={'facecolor': 'white', 'alpha': 1., 'pad': 10}, transform=ax.transAxes)
      #plt.xlim(1e-5, 3e-4)
      #plt.ylim(1e-1, 1e4)
      #plt.yscale('linear')
      #if self.fe == None and self.fu == None and self.ft == None:
      #  name_lin = 'limit_m_{}_lin'.format(mass.replace('.', 'p')) 
      #else:
      #  name_lin = 'limit_m_{}_scenario_{}_{}_{}_lin'.format(mass.replace('.', 'p'), self.fe, self.fu, self.ft) 
      #plt.savefig('{}/{}.pdf'.format(plotDir, name_lin))
      #plt.savefig('{}/{}.png'.format(plotDir, name_lin))
      plt.yscale('log')
      plt.xscale('log')
      if not self.do_coupling_scenario:
        name_log = 'limit_{}_m_{}_log'.format(self.scenario, mass.replace('.', 'p')) 
      else:
        name_log = 'limit_{}_m_{}_scenario_{}_{}_{}_log'.format(self.scenario, mass.replace('.', 'p'), self.fe, self.fu, self.ft) 
      print '--> {}/{}.png created'.format(plotDir, name_log)
      plt.savefig('{}/{}.pdf'.format(plotDir, name_log))
      plt.savefig('{}/{}.png'.format(plotDir, name_log))
      

      # save the crossing for 2D limits
      limits2D[mass] = OrderedDict()

      # find the intersections    
      x_minus_two = self.get_intersection(v2s, minus_two)
      x_minus_one = self.get_intersection(v2s, minus_one)
      x_central = self.get_intersection(v2s, central)
      x_plus_one = self.get_intersection(v2s, plus_one)
      x_plus_two = self.get_intersection(v2s, plus_two)
      print 'central {}'.format(x_central)

      #if float(mass) < 2.8 or float(mass) > 3.5:
      if x_plus_one == -99:
        print '\nWARNING - could not find crossing for +1sigma'
        crossings = np.linspace(1, 3, 50)
        for crossing in crossings:
          x_central_tmp = self.get_intersection(v2s, central, crossing)
          x_plus_one_tmp = self.get_intersection(v2s, plus_one, crossing)
          if x_plus_one_tmp != -99: 
            break
        x_diff_plus_one = x_plus_one_tmp - x_central_tmp
        x_plus_one = x_central + x_diff_plus_one

      if x_plus_two == -99:
        print '\nWARNING - could not find crossing for +2sigma'
        crossings = np.linspace(1, 3, 50)
        for crossing in crossings:
          x_central_tmp = self.get_intersection(v2s, central, crossing)
          x_plus_two_tmp = self.get_intersection(v2s, plus_two, crossing)
          if x_plus_two_tmp != -99: 
            break
        x_diff_plus_two = x_plus_two_tmp - x_central_tmp
        x_plus_two = x_central + x_diff_plus_two

      #if x_minus_one == -99:
      #  print '\nWARNING - could not find crossing for +1sigma'
      #  crossings = np.linspace(1, 2, 20)
      #  for crossing in crossings:
      #    x_central_tmp = self.get_intersection(v2s, central, crossing)
      #    x_minus_one_tmp = self.get_intersection(v2s, minus_one, crossing)
      #    if x_minus_one_tmp != -99: 
      #      break
      #  x_diff_minus_one = x_minus_one_tmp - x_central_tmp
      #  x_minus_one = x_central + x_diff_minus_one

      #if x_minus_two == -99:
      #  print '\nWARNING - could not find crossing for +2sigma'
      #  crossings = np.linspace(1, 2, 20)
      #  for crossing in crossings:
      #    x_central_tmp = self.get_intersection(v2s, central, crossing)
      #    x_minus_two_tmp = self.get_intersection(v2s, minus_two, crossing)
      #    if x_minus_two_tmp != -99: 
      #      break
      #  x_diff_minus_two = x_minus_two_tmp - x_central_tmp
      #  x_minus_two = x_central + x_diff_minus_two

      #if x_minus_two == -99:
      #  x_minus_two = 4e-3 #TODO remove, this is a temporary fix for mass 4.5

      if not self.do_blind:
        x_obs = self.get_intersection(v2s, obs)

      limits2D[mass]['exp_minus_two'] = x_minus_two
      limits2D[mass]['exp_minus_one'] = x_minus_one
      limits2D[mass]['exp_central'  ] = x_central  
      limits2D[mass]['exp_plus_one' ] = x_plus_one 
      limits2D[mass]['exp_plus_two' ] = x_plus_two 
      if not self.do_blind:
        limits2D[mass]['obs'] = x_obs 

      #print '({}, {}, {}): {}'.format(self.fe, self.fu, self.ft, x_central)

    print '\nwill plot 2D limits' 
    with open('{}/results.pck'.format(plotDir), 'w') as ff:
        pickle.dump(limits2D, ff)

    masses_obs       = []
    masses_central   = []
    #masses_one_sigma = []
    #masses_two_sigma = []
    masses_plus_one_sigma = []
    masses_minus_one_sigma = []
    masses_plus_two_sigma = []
    masses_minus_two_sigma = []

    boundary_plus_two = []
    boundary_plus_one = []
    boundary_minus_two = []
    boundary_minus_one = []

    minus_two = []
    minus_one = []
    central   = []
    plus_one  = []
    plus_two  = []
    obs       = []

    # go through the different mass points first left to right to catch the lower exclusion bound
    # then right to left to catch the upper exclusion bound
    for mass in sorted(limits2D.keys(), key=self.sortList):

        if self.do_coupling_scenario:
          exclusion_coupling_filename = '{}/exclusion_m_{}_{}_{}_{}.txt'.format(plotDir, str(mass).replace('.', 'p'), self.fe, self.fu, self.ft)
          exclusion_coupling_file = open(exclusion_coupling_filename, 'w+')
        
        #if len(limits2D[mass]['exp_central'])>0 and len(limits2D[mass]['exp_minus_one'])>0 and len(limits2D[mass]['exp_plus_one' ])>0 and len(limits2D[mass]['exp_minus_two'])>0 and len(limits2D[mass]['exp_plus_two' ])>0:
        #central.append(limits2D[mass]['exp_central'])
        #minus_one.append(limits2D[mass]['exp_minus_one'])
        #plus_one.append(limits2D[mass]['exp_plus_one' ])
        #minus_two.append(limits2D[mass]['exp_minus_two'])
        #plus_two.append(limits2D[mass]['exp_plus_two' ])

        #masses_central.append(float(mass))
        #masses_one_sigma.append(float(mass))
        #masses_two_sigma.append(float(mass))

        central_value = limits2D[mass]['exp_central']
        #central.append(limits2D[mass]['exp_central'])
        central.append(central_value)
        masses_central.append(float(mass))
        if limits2D[mass]['exp_plus_one' ] != -99. and limits2D[mass]['exp_plus_two' ] != -99.:
          plus_two_value = limits2D[mass]['exp_plus_two' ]
          plus_one_value = limits2D[mass]['exp_plus_one' ]

          minus_two_value = limits2D[mass]['exp_minus_two']
          minus_one_value = limits2D[mass]['exp_minus_one']
          #if minus_one_value == -99.:
          #  diff = plus_one_value - central_value
          #  minus_one_value = abs(central_value - diff)
          #if minus_two_value == -99.:
          #  diff = plus_two_value - central_value
          #  minus_two_value = abs(central_value - diff)

          if plus_two_value != -99.:
            plus_two.append(plus_two_value)
            boundary_plus_two.append(central_value)
            masses_plus_two_sigma.append(float(mass))

          if minus_two_value != -99.:
            minus_two.append(minus_two_value)
            boundary_minus_two.append(central_value)
            masses_minus_two_sigma.append(float(mass))

          if plus_one_value != -99.:
            plus_one.append(plus_one_value)
            boundary_plus_one.append(central_value)
            masses_plus_one_sigma.append(float(mass))

          if minus_one_value != -99.:
            minus_one.append(minus_one_value)
            boundary_minus_one.append(central_value)
            masses_minus_one_sigma.append(float(mass))

          #minus_one.append(minus_one_value)
          #plus_one.append(plus_one_value)

          #masses_one_sigma.append(float(mass))

        if not self.do_blind:
          obs.append(limits2D[mass]['obs'])
          masses_obs.append(float(mass))

        if self.do_coupling_scenario:
          exclusion_coupling_file.write('\n{} {} {} {}'.format(self.fe.replace('p', '.'), self.fu.replace('p', '.'), self.ft.replace('p', '.'), limits2D[mass]['exp_central']))
          exclusion_coupling_file.close()
          print '--> {} created'.format(exclusion_coupling_filename) 

       
    '''
    for mass in sorted(limits2D.keys(), key=self.sortList, reverse=True):

        if not self.do_blind:
          if len(limits2D[mass]['obs'])>1: 
              obs.append( max(limits2D[mass]['obs']) )
              masses_obs.append(mass)
        
        #if len(limits2D[mass]['exp_central'])>1 and len(limits2D[mass]['exp_minus_one'])>1 and len(limits2D[mass]['exp_plus_one' ])>1 and len(limits2D[mass]['exp_minus_two'])>1 and len(limits2D[mass]['exp_plus_two' ])>1: 
        central.append( max(limits2D[mass]['exp_central'  ]) )
        masses_central.append(float(mass))

        minus_one       .append( max(limits2D[mass]['exp_minus_one']) )
        plus_one        .append( max(limits2D[mass]['exp_plus_one' ]) )
        masses_one_sigma.append(float(mass))

        minus_two       .append( max(limits2D[mass]['exp_minus_two']) )
        plus_two        .append( max(limits2D[mass]['exp_plus_two' ]) )
        masses_two_sigma.append(float(mass))
    '''


    # plot the 2D limits
    print masses_central
    #print 'the_central = np.array({})'.format(central)
    #print 'the_minus_two = np.array({})'.format(minus_two)
    #print 'the_minus_one = np.array({})'.format(minus_one)
    #print 'the_plus_one = np.array({})'.format(plus_one)
    #print 'the_plus_two = np.array({})'.format(plus_two)
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=self.fe.replace('p', '.'), mu=r'\mu', fu=self.fu.replace('p', '.'), tau=r'\tau', ft=self.ft.replace('p', '.'))
    ax.text(0.1, 0.93, 'CMS', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30, fontweight='bold')
    ax.text(0.17, 0.84, 'Preliminary', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=25, fontstyle='italic')
    ax.text(0.26, 0.75, coupling_scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.25, 0.66, 'Lepton universality tests', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='blue', fontsize=18)
    ax.text(0.84, 0.93, self.scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='black', fontsize=23, fontweight='bold')
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--')
    #f1 = plt.fill_between(masses_two_sigma, minus_two, plus_two, color='gold'       , label=r'95% expected')
    #f2 = plt.fill_between(masses_one_sigma, minus_one, plus_one, color='forestgreen', label=r'68% expected')
    f1 = plt.fill_between(masses_minus_two_sigma, minus_two, boundary_minus_two, color='gold'     , label=r'95% expected')
    f2 = plt.fill_between(masses_minus_one_sigma, minus_one, boundary_minus_one, color='forestgreen', label=r'68% expected')
    f3 = plt.fill_between(masses_plus_two_sigma, boundary_plus_two, plus_two, color='gold'       , label=r'95% expected')
    f4 = plt.fill_between(masses_plus_one_sigma, boundary_plus_one, plus_one, color='forestgreen', label=r'68% expected')
    p1, = plt.plot(masses_central, central, color='red', label='Median expected', linewidth=2)

    #p2, = plt.plot(db.masses_delphidisplaced, db.exp_delphidisplaced, color='darkorange', label='Delphi displaced', linewidth=1.3, linestyle='dashed')
    p2, = plt.plot(db.masses_atlas_lower, db.exp_atlas_lower, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed')
    p2_2, = plt.plot(db.masses_atlas_upper, db.exp_atlas_upper, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed')
    p3, = plt.plot(db.masses_cmsdisplacedmuon, db.exp_cmsdisplacedmuon, color='blueviolet', label='CMS displaced', linewidth=1.3, linestyle='dashed')
    #p4, = plt.plot(db.masses_lhcb, db.exp_lhcb, color='darkred', label='LHCb', linewidth=1.3, linestyle='dashed')
    p4, = plt.plot(db.masses_lhcb_peskin, db.exp_lhcb_peskin, color='darkred', label='LHCb', linewidth=1.3, linestyle='dashed')
    p5, = plt.plot(db.masses_belle, db.exp_belle, color='magenta', label='Belle', linewidth=1.3, linestyle='dashed')

    #p3, = plt.plot(db.masses_delphiprompt, db.exp_delphiprompt, color='blueviolet', label='Delphi prompt', linewidth=1.3, linestyle='dashed')
    #if 'mmm' in self.channels or 'mem' in self.channels:
    #p4, = plt.plot(db.masses_atlasdisplacedmuonLNV, db.exp_atlasdisplacedmuonLNV, color='firebrick', label='Atlas displaced muon LNV', linewidth=1.3, linestyle='dashed')
    #p5, = plt.plot(db.masses_atlasdisplacedmuonLNC, db.exp_atlasdisplacedmuonLNC, color='darkorange', label='Atlas displaced muon LNC', linewidth=1.3, linestyle='dashed')
    #p7, = plt.plot(db.masses_cmspromptmuon, db.exp_cmspromptmuon, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')
    #else: 
    #  p7, = plt.plot(db.masses_cmspromptelectron, db.exp_cmspromptelectron, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')

    if not self.do_blind:
      p8, = plt.plot(masses_obs, obs, color='black', label='Observed', linewidth=2)

    #if not self.do_blind:
    #  first_legend = plt.legend(handles=[p1, p8, f1, f2], loc='lower right', fontsize=18)
    #else:
    #  first_legend = plt.legend(handles=[p1, f2, f1], loc='lower right', fontsize=18)
    #ax = plt.gca().add_artist(first_legend)

    second_legend = plt.legend(handles=[p2, p3, p4, p5], loc='lower left', fontsize=18)
    #ax = plt.gca().add_artist(second_legend)
    #if 'mmm' in self.channels or 'mem' in self.channels:
    #second_legend = plt.legend(handles=[p2, p3, p4, p5, p7], loc='lower left')
    #else: 
    #  second_legend = plt.legend(handles=[p2, p3, p7], loc='upper left')

    plt.title(lumi + ' (13 TeV)', loc='right', fontsize=23)
    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(3e-7, 1e-0)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    #plt.xlim(0, max(masses_central))
    plt.xlim(min(masses_central), max(masses_central))
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = '2d_hnl_limit_{}'.format(self.scenario) 
    else:
      name_2d = '2d_hnl_limit_scenario_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


    ## smoothing the structures
    #plt.clf()
    #smooth_index = 4
    #f, ax = plt.subplots(figsize=(9, 8))
    #if not self.do_coupling_scenario:
    #  self.fe = '0.0'
    #  self.fu = '1.0'
    #  self.ft = '0.0'
    #coupling_scenario = r'(f$_{e}$={fe}, f$_{mu}$={fu}, f$_{tau}$={ft})'.format(e='e', fe=self.fe.replace('p', '.'), mu=r'\mu', fu=self.fu.replace('p', '.'), tau=r'\tau', ft=self.ft.replace('p', '.'))
    #ax.text(0.1, 0.93, 'CMS', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30, fontweight='bold')
    #ax.text(0.17, 0.84, 'Preliminary', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=25, fontstyle='italic')
    #ax.text(0.26, 0.73, coupling_scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
    #ax.text(0.25, 0.63, 'Lepton universality tests', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='blue', fontsize=18)
    #ax.text(0.84, 0.93, self.scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='black', fontsize=23, fontweight='bold')
    #plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--')
    #f1 = plt.fill_between(masses_two_sigma, self.smooth(minus_two, smooth_index), self.smooth(plus_two, smooth_index), color='gold', label=r'95% expected')
    #f2 = plt.fill_between(masses_one_sigma, self.smooth(minus_one, smooth_index), self.smooth(plus_one, smooth_index), color='forestgreen', label=r'68% expected')
    ##f1 = plt.fill_between(masses_two_sigma, self.smooth(minus_two, smooth_index), self.smooth(plus_two, smooth_index-1), color='gold', label=r'95% expected')
    ##f2 = plt.fill_between(masses_one_sigma, self.smooth(minus_one, smooth_index), self.smooth(plus_one, smooth_index-1), color='forestgreen', label=r'68% expected')
    #p1, = plt.plot(masses_central, self.smooth(central, smooth_index), color='red', label='Median expected', linewidth=2)
    ##p2, = plt.plot(db.masses_delphidisplaced, db.exp_delphidisplaced, color='black', label='Delphi displaced', linewidth=1.3, linestyle='dashed')
    ##p3, = plt.plot(db.masses_delphiprompt, db.exp_delphiprompt, color='blueviolet', label='Delphi prompt', linewidth=1.3, linestyle='dashed')
    ##if 'mmm' in self.channels or 'mem' in self.channels:
    ##  p4, = plt.plot(db.masses_atlasdisplacedmuonLNV, db.exp_atlasdisplacedmuonLNV, color='firebrick', label='Atlas displaced muon LNV', linewidth=1.3, linestyle='dashed')
    ##  p5, = plt.plot(db.masses_atlasdisplacedmuonLNC, db.exp_atlasdisplacedmuonLNC, color='darkorange', label='Atlas displaced muon LNC', linewidth=1.3, linestyle='dashed')
    ##  p7, = plt.plot(db.masses_cmspromptmuon, db.exp_cmspromptmuon, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')
    ##else: 
    ##  p7, = plt.plot(db.masses_cmspromptelectron, db.exp_cmspromptelectron, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')

    #if not self.do_blind:
    #  p8, = plt.plot(masses_obs, obs, color='black', label='observed', linewidth=2)

    #if not self.do_blind:
    #  first_legend = plt.legend(handles=[p1, p8, f1, f2], loc='lower right', fontsize=20)
    #else:
    #  first_legend = plt.legend(handles=[p1, f2, f1], loc='lower right', fontsize=20)
    ##ax = plt.gca().add_artist(first_legend)
    ##if 'mmm' in self.channels or 'mem' in self.channels:
    ##  second_legend = plt.legend(handles=[p2, p3, p4, p5, p7], loc='lower left')
    ##else: 
    ##  second_legend = plt.legend(handles=[p2, p3, p7], loc='upper left')

    #plt.title(lumi + ' (13 TeV)', loc='right', fontsize=23)
    #plt.ylabel(r'$|V|^2$', fontsize=23)
    #plt.yticks(fontsize=17)
    ##plt.ylim(1e-10, 1e-0)
    #plt.ylim(1e-5, 1e-0)
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    ##plt.xlim(0, max(masses_central))
    #plt.xlim(min(masses_central), max(masses_central))
    #plt.xticks(fontsize=17)
    #plt.yscale('log')
    #plt.xscale('linear')
    #plt.grid(True)
    #if not self.do_coupling_scenario:
    #  name_2d = '2d_hnl_limit_{}_smoothed'.format(self.scenario) 
    #else:
    #  name_2d = '2d_hnl_limit_scenario_{}_{}_{}_{}_smoothed'.format(self.scenario, self.fe, self.fu, self.ft) 
    #plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    #plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    #print '--> {}/{}.png created'.format(plotDir, name_2d)


    #print '\n-> Plots saved in {}'.format(plotDir)


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


