import sys
import os
from os import path
import glob
import re
import pickle
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
  parser.add_argument('--homedir'           , type=str, dest='homedir'           , help='name of the homedir'                           , default=None)
  parser.add_argument('--outdirlabel'       , type=str, dest='outdirlabel'       , help='name of the outdir'                           , default=None)
  parser.add_argument('--subdirlabel'       , type=str, dest='subdirlabel'       , help='name of the subdir'                           , default=None)
  parser.add_argument('--signal_type'       , type=str, dest='signal_type'       , help='signal under consideration'                   , default='majorana', choices=['majorana', 'dirac'])
  parser.add_argument('--mass_whitelist'    , type=str, dest='mass_whitelist'    , help='allowed values for masses'                    , default='None')
  parser.add_argument('--mass_blacklist'    , type=str, dest='mass_blacklist'    , help='values for masses to skip'                    , default='None')
  parser.add_argument('--coupling_whitelist', type=str, dest='coupling_whitelist', help='allowed values for couplings'                 , default='None')
  parser.add_argument('--coupling_blacklist', type=str, dest='coupling_blacklist', help='values for couplings to skip'                 , default='None')
  parser.add_argument('--fe'                , type=str, dest='fe'                , help='electron coupling fraction'                   , default='1.0')
  parser.add_argument('--fu'                , type=str, dest='fu'                , help='muon coupling fraction'                       , default='1.0')
  parser.add_argument('--ft'                , type=str, dest='ft'                , help='tau coupling fraction'                        , default='1.0')
  parser.add_argument('--run_blind'                   , dest='run_blind'         , help='run blinded or unblinded', action='store_true', default=False)
  return parser.parse_args()


class LimitPlotter(object):
  def __init__(self, signal_type, homedir, outdirlabel, subdirlabel, mass_whitelist, mass_blacklist, coupling_whitelist, coupling_blacklist, run_blind, fe, fu, ft):
    self.signal_type = signal_type
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.mass_whitelist = mass_whitelist
    self.mass_blacklist = mass_blacklist
    self.coupling_whitelist = coupling_whitelist
    self.coupling_blacklist = coupling_blacklist
    self.run_blind = run_blind
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


  def process(self):
    signal_type = self.signal_type 
    flavour = r'$|V|^2_{\mu}$'
    lumi =  '41.6 fb'+r'$^{-1}$'

    # get the files 
    if not self.do_coupling_scenario:
      pathToResults = './results/'
    else:
      pathToResults = './results_{}_{}_{}/'.format(self.fe, self.fu, self.ft)

    fileName = 'result*.txt'

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
     
        # for each mass, get the list of the couplings from the file name
        coupling = limitFile[limitFile.find('v2_')+3:limitFile.find('.txt')]
        val_coupling = float(coupling)
      
        # get white/black listed coupling
        if self.coupling_whitelist!='None':
          if str(val_coupling) not in self.coupling_whitelist.split(','): continue 
        
        if self.coupling_blacklist!='None':
          if str(val_coupling) in self.coupling_blacklist.split(','): continue 
        
        try:
          thefile = open('{}/result_m_{}_v2_{}.txt'.format(pathToResults, mass, coupling), 'r')
          
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
          if not self.run_blind:
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
          if not self.run_blind: obs.append(float(val_obs))

        except:
          print 'Cannot open {}result_m_{}_v2_{}.txt'.format(pathToResults, mass, coupling)


      print '-> will plot 1D limit for mass {}'.format(mass)
      
      if not self.run_blind:
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
      if not self.run_blind:
          obs = [jj[6] for jj in graph]
      
      print '\n-> Plots will be saved in {}'.format(plotDir)
          
      plt.clf()
      ax = plt.axes()
      print '   couplings: {}'.format(v2s)
      coupling_scenario = r'(f$_{e}$={fe}, f$_{mu}$={fu}, f$_{tau}$={ft})'.format(e='e', fe=self.fe.replace('p', '.'), mu=r'\mu', fu=self.fu.replace('p', '.'), tau=r'\tau', ft=self.ft.replace('p', '.'))
      plt.fill_between(v2s, minus_two, plus_two, color='gold', label=r'$\pm 2 \sigma$')
      plt.fill_between(v2s, minus_one, plus_one, color='forestgreen' , label=r'$\pm 1 \sigma$')
      plt.plot(v2s, central, color='red', label='central expected', linewidth=2)
      if not self.run_blind:
          plt.plot(v2s, obs, color='black', label='observed')    
      
      plt.axhline(y=1, color='black', linestyle='-')
      #plt.xlabel(flavour)
      plt.xlabel(r'$|V|^2$')
      plt.ylabel('exclusion limit 95% CL')

      #plt.title('HNL m = %s GeV %s' %(mass, signal_type))
      if self.do_coupling_scenario:
        plt.title(r'Majorana HN, m = %s GeV, %s' %(mass, coupling_scenario))
      else:
        plt.title(r'Majorana HN, m = %s GeV' %(mass,))
      plt.legend()
      plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
      #ax.text(0.5, 0.5, coupling_scenario, style='normal', bbox={'facecolor': 'white', 'alpha': 1., 'pad': 10}, transform=ax.transAxes)
      #plt.xlim(1e-5, 1e-3)
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
        name_log = 'limit_m_{}_log'.format(mass.replace('.', 'p')) 
      else:
        name_log = 'limit_m_{}_scenario_{}_{}_{}_log'.format(mass.replace('.', 'p'), self.fe, self.fu, self.ft) 
      plt.savefig('{}/{}.pdf'.format(plotDir, name_log))
      plt.savefig('{}/{}.png'.format(plotDir, name_log))
      

      # save the crossing for 2D limits
      limits2D[mass] = OrderedDict()

      # find the intersections    
      x_minus_two, y = intersection(np.array(v2s), np.array(minus_two), np.array(v2s), np.ones(len(v2s)))    
      x_minus_one, y = intersection(np.array(v2s), np.array(minus_one), np.array(v2s), np.ones(len(v2s)))    
      x_central  , y = intersection(np.array(v2s), np.array(central)  , np.array(v2s), np.ones(len(v2s)))    
      x_plus_one , y = intersection(np.array(v2s), np.array(plus_one) , np.array(v2s), np.ones(len(v2s)))    
      x_plus_two , y = intersection(np.array(v2s), np.array(plus_two) , np.array(v2s), np.ones(len(v2s)))    
      if not self.run_blind:
          x_obs , y = intersection(np.array(v2s), np.array(obs) , np.array(v2s), np.ones(len(v2s)))    

      limits2D[mass]['exp_minus_two'] = x_minus_two
      limits2D[mass]['exp_minus_one'] = x_minus_one
      limits2D[mass]['exp_central'  ] = x_central  
      limits2D[mass]['exp_plus_one' ] = x_plus_one 
      limits2D[mass]['exp_plus_two' ] = x_plus_two 
      if not self.run_blind:
          limits2D[mass]['obs'] = x_obs 

      #print '({}, {}, {}): {}'.format(self.fe, self.fu, self.ft, x_central)

    print '\nwill plot 2D limits' 
    with open('{}/results.pck'.format(plotDir), 'w') as ff:
        pickle.dump(limits2D, ff)

    masses_obs       = []
    masses_central   = []
    masses_one_sigma = []
    masses_two_sigma = []

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
          exclusion_coupling_filename = '{}/exclusion_m_{}.txt'.format(plotDir, str(mass).replace('.', 'p'))
          if not path.exists(exclusion_coupling_filename):
            exclusion_coupling_file = open(exclusion_coupling_filename, 'w+')
          else:
            exclusion_coupling_file = open(exclusion_coupling_filename, 'a+')
        
        if not self.run_blind:
          if len(limits2D[mass]['obs'])>0: 
              obs.append( min(limits2D[mass]['obs']) )
              masses_obs.append(mass)
        
        #if len(limits2D[mass]['exp_central'])>0 and len(limits2D[mass]['exp_minus_one'])>0 and len(limits2D[mass]['exp_plus_one' ])>0 and len(limits2D[mass]['exp_minus_two'])>0 and len(limits2D[mass]['exp_plus_two' ])>0:
        central.append( min(limits2D[mass]['exp_central']) )
        masses_central.append(float(mass))
        #NOTE make sure that same coupling line is overwritten
        if self.do_coupling_scenario:
          #print '{} {} {} {}'.format(self.fe.replace('p', '.'), self.fu.replace('p', '.'), self.ft.replace('p', '.'), min(limits2D[mass]['exp_central']))
          exclusion_coupling_file.write('\n{} {} {} {}'.format(self.fe.replace('p', '.'), self.fu.replace('p', '.'), self.ft.replace('p', '.'), min(limits2D[mass]['exp_central']))) 

        minus_one.append( min(limits2D[mass]['exp_minus_one']) )
        plus_one.append( min(limits2D[mass]['exp_plus_one' ]) )
        masses_one_sigma.append(float(mass))

        minus_two.append( min(limits2D[mass]['exp_minus_two']) )
        plus_two.append( min(limits2D[mass]['exp_plus_two' ]) )
        masses_two_sigma.append(float(mass))

        if self.do_coupling_scenario:
          exclusion_coupling_file.close()

       
    '''
    for mass in sorted(limits2D.keys(), key=self.sortList, reverse=True):

        if not self.run_blind:
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
    plt.clf()
    ax = plt.axes()
    f1 = plt.fill_between(masses_two_sigma, minus_two, plus_two, color='gold'       , label=r'$\pm 2 \sigma$')
    f2 = plt.fill_between(masses_one_sigma, minus_one, plus_one, color='forestgreen', label=r'$\pm 1 \sigma$')
    p1, = plt.plot(masses_central, central, color='red', label='central expected', linewidth=2)
    #p2, = plt.plot(db.masses_delphidisplaced, db.exp_delphidisplaced, color='black', label='Delphi displaced', linewidth=1.3, linestyle='dashed')
    #p3, = plt.plot(db.masses_delphiprompt, db.exp_delphiprompt, color='blueviolet', label='Delphi prompt', linewidth=1.3, linestyle='dashed')
    #if 'mmm' in self.channels or 'mem' in self.channels:
    #  p4, = plt.plot(db.masses_atlasdisplacedmuonLNV, db.exp_atlasdisplacedmuonLNV, color='firebrick', label='Atlas displaced muon LNV', linewidth=1.3, linestyle='dashed')
    #  p5, = plt.plot(db.masses_atlasdisplacedmuonLNC, db.exp_atlasdisplacedmuonLNC, color='darkorange', label='Atlas displaced muon LNC', linewidth=1.3, linestyle='dashed')
    #  p7, = plt.plot(db.masses_cmspromptmuon, db.exp_cmspromptmuon, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')
    #else: 
    #  p7, = plt.plot(db.masses_cmspromptelectron, db.exp_cmspromptelectron, color='blue', label='CMS prompt muon', linewidth=1.3, linestyle='dashed')

    if not self.run_blind:
      p8, = plt.plot(masses_obs, obs, color='black', label='observed', linewidth=2)

    if not self.run_blind:
      first_legend = plt.legend(handles=[p1, p8, f1, f2], loc='upper right')
    else:
      first_legend = plt.legend(handles=[p1, f1, f2], loc='upper right')
    ax = plt.gca().add_artist(first_legend)
    #if 'mmm' in self.channels or 'mem' in self.channels:
    #  second_legend = plt.legend(handles=[p2, p3, p4, p5, p7], loc='lower left')
    #else: 
    #  second_legend = plt.legend(handles=[p2, p3, p7], loc='upper left')

    plt.title('CMS Preliminary', loc='left')
    plt.title(lumi + ' (13 TeV)', loc='right')
    plt.ylabel(flavour)
    #plt.ylim(1e-10, 1e-0)
    plt.ylim(1e-6, 1e-0)
    plt.xlabel('mass (GeV)')
    #plt.xlim(0, max(masses_central))
    plt.xlim(min(masses_central), max(masses_central))
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = '2d_hnl_limit'.format(mass.replace('.', 'p')) 
    else:
      name_2d = '2d_hnl_limit_scenario_{}_{}_{}'.format(mass.replace('.', 'p'), self.fe, self.fu, self.ft) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))


    #print '\n-> Plots saved in {}'.format(plotDir)


if __name__ == "__main__":

  opt=getOptions()

  plotter = LimitPlotter(
      signal_type = opt.signal_type,
      homedir = opt.homedir,
      outdirlabel = opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      mass_whitelist = opt.mass_whitelist,
      mass_blacklist = opt.mass_blacklist,
      coupling_whitelist = opt.coupling_whitelist,
      coupling_blacklist = opt.coupling_blacklist,
      run_blind = True if opt.run_blind else False,
      fe = opt.fe,
      fu = opt.fu,
      ft = opt.ft,
      )

  plotter.process()


