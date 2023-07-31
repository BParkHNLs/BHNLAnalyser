# https://sukhbinder.wordpress.com/2017/06/13/intersection-of-two-curves-in-pure-numpy/
# handle python subprocess, it might need to be updated when switching to python3
# https://stackoverflow.com/questions/4760215/running-shell-command-and-capturing-the-output

import os
from os import path
import re
from glob import glob
from itertools import product
from collections import OrderedDict
from decimal import Decimal

import sys
sys.path.append('../objects')
from categories import categories

'''

Script to combine the datacards between the different categories

'''

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to combine the datacards among categories', add_help=True)
  parser.add_argument('--homedir'           , type=str, dest='homedir'           , help='name of the homedir'                          , default=None)
  parser.add_argument('--outdirlabel'       , type=str, dest='outdirlabel'       , help='name of the outdir'                           , default=None)
  parser.add_argument('--subdirlabel'       , type=str, dest='subdirlabel'       , help='name of the subdir'                           , default=None)
  parser.add_argument('--mass'              , type=str, dest='mass'              , help='mass'                                         , default='1.0')
  #parser.add_argument('--signal_type'       , type=str, dest='signal_type'       , help='signal under consideration'                   , default='majorana', choices=['majorana', 'dirac'])
  parser.add_argument('--scenario'          , type=str, dest='scenario'          , help='signal under consideration'                   , default='Majorana', choices=['Majorana', 'Dirac'])
  parser.add_argument('--categories_label'  , type=str, dest='categories_label'  , help='label of the list of categories'              , default='standard')
  parser.add_argument('--mass_whitelist'    , type=str, dest='mass_whitelist'    , help='allowed values for masses'                    , default=None)
  parser.add_argument('--mass_blacklist'    , type=str, dest='mass_blacklist'    , help='values for masses to skip'                    , default=None)
  parser.add_argument('--coupling_whitelist', type=str, dest='coupling_whitelist', help='allowed values for couplings'                 , default=None)
  parser.add_argument('--coupling_blacklist', type=str, dest='coupling_blacklist', help='values for couplings to skip'                 , default=None)
  parser.add_argument('--wildcard'          , type=str, dest='wildcard'          , help='datacard generic string'                      , default='datacard_bhnl*.txt')
  parser.add_argument('--do_blind'                    , dest='do_blind'          , help='run blinded or unblinded', action='store_true', default=False)
  return parser.parse_args()


class DatacardCombiner(object):
  def __init__(self, opt):
    self.homedir = opt.homedir
    self.outdirlabel = opt.outdirlabel
    self.subdirlabel = opt.subdirlabel
    #signal_type = opt.signal_type #TODO for the moment, does not do anything 
    self.scenario = opt.scenario
    self.the_mass = opt.mass
    self.categories = categories[opt.categories_label]
    self.path_to_datacards = './'
    self.datacard_wildcard = opt.wildcard
    self.mass_whitelist = opt.mass_whitelist
    self.mass_blacklist = opt.mass_blacklist
    self.coupling_whitelist = opt.coupling_whitelist
    self.coupling_blacklist = opt.coupling_blacklist
    

  def getRateList(self, datacard_name):
    '''
      Get the list of the signal and background rates from the initial datacard
    '''
    f = open(datacard_name)
    lines = f.readlines()
    for line in lines:
      if 'rate ' not in line: continue
      rate_line = line
      break
    idx = 0
    rate_list = []
    while idx < len(rate_line) and idx != -1:
      if rate_line[idx].isdigit():
        idx1 = idx
        idx2 = rate_line.find(' ', idx+1)
        rate_list.append(float(rate_line[idx1:idx2]))
        idx = idx2
      else:
        idx = idx + 1

    return rate_list


  def getAlpha(self, line):
    '''
      Get the alpha of the statistical uncertainty from line
    '''
    # search where to start fetching for alpha (after the number of events)
    idx0 = line.find('gmN')+3
    idx = idx0
    while idx < len(line) and idx != -1:
      if line[idx].isdigit():
        idx1 = idx
        idx2 = line.find(' ', idx1+1)
        idx_begin = idx2
        break
      else:
        idx = idx + 1

    for idx in range(idx_begin, len(line)):
      if line[idx].isdigit():
        idx1 = idx
        idx2 = line.find(' ', idx1+1)
        alpha = line[idx1:idx2]
        break
      else:
        idx = idx + 1

    return alpha


  def updateRateList(self, rate_list):
    '''
      Update the signal rates according to given coupling scenario
    '''
    updated_rate_list = []
    for rate in rate_list:
      if rate != 1.:
        updated_rate = rate * 2.0
      else:
        # do not modify background rate, keep it to 1
        updated_rate = 1.

      updated_rate_list.append(updated_rate)

    return updated_rate_list


  def updateAlpha(self, alpha):
    '''
      Update the value of alpha
    '''
    updated_alpha = float(alpha) * 2.0
    updated_alpha = str(updated_alpha)

    return updated_alpha


  def updateDatacard(self, datacard_name, updated_rate_list):
    '''
      Update the datacard with the updated signal rates
    '''
    updated_datacard_name = datacard_name + '_tmp'
    updated_datacard = open(updated_datacard_name, 'w+')
    datacard = open(datacard_name)
    lines = datacard.readlines()
    for line in lines:
      if 'rate ' in line:
        rate_line = 'rate                                                           '            
        for updated_rate in updated_rate_list:
          rate_line += '{}  '.format(updated_rate)    
        updated_datacard.write(rate_line + '\n')

      elif 'gmN' in line:
        alpha = self.getAlpha(line)
        updated_alpha = self.updateAlpha(alpha)
        updated_line = line.replace(alpha, updated_alpha)
        updated_datacard.write(updated_line + '\n')

      else:
        updated_datacard.write(line)

    updated_datacard.close()

    print 'mv {} {}'.format(updated_datacard_name, datacard_name)
    os.system('mv {} {}'.format(updated_datacard_name, datacard_name))


  def getCategoryLabel(self, datacard):
    return datacard[datacard.find('cat_')+4:datacard.find('.txt')]


  def process(self):
    # create directories
    outputdir = '{}/outputs/{}/datacards_combined/{}'.format(self.homedir, self.outdirlabel, self.subdirlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))    

    print 'loading cards...'
    all_datacards = []
    all_datacards = glob('/'.join([self.path_to_datacards, self.datacard_wildcard]))
    print '... datacards loaded'

    categories_to_combine = []
    for category in self.categories:
      if 'incl' in category.label: continue
      categories_to_combine.append(category.label)

    # nested dictionary with mass and coupling as keys
    digested_datacards = OrderedDict()
    channel_labels = OrderedDict()

    # store results for 2D limits
    limits2D = OrderedDict()

    the_set_datacards = all_datacards

    #for card in all_datacards:
    #  print card

    for idc_ref in the_set_datacards:
        name = idc_ref.split('/')[-1]
        signal_mass = name[name.find('m_')+2:name.find('_', name.find('m_')+2)]
        signal_coupling = name[name.find('v2_')+3:name.find('_cat')]
       
        # get white/black listed mass/couplings
        if self.mass_whitelist != None:
          if str(signal_mass) not in self.mass_whitelist.split(','): continue
        
        if self.mass_blacklist != None:
          if str(signal_mass) in self.mass_blacklist.split(','): continue
        
        if self.coupling_whitelist != None:
          if str(signal_coupling) not in self.coupling_whitelist.split(','): continue 
        
        if self.coupling_blacklist != None:
          if str(signal_coupling) in self.coupling_blacklist.split(','): continue 

        # will fetch the datacards
        if signal_mass not in digested_datacards.keys():
            digested_datacards[signal_mass] = OrderedDict()
            channel_labels[signal_mass] = OrderedDict()
        
        if signal_coupling not in digested_datacards[signal_mass].keys():
            digested_datacards[signal_mass][signal_coupling] = []
            channel_labels[signal_mass][signal_coupling] = []

        for cat in categories_to_combine:
          #print 'cat ',cat
          if self.scenario == 'Dirac' and '_SS' in cat: continue
          if cat+'.txt' in idc_ref:
            #print '{} {} {}'.format(signal_mass, signal_coupling, idc_ref)
            digested_datacards[signal_mass][signal_coupling].append(idc_ref) 
            channel_labels[signal_mass][signal_coupling].append(self.getCategoryLabel(idc_ref)) 
        
        
    for mass, couplings in digested_datacards.iteritems():
        if mass != self.the_mass: continue
        print 'mass =', mass
        
        v2s       = []
        obs       = []
        minus_two = []
        minus_one = []
        central   = []
        plus_one  = []
        plus_two  = []

        for coupling in couplings.keys():
            print '\tcoupling =', coupling

            datacards_to_combine = digested_datacards[mass][coupling]
            labels = channel_labels[mass][coupling]
                                                                  
            # combine the cards    
            command = 'combineCards.py'
            for idc, datacard in enumerate(datacards_to_combine):
              #print datacard
              #print '{} \t {}'.format(cat, datacard)

              # fetch ctau to add to combined datacard name
              ctau = datacard[datacard.find('ctau_')+5:datacard.find('_', datacard.find('ctau_')+5)]

              if any([v in datacard for v in categories_to_combine]):
                command += ' {}={}'.format(labels[idc], datacard)

            combined_datacard = '{o}/datacard_combined_{sc}_m_{m}_ctau_{ctau}_v2_{v2}.txt'.format(o=outputdir, m=str(mass), ctau=ctau, v2=coupling, sc=self.scenario)
            command += (' > {}'.format(combined_datacard)) 

            #print command
            os.system(command)
            
            #print ('\t\t -> combined datacards between the categories in {o}/datacard_combined_m_{m}_ctau_{ctau}_v2_{v2}_{sc}.txt'.format(o=outputdir, m=str(mass), ctau=ctau, v2=coupling, sc=self.scenario))
            print ('\t\t -> combined datacards between the categories in {}'.format(combined_datacard))

            # in the dirac scenario, correct the signal rate
            if self.scenario == 'Dirac':
              rate_list = self.getRateList(datacard_name=combined_datacard)
              updated_rate_list = self.updateRateList(rate_list=rate_list)
              self.updateDatacard(datacard_name=combined_datacard, updated_rate_list=updated_rate_list)



if __name__ == "__main__":
  # getting the parsed info
  opt = getOptions()

  DatacardCombiner(opt=opt).process()



