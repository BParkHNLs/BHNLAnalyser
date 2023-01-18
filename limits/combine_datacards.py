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
  parser = ArgumentParser(description='Script to combine the datacards among displacement bins, years, and flavour channels', add_help=True)
  parser.add_argument('--homedir'           , type=str, dest='homedir'           , help='name of the homedir'                          , default=None)
  parser.add_argument('--outdirlabel'       , type=str, dest='outdirlabel'       , help='name of the outdir'                           , default=None)
  parser.add_argument('--subdirlabel'       , type=str, dest='subdirlabel'       , help='name of the subdir'                           , default=None)
  parser.add_argument('--mass'              , type=str, dest='mass'              , help='mass'                                         , default='1.0')
  parser.add_argument('--signal_type'       , type=str, dest='signal_type'       , help='signal under consideration'                   , default='majorana', choices=['majorana', 'dirac'])
  parser.add_argument('--categories_label'  , type=str, dest='categories_label'  , help='label of the list of categories'              , default='standard')
  parser.add_argument('--mass_whitelist'    , type=str, dest='mass_whitelist'    , help='allowed values for masses'                    , default=None)
  parser.add_argument('--mass_blacklist'    , type=str, dest='mass_blacklist'    , help='values for masses to skip'                    , default=None)
  parser.add_argument('--coupling_whitelist', type=str, dest='coupling_whitelist', help='allowed values for couplings'                 , default=None)
  parser.add_argument('--coupling_blacklist', type=str, dest='coupling_blacklist', help='values for couplings to skip'                 , default=None)
  parser.add_argument('--wildcard'          , type=str, dest='wildcard'          , help='datacard generic string'                      , default='datacard_bhnl*.txt')
  parser.add_argument('--run_blind'                   , dest='run_blind'         , help='run blinded or unblinded', action='store_true', default=False)
  return parser.parse_args()


# getting the parsed info
opt = getOptions()

homedir = opt.homedir
outdirlabel = opt.outdirlabel
subdirlabel = opt.subdirlabel
signal_type = opt.signal_type #TODO for the moment, does not do anything 
the_mass = opt.mass
categories = categories[opt.categories_label]
path_to_datacards = subdirlabel
datacard_wildcard = opt.wildcard

# create directories
outputdir = '{}/outputs/{}/datacards_combined/{}'.format(homedir, outdirlabel, subdirlabel)
if not path.exists(outputdir):
  os.system('mkdir -p {}'.format(outputdir))    

print 'loading cards...'
all_datacards = []
all_datacards = glob('/'.join([path_to_datacards, datacard_wildcard]))
print '... datacards loaded'

categories_to_combine = []
for category in categories:
  if 'incl' in category.label: continue
  categories_to_combine.append(category.label)

# nested dictionary with mass and coupling as keys
digested_datacards = OrderedDict()

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
    if opt.mass_whitelist != None:
      if str(signal_mass) not in opt.mass_whitelist.split(','): continue
    
    if opt.mass_blacklist != None:
      if str(signal_mass) in opt.mass_blacklist.split(','): continue
    
    if opt.coupling_whitelist != None:
      if str(signal_coupling) not in opt.coupling_whitelist.split(','): continue 
    
    if opt.coupling_blacklist != None:
      if str(signal_coupling) in opt.coupling_blacklist.split(','): continue 

    # will fetch the datacards
    if signal_mass not in digested_datacards.keys():
        digested_datacards[signal_mass] = OrderedDict()
    
    if signal_coupling not in digested_datacards[signal_mass].keys():
        digested_datacards[signal_mass][signal_coupling] = []

    for cat in categories_to_combine:
      #print 'cat ',cat
      if cat+'.txt' in idc_ref:
        #print '{} {} {}'.format(signal_mass, signal_coupling, idc_ref)
        digested_datacards[signal_mass][signal_coupling].append(idc_ref) 
    
    
for mass, couplings in digested_datacards.iteritems():
    if mass != the_mass: continue
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
                                                              
        # combine the cards    
        command = 'combineCards.py'
        for idc in datacards_to_combine:
          #print idc
          #print '{} \t {}'.format(cat, idc)

          # fetch ctau to add to combined datacard name
          ctau = idc[idc.find('ctau_')+5:idc.find('_', idc.find('ctau_')+5)]

          if any([v in idc for v in categories_to_combine]):
            command += ' {}'.format(idc)

        command += (' > {o}/datacard_combined_m_{m}_ctau_{ctau}_v2_{v2}.txt'.format(o=outputdir, m=str(mass), ctau=ctau, v2=coupling)) 

        #print command
        os.system(command)
        
        print ('\t\t -> combined datacards between the categories in {o}/datacard_combined_m_{m}_ctau_{ctau}_v2_{v2}.txt'.format(o=outputdir, m=str(mass), ctau=ctau, v2=coupling))
