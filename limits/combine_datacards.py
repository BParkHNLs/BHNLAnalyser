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

'''

Script to combine the datacards among displacement bins, years and channels

'''

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to combine the datacards among displacement bins, years, and flavour channels', add_help=True)
  parser.add_argument('--version', type=str, dest='version', help='version label', default='L1')
  parser.add_argument('--signal', type=str, dest='signal', help='signal under consideration', default='majorana', choices=['majorana', 'dirac'])
  parser.add_argument('--run_blind', dest='run_blind', help='run blinded or unblinded', action='store_true', default=False)
  parser.add_argument('--years', type=str, dest='years', help='years to combine', default='2016,2017,2018')
  parser.add_argument('--channels', type=str, dest='channels', help='channels to combine', default='mmm,mem_os,mem_ss')
  parser.add_argument('--mass_whitelist', type=str, dest='mass_whitelist', help='allowed values for masses', default=None)
  parser.add_argument('--mass_blacklist', type=str, dest='mass_blacklist', help='values for masses to skip', default=None)
  parser.add_argument('--coupling_whitelist', type=str, dest='coupling_whitelist', help='allowed values for couplings', default=None)
  parser.add_argument('--coupling_blacklist', type=str, dest='coupling_blacklist', help='values for couplings to skip', default=None)
  parser.add_argument('--pathDC', type=str, dest='pathDC', help='path to datacards to be analysed', default='./datacards')
  parser.add_argument('--wildcard', type=str, dest='wildcard', help='datacard generic string', default='datacard*hnl_m_12*.txt')
  parser.add_argument('--submit_batch', dest='submit_batch', help='submit on the batch?', action='store_true', default=False)
  return parser.parse_args()


# getting the parsed info
opt = getOptions()

version = opt.version

signal_type = opt.signal 

#years = opt.years.split(',') 

#channels = opt.channels.split(',')

path_to_datacards = opt.pathDC

#path_year = {}
#for year in years:
#  path_year[year] = '{p}/{y}'.format(p=path_to_datacards, y=year)

#path_channel = {}
#path_channel['mmm'] = 'mmm'
#path_channel['mem_os'] = 'mem_os'
#path_channel['mem_ss'] = 'mem_ss'
#path_channel['eee'] = 'eee'
#path_channel['eem_os'] = 'eem_os'
#path_channel['eem_ss'] = 'eem_ss'

datacard_wildcard = opt.wildcard

# create directories
if not path.exists('./datacards_combined/{}'.format(version)):
  os.system('mkdir -p ./datacards_combined/{}'.format(version))    
if not path.exists('./logs/{}'.format(version)):
  os.system('mkdir ./logs/{}'.format(version))

print 'loading cards...'
#all_datacards = {}
#for channel in channels:
#  all_datacards[channel] = {}
#  for year in years:
#    all_datacards[channel][year] = []

all_datacards = []

#for channel in channels:
#  for year in years:
    #all_datacards[channel][year] = glob('/'.join([path_year[year], path_channel[channels[0]], datacard_wildcard]))
all_datacards = glob('/'.join([path_to_datacards, datacard_wildcard]))

print '... datacards loaded'


#categories_to_combine = {}
#for channel in channels:
#  categories_to_combine[channel] = {}
#  for year in years:
#    categories_to_combine[channel][year] = []

# make sure that the displacements bins are correct
#for channel in channels:
  #for year in years:
    #categories_to_combine[channel][year] = OrderedDict(zip(['lxy_lt_0p5', 'lxy_0p5_to_1p5','lxy_1p5_to_4p0', 'lxy_mt_4p0'], ['{c}_{y}_disp1'.format(c=channel, y=year), '{c}_{y}_disp2'.format(c=channel, y=year), '{c}_{y}_disp3'.format(c=channel, y=year), '{c}_{y}_disp4'.format(c=channel, y=year)]))
  
categories_to_combine = ['lxy0to1_OS', 'lxy1to5_OS', 'lxygt5_OS', 'lxy0to1_SS', 'lxy1to5_SS', 'lxygt5_SS']

# nested dictionary with mass and coupling as keys
digested_datacards = OrderedDict()

# store results for 2D limits
limits2D = OrderedDict()

# will store paths to the datacards
#idc_datacard= OrderedDict()
#for channel in channels:
#  idc_datacard[channel] = {}
#  for year in years: 
#    idc_datacard[channel][year] = []

# makes the assumption that all years/channels have same sampling in mass/coupling/displacement
#the_set_datacards = all_datacards[channels[0]][years[0]]
the_set_datacards = all_datacards

#for card in all_datacards:
#  print card

for idc_ref in the_set_datacards:
    name = idc_ref.split('/')[-1]
#    for channel in channels:
#      for year in years:
#        idc_datacard[channel][year] = '/'.join([path_year[year], path_channel[channel], idc])

    # string mangling
    #print 'idc_ref ',idc_ref
    #print 'idc ',idc
    #name = idc.split('.')[0]
    #print 'name ',name
    #signal_name = re.findall(r'hnl_m_\d+_v2_\d+p\d+Em\d+', name)[0]
    #signal_mass = float(re.findall(r'\d+', re.findall(r'hnl_m_\d+_', signal_name)[0])[0])
    #signal_coupling_raw = re.findall(r'\d+', re.findall(r'_\d+p\d+Em\d+', signal_name)[0])
    #signal_coupling = float('%s.%se-%s' %(signal_coupling_raw[0], signal_coupling_raw[1], signal_coupling_raw[2]))
    signal_mass = name[name.find('m_')+2:name.find('_', name.find('m_')+2)]
    #print signal_mass
    signal_coupling = name[name.find('v2_')+3:name.find('.txt')]
    #print signal_coupling
   
    # get white/black listed mass/couplings
    if opt.mass_whitelist!='None':
      if str(signal_mass) not in opt.mass_whitelist.split(','): continue
    
    if opt.mass_blacklist!='None':
      if str(signal_mass) in opt.mass_blacklist.split(','): continue
    
    if opt.coupling_whitelist!='None':
      if str(signal_coupling) not in opt.coupling_whitelist.split(','): continue 
    
    if opt.coupling_blacklist!='None':
      if str(signal_coupling) in opt.coupling_blacklist.split(','): continue 
    

    # will fetch the datacards per year/channel
    if signal_mass not in digested_datacards.keys():
        digested_datacards[signal_mass] = OrderedDict()
    
    if signal_coupling not in digested_datacards[signal_mass].keys():
        #digested_datacards[signal_mass][signal_coupling] = OrderedDict()
        digested_datacards[signal_mass][signal_coupling] = []

    #for channel in channels:
    #  for year in years:
    #    if '{c}_{y}'.format(c=channel, y=year) not in digested_datacards[signal_mass][signal_coupling].keys():
    #      digested_datacards[signal_mass][signal_coupling]['{c}_{y}'.format(c=channel, y=year)] = []
    
        #if any([v in idc_datacard[channel][year] for v in categories_to_combine[channel][year].keys()]):
          #digested_datacards[signal_mass][signal_coupling]['{c}_{y}'.format(c=channel, y=year)].append(idc_datacard[channel][year]) 

    #digested_datacards[signal_mass][signal_coupling] = []
    #if any([v in idc_datacard[channel][year] for v in categories_to_combine[channel][year].keys()]):
    #if any([v in the_set_datacards for v in categories_to_combine]):
    for cat in categories_to_combine:
      if cat in idc_ref:
        print '{} {} {}'.format(signal_mass, signal_coupling, idc_ref)
        digested_datacards[signal_mass][signal_coupling].append(idc_ref) 
    
    
for mass, couplings in digested_datacards.iteritems():

    print 'mass =', mass
    
    v2s       = []
    obs       = []
    minus_two = []
    minus_one = []
    central   = []
    plus_one  = []
    plus_two  = []
    
    #datacards_to_combine = {}
    #for channel in channels:
    #  datacards_to_combine[channel] = {}
    #  for year in years:
    #    datacards_to_combine[channel][year] = []


    for coupling in couplings.keys():
        print '\tcoupling =', coupling
      
        # needed in case not all the years/channels have the same signal grid points
        #bad_channels = []

        #for channel in channels:
        #  for year in years:
            #datacards_to_combine[channel][year] = couplings[coupling]['{c}_{y}'.format(c=channel, y=year)]
        datacards_to_combine = digested_datacards[mass][coupling] #couplings[coupling]
            # check whether datacards for this given mass/coupling exist for each channel/year
            #for cat, idx in enumerate(datacards_to_combine[channel][year]):
            #  try:
            #    datacardtest = open(datacards_to_combine[channel][year][cat], 'r')
            #  except:
            #    print "WARNING: {} doesn't exist".format(datacards_to_combine[channel][year][cat])
            #    print "--> the grid point {m}-{c} will be ignored for {y}_{ch}".format(m=mass, c=coupling, y=year, ch=channel)
            #    bad_channels.append([year, channel])
                                                              
        # gonna combine the cards    
        command = 'combineCards.py'
        #for channel in channels:
        #  for year in years:
        #    if [year, channel] not in bad_channels: 
        #for cat in datacards_to_combine:
        #  print cat
        #for cat, idc in product(categories_to_combine, datacards_to_combine):
          #print 'idc ',idc
          #if cat in idc:
        #for icat, cat in enumerate(categories_to_combine):
        for idc in datacards_to_combine:
          #print idc
          #print '{} \t {}'.format(cat, idc)
          if any([v in idc for v in categories_to_combine]):
            #if cat in idc:
            #print 'bon'
            command += ' {}'.format(idc)

        #command += (' > datacards_combined/{v}/datacard_combined_{m}_{c}.txt' %(v=opt.version, m=str(mass), v=Decimal(coupling))).replace('-', 'm') 
        command += (' > datacards_combined/{v}/datacard_combined_m_{m}_v2_{c}.txt'.format(v=opt.version, m=str(mass), c=coupling)) 

        print command
        os.system(command)
        
        print ('\t\t -> combined datacards between the categories in datacards_combined/{v}/datacard_combined_m_{m}_v2_{c}.txt'.format(v=opt.version, m=str(mass), c=coupling))
