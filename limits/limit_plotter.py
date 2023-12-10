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
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt
from intersection import intersection
from utils import getMassList
sys.path.append('../scripts')
from decays import HNLDecays 

# use interpolation/spline 2d instead of 3d?

# create pair mass-value (median etc)
# create separate vectors for 1st/2nd/3rd intersections
# design strategy in to insert the vectors into one another to create one correctly sorted vector

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


class Points(object):
  def __init__(self, masses_1=None, masses_2=None, masses_3=None, masses_tot=None, values_1=None, values_2=None, values_3=None, values_tot=None):
    self.masses_1 = [masses_1]
    self.masses_2 = [masses_2]
    self.masses_3 = [masses_3]
    self.masses_tot = masses_tot
    self.values_1 = [values_1]
    self.values_2 = [values_2]
    self.values_3 = [values_3]
    self.values_tot = values_tot

    #if len(masses_1) == 1: self.masses_1 = [self.masses_1]
    #if len(values_1) == 1: self.values_1 = [self.values_1]


  def split_list(self, masses_lists, values_lists):
    masses_stuecke = []
    values_stuecke = []

    for i, masses_list in enumerate(masses_lists):
      start_idx = 0
      for j, mass in enumerate(masses_list):
        ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
        if j+1 <= len(masses_list)-1:
          reverse = False
          if masses_list[j+1] < masses_list[j]: reverse = True
          if reverse and masses_list[j+1] != self.mass_list[ref_idx-1] and masses_list[j+1] != masses_list[j]: 
            masses_stuecke.append(masses_list[start_idx:j+1])
            values_stuecke.append(values_lists[i][start_idx:j+1])
            start_idx = j+1
          elif not reverse and masses_list[j+1] != self.mass_list[ref_idx+1] and masses_list[j+1] != masses_list[j]: 
            masses_stuecke.append(masses_list[start_idx:j+1])
            values_stuecke.append(values_lists[i][start_idx:j+1])
            start_idx = j+1
        elif j == len(masses_list)-1:
            masses_stuecke.append(masses_list[start_idx:])
            values_stuecke.append(values_lists[i][start_idx:])

    #start_idx = 0
    #for i, mass in enumerate(masses):
    #  ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
    #  if i+1 <= len(masses)-1:
    #    reverse = False
    #    if masses[i+1] < masses[i]: reverse = True
    #    if reverse and masses[i+1] != self.mass_list[ref_idx-1] and masses[i+1] != masses[i]: 
    #      masses_stuecke.append(masses[start_idx:i+1])
    #      values_stuecke.append(values[start_idx:i+1])
    #      start_idx = i+1
    #    elif not reverse and masses[i+1] != self.mass_list[ref_idx+1] and masses[i+1] != masses[i]: 
    #      masses_stuecke.append(masses[start_idx:i+1])
    #      values_stuecke.append(values[start_idx:i+1])
    #      start_idx = i+1
    #  elif i == len(masses)-1:
    #      masses_stuecke.append(masses[start_idx:])
    #      values_stuecke.append(values[start_idx:])


    return masses_stuecke, values_stuecke



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
    self.plot_1D = True
    self.plot_scatter_figure = True
    self.plot_countour_figure = False
    self.no_exclusion_value = 2e-5
    #self.mass_list = [2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9]
    #self.mass_list = [2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9]
    self.mass_list = [1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.53, 1.56, 1.59, 1.62, 1.65, 1.68, 1.71, 1.74, 1.77, 1.8, 1.83, 1.86, 1.89, 1.92, 1.95, 1.98, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0]


    #masses = [2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.15, 3.2, 3.25, 3.3, 3.35, 3.35, 3.3, 3.25, 3.2, 3.15, 3.05, 3.0, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0]
    #start_idx = 0
    #for i, mass in enumerate(masses):
    #  ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
    #  if i+1 <= len(masses)-1:
    #    reverse = False
    #    if masses[i+1] < masses[i]: reverse = True
    #    if reverse and masses[i+1] != self.mass_list[ref_idx-1] and masses[i+1] != masses[i]: 
    #      print masses[start_idx:i+1]
    #      start_idx = i+1
    #    elif not reverse and masses[i+1] != self.mass_list[ref_idx+1] and masses[i+1] != masses[i]: 
    #      print masses[start_idx:i+1]
    #      start_idx = i+1
    #  elif i == len(masses)-1:
    #      print masses[start_idx:]

    #masses_1 =[2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.75, 3.9, 3.95]
    #masses_2 = [3.35, 3.3, 3.25, 3.2, 3.15, 3.05, 3.0, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0]
    #masses_3 = [3.0, 3.15, 3.2, 3.3, 3.35]
    #central_1 =[8.285167094634209e-05, 7.269881716046762e-05, 7.38980865550045e-05, 7.27598789444922e-05, 7.865094278461616e-05, 8.129976416183363e-05, 8.534900956985149e-05, 8.906858615782552e-05, 9.3181555785507e-05, 9.675002568920538e-05, 0.00010381556566877787, 9.01121723827393e-05, 0.00012486944563202377, 0.00013134068836221402, 0.0001356780763214012, 0.00015178203496595726, 0.00016902955456942539, 0.00019076005662924152, 0.00021250057410330143, 0.0002326229573804816, 0.00025338203268197846, 0.0002642123664431638, 0.00034965762317725553, 0.0003838163342779868, 0.0004115239066131318, 0.0005335779124314538, 0.0006451678307988079, 0.02583586690131893, 0.024326866511875926, 0.06274351343843303, 0.03582253543989143, 0.04452421835221904, 0.029593104479646564, 0.0751931013714551, 0.03566699029165148, 0.050321346067988576]
    #central_2 = [0.0010226986854741256, 0.0016929902718423521, 0.0022265502443276756, 0.0026341268718783847, 0.0029795903342881916, 0.004624298630692722, 0.00517875728336955, 0.006463735209292059, 0.006166783146800847, 0.006695417907491095, 0.007454099322105787, 0.008823399143590123, 0.010797657032069146, 0.1068774087557292, 0.032371458406655854, 0.015141675775902693, 0.02009883383463077, 0.03492439500215968, 0.05247430737665606, 0.05798308921401224, 0.06495971427465606, 0.11853927177971275, 0.13761696086357772, 0.11658101683467371, 0.4005247365389599]
    #central_3 = [0.06818633903139423, 0.058778376871920224, 0.0869582314084902, 0.032383098539866075, 0.06456445774287749]

    #print masses_1
    #print masses_2
    #print masses_3

    #masses_tot = []
    #values_tot = []
    #masses_tot_part1 = []
    #values_tot_part1 = []

    #idx_start = 0

    #for i, mass_1 in enumerate(masses_1):
    #  if mass_1 <= masses_2[0]:
    #    masses_tot_part1.append(mass_1)
    #    values_tot_part1.append(central_1[i])
    #  else:
    #    idx_start = i
    #    if central_2[0] < central_1[i+1]:
    #      for j, mass_2 in enumerate(masses_2):
    #        masses_tot_part1.append(mass_2)
    #        values_tot_part1.append(central_2[j])
    #    break
    #masses_tot.append(masses_tot_part1)
    #values_tot.append(values_tot_part1)

    #if idx_start < len(masses_1):
    #  masses_tot_part2 = []
    #  values_tot_part2 = []
    #  for i, mass_1 in enumerate(masses_1):
    #    if i < idx_start: continue
    #    masses_tot_part2.append(mass_1)
    #    values_tot_part2.append(mass_1)

    #masses_tot.append(masses_tot_part2)
    #values_tot.append(values_tot_part2)

    #print masses_tot

    #masses = [[2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.15, 3.2, 3.25, 3.3, 3.35, 3.35, 3.3, 3.25, 3.2, 3.15, 3.05, 3.0, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0], [3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.75, 3.9, 3.95]]
    #start_idx = 0
    #for i, masses_list in enumerate(masses_lists):
    #  for j, mass in enumerate(masses_list):
    #    ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
    #    if j+1 <= len(masses_list)-1:
    #      reverse = False
    #      if masses_list[j+1] < masses_list[j]: reverse = True
    #      if reverse and masses_list[j+1] != self.mass_list[ref_idx-1] and masses_list[j+1] != masses_list[j]: 
    #        masses_stuecke.append(masses_list[start_idx:j+1])
    #        values_stuecke.append(values_lists[i][start_idx:j+1])
    #        start_idx = j+1
    #      elif not reverse and masses_list[j+1] != self.mass_list[ref_idx+1] and masses_list[j+1] != masses_list[j]: 
    #        masses_stuecke.append(masses_list[start_idx:j+1])
    #        values_stuecke.append(values_lists[i][start_idx:j+1])
    #        start_idx = j+1
    #    elif j == len(masses_list)-1:
    #        masses_stuecke.append(masses_list[start_idx:])
    #        values_stuecke.append(values_lists[i][start_idx:])


    #masses = [2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0]
    #masses_central = [[2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05], [3.15, 3.2, 3.25, 3.3, 3.35, 3.35, 3.3, 3.25, 3.2, 3.15], [3.05, 3.0, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5], [2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0]]
    #print masses
    #print masses.index(2.0)
    #print masses_central

    #masses_central_all = []
    #for masses in masses_central:
    #  masses_central_all += masses
    #print masses_central_all
    #masses_central_all_np = np.array(masses_central_all)
    #indices = list(np.where(masses_central_all_np == 2.0)[0])
    #print indices

    #print masses_central_all.index(2.0) 


  def sortList(self, input):
    return float(input)


  def getCouplingTarget(self, mass, coupling):
    val_coupling = float(coupling)

    if self.do_coupling_scenario:
      val_fe = float(self.fe.replace('p', '.'))
      val_fu = float(self.fu.replace('p', '.'))
      val_ft = float(self.ft.replace('p', '.'))
      if val_fe == 0.3: val_fe = 1./3.
      if val_fu == 0.3: val_fu = 1./3.
      if val_ft == 0.3: val_ft = 1./3.
      #TODO check that 0p3 is 1/3

      decay_width_ini = HNLDecays(mass=float(mass), fe=0., fu=1., ft=0.).decay_rate['tot']
      decay_width_new = HNLDecays(mass=float(mass), fe=val_fe, fu=val_fu, ft=val_ft).decay_rate['tot']
      corr = decay_width_ini / decay_width_new
      #print corr
      val_coupling = corr * val_coupling 

    if self.scenario == 'Dirac':
      val_coupling = 2.0 * val_coupling
      #print 'dirac'

    return val_coupling


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
      intersection = self.no_exclusion_value
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


  def get_intersection_new(self, couplings, coupling_start, values, up_direction=False, crossing=1):
    '''
      Function that returns the coupling at which the limit intersects with 1 (crossing) in the log-log plane
      Argument 'up_direction': True means that the crossing is searched for increasing limit values
                               False means that the crossing is searched for decreasing limit values
    '''
    # first, search for couplings whose associated limit is the closest to up and down 1
    coupling_up = 0
    coupling_down = 1e9
    value_up = 0
    value_down = 0
    for icoupling, coupling in enumerate(couplings):
      if coupling < coupling_start: continue
      if icoupling+1 >= len(couplings): continue
      #print '{} {} {} {}'.format(values[icoupling], values[icoupling+1], crossing, values[icoupling] > crossing and values[icoupling+1] < crossing)
      if not up_direction and values[icoupling] >= crossing and values[icoupling+1] <= crossing:
        coupling_up = couplings[icoupling+1]
        value_up = values[icoupling+1]
        coupling_down = couplings[icoupling]
        value_down = values[icoupling]
        break
      elif up_direction and values[icoupling] <= crossing and values[icoupling+1] >= crossing:
        coupling_up = couplings[icoupling+1]
        value_up = values[icoupling+1]
        coupling_down = couplings[icoupling]
        value_down = values[icoupling]
        break

    #print 'coupling up: {} {}'.format(coupling_up, value_up)
    #print 'coupling down: {} {}'.format(coupling_down, value_down)

    # do not proceed if crossing is not found
    if coupling_up == 0 or coupling_down == 1e9:
      intersection = self.no_exclusion_value
    else:
      # then, search for the intersection with 1 in the log-log plane
      # powerlaw y = kx^m behaves as a linear law in the log-log plane: log(y) = mlog(x) + log(k)
      # get the slope m
      m = (math.log(value_up) - math.log(value_down)) / (math.log(coupling_up) - math.log(coupling_down))

      # and the coordinate k
      k = value_up / math.pow(coupling_up, m)

      # to finally compute the coupling where there is the intersection
      intersection = math.pow(float(crossing)/float(k), 1./float(m)) 

    #coupling_right = coupling_up if coupling_up > coupling_down else coupling_down
    coupling_right = intersection

    # do not consider coupling values above 1e-2 (with some margin)
    #if intersection > 5.e-2: intersection = self.no_exclusion_value
    #if intersection > 1.e-2: intersection = self.no_exclusion_value

    return intersection, coupling_right


  def get_intersection_list(self, couplings, values, crossing=1):
    '''
      Get the list of all the intersection points
    '''
    intersections = []
    intersection = 1e9 # default starting value
    coupling_start = 0.
    up_direction = False

    while intersection != self.no_exclusion_value:
      intersection, coupling_start = self.get_intersection_new(couplings=couplings, coupling_start=coupling_start, values=values, up_direction=up_direction, crossing=crossing)
      if intersection != self.no_exclusion_value:
        intersections.append(intersection)
      # switch direction
      if up_direction: up_direction = False
      elif not up_direction: up_direction = True

    return intersections

  def get_turning_point(self, masses, values, quantity):
    '''
      Returns the mass at which there is the turning point
    '''
    #exclusion_list_size = []
    #is_turning_mass = False
    turning_mass = 10. # defaults to the full mass range
    for mass in sorted(masses, key=self.sortList, reverse=True):
      if float(mass) < 3.: continue # turning point should find itself in the double exclusion region
      #exclusion_points = values[mass][quantity]
      exclusion_list_size = len(values[mass][quantity])
      if exclusion_list_size > 1:
        turning_mass = mass
        break

    return turning_mass

    #for imass, mass in enumerate(masses):
    #  exclusion_points = values[mass][quantity]
    #  if len(exclusion_points) > 


  def get_lower_limit(self, masses, values, quantity):
    mass_list = []
    exclusion_list = []

    mass_turning = self.get_turning_point(masses=masses, values=values, quantity=quantity)

    for mass in sorted(masses, key=self.sortList):
      if float(mass) > float(mass_turning): continue
      value_list = values[mass][quantity]
      if len(value_list) > 0:
        exclusion_list.append(value_list[0])
        mass_list.append(float(mass))

    return mass_list, exclusion_list
      

  def get_upper_limit(self, masses, values, quantity):
    mass_list = []
    exclusion_list = []

    counter = 0
    for mass in sorted(masses, key=self.sortList, reverse=True):
      value_list = values[mass][quantity]
      if len(value_list) > 1:
        # plot the upper limit only in the double exclusion regime (3 -4 GeV)
        #if counter == 0 and float(mass) < 3: 
        #  break
        counter += 1
        exclusion_list.append(value_list[1])
        mass_list.append(float(mass))

    return mass_list, exclusion_list


  def get_total_limit(self, masses_lower, masses_upper, limits_lower, limits_upper):
    '''
      Combine the upper and lower limits
    '''
    #TODO address case where points are missing
    masses_all = masses_lower + masses_upper
    limits_all = limits_lower + limits_upper

    return masses_all, limits_all


  def get_missing_points(self, masses, values, quantity):
    mass_list = []
    exclusion_list = []

    for mass in sorted(masses, key=self.sortList):
      value_list = values[mass][quantity]
      if len(value_list) == 0:
        exclusion_list.append(2e-5)
        mass_list.append(float(mass))

    return mass_list, exclusion_list


  #def get_intersection_up(self, couplings, coupling_start, values, crossing=1):
  #  '''
  #    Function that returns the coupling at which the limit intersects with 1 (crossing) in the log-log plane
  #  '''
  #  # first, search for couplings whose associated limit is the closest to up and down 1
  #  coupling_up = 0
  #  coupling_down = 1e9
  #  value_up = 0
  #  value_down = 0
  #  for icoupling, coupling in enumerate(couplings):
  #    if coupling < coupling_start: continue
  #    if icoupling+1 >= len(couplings): continue
  #    #print '{} {} {} {}'.format(values[icoupling], values[icoupling+1], crossing, values[icoupling] > crossing and values[icoupling+1] < crossing)
  #    if values[icoupling] <= crossing and values[icoupling+1] >= crossing:
  #      coupling_up = couplings[icoupling+1]
  #      value_up = values[icoupling+1]
  #      coupling_down = couplings[icoupling]
  #      value_down = values[icoupling]
  #      break

  #  #print 'coupling up: {} {}'.format(coupling_up, value_up)
  #  #print 'coupling down: {} {}'.format(coupling_down, value_down)

  #  # do not proceed if crossing is not found
  #  if coupling_up == 0 or coupling_down == 1e9:
  #    intersection = -99
  #  else:
  #    # then, search for the intersection with 1 in the log-log plane
  #    # powerlaw y = kx^m behaves as a linear law in the log-log plane: log(y) = mlog(x) + log(k)
  #    # get the slope m
  #    m = (math.log(value_up) - math.log(value_down)) / (math.log(coupling_up) - math.log(coupling_down))

  #    # and the coordinate k
  #    k = value_up / math.pow(coupling_up, m)

  #    # to finally compute the coupling where there is the intersection
  #    intersection = math.pow(float(crossing)/float(k), 1./float(m)) 

  #  #coupling_right = coupling_up if coupling_up > coupling_down else coupling_down
  #  coupling_right = intersection

  #  return intersection, coupling_right



  def smooth(self, y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


  def get_lists(self, limits2D, list_name):
    quantity_1 = []
    quantity_2 = []
    quantity_3 = []
    quantity_missing = []
    masses_1 = []
    masses_2 = []
    masses_3 = []
    masses_missing = []

    cutoff = 1e-1#2e-2

    for mass in sorted(limits2D.keys(), key=self.sortList):
      if len(limits2D[mass][list_name]) > 0:
        if limits2D[mass][list_name][0] > cutoff: continue 
        quantity_1.append(limits2D[mass][list_name][0])
        masses_1.append(float(mass))
      else:
        quantity_missing.append(self.no_exclusion_value)
        masses_missing.append(float(mass))

      if len(limits2D[mass][list_name]) > 1:
        if limits2D[mass][list_name][1] > cutoff: continue 
        quantity_2.append(limits2D[mass][list_name][1])
        masses_2.append(float(mass))

      if len(limits2D[mass][list_name]) > 2:
        # do not consider third exclusion above 1e-2
        #if limits2D[mass][list_name][2] <= 1e-2: 
        if limits2D[mass][list_name][2] > cutoff: continue 
        quantity_3.append(limits2D[mass][list_name][2])
        masses_3.append(float(mass))

    # reverse intermediate exclusion
    masses_2.reverse()
    quantity_2.reverse()

    return masses_1, masses_2, masses_3, masses_missing, quantity_1, quantity_2, quantity_3, quantity_missing


  def concatenate_lists(self, masses_1, masses_2, masses_3, values_1, values_2, values_3):
    masses_tot = []
    values_tot = []
    masses_tot_part1 = []
    values_tot_part1 = []

    idx_start = 0

    for i, mass_1 in enumerate(masses_1):
      #if mass_1 > 3e-2: continue # remove point above cutoff
      if mass_1 < masses_2[0]:
        masses_tot_part1.append(mass_1)
        values_tot_part1.append(values_1[i])
      #else:
      elif mass_1 == masses_2[0]:
        masses_tot_part1.append(mass_1)
        values_tot_part1.append(values_1[i])
        idx_start = i
        if i+1 < len(masses_1)-1 and values_2[0] < values_1[i+1]:
          for j, mass_2 in enumerate(masses_2):
            masses_tot_part1.append(mass_2)
            values_tot_part1.append(values_2[j])
        elif i+1 >= len(masses_1)-1: 
          for j, mass_2 in enumerate(masses_2):
            masses_tot_part1.append(mass_2)
            values_tot_part1.append(values_2[j])
        else:
          print 'non apply'
        break
    masses_tot.append(masses_tot_part1)
    values_tot.append(values_tot_part1)

    if idx_start+1 < len(masses_1):
      masses_tot_part2 = []
      values_tot_part2 = []
      for i, mass_1 in enumerate(masses_1):
        #if mass_1 > 3e-2: continue # remove point above cutoff
        if i < idx_start+1: continue
        masses_tot_part2.append(mass_1)
        values_tot_part2.append(values_1[i])

      masses_tot.append(masses_tot_part2)
      values_tot.append(values_tot_part2)

    points = Points(
        masses_1 = masses_1, 
        masses_2 = masses_2, 
        masses_3 = masses_3, 
        masses_tot = masses_tot, 
        values_1 = values_1, 
        values_2 = values_2, 
        values_3 = values_3, 
        values_tot = values_tot, 
        )

    return points
    #return masses_tot, values_tot


  def split_list(self, masses_lists, values_lists):
    masses_stuecke = []
    values_stuecke = []

    for i, masses_list in enumerate(masses_lists):
      start_idx = 0
      for j, mass in enumerate(masses_list):
        ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
        if j+1 <= len(masses_list)-1:
          reverse = False
          if masses_list[j+1] < masses_list[j]: reverse = True
          if reverse and masses_list[j+1] != self.mass_list[ref_idx-1] and masses_list[j+1] != masses_list[j]: 
            masses_stuecke.append(masses_list[start_idx:j+1])
            values_stuecke.append(values_lists[i][start_idx:j+1])
            start_idx = j+1
          elif not reverse and masses_list[j+1] != self.mass_list[ref_idx+1] and masses_list[j+1] != masses_list[j]: 
            masses_stuecke.append(masses_list[start_idx:j+1])
            values_stuecke.append(values_lists[i][start_idx:j+1])
            start_idx = j+1
        elif j == len(masses_list)-1:
            masses_stuecke.append(masses_list[start_idx:])
            values_stuecke.append(values_lists[i][start_idx:])

    #start_idx = 0
    #for i, mass in enumerate(masses):
    #  ref_idx = self.mass_list.index(mass) # index of mass in reference mass_list
    #  if i+1 <= len(masses)-1:
    #    reverse = False
    #    if masses[i+1] < masses[i]: reverse = True
    #    if reverse and masses[i+1] != self.mass_list[ref_idx-1] and masses[i+1] != masses[i]: 
    #      masses_stuecke.append(masses[start_idx:i+1])
    #      values_stuecke.append(values[start_idx:i+1])
    #      start_idx = i+1
    #    elif not reverse and masses[i+1] != self.mass_list[ref_idx+1] and masses[i+1] != masses[i]: 
    #      masses_stuecke.append(masses[start_idx:i+1])
    #      values_stuecke.append(values[start_idx:i+1])
    #      start_idx = i+1
    #  elif i == len(masses)-1:
    #      masses_stuecke.append(masses[start_idx:])
    #      values_stuecke.append(values[start_idx:])


    return masses_stuecke, values_stuecke


  def get_islands_coordinates(self, masses_2_stuecke, masses_3_stuecke, values_2_stuecke, values_3_stuecke):
    islands_coordinates = []
    for i, masses_2 in enumerate(masses_2_stuecke):
      if len(masses_2) == 1: # single point
        #print 'single'
        #print masses_2
        if masses_2 in masses_3_stuecke: # single point in both the middle and upper limit
          ref_idx = self.mass_list.index(masses_2[0]) # index of mass in reference mass_list
          #center = (masses_2[0], (values_3_stuecke[i][0] - values_2_stuecke[i][0]) / 2.)
          #center_x = masses_2[0]
          #center_y = (values_3_stuecke[i][0] - values_2_stuecke[i][0]) / 2.
          #width = self.mass_list[ref_idx+1] - self.mass_list[ref_idx-1]
          #height = values_3_stuecke[i][0] - values_2_stuecke[i][0]
          #island_coordinates = [center_x, center_y, width, height]
          #islands_coordinates.append(island_coordinates)

          x =  self.mass_list[ref_idx-1]
          #x =  (self.mass_list[ref_idx-1] - masses_2[0]) / 2.
          y = values_2_stuecke[i][0]
          #width = ((self.mass_list[ref_idx+1] - masses_2[0]) / 2.) - ((masses_2[0]- self.mass_list[ref_idx-1]) / 2.)
          width = self.mass_list[ref_idx+1] - self.mass_list[ref_idx-1]
          height =  values_3_stuecke[i][0] - values_2_stuecke[i][0]
          #island_coordinates = [center_x, center_y, width, height]
          island_coordinates = [x, y, width, height]
          islands_coordinates.append(island_coordinates)
        else:
          island_coordinates = None

    return islands_coordinates


  def find_boundaries(self, masses_tot, masses_stuecke_central, values_tot, values_stuecke_central):
    #print '\n'
    #print masses_tot
    #print masses_stuecke_central

    ## aggregate all the chunks
    #masses_central_all = []
    #for masses in masses_stuecke_central:
    #  masses_central_all += masses

    #values_central_all = []
    #for values in values_stuecke_central:
    #  values_central_all += values

    ## needed to search for multiple indices
    #masses_central_all_np = np.array(masses_central_all)

    #boundaries = []
    #for i, mass in enumerate(masses_tot):
    #  print '\n'
    #  print mass
    #  print values_tot[i]
    #  # find indices in central list
    #  indices = list(np.where(masses_central_all_np == mass)[0])
    #  if len(indices) == 1:
    #    print 'len 1'
    #    value = values_central_all[indices[0]]
    #    boundaries.append(value)
    #  elif len(indices) == 2:
    #    print 'len 2'
    #    value_1 = values_central_all[indices[0]]
    #    value_2 = values_central_all[indices[1]]
    #    print 'val 1: {} val 2: {}'.format(value_1, value_2)
    #    print 'diff 1:{} diff 2: {}'.format(abs(values_tot[i] - value_1)/value_1, abs(values_tot[i] - value_2)/value_2)
    #    if abs(values_tot[i] - value_1) /value_1 < abs(values_tot[i] - value_2) / value_2:
    #      value = value_1
    #    else:
    #      value = value_2
    #    print 'value: {}'.format(value)
    #    boundaries.append(value)
    #  else:
    #    print 'no central value found'

    boundaries = []
    for i, mass in enumerate(masses_tot):
      central_exists = False
      for i, mass_stueck_central in enumerate(masses_stuecke_central):
        for j, mass_central in enumerate(mass_stueck_central):
          if mass == mass_central:
            central_exists = True
            value = values_stuecke_central[i][j]
            if value not in boundaries:
              boundaries.append(values_stuecke_central[i][j])
              #print 'append'
      #if not central_exists: boundaries.append(values_tot[i])
      if not central_exists: boundaries.append(2e-2) # or average of adjacent points

    #print boundaries

    return boundaries


  def plot_scatter(self, masses_1, masses_2, masses_3, masses_missing_1, values_1, values_2, values_3, values_missing_1, label, plotDir):
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}, {}'.format(label, coupling_scenario, self.scenario), loc='right', fontsize=23)

    plt.plot(masses_1, values_1, marker='X', color='red', label='Median expected', linewidth=0)
    plt.plot(masses_2, values_2, marker='X', color='blue', label='Median expected', linewidth=0)
    plt.plot(masses_3, values_3, marker='X', color='darkgreen', label='Median expected', linewidth=0)
    plt.plot(masses_missing_1, values_missing_1, marker='o', color='black', label='Median expected', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    masses_tot = masses_1 + masses_2 + masses_3
    #plt.xlim(min(masses_tot), max(masses_tot))
    plt.xlim(2., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'scatter_{}'.format(self.scenario) 
    else:
      name_2d = 'scatter_scenario_{}_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft, label) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


  def plot_scatter_all(self, points_all, plotDir):
    plt.clf()
    f, ax = plt.subplots(figsize=(15, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}'.format(coupling_scenario, self.scenario), loc='right', fontsize=23)

    masses_central_1_chunks, values_central_1_chunks = self.split_list(masses_lists=points_all['central'].masses_1, values_lists=points_all['central'].values_1)
    masses_central_2_chunks, values_central_2_chunks = self.split_list(masses_lists=points_all['central'].masses_2, values_lists=points_all['central'].values_2)
    masses_central_3_chunks, values_central_3_chunks = self.split_list(masses_lists=points_all['central'].masses_3, values_lists=points_all['central'].values_3)
    for i, chunk in enumerate(masses_central_1_chunks):
      plt.plot(masses_central_1_chunks[i], values_central_1_chunks[i], marker='X', color='red', label='Median expected', linewidth=0)
    for i, chunk in enumerate(masses_central_2_chunks):
      plt.plot(masses_central_2_chunks[i], values_central_2_chunks[i], marker='X', color='red', label='Median expected', linewidth=0)
    for i, chunk in enumerate(masses_central_3_chunks):
      plt.plot(masses_central_3_chunks[i], values_central_3_chunks[i], marker='X', color='red', label='Median expected', linewidth=0)

    masses_plus_two_1_chunks, values_plus_two_1_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_1, values_lists=points_all['plus_two'].values_1)
    masses_plus_two_2_chunks, values_plus_two_2_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_2, values_lists=points_all['plus_two'].values_2)
    masses_plus_two_3_chunks, values_plus_two_3_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_3, values_lists=points_all['plus_two'].values_3)
    for i, chunk in enumerate(masses_plus_two_1_chunks):
      plt.plot(masses_plus_two_1_chunks[i], values_plus_two_1_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
    for i, chunk in enumerate(masses_plus_two_2_chunks):
      plt.plot(masses_plus_two_2_chunks[i], values_plus_two_2_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
    for i, chunk in enumerate(masses_plus_two_3_chunks):
      plt.plot(masses_plus_two_3_chunks[i], values_plus_two_3_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)

    masses_plus_one_1_chunks, values_plus_one_1_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_1, values_lists=points_all['plus_one'].values_1)
    masses_plus_one_2_chunks, values_plus_one_2_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_2, values_lists=points_all['plus_one'].values_2)
    masses_plus_one_3_chunks, values_plus_one_3_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_3, values_lists=points_all['plus_one'].values_3)
    for i, chunk in enumerate(masses_plus_one_1_chunks):
      plt.plot(masses_plus_one_1_chunks[i], values_plus_one_1_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
    for i, chunk in enumerate(masses_plus_one_2_chunks):
      plt.plot(masses_plus_one_2_chunks[i], values_plus_one_2_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
    for i, chunk in enumerate(masses_plus_one_3_chunks):
      plt.plot(masses_plus_one_3_chunks[i], values_plus_one_3_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
      
    masses_minus_two_1_chunks, values_minus_two_1_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_1, values_lists=points_all['minus_two'].values_1)
    masses_minus_two_2_chunks, values_minus_two_2_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_2, values_lists=points_all['minus_two'].values_2)
    masses_minus_two_3_chunks, values_minus_two_3_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_3, values_lists=points_all['minus_two'].values_3)
    for i, chunk in enumerate(masses_minus_two_1_chunks):
      plt.plot(masses_minus_two_1_chunks[i], values_minus_two_1_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
    for i, chunk in enumerate(masses_minus_two_2_chunks):
      plt.plot(masses_minus_two_2_chunks[i], values_minus_two_2_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
    for i, chunk in enumerate(masses_minus_two_3_chunks):
      plt.plot(masses_minus_two_3_chunks[i], values_minus_two_3_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)

    masses_minus_one_1_chunks, values_minus_one_1_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_1, values_lists=points_all['minus_one'].values_1)
    masses_minus_one_2_chunks, values_minus_one_2_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_2, values_lists=points_all['minus_one'].values_2)
    masses_minus_one_3_chunks, values_minus_one_3_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_3, values_lists=points_all['minus_one'].values_3)
    for i, chunk in enumerate(masses_minus_one_1_chunks):
      plt.plot(masses_minus_one_1_chunks[i], values_minus_one_1_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
    for i, chunk in enumerate(masses_minus_one_2_chunks):
      plt.plot(masses_minus_one_2_chunks[i], values_minus_one_2_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
    for i, chunk in enumerate(masses_minus_one_3_chunks):
      plt.plot(masses_minus_one_3_chunks[i], values_minus_one_3_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)

    if not self.do_blind:
      masses_obs_1_chunks, values_obs_1_chunks = self.split_list(masses_lists=points_all['obs'].masses_1, values_lists=points_all['obs'].values_1)
      masses_obs_2_chunks, values_obs_2_chunks = self.split_list(masses_lists=points_all['obs'].masses_2, values_lists=points_all['obs'].values_2)
      masses_obs_3_chunks, values_obs_3_chunks = self.split_list(masses_lists=points_all['obs'].masses_3, values_lists=points_all['obs'].values_3)
      for i, chunk in enumerate(masses_obs_1_chunks):
        plt.plot(masses_obs_1_chunks[i], values_obs_1_chunks[i], marker='X', color='black', label='Observed', linewidth=0)
      for i, chunk in enumerate(masses_obs_2_chunks):
        plt.plot(masses_obs_2_chunks[i], values_obs_2_chunks[i], marker='X', color='black', label='Observed', linewidth=0)
      for i, chunk in enumerate(masses_obs_3_chunks):
        plt.plot(masses_obs_3_chunks[i], values_obs_3_chunks[i], marker='X', color='black', label='Observed', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    plt.xlim(2., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'countour_{}'.format(self.scenario) 
    else:
      name_2d = 'scatter_all_scenario_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


  def plot_countour(self, points, masses_missing, values_missing, label, plotDir):
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}, {}'.format(label, coupling_scenario, self.scenario), loc='right', fontsize=23)

    #islands_coordinates = self.get_islands_coordinates(masses_2_stuecke, masses_3_stuecke, values_2_stuecke, values_3_stuecke)
    #for island_coordinates in islands_coordinates:
    #  #a = island_coordinates[0]
    #  #b = island_coordinates[1]
    #  #c = island_coordinates[2]
    #  #d = island_coordinates[3]
    #  #center_x = island_coordinates[0]
    #  #center_y = island_coordinates[1]
    #  #width = island_coordinates[2]
    #  #height = island_coordinates[3]
    #  x = island_coordinates[0]
    #  y = island_coordinates[1]
    #  width = island_coordinates[2]
    #  height = island_coordinates[3]
    #  #print center_x, center_y, width, height

    #  #from matplotlib.transforms import ScaledTranslation
    #  ## Ellipse centre coordinates
    #  ##x, y = 3, 8

    #  ## use the axis scale tform to figure out how far to translate 
    #  #ell_offset = ScaledTranslation(center_x, center_y, ax.transScale)

    #  ## construct the composite tform
    #  #ell_tform = ell_offset + ax.transLimits + ax.transAxes

    #  plt.gca().add_patch(Rectangle((x, y), width, height, edgecolor='red', facecolor='white', fill=False, linewidth=2, zorder=3)) 
    #  #plt.gca().add_patch(FancyBboxPatch((x, y), width, height, edgecolor='red', boxstyle='round', fill=False, mutation_scale=0.1, linewidth=2, zorder=3)) 
    #  #plt.gca().add_patch(Ellipse((center_x, center_y), width, height, edgecolor='red', facecolor='white')) 
    #  #plt.gca().add_patch(Ellipse((center_x, center_y), width, height, edgecolor='red', transform=ell_tform)) 
    #  #ax.add_patch(Ellipse((center_x, center_y), width, 0.1, edgecolor='red', transform=ell_tform)) 
    #  #ellipse = plt.Ellipse((center_x, center_y), width, 0.1, edgecolor='red', transform=ell_tform)
    #  #ax.add_artist(ellipse)
    #  #plt.gca().add_patch(Ellipse((3.5, 1e-3), 1, 0.1, edgecolor='red', facecolor='white')) 

    masses_tot_chunks, values_tot_chunks = self.split_list(masses_lists=points.masses_tot, values_lists=points.values_tot)
    for i, chunk in enumerate(masses_tot_chunks):
      plt.plot(masses_tot_chunks[i], values_tot_chunks[i], marker='X', color='red', label='Median expected', linewidth=0)
      plt.plot(masses_tot_chunks[i], values_tot_chunks[i], color='red', label='Median expected', linewidth=2)

    plt.plot(masses_missing, values_missing, marker='o', color='black', label='Median expected', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    #masses_tot = masses_1 + masses_2 + masses_3
    #plt.xlim(min(masses_tot), max(masses_tot))
    plt.xlim(1., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'countour_{}'.format(self.scenario) 
    else:
      name_2d = 'countour_scenario_{}_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft, label) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


  #def plot_countour(self, masses_1, masses_2, masses_3, masses_missing_1, values_1, values_2, values_3, values_missing_1, label, plotDir):
  def plot_countour_1(self, masses_stuecke, masses_missing, values_stuecke, values_missing, label, plotDir):
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}, {}'.format(label, coupling_scenario, self.scenario), loc='right', fontsize=23)

    #masses_1_stuecke, values_1_stuecke = self.split_list(masses=masses_1, values=values_1)
    #masses_2_stuecke, values_2_stuecke = self.split_list(masses=masses_2, values=values_2)
    #masses_3_stuecke, values_3_stuecke = self.split_list(masses=masses_3, values=values_3)

    #islands_coordinates = self.get_islands_coordinates(masses_2_stuecke, masses_3_stuecke, values_2_stuecke, values_3_stuecke)
    #for island_coordinates in islands_coordinates:
    #  #a = island_coordinates[0]
    #  #b = island_coordinates[1]
    #  #c = island_coordinates[2]
    #  #d = island_coordinates[3]
    #  #center_x = island_coordinates[0]
    #  #center_y = island_coordinates[1]
    #  #width = island_coordinates[2]
    #  #height = island_coordinates[3]
    #  x = island_coordinates[0]
    #  y = island_coordinates[1]
    #  width = island_coordinates[2]
    #  height = island_coordinates[3]
    #  #print center_x, center_y, width, height

    #  #from matplotlib.transforms import ScaledTranslation
    #  ## Ellipse centre coordinates
    #  ##x, y = 3, 8

    #  ## use the axis scale tform to figure out how far to translate 
    #  #ell_offset = ScaledTranslation(center_x, center_y, ax.transScale)

    #  ## construct the composite tform
    #  #ell_tform = ell_offset + ax.transLimits + ax.transAxes

    #  plt.gca().add_patch(Rectangle((x, y), width, height, edgecolor='red', facecolor='white', fill=False, linewidth=2, zorder=3)) 
    #  #plt.gca().add_patch(FancyBboxPatch((x, y), width, height, edgecolor='red', boxstyle='round', fill=False, mutation_scale=0.1, linewidth=2, zorder=3)) 
    #  #plt.gca().add_patch(Ellipse((center_x, center_y), width, height, edgecolor='red', facecolor='white')) 
    #  #plt.gca().add_patch(Ellipse((center_x, center_y), width, height, edgecolor='red', transform=ell_tform)) 
    #  #ax.add_patch(Ellipse((center_x, center_y), width, 0.1, edgecolor='red', transform=ell_tform)) 
    #  #ellipse = plt.Ellipse((center_x, center_y), width, 0.1, edgecolor='red', transform=ell_tform)
    #  #ax.add_artist(ellipse)
    #  #plt.gca().add_patch(Ellipse((3.5, 1e-3), 1, 0.1, edgecolor='red', facecolor='white')) 

    for i, masses_1 in enumerate(masses_stuecke):
      plt.plot(masses_stuecke[i], values_stuecke[i], marker='X', color='red', label='Median expected', linewidth=0)
      plt.plot(masses_stuecke[i], values_stuecke[i], color='red', label='Median expected', linewidth=2)

    #plt.plot(masses_2, values_2, marker='X', color='blue', label='Median expected', linewidth=0)
    #for i, masses_2 in enumerate(masses_2_stuecke):
    #  plt.plot(masses_2_stuecke[i], values_2_stuecke[i], color='blue', label='Median expected', linewidth=2)

    #plt.plot(masses_3, values_3, marker='X', color='darkgreen', label='Median expected', linewidth=0)
    #for i, masses_3 in enumerate(masses_3_stuecke):
    #  plt.plot(masses_3_stuecke[i], values_3_stuecke[i], color='darkgreen', label='Median expected', linewidth=2)

    plt.plot(masses_missing, values_missing, marker='o', color='black', label='Median expected', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    #masses_tot = masses_1 + masses_2 + masses_3
    #plt.xlim(min(masses_tot), max(masses_tot))
    plt.xlim(1., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'countour_{}'.format(self.scenario) 
    else:
      name_2d = 'countour_scenario_{}_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft, label) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


  def plot_2dlimit(self, points_all, plotDir):
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}'.format(coupling_scenario, self.scenario), loc='right', fontsize=23)

    #print '\n'
    #print masses_stuecke_all_1['plus_two'][1]
    #print values_stuecke_all_1['plus_two'][1]
    #print '\n'
    #print masses_stuecke_all_2['plus_two'][0]
    #print values_stuecke_all_2['plus_two'][0]
    #print '\n'
    #print masses_stuecke_all_3['plus_two'][0]
    #print values_stuecke_all_3['plus_two'][0]

    #[3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65]
    #[3.65, 3.6, 3.55, 3.5]
    #[3.5, 3.55, 3.6, 3.65]
    #[0.00012034689271124365, 0.00012425505596183535, 0.00012387651316755617, 0.0001357589514113421, 0.00012398053056366553, 0.0001621255964163732, 0.00019170782693026195, 0.00019805437170792586, 0.0006330616450559552, 0.0002532928957738135, 0.0003514183210573947]
    #[0.0011121770802836684, 0.0019333104686060747, 0.0014249369293746052, 0.002675377363665742]
    #[0.00551420027591585, 0.010041985319261523, 0.0046999341557779975, 0.004726644830799094]

    #masses_tot = masses_stuecke_all_1['plus_two'][1] + masses_stuecke_all_2['plus_two'][0] + masses_stuecke_all_3['plus_two'][0]
    #values_tot = values_stuecke_all_1['plus_two'][1] + values_stuecke_all_2['plus_two'][0] + values_stuecke_all_3['plus_two'][0]
    #masses_tot = [3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.65, 3.6, 3.55, 3.5, 3.5, 3.55, 3.6, 3.65, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9]
    #values_tot = [0.00012034689271124365, 0.00012425505596183535, 0.00012387651316755617, 0.0001357589514113421, 0.00012398053056366553, 0.0001621255964163732, 0.00019170782693026195, 0.00019805437170792586, 0.0006330616450559552, 0.0002532928957738135, 0.0003514183210573947, 0.0011121770802836684, 0.0019333104686060747, 0.0014249369293746052, 0.002675377363665742, 0.00551420027591585, 0.010041985319261523, 0.0046999341557779975, 0.004726644830799094, 0.010422340034747238, 0.009992252977129756, 0.009730473592131002, 0.00823294826799931, 0.009853815531909501, 0.010269952092311059, 0.010586853974853856, 0.010166377180102508, 0.009314321439711618, 0.00849178867340154, 0.011485044420968676, 0.00921101394854829, 0.012493053272891188, 0.014662758275220093, 0.02478249047358668]
    #boundaries = self.find_boundaries(masses_tot, masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #plt.fill_between(masses_tot, values_tot, boundaries, color='gold', label=r'95% expected')

    # readd masses_central_1 and _2 to fill colour inbetween

    masses_central_tot_chunks, values_central_tot_chunks = self.split_list(masses_lists=points_all['central'].masses_tot, values_lists=points_all['central'].values_tot)
    masses_central_1_chunks, values_central_1_chunks = self.split_list(masses_lists=points_all['central'].masses_1, values_lists=points_all['central'].values_1)
    masses_central_2_chunks, values_central_2_chunks = self.split_list(masses_lists=points_all['central'].masses_2, values_lists=points_all['central'].values_2)
    for i, chunk in enumerate(masses_central_tot_chunks):
      plt.plot(masses_central_tot_chunks[i], values_central_tot_chunks[i], marker='X', color='red', label='Median expected', linewidth=0)
      plt.plot(masses_central_tot_chunks[i], values_central_tot_chunks[i], color='red', label='Median expected', linewidth=2)

    masses_plus_two_tot_chunks, values_plus_two_tot_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_tot, values_lists=points_all['plus_two'].values_tot)
    masses_plus_two_1_chunks, values_plus_two_1_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_1, values_lists=points_all['plus_two'].values_1)
    masses_plus_two_2_chunks, values_plus_two_2_chunks = self.split_list(masses_lists=points_all['plus_two'].masses_2, values_lists=points_all['plus_two'].values_2)
    for i, chunk in enumerate(masses_plus_two_tot_chunks):
      plt.plot(masses_plus_two_tot_chunks[i], values_plus_two_tot_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
      plt.plot(masses_plus_two_tot_chunks[i], values_plus_two_tot_chunks[i], color='gold', label=r'95% expected', linewidth=2)
      #boundaries_plus_two_1 = self.find_boundaries(masses_plus_two_1_chunks[i], masses_central_1_chunks, values_plus_two_1_chunks[i], values_central_1_chunks)
      #plt.fill_between(masses_plus_two_1_chunks[i], values_plus_two_1_chunks[i], boundaries_plus_two_1, color='gold', label=r'95% expected')
      #boundaries_plus_two_2 = self.find_boundaries(masses_plus_two_2_chunks[i], masses_central_2_chunks, values_plus_two_2_chunks[i], values_central_2_chunks)
      #plt.fill_between(masses_plus_two_2_chunks[i], values_plus_two_2_chunks[i], boundaries_plus_two_2, color='gold', label=r'95% expected')

    masses_plus_one_tot_chunks, values_plus_one_tot_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_tot, values_lists=points_all['plus_one'].values_tot)
    masses_plus_one_1_chunks, values_plus_one_1_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_1, values_lists=points_all['plus_one'].values_1)
    masses_plus_one_2_chunks, values_plus_one_2_chunks = self.split_list(masses_lists=points_all['plus_one'].masses_2, values_lists=points_all['plus_one'].values_2)
    for i, chunk in enumerate(masses_plus_one_tot_chunks):
      plt.plot(masses_plus_one_tot_chunks[i], values_plus_one_tot_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
      plt.plot(masses_plus_one_tot_chunks[i], values_plus_one_tot_chunks[i], color='forestgreen', label=r'68% expected', linewidth=2)
      #boundaries_plus_one_1 = self.find_boundaries(masses_plus_one_1_chunks[i], masses_central_1_chunks, values_plus_one_1_chunks[i], values_central_1_chunks)
      #plt.fill_between(masses_plus_one_1_chunks[i], values_plus_one_1_chunks[i], boundaries_plus_one_1, color='forestgreen', label=r'68% expected')
      #boundaries_plus_one_2 = self.find_boundaries(masses_plus_one_2_chunks[i], masses_central_2_chunks, values_plus_one_2_chunks[i], values_central_2_chunks)
      #plt.fill_between(masses_plus_one_2_chunks[i], values_plus_one_2_chunks[i], boundaries_plus_one_2, color='forestgreen', label=r'68% expected')
      
    masses_minus_two_tot_chunks, values_minus_two_tot_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_tot, values_lists=points_all['minus_two'].values_tot)
    masses_minus_two_1_chunks, values_minus_two_1_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_1, values_lists=points_all['minus_two'].values_1)
    masses_minus_two_2_chunks, values_minus_two_2_chunks = self.split_list(masses_lists=points_all['minus_two'].masses_2, values_lists=points_all['minus_two'].values_2)
    for i, chunk in enumerate(masses_minus_two_tot_chunks):
      plt.plot(masses_minus_two_tot_chunks[i], values_minus_two_tot_chunks[i], marker='X', color='gold', label=r'95% expected', linewidth=0)
      plt.plot(masses_minus_two_tot_chunks[i], values_minus_two_tot_chunks[i], color='gold', label=r'95% expected', linewidth=2)
      #boundaries_minus_two_1 = self.find_boundaries(masses_minus_two_1_chunks[i], masses_central_1_chunks, values_minus_two_1_chunks[i], values_central_1_chunks)
      #plt.fill_between(masses_minus_two_1_chunks[i], values_minus_two_1_chunks[i], boundaries_minus_two_1, color='gold', label=r'95% expected')
      #boundaries_minus_two_2 = self.find_boundaries(masses_minus_two_2_chunks[i], masses_central_2_chunks, values_minus_two_2_chunks[i], values_central_2_chunks)
      #plt.fill_between(masses_minus_two_2_chunks[i], values_minus_two_2_chunks[i], boundaries_minus_two_2, color='gold', label=r'95% expected')

    masses_minus_one_tot_chunks, values_minus_one_tot_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_tot, values_lists=points_all['minus_one'].values_tot)
    masses_minus_one_1_chunks, values_minus_one_1_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_1, values_lists=points_all['minus_one'].values_1)
    masses_minus_one_2_chunks, values_minus_one_2_chunks = self.split_list(masses_lists=points_all['minus_one'].masses_2, values_lists=points_all['minus_one'].values_2)
    for i, chunk in enumerate(masses_minus_one_tot_chunks):
      plt.plot(masses_minus_one_tot_chunks[i], values_minus_one_tot_chunks[i], marker='X', color='forestgreen', label=r'68% expected', linewidth=0)
      plt.plot(masses_minus_one_tot_chunks[i], values_minus_one_tot_chunks[i], color='forestgreen', label=r'68% expected', linewidth=2)

    if not self.do_blind:
      masses_obs_tot_chunks, values_obs_tot_chunks = self.split_list(masses_lists=points_all['obs'].masses_tot, values_lists=points_all['obs'].values_tot)
      for i, chunk in enumerate(masses_obs_tot_chunks):
        plt.plot(masses_obs_tot_chunks[i], values_obs_tot_chunks[i], marker='X', color='black', label='Observed', linewidth=0)
        plt.plot(masses_obs_tot_chunks[i], values_obs_tot_chunks[i], color='black', label='Observed', linewidth=2)

    #for i, masses in enumerate(masses_stuecke_all['plus_two']):
    #  plt.plot(masses_stuecke_all['plus_two'][i], values_stuecke_all['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_plus_two = self.find_boundaries(masses_stuecke_all['plus_two'][i], masses_stuecke_all['central'], values_stuecke_all['plus_two'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['plus_two'][i], values_stuecke_all['plus_two'][i], boundaries_plus_two, color='gold', label=r'95% expected')

    ##for i, masses in enumerate(masses_stuecke_all['plus_one']):
    ##  plt.plot(masses_stuecke_all['plus_one'][i], values_stuecke_all['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    ##  boundaries_plus_one = self.find_boundaries(masses_stuecke_all['plus_one'][i], masses_stuecke_all['central'], values_stuecke_all['plus_one'][i], values_stuecke_all['central'])
    ##  plt.fill_between(masses_stuecke_all['plus_one'][i], values_stuecke_all['plus_one'][i], boundaries_plus_one, color='forestgreen', label=r'68% expected')

    #for i, masses in enumerate(masses_stuecke_all['central']):
    #  plt.plot(masses_stuecke_all['central'][i], values_stuecke_all['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
    #  plt.plot(masses_stuecke_all['central'][i], values_stuecke_all['central'][i], color='red', label='Median expected', linewidth=2)

    #for i, masses in enumerate(masses_stuecke_all['minus_one']):
    #  plt.plot(masses_stuecke_all['minus_one'][i], values_stuecke_all['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_minus_one = self.find_boundaries(masses_stuecke_all['minus_one'][i], masses_stuecke_all['central'], values_stuecke_all['minus_one'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_one'][i], values_stuecke_all['minus_one'][i], boundaries_minus_one, color='forestgreen', label=r'68% expected')

    #for i, masses in enumerate(masses_stuecke_all['minus_two']):
    #  plt.plot(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all['minus_two'][i], masses_stuecke_all['central'], values_stuecke_all['minus_two'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')

    #for i, masses in enumerate(masses_stuecke_all['minus_two']):
    #  plt.plot(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all['minus_two'][i], masses_stuecke_all['central'], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')

    #for i, masses in enumerate(masses_stuecke_all_1['plus_two']):
    #  plt.plot(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_plus_two = self.find_boundaries(masses_stuecke_all_1['plus_two'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  #plt.fill_between(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], values_stuecke_all_1['central'][i], color='gold', label=r'95% expected')
    #  plt.fill_between(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], boundaries_plus_two, color='gold', label=r'95% expected')
    #for i, masses in enumerate(masses_stuecke_all_2['plus_two']):
    #  plt.plot(masses_stuecke_all_2['plus_two'][i], values_stuecke_all_2['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #for i, masses in enumerate(masses_stuecke_all_3['plus_two']):
    #  plt.plot(masses_stuecke_all_3['plus_two'][i], values_stuecke_all_3['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['plus_one']):
    #  plt.plot(masses_stuecke_all_1['plus_one'][i], values_stuecke_all_1['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_plus_one = self.find_boundaries(masses_stuecke_all_1['plus_one'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['plus_one'][i], values_stuecke_all_1['plus_one'][i], boundaries_plus_one, color='forestgreen', label=r'68% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['plus_one']):
    #  plt.plot(masses_stuecke_all_2['plus_one'][i], values_stuecke_all_2['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['plus_one']):
    #  plt.plot(masses_stuecke_all_3['plus_one'][i], values_stuecke_all_3['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['central']):
    #  plt.plot(masses_stuecke_all_1['central'][i], values_stuecke_all_1['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
    #  plt.plot(masses_stuecke_all_1['central'][i], values_stuecke_all_1['central'][i], color='red', label='Median expected', linewidth=2)
    #for i, masses_2 in enumerate(masses_stuecke_all_2['central']):
    #  plt.plot(masses_stuecke_all_2['central'][i], values_stuecke_all_2['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['central']):
    #  plt.plot(masses_stuecke_all_3['central'][i], values_stuecke_all_3['central'][i], color='red', label='Median expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['minus_two']):
    #  plt.plot(masses_stuecke_all_1['minus_two'][i], values_stuecke_all_1['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all_1['minus_two'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['minus_two'][i], values_stuecke_all_1['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['minus_two']):
    #  plt.plot(masses_stuecke_all_2['minus_two'][i], values_stuecke_all_2['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['minus_two']):
    #  plt.plot(masses_stuecke_all_3['minus_two'][i], values_stuecke_all_3['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['minus_one']):
    #  plt.plot(masses_stuecke_all_1['minus_one'][i], values_stuecke_all_1['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_minus_one = self.find_boundaries(masses_stuecke_all_1['minus_one'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['minus_one'][i], values_stuecke_all_1['minus_one'][i], boundaries_minus_one, color='forestgreen', label=r'68% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['minus_one']):
    #  plt.plot(masses_stuecke_all_2['minus_one'][i], values_stuecke_all_2['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['minus_one']):
    #  plt.plot(masses_stuecke_all_3['minus_one'][i], values_stuecke_all_3['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)

    #if not self.do_blind:
    #  for i, masses_1 in enumerate(masses_stuecke_all_1['obs']):
    #    plt.plot(masses_stuecke_all_1['obs'][i], values_stuecke_all_1['obs'][i], color='black', label='Observed', marker='X', linewidth=0)
    #    plt.plot(masses_stuecke_all_1['obs'][i], values_stuecke_all_1['obs'][i], color='black', label='Observed', linewidth=2)
    #  for i, masses_2 in enumerate(masses_stuecke_all_2['obs']):
    #    plt.plot(masses_stuecke_all_2['obs'][i], values_stuecke_all_2['obs'][i], color='black', label='Observed', marker='X', linewidth=0)
    #  for i, masses_3 in enumerate(masses_stuecke_all_3['obs']):
    #    plt.plot(masses_stuecke_all_3['obs'][i], values_stuecke_all_3['obs'][i], color='black', label='Observed', marker='X', linewidth=0)

    #for key in masses_stuecke_all_2.keys():
    #  islands_coordinates = self.get_islands_coordinates(masses_stuecke_all_2[key], masses_stuecke_all_3[key], values_stuecke_all_2[key], values_stuecke_all_3[key])
    #  for island_coordinates in islands_coordinates:
    #    if island_coordinates == None: continue
    #    x = island_coordinates[0]
    #    y = island_coordinates[1]
    #    width = island_coordinates[2]
    #    height = island_coordinates[3]
    #    if key in ['minus_one', 'plus_one']: 
    #      edgecolor = 'gold' #'forestgreen'
    #      facecolor = 'gold' #'forestgreen'
    #      linewidth = 0
    #      fill = True
    #    if key in ['minus_two', 'plus_two']: 
    #      edgecolor = 'gold'
    #      facecolor = 'gold'
    #      linewidth = 0
    #      fill = True
    #    if key in ['obs']: 
    #      edgecolor = 'black'
    #      facecolor = 'white'
    #      linewidth = 2
    #      fill = False
    #    plt.gca().add_patch(Rectangle((x, y), width, height, edgecolor=edgecolor, facecolor=facecolor, fill=fill, linewidth=linewidth)) 

    ##plt.plot(masses_2, values_2, marker='X', color='blue', label='Median expected', linewidth=0)
    ##for i, masses_2 in enumerate(masses_2_stuecke):
    ##  plt.plot(masses_2_stuecke[i], values_2_stuecke[i], color='blue', label='Median expected', linewidth=2)

    ##plt.plot(masses_3, values_3, marker='X', color='darkgreen', label='Median expected', linewidth=0)
    ##for i, masses_3 in enumerate(masses_3_stuecke):
    ##  plt.plot(masses_3_stuecke[i], values_3_stuecke[i], color='darkgreen', label='Median expected', linewidth=2)

    ##plt.plot(masses_missing_1, values_missing_1, marker='o', color='black', label='Median expected', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    #masses_tot = masses_1 + masses_2 + masses_3
    #plt.xlim(min(masses_tot), max(masses_tot))
    plt.xlim(1., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'countour_{}'.format(self.scenario) 
    else:
      name_2d = 'scatter_2dcountour_scenario_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)

  def plot_2dlimit_1(self, masses_stuecke_all, values_stuecke_all, plotDir):
    # have dict masses, values with key central etc.?
    plt.clf()
    f, ax = plt.subplots(figsize=(9, 8))
    y_range_min = 1e-5
    y_range_max = 1e-1
    if not self.do_coupling_scenario:
      self.fe = '0.0'
      self.fu = '1.0'
      self.ft = '0.0'
    if self.fe == '0p5': fe_label = '1/2'
    elif self.fe == '0p3': fe_label = '1/3'
    else: fe_label = self.fe.replace('p', '.')
    if self.fu == '0p5': fu_label = '1/2'
    elif self.fu == '0p3': fu_label = '1/3'
    else: fu_label = self.fu.replace('p', '.')
    if self.ft == '0p5': ft_label = '1/2'
    elif self.ft == '0p3': ft_label = '1/3'
    else: ft_label = self.ft.replace('p', '.')
    coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)
    plt.title('{}, {}'.format(coupling_scenario, self.scenario), loc='right', fontsize=23)

    #print '\n'
    #print masses_stuecke_all_1['plus_two'][1]
    #print values_stuecke_all_1['plus_two'][1]
    #print '\n'
    #print masses_stuecke_all_2['plus_two'][0]
    #print values_stuecke_all_2['plus_two'][0]
    #print '\n'
    #print masses_stuecke_all_3['plus_two'][0]
    #print values_stuecke_all_3['plus_two'][0]

    #[3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65]
    #[3.65, 3.6, 3.55, 3.5]
    #[3.5, 3.55, 3.6, 3.65]
    #[0.00012034689271124365, 0.00012425505596183535, 0.00012387651316755617, 0.0001357589514113421, 0.00012398053056366553, 0.0001621255964163732, 0.00019170782693026195, 0.00019805437170792586, 0.0006330616450559552, 0.0002532928957738135, 0.0003514183210573947]
    #[0.0011121770802836684, 0.0019333104686060747, 0.0014249369293746052, 0.002675377363665742]
    #[0.00551420027591585, 0.010041985319261523, 0.0046999341557779975, 0.004726644830799094]

    #masses_tot = masses_stuecke_all_1['plus_two'][1] + masses_stuecke_all_2['plus_two'][0] + masses_stuecke_all_3['plus_two'][0]
    #values_tot = values_stuecke_all_1['plus_two'][1] + values_stuecke_all_2['plus_two'][0] + values_stuecke_all_3['plus_two'][0]
    #masses_tot = [3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.65, 3.6, 3.55, 3.5, 3.5, 3.55, 3.6, 3.65, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9]
    #values_tot = [0.00012034689271124365, 0.00012425505596183535, 0.00012387651316755617, 0.0001357589514113421, 0.00012398053056366553, 0.0001621255964163732, 0.00019170782693026195, 0.00019805437170792586, 0.0006330616450559552, 0.0002532928957738135, 0.0003514183210573947, 0.0011121770802836684, 0.0019333104686060747, 0.0014249369293746052, 0.002675377363665742, 0.00551420027591585, 0.010041985319261523, 0.0046999341557779975, 0.004726644830799094, 0.010422340034747238, 0.009992252977129756, 0.009730473592131002, 0.00823294826799931, 0.009853815531909501, 0.010269952092311059, 0.010586853974853856, 0.010166377180102508, 0.009314321439711618, 0.00849178867340154, 0.011485044420968676, 0.00921101394854829, 0.012493053272891188, 0.014662758275220093, 0.02478249047358668]
    #boundaries = self.find_boundaries(masses_tot, masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #plt.fill_between(masses_tot, values_tot, boundaries, color='gold', label=r'95% expected')

    # readd masses_central_1 and _2 to fill colour inbetween

    for i, masses in enumerate(masses_stuecke_all['plus_two']):
      plt.plot(masses_stuecke_all['plus_two'][i], values_stuecke_all['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
      boundaries_plus_two = self.find_boundaries(masses_stuecke_all['plus_two'][i], masses_stuecke_all['central'], values_stuecke_all['plus_two'][i], values_stuecke_all['central'])
      plt.fill_between(masses_stuecke_all['plus_two'][i], values_stuecke_all['plus_two'][i], boundaries_plus_two, color='gold', label=r'95% expected')

    #for i, masses in enumerate(masses_stuecke_all['plus_one']):
    #  plt.plot(masses_stuecke_all['plus_one'][i], values_stuecke_all['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_plus_one = self.find_boundaries(masses_stuecke_all['plus_one'][i], masses_stuecke_all['central'], values_stuecke_all['plus_one'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['plus_one'][i], values_stuecke_all['plus_one'][i], boundaries_plus_one, color='forestgreen', label=r'68% expected')

    for i, masses in enumerate(masses_stuecke_all['central']):
      plt.plot(masses_stuecke_all['central'][i], values_stuecke_all['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
      plt.plot(masses_stuecke_all['central'][i], values_stuecke_all['central'][i], color='red', label='Median expected', linewidth=2)

    #for i, masses in enumerate(masses_stuecke_all['minus_one']):
    #  plt.plot(masses_stuecke_all['minus_one'][i], values_stuecke_all['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_minus_one = self.find_boundaries(masses_stuecke_all['minus_one'][i], masses_stuecke_all['central'], values_stuecke_all['minus_one'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_one'][i], values_stuecke_all['minus_one'][i], boundaries_minus_one, color='forestgreen', label=r'68% expected')

    #for i, masses in enumerate(masses_stuecke_all['minus_two']):
    #  plt.plot(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all['minus_two'][i], masses_stuecke_all['central'], values_stuecke_all['minus_two'][i], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')

    #for i, masses in enumerate(masses_stuecke_all['minus_two']):
    #  plt.plot(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all['minus_two'][i], masses_stuecke_all['central'], values_stuecke_all['central'])
    #  plt.fill_between(masses_stuecke_all['minus_two'][i], values_stuecke_all['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')

    #for i, masses in enumerate(masses_stuecke_all_1['plus_two']):
    #  plt.plot(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_plus_two = self.find_boundaries(masses_stuecke_all_1['plus_two'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  #plt.fill_between(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], values_stuecke_all_1['central'][i], color='gold', label=r'95% expected')
    #  plt.fill_between(masses_stuecke_all_1['plus_two'][i], values_stuecke_all_1['plus_two'][i], boundaries_plus_two, color='gold', label=r'95% expected')
    #for i, masses in enumerate(masses_stuecke_all_2['plus_two']):
    #  plt.plot(masses_stuecke_all_2['plus_two'][i], values_stuecke_all_2['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #for i, masses in enumerate(masses_stuecke_all_3['plus_two']):
    #  plt.plot(masses_stuecke_all_3['plus_two'][i], values_stuecke_all_3['plus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['plus_one']):
    #  plt.plot(masses_stuecke_all_1['plus_one'][i], values_stuecke_all_1['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_plus_one = self.find_boundaries(masses_stuecke_all_1['plus_one'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['plus_one'][i], values_stuecke_all_1['plus_one'][i], boundaries_plus_one, color='forestgreen', label=r'68% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['plus_one']):
    #  plt.plot(masses_stuecke_all_2['plus_one'][i], values_stuecke_all_2['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['plus_one']):
    #  plt.plot(masses_stuecke_all_3['plus_one'][i], values_stuecke_all_3['plus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['central']):
    #  plt.plot(masses_stuecke_all_1['central'][i], values_stuecke_all_1['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
    #  plt.plot(masses_stuecke_all_1['central'][i], values_stuecke_all_1['central'][i], color='red', label='Median expected', linewidth=2)
    #for i, masses_2 in enumerate(masses_stuecke_all_2['central']):
    #  plt.plot(masses_stuecke_all_2['central'][i], values_stuecke_all_2['central'][i], color='red', label='Median expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['central']):
    #  plt.plot(masses_stuecke_all_3['central'][i], values_stuecke_all_3['central'][i], color='red', label='Median expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['minus_two']):
    #  plt.plot(masses_stuecke_all_1['minus_two'][i], values_stuecke_all_1['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #  boundaries_minus_two = self.find_boundaries(masses_stuecke_all_1['minus_two'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['minus_two'][i], values_stuecke_all_1['minus_two'][i], boundaries_minus_two, color='gold', label=r'95% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['minus_two']):
    #  plt.plot(masses_stuecke_all_2['minus_two'][i], values_stuecke_all_2['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['minus_two']):
    #  plt.plot(masses_stuecke_all_3['minus_two'][i], values_stuecke_all_3['minus_two'][i], color='gold', label=r'95% expected', marker='X', linewidth=0)

    #for i, masses_1 in enumerate(masses_stuecke_all_1['minus_one']):
    #  plt.plot(masses_stuecke_all_1['minus_one'][i], values_stuecke_all_1['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #  boundaries_minus_one = self.find_boundaries(masses_stuecke_all_1['minus_one'][i], masses_stuecke_all_1['central'], values_stuecke_all_1['central'])
    #  plt.fill_between(masses_stuecke_all_1['minus_one'][i], values_stuecke_all_1['minus_one'][i], boundaries_minus_one, color='forestgreen', label=r'68% expected')
    #for i, masses_2 in enumerate(masses_stuecke_all_2['minus_one']):
    #  plt.plot(masses_stuecke_all_2['minus_one'][i], values_stuecke_all_2['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)
    #for i, masses_3 in enumerate(masses_stuecke_all_3['minus_one']):
    #  plt.plot(masses_stuecke_all_3['minus_one'][i], values_stuecke_all_3['minus_one'][i], color='forestgreen', label=r'68% expected', marker='X', linewidth=0)

    #if not self.do_blind:
    #  for i, masses_1 in enumerate(masses_stuecke_all_1['obs']):
    #    plt.plot(masses_stuecke_all_1['obs'][i], values_stuecke_all_1['obs'][i], color='black', label='Observed', marker='X', linewidth=0)
    #    plt.plot(masses_stuecke_all_1['obs'][i], values_stuecke_all_1['obs'][i], color='black', label='Observed', linewidth=2)
    #  for i, masses_2 in enumerate(masses_stuecke_all_2['obs']):
    #    plt.plot(masses_stuecke_all_2['obs'][i], values_stuecke_all_2['obs'][i], color='black', label='Observed', marker='X', linewidth=0)
    #  for i, masses_3 in enumerate(masses_stuecke_all_3['obs']):
    #    plt.plot(masses_stuecke_all_3['obs'][i], values_stuecke_all_3['obs'][i], color='black', label='Observed', marker='X', linewidth=0)

    #for key in masses_stuecke_all_2.keys():
    #  islands_coordinates = self.get_islands_coordinates(masses_stuecke_all_2[key], masses_stuecke_all_3[key], values_stuecke_all_2[key], values_stuecke_all_3[key])
    #  for island_coordinates in islands_coordinates:
    #    if island_coordinates == None: continue
    #    x = island_coordinates[0]
    #    y = island_coordinates[1]
    #    width = island_coordinates[2]
    #    height = island_coordinates[3]
    #    if key in ['minus_one', 'plus_one']: 
    #      edgecolor = 'gold' #'forestgreen'
    #      facecolor = 'gold' #'forestgreen'
    #      linewidth = 0
    #      fill = True
    #    if key in ['minus_two', 'plus_two']: 
    #      edgecolor = 'gold'
    #      facecolor = 'gold'
    #      linewidth = 0
    #      fill = True
    #    if key in ['obs']: 
    #      edgecolor = 'black'
    #      facecolor = 'white'
    #      linewidth = 2
    #      fill = False
    #    plt.gca().add_patch(Rectangle((x, y), width, height, edgecolor=edgecolor, facecolor=facecolor, fill=fill, linewidth=linewidth)) 

    ##plt.plot(masses_2, values_2, marker='X', color='blue', label='Median expected', linewidth=0)
    ##for i, masses_2 in enumerate(masses_2_stuecke):
    ##  plt.plot(masses_2_stuecke[i], values_2_stuecke[i], color='blue', label='Median expected', linewidth=2)

    ##plt.plot(masses_3, values_3, marker='X', color='darkgreen', label='Median expected', linewidth=0)
    ##for i, masses_3 in enumerate(masses_3_stuecke):
    ##  plt.plot(masses_3_stuecke[i], values_3_stuecke[i], color='darkgreen', label='Median expected', linewidth=2)

    ##plt.plot(masses_missing_1, values_missing_1, marker='o', color='black', label='Median expected', linewidth=0)

    plt.ylabel(r'$|V|^2$', fontsize=23)
    plt.yticks(fontsize=17)
    plt.ylim(y_range_min, y_range_max)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    #masses_tot = masses_1 + masses_2 + masses_3
    #plt.xlim(min(masses_tot), max(masses_tot))
    plt.xlim(1., 6.)
    plt.xticks(fontsize=17)
    plt.yscale('log')
    plt.xscale('linear')
    plt.grid(True)
    if not self.do_coupling_scenario:
      name_2d = 'countour_{}'.format(self.scenario) 
    else:
      name_2d = 'scatter_2dcountour_scenario_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    print '--> {}/{}.png created'.format(plotDir, name_2d)


  def process(self):
    #lumi =  '5.3 fb'+r'$^{-1}$'+' projected to 41.6 fb'+r'$^{-1}$'
    #lumi =  '40.0 fb'+r'$^{-1}$'
    lumi =  '41.6 fb'+r'$^{-1}$'

    # get the files 
    if not self.do_coupling_scenario:
      pathToResults = '{}/outputs/{}/limits/{}/results/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    else:
      pathToResults = '{}/outputs/{}/limits/{}/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.fe, self.fu, self.ft) 
      #FIXME
      #pathToResults = '{}/outputs/{}/limits/{}/electron/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.fe, self.fu, self.ft) 

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

      #if float(mass) < 2. or float(mass) > 5.: continue
      #if float(mass) > 3.: continue

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
        val_coupling = self.getCouplingTarget(mass=mass, coupling=coupling)
      
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
            # remove files that didn't finish without an error
            #if 'ERROR' in line: 
            #  print 'there was an error'
            #  break

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

      if self.plot_1D:
        print '-> will plot 1D limit for mass {}'.format(mass)
        
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
          plt.title(r'%s N, m = %s GeV, %s' %(self.scenario, mass, coupling_scenario))
        else:
          if str(mass) == '1': 
            mass = '1.0'
          if str(mass) == '2': 
            mass = '2.0'
          if str(mass) == '3': 
            mass = '3.0'
          plt.title('{} N, m = {} GeV'.format(self.scenario, mass))
        plt.legend()
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
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
      x_minus_two_list = self.get_intersection_list(v2s, minus_two)
      x_minus_one = self.get_intersection(v2s, minus_one)
      x_minus_one_list = self.get_intersection_list(v2s, minus_one)
      x_central = self.get_intersection(v2s, central) #TODO remove
      x_central_list = self.get_intersection_list(v2s, central)
      print self.get_intersection_list(v2s, central)
      x_plus_one = self.get_intersection(v2s, plus_one)
      x_plus_one_list = self.get_intersection_list(v2s, plus_one)
      x_plus_two = self.get_intersection(v2s, plus_two)
      x_plus_two_list = self.get_intersection_list(v2s, plus_two)
      #print 'central {}'.format(x_central)

      if not self.do_blind:
        x_obs = self.get_intersection(v2s, obs)
        x_obs_list = self.get_intersection_list(v2s, obs)

      #if float(mass) < 2.8 or float(mass) > 3.5:
      if x_plus_one == self.no_exclusion_value:
        print '\nWARNING - could not find crossing for +1sigma'
        crossings = np.linspace(1, 3, 50)
        for crossing in crossings:
          x_central_tmp = self.get_intersection(v2s, central, crossing)
          x_plus_one_tmp = self.get_intersection(v2s, plus_one, crossing)
          if x_plus_one_tmp != self.no_exclusion_value: 
            break
        x_diff_plus_one = x_plus_one_tmp - x_central_tmp
        x_plus_one = x_central + x_diff_plus_one

      if x_plus_two == self.no_exclusion_value:
        print '\nWARNING - could not find crossing for +2sigma'
        crossings = np.linspace(1, 3, 50)
        for crossing in crossings:
          x_central_tmp = self.get_intersection(v2s, central, crossing)
          x_plus_two_tmp = self.get_intersection(v2s, plus_two, crossing)
          if x_plus_two_tmp != self.no_exclusion_value: 
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

      #if not self.do_blind:
      #  x_obs = self.get_intersection(v2s, obs)

      limits2D[mass]['exp_minus_two'] = x_minus_two
      limits2D[mass]['exp_minus_two_list'] = x_minus_two_list
      limits2D[mass]['exp_minus_one'] = x_minus_one
      limits2D[mass]['exp_minus_one_list'] = x_minus_one_list
      #limits2D[mass]['exp_central'  ] = x_central  
      limits2D[mass]['exp_central_list'  ] = x_central_list  
      limits2D[mass]['exp_plus_one' ] = x_plus_one 
      limits2D[mass]['exp_plus_one_list' ] = x_plus_one_list 
      limits2D[mass]['exp_plus_two' ] = x_plus_two 
      limits2D[mass]['exp_plus_two_list' ] = x_plus_two_list 
      if not self.do_blind:
        limits2D[mass]['obs'] = x_obs 
        limits2D[mass]['obs_list'] = x_obs_list 

    #mass_minus_one_turning = self.get_turning_point(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_one_list')

    print '\nwill plot 2D limits' 
    with open('{}/results.pck'.format(plotDir), 'w') as ff:
        pickle.dump(limits2D, ff)

    #masses_obs = []
    #masses_obs_2 = []
    #masses_obs_3 = []
    #masses_obs_missing = []
    #masses_central = []
    #masses_central_2 = []
    #masses_central_3 = []
    #masses_central_missing = []
    #masses_minus_one_ = []
    #masses_minus_one_2 = []
    #masses_minus_one_3 = []
    #masses_minus_one_missing = []

    #masses_plus_one_ = []
    #masses_plus_two_ = []
    #masses_minus_two_ = []

    #boundary_plus_two = []
    #boundary_plus_one = []
    #boundary_minus_two = []
    #boundary_minus_one = []

    #minus_two = []
    #minus_one = []
    #minus_one_2 = []
    #minus_one_3 = []
    #minus_one_missing = []
    #central   = []
    #central_2   = []
    #central_3   = []
    #central_missing   = []
    #central_upper   = []
    ##central_3   = []
    #plus_one  = []
    #plus_two  = []
    #obs       = []
    #obs_2       = []
    #obs_3       = []
    #obs_missing       = []

    masses_minus_two_1, masses_minus_two_2, masses_minus_two_3, masses_minus_two_missing, minus_two_1, minus_two_2, minus_two_3, minus_two_missing = self.get_lists(limits2D, 'exp_minus_two_list')
    masses_minus_one_1, masses_minus_one_2, masses_minus_one_3, masses_minus_one_missing, minus_one_1, minus_one_2, minus_one_3, minus_one_missing = self.get_lists(limits2D, 'exp_minus_one_list')
    masses_central_1, masses_central_2, masses_central_3, masses_central_missing, central_1, central_2, central_3, central_missing = self.get_lists(limits2D, 'exp_central_list')
    masses_plus_one_1, masses_plus_one_2, masses_plus_one_3, masses_plus_one_missing, plus_one_1, plus_one_2, plus_one_3, plus_one_missing = self.get_lists(limits2D, 'exp_plus_one_list')
    masses_plus_two_1, masses_plus_two_2, masses_plus_two_3, masses_plus_two_missing, plus_two_1, plus_two_2, plus_two_3, plus_two_missing = self.get_lists(limits2D, 'exp_plus_two_list')
    if not self.do_blind:
      masses_obs_1, masses_obs_2, masses_obs_3, masses_obs_missing, obs_1, obs_2, obs_3, obs_missing = self.get_lists(limits2D, 'obs_list')

    # plot scatter figures
    if self.plot_scatter_figure:
      self.plot_scatter(masses_1=masses_minus_two_1, masses_2=masses_minus_two_2, masses_3=masses_minus_two_3, masses_missing_1=masses_minus_two_missing, values_1=minus_two_1, values_2=minus_two_2, values_3=minus_two_3, values_missing_1=minus_two_missing, label='minus_two', plotDir=plotDir)
      self.plot_scatter(masses_1=masses_minus_one_1, masses_2=masses_minus_one_2, masses_3=masses_minus_one_3, masses_missing_1=masses_minus_one_missing, values_1=minus_one_1, values_2=minus_one_2, values_3=minus_one_3, values_missing_1=minus_one_missing, label='minus_one', plotDir=plotDir)
      self.plot_scatter(masses_1=masses_central_1, masses_2=masses_central_2, masses_3=masses_central_3, masses_missing_1=masses_central_missing, values_1=central_1, values_2=central_2, values_3=central_3, values_missing_1=central_missing, label='central', plotDir=plotDir)
      self.plot_scatter(masses_1=masses_plus_one_1, masses_2=masses_plus_one_2, masses_3=masses_plus_one_3, masses_missing_1=masses_plus_one_missing, values_1=plus_one_1, values_2=plus_one_2, values_3=plus_one_3, values_missing_1=plus_one_missing, label='plus_one', plotDir=plotDir)
      self.plot_scatter(masses_1=masses_plus_two_1, masses_2=masses_plus_two_2, masses_3=masses_plus_two_3, masses_missing_1=masses_plus_two_missing, values_1=plus_two_1, values_2=plus_two_2, values_3=plus_two_3, values_missing_1=plus_two_missing, label='plus_two', plotDir=plotDir)
      if not self.do_blind:
        self.plot_scatter(masses_1=masses_obs_1, masses_2=masses_obs_2, masses_3=masses_obs_3, masses_missing_1=masses_obs_missing, values_1=obs_1, values_2=obs_2, values_3=obs_3, values_missing_1=obs_missing, label='obs', plotDir=plotDir)

    # concatenate lists
    #masses_minus_two_tot, values_minus_two_tot = self.concatenate_lists(masses_1=masses_minus_two_1, masses_2=masses_minus_two_2, masses_3=masses_minus_two_3, values_1=minus_two_1, values_2=minus_two_2, values_3=minus_two_3)
    #masses_minus_one_tot, values_minus_one_tot = self.concatenate_lists(masses_1=masses_minus_one_1, masses_2=masses_minus_one_2, masses_3=masses_minus_one_3, values_1=minus_one_1, values_2=minus_one_2, values_3=minus_one_3)

    #masses_central_tot, values_central_tot = self.concatenate_lists(masses_1=masses_central_1, masses_2=masses_central_2, masses_3=masses_central_3, values_1=central_1, values_2=central_2, values_3=central_3)
    #masses_plus_one_tot, values_plus_one_tot = self.concatenate_lists(masses_1=masses_plus_one_1, masses_2=masses_plus_one_2, masses_3=masses_plus_one_3, values_1=plus_one_1, values_2=plus_one_2, values_3=plus_one_3)
    #masses_plus_two_tot, values_plus_two_tot = self.concatenate_lists(masses_1=masses_plus_two_1, masses_2=masses_plus_two_2, masses_3=masses_plus_two_3, values_1=plus_two_1, values_2=plus_two_2, values_3=plus_two_3)
    #if not self.do_blind:
    #  masses_obs_tot, values_obs_tot = self.concatenate_lists(masses_1=masses_obs_1, masses_2=masses_obs_2, masses_3=masses_obs_3, values_1=obs_1, values_2=obs_2, values_3=obs_3)

    points_minus_two = self.concatenate_lists(masses_1=masses_minus_two_1, masses_2=masses_minus_two_2, masses_3=masses_minus_two_3, values_1=minus_two_1, values_2=minus_two_2, values_3=minus_two_3)
    points_minus_one = self.concatenate_lists(masses_1=masses_minus_one_1, masses_2=masses_minus_one_2, masses_3=masses_minus_one_3, values_1=minus_one_1, values_2=minus_one_2, values_3=minus_one_3)
    points_central = self.concatenate_lists(masses_1=masses_central_1, masses_2=masses_central_2, masses_3=masses_central_3, values_1=central_1, values_2=central_2, values_3=central_3)
    points_plus_one  = self.concatenate_lists(masses_1=masses_plus_one_1, masses_2=masses_plus_one_2, masses_3=masses_plus_one_3, values_1=plus_one_1, values_2=plus_one_2, values_3=plus_one_3)
    points_plus_two = self.concatenate_lists(masses_1=masses_plus_two_1, masses_2=masses_plus_two_2, masses_3=masses_plus_two_3, values_1=plus_two_1, values_2=plus_two_2, values_3=plus_two_3)
    if not self.do_blind:
      points_obs = self.concatenate_lists(masses_1=masses_obs_1, masses_2=masses_obs_2, masses_3=masses_obs_3, values_1=obs_1, values_2=obs_2, values_3=obs_3)

    #masses_minus_two_stuecke, values_minus_two_stuecke = self.split_list(masses_lists=masses_minus_two_tot, values_lists=values_minus_two_tot)
    #masses_minus_one_stuecke, values_minus_one_stuecke = self.split_list(masses_lists=masses_minus_one_tot, values_lists=values_minus_one_tot)
    #masses_central_stuecke, values_central_stuecke = self.split_list(masses_lists=masses_central_tot, values_lists=values_central_tot)
    #masses_plus_one_stuecke, values_plus_one_stuecke = self.split_list(masses_lists=masses_plus_one_tot, values_lists=values_plus_one_tot)
    #masses_plus_two_stuecke, values_plus_two_stuecke = self.split_list(masses_lists=masses_plus_two_tot, values_lists=values_plus_two_tot)
    #masses_obs_stuecke, values_obs_stuecke = self.split_list(masses_lists=masses_obs_tot, values_lists=values_obs_tot)
      
    # plot contour
    #if self.plot_countour_figure:
    #  self.plot_countour(masses_stuecke=masses_minus_two_stuecke, masses_missing=masses_minus_two_missing, values_stuecke=values_minus_two_stuecke, values_missing=minus_two_missing, label='minus_two', plotDir=plotDir)
    #  self.plot_countour(masses_stuecke=masses_minus_one_stuecke, masses_missing=masses_minus_one_missing, values_stuecke=values_minus_one_stuecke, values_missing=minus_one_missing, label='minus_one', plotDir=plotDir)
    #  self.plot_countour(masses_stuecke=masses_central_stuecke, masses_missing=masses_central_missing, values_stuecke=values_central_stuecke, values_missing=central_missing, label='central', plotDir=plotDir)
    #  self.plot_countour(masses_stuecke=masses_plus_one_stuecke, masses_missing=masses_plus_one_missing, values_stuecke=values_plus_one_stuecke, values_missing=plus_one_missing, label='plus_one', plotDir=plotDir)
    #  self.plot_countour(masses_stuecke=masses_plus_two_stuecke, masses_missing=masses_plus_two_missing, values_stuecke=values_plus_two_stuecke, values_missing=plus_two_missing, label='plus_two', plotDir=plotDir)
    #  if not self.do_blind:
    #    self.plot_countour(masses_stuecke=masses_obs_stuecke, masses_missing=masses_obs_missing, values_stuecke=values_obs_stuecke, values_missing=obs_missing, label='obs', plotDir=plotDir)

    if self.plot_countour_figure:
      self.plot_countour(points=points_minus_two, masses_missing=masses_minus_two_missing, values_missing=minus_two_missing, label='minus_two', plotDir=plotDir)
      self.plot_countour(points=points_minus_one, masses_missing=masses_minus_one_missing, values_missing=minus_one_missing, label='minus_one', plotDir=plotDir)
      self.plot_countour(points=points_central, masses_missing=masses_central_missing, values_missing=central_missing, label='central', plotDir=plotDir)
      self.plot_countour(points=points_plus_one, masses_missing=masses_plus_one_missing, values_missing=plus_one_missing, label='plus_one', plotDir=plotDir)
      self.plot_countour(points=points_plus_two, masses_missing=masses_plus_two_missing, values_missing=plus_two_missing, label='plus_two', plotDir=plotDir)
      if not self.do_blind:
        self.plot_countour(points=points_obs, masses_missing=masses_obs_missing, values_missing=obs_missing, label='obs', plotDir=plotDir)

    #masses_stuecke_all = OrderedDict()
    #masses_stuecke_all['minus_two'] = masses_minus_two_stuecke
    #masses_stuecke_all['minus_one'] = masses_minus_one_stuecke
    #masses_stuecke_all['central'] = masses_central_stuecke
    #masses_stuecke_all['plus_one'] = masses_plus_one_stuecke
    #masses_stuecke_all['plus_two'] = masses_plus_two_stuecke
    #if not self.do_blind:
    #  masses_stuecke_all['obs'] = masses_obs_stuecke
    #  
    #values_stuecke_all = OrderedDict()
    #values_stuecke_all['minus_two'] = values_minus_two_stuecke
    #values_stuecke_all['minus_one'] = values_minus_one_stuecke
    #values_stuecke_all['central'] = values_central_stuecke
    #values_stuecke_all['plus_one'] = values_plus_one_stuecke
    #values_stuecke_all['plus_two'] = values_plus_two_stuecke
    #if not self.do_blind:
    #  values_stuecke_all['obs'] = values_obs_stuecke

    points_all = OrderedDict()
    points_all['minus_two'] = points_minus_two
    points_all['minus_one'] = points_minus_one
    points_all['central'] = points_central
    points_all['plus_one'] = points_plus_one
    points_all['plus_two'] = points_plus_two
    if not self.do_blind:
      points_all['obs'] = points_obs

    #self.plot_2dlimit(points_all=points_all, plotDir=plotDir)
    self.plot_scatter_all(points_all=points_all, plotDir=plotDir)

    #self.plot_2dlimit(masses_stuecke_all=masses_stuecke_all, values_stuecke_all=values_stuecke_all, plotDir=plotDir)
    #self.plot_2dlimit(masses_stuecke_all_1=masses_stuecke_all_1, masses_stuecke_all_2=masses_stuecke_all_2, masses_stuecke_all_3=masses_stuecke_all_3, values_stuecke_all_1=values_stuecke_all_1, values_stuecke_all_2=values_stuecke_all_2, values_stuecke_all_3=values_stuecke_all_3, plotDir=plotDir)

    #masses_minus_two_1_stuecke, values_minus_two_1_stuecke = self.split_list(masses=masses_minus_two_1, values=minus_two_1)
    #masses_minus_two_2_stuecke, values_minus_two_2_stuecke = self.split_list(masses=masses_minus_two_2, values=minus_two_2)
    #masses_minus_two_3_stuecke, values_minus_two_3_stuecke = self.split_list(masses=masses_minus_two_3, values=minus_two_3)

    #masses_minus_one_1_stuecke, values_minus_one_1_stuecke = self.split_list(masses=masses_minus_one_1, values=minus_one_1)
    #masses_minus_one_2_stuecke, values_minus_one_2_stuecke = self.split_list(masses=masses_minus_one_2, values=minus_one_2)
    #masses_minus_one_3_stuecke, values_minus_one_3_stuecke = self.split_list(masses=masses_minus_one_3, values=minus_one_3)

    #masses_central_1_stuecke, values_central_1_stuecke = self.split_list(masses=masses_central_1, values=central_1)
    #masses_central_2_stuecke, values_central_2_stuecke = self.split_list(masses=masses_central_2, values=central_2)
    #masses_central_3_stuecke, values_central_3_stuecke = self.split_list(masses=masses_central_3, values=central_3)

    #masses_plus_one_1_stuecke, values_plus_one_1_stuecke = self.split_list(masses=masses_plus_one_1, values=plus_one_1)
    #masses_plus_one_2_stuecke, values_plus_one_2_stuecke = self.split_list(masses=masses_plus_one_2, values=plus_one_2)
    #masses_plus_one_3_stuecke, values_plus_one_3_stuecke = self.split_list(masses=masses_plus_one_3, values=plus_one_3)

    #masses_plus_two_1_stuecke, values_plus_two_1_stuecke = self.split_list(masses=masses_plus_two_1, values=plus_two_1)
    #masses_plus_two_2_stuecke, values_plus_two_2_stuecke = self.split_list(masses=masses_plus_two_2, values=plus_two_2)
    #masses_plus_two_3_stuecke, values_plus_two_3_stuecke = self.split_list(masses=masses_plus_two_3, values=plus_two_3)

    #if not self.do_blind:
    #  masses_obs_1_stuecke, values_obs_1_stuecke = self.split_list(masses=masses_obs_1, values=obs_1)
    #  masses_obs_2_stuecke, values_obs_2_stuecke = self.split_list(masses=masses_obs_2, values=obs_2)
    #  masses_obs_3_stuecke, values_obs_3_stuecke = self.split_list(masses=masses_obs_3, values=obs_3)

    #masses_stuecke_all_1 = OrderedDict()
    #masses_stuecke_all_1['minus_two'] = masses_minus_two_1_stuecke
    #masses_stuecke_all_1['minus_one'] = masses_minus_one_1_stuecke
    #masses_stuecke_all_1['central'] = masses_central_1_stuecke
    #masses_stuecke_all_1['plus_one'] = masses_plus_one_1_stuecke
    #masses_stuecke_all_1['plus_two'] = masses_plus_two_1_stuecke
    #if not self.do_blind:
    #  masses_stuecke_all_1['obs'] = masses_obs_1_stuecke

    #masses_stuecke_all_2 = OrderedDict()
    #masses_stuecke_all_2['minus_two'] = masses_minus_two_2_stuecke
    #masses_stuecke_all_2['minus_one'] = masses_minus_one_2_stuecke
    #masses_stuecke_all_2['central'] = masses_central_2_stuecke
    #masses_stuecke_all_2['plus_one'] = masses_plus_one_2_stuecke
    #masses_stuecke_all_2['plus_two'] = masses_plus_two_2_stuecke
    #if not self.do_blind:
    #  masses_stuecke_all_2['obs'] = masses_obs_2_stuecke

    #masses_stuecke_all_3 = OrderedDict()
    #masses_stuecke_all_3['minus_two'] = masses_minus_two_3_stuecke
    #masses_stuecke_all_3['minus_one'] = masses_minus_one_3_stuecke
    #masses_stuecke_all_3['central'] = masses_central_3_stuecke
    #masses_stuecke_all_3['plus_one'] = masses_plus_one_3_stuecke
    #masses_stuecke_all_3['plus_two'] = masses_plus_two_3_stuecke
    #if not self.do_blind:
    #  masses_stuecke_all_3['obs'] = masses_obs_3_stuecke

    #values_stuecke_all_1 = OrderedDict()
    #values_stuecke_all_1['minus_two'] = values_minus_two_1_stuecke
    #values_stuecke_all_1['minus_one'] = values_minus_one_1_stuecke
    #values_stuecke_all_1['central'] = values_central_1_stuecke
    #values_stuecke_all_1['plus_one'] = values_plus_one_1_stuecke
    #values_stuecke_all_1['plus_two'] = values_plus_two_1_stuecke
    #if not self.do_blind:
    #  values_stuecke_all_1['obs'] = values_obs_1_stuecke

    #values_stuecke_all_2 = OrderedDict()
    #values_stuecke_all_2['minus_two'] = values_minus_two_2_stuecke
    #values_stuecke_all_2['minus_one'] = values_minus_one_2_stuecke
    #values_stuecke_all_2['central'] = values_central_2_stuecke
    #values_stuecke_all_2['plus_one'] = values_plus_one_2_stuecke
    #values_stuecke_all_2['plus_two'] = values_plus_two_2_stuecke
    #if not self.do_blind:
    #  values_stuecke_all_2['obs'] = values_obs_2_stuecke

    #values_stuecke_all_3 = OrderedDict()
    #values_stuecke_all_3['minus_two'] = values_minus_two_3_stuecke
    #values_stuecke_all_3['minus_one'] = values_minus_one_3_stuecke
    #values_stuecke_all_3['central'] = values_central_3_stuecke
    #values_stuecke_all_3['plus_one'] = values_plus_one_3_stuecke
    #values_stuecke_all_3['plus_two'] = values_plus_two_3_stuecke
    #if not self.do_blind:
    #  values_stuecke_all_3['obs'] = values_obs_3_stuecke

    #self.plot_2dlimit(masses_stuecke_all_1=masses_stuecke_all_1, masses_stuecke_all_2=masses_stuecke_all_2, masses_stuecke_all_3=masses_stuecke_all_3, values_stuecke_all_1=values_stuecke_all_1, values_stuecke_all_2=values_stuecke_all_2, values_stuecke_all_3=values_stuecke_all_3, plotDir=plotDir)



    ## go through the different mass points first left to right to catch the lower exclusion bound
    ## then right to left to catch the upper exclusion bound

    #masses_central, central = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_central_list')
    #masses_central_upper, central_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_central_list')
    #masses_central_all, central_all = self.get_total_limit(masses_lower=masses_central, masses_upper=masses_central_upper, limits_lower=central, limits_upper=central_upper)

    ##masses_central_all = masses_central + masses_central_upper
    ##central_all = central + central_upper

    ##mass_central_turning = self.get_turning_point(masses=limits2D.keys(), values=limits2D, quantity='exp_central_list')
    ##for mass in sorted(limits2D.keys(), key=self.sortList):
    ##  if float(mass) > mass_central_turning: continue 
    ##  #central_value = limits2D[mass]['exp_central']
    ##  if len(limits2D[mass]['exp_central_list']) > 0:
    ##  #central_value = limits2D[mass]['exp_central_list'][0]
    ##  #if central_value != -99.:
    ##    #central.append(central_value)
    ##    central.append(limits2D[mass]['exp_central_list'][0])
    ##    masses_central.append(float(mass))

    ##  #if len(limits2D[mass]['exp_central_list']) > 1:
    ##  #  central_2.append(limits2D[mass]['exp_central_list'][1])
    ##  #  masses_central_2.append(float(mass))

    ##  #if len(limits2D[mass]['exp_central_list']) > 2:
    ##  #  central_3.append(limits2D[mass]['exp_central_list'][2])
    ##  #  masses_central_3.append(float(mass))

    #if not self.do_blind:
    #  masses_obs, obs = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='obs_list')
    #  masses_obs_upper, obs_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='obs_list')
    #  masses_obs_missing, obs_missing = self.get_missing_points(masses=limits2D.keys(), values=limits2D, quantity='obs_list')
    #  masses_obs_all, obs_all = self.get_total_limit(masses_lower=masses_obs, masses_upper=masses_obs_upper, limits_lower=obs, limits_upper=obs_upper)

    #masses_minus_one, minus_one = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_one_list')
    #masses_minus_one_upper, minus_one_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_one_list')
    #masses_minus_one_missing, minus_one_missing = self.get_missing_points(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_one_list')
    #masses_minus_one_all, minus_one_all = self.get_total_limit(masses_lower=masses_minus_one, masses_upper=masses_minus_one_upper, limits_lower=minus_one, limits_upper=minus_one_upper)

    #masses_plus_one, plus_one = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_one_list')
    #masses_plus_one_upper, plus_one_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_one_list')
    #masses_plus_one_missing, plus_one_missing = self.get_missing_points(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_one_list')
    #masses_plus_one_all, plus_one_all = self.get_total_limit(masses_lower=masses_plus_one, masses_upper=masses_plus_one_upper, limits_lower=plus_one, limits_upper=plus_one_upper)

    ##masses_plus_one_3 = []
    ##plus_one_3 = []
    ##for mass in sorted(limits2D.keys(), key=self.sortList):
    ##  if len(limits2D[mass]['exp_plus_one_list']) > 2:
    ##    plus_one_3.append(limits2D[mass]['exp_plus_one_list'][2])
    ##    masses_plus_one_3.append(float(mass))

    #masses_minus_two, minus_two = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_two_list')
    #masses_minus_two_upper, minus_two_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_two_list')
    #masses_minus_two_missing, minus_two_missing = self.get_missing_points(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_two_list')
    #masses_minus_two_all, minus_two_all = self.get_total_limit(masses_lower=masses_minus_two, masses_upper=masses_minus_two_upper, limits_lower=minus_two, limits_upper=minus_two_upper)

    ##masses_minus_two_3 = []
    ##minus_two_3 = []
    ##for mass in sorted(limits2D.keys(), key=self.sortList):
    ##  if len(limits2D[mass]['exp_minus_two_list']) > 2:
    ##    minus_two_3.append(limits2D[mass]['exp_minus_two_list'][2])
    ##    masses_minus_two_3.append(float(mass))

    #masses_plus_two, plus_two = self.get_lower_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_two_list')
    #masses_plus_two_upper, plus_two_upper = self.get_upper_limit(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_two_list')
    #masses_plus_two_missing, plus_two_missing = self.get_missing_points(masses=limits2D.keys(), values=limits2D, quantity='exp_plus_two_list')
    #masses_plus_two_all, plus_two_all = self.get_total_limit(masses_lower=masses_plus_two, masses_upper=masses_plus_two_upper, limits_lower=plus_two, limits_upper=plus_two_upper)

    #masses_plus_two_3 = []
    #plus_two_3 = []
    #for mass in sorted(limits2D.keys(), key=self.sortList):
    #  if len(limits2D[mass]['exp_plus_two_list']) > 2:
    #    plus_two_3.append(limits2D[mass]['exp_plus_two_list'][2])
    #    masses_plus_two_3.append(float(mass))

    ##mass_minus_one_turning = self.get_turning_point(masses=limits2D.keys(), values=limits2D, quantity='exp_minus_one_list')
    ##for mass in sorted(limits2D.keys(), key=self.sortList):
    ##  if float(mass) > mass_minus_one_turning: continue
    ##  if len(limits2D[mass]['exp_minus_one_list']) > 0:
    ##    minus_one.append(limits2D[mass]['exp_minus_one_list'][0])
    ##    masses_minus_one_sigma.append(float(mass))

    #  #if len(limits2D[mass]['exp_minus_one_list']) > 1:
    #  #  minus_one_2.append(limits2D[mass]['exp_minus_one_list'][1])
    #  #  masses_minus_one_sigma_2.append(float(mass))

    #  #if len(limits2D[mass]['exp_minus_one_list']) > 2:
    #  #  minus_one_3.append(limits2D[mass]['exp_minus_one_list'][2])
    #  #  masses_minus_one_sigma_3.append(float(mass))

    #  #if len(limits2D[mass]['exp_minus_one_list']) == 0:
    #  #  minus_one_missing.append(2e-5)
    #  #  masses_minus_one_sigma_missing.append(float(mass))

    #    #if limits2D[mass]['exp_plus_one' ] != self.no_exclusion_value and limits2D[mass]['exp_plus_two' ] != self.no_exclusion_value:
    #    #  plus_two_value = limits2D[mass]['exp_plus_two' ]
    #    #  plus_one_value = limits2D[mass]['exp_plus_one' ]

    #    #  minus_two_value = limits2D[mass]['exp_minus_two']
    #    #  minus_one_value = limits2D[mass]['exp_minus_one']
    #    
    #      #central_value = limits2D[mass]['exp_central_list'][0] #FIXME

    #      #if plus_two_value != -99.:
    #      #  plus_two.append(plus_two_value)
    #      #  
    #      #  boundary_plus_two.append(central_value)
    #      #  masses_plus_two_sigma.append(float(mass))

    #      #if minus_two_value != -99.:
    #      #  minus_two.append(minus_two_value)
    #      #  boundary_minus_two.append(central_value)
    #      #  masses_minus_two_sigma.append(float(mass))

    #      #if plus_one_value != -99.:
    #      #  plus_one.append(plus_one_value)
    #      #  boundary_plus_one.append(central_value)
    #      #  masses_plus_one_sigma.append(float(mass))

    #      #if minus_one_value != -99.:
    #      #  minus_one.append(minus_one_value)
    #      #  boundary_minus_one.append(central_value)
    #      #  masses_minus_one_sigma.append(float(mass))

    ##counter = 0
    ##for mass in sorted(limits2D.keys(), key=self.sortList, reverse=True):
    ##  if len(limits2D[mass]['exp_central_list']) > 1:
    ##    # plot the upper limit only in the double exclusion regime (3 -4 GeV)
    ##    if counter == 0 and float(mass) < 3: 
    ##      break
    ##    counter += 1
    ##    central_upper.append(limits2D[mass]['exp_central_list'][1])
    ##    masses_central_upper.append(float(mass))

    ##masses_central_all = masses_central + masses_central_upper
    ##central_all = central + central_upper

    ##for mass in sorted(limits2D.keys(), key=self.sortList, reverse=True):
    ##  if len(limits2D[mass]['exp_minus_one_list']) > 1:
    ##    minus_one_2.append(limits2D[mass]['exp_minus_one_list'][1])
    ##    masses_minus_one_sigma_2.append(float(mass))

    ##masses_minus_one_sigma_all = masses_minus_one_sigma + masses_minus_one_sigma_2
    ##minus_one_all = minus_one + minus_one_2

    ##if not self.do_blind:
    ##  for mass in sorted(limits2D.keys(), key=self.sortList, reverse=True):
    ##    if len(limits2D[mass]['obs_list']) > 1:
    ##      obs_2.append(limits2D[mass]['obs_list'][1])
    ##      masses_obs_2.append(float(mass))

    ##masses_obs_all = masses_obs + masses_obs_2
    ##obs_all = obs + obs_2

    #masses_tot = masses_central_all + masses_obs_all + masses_minus_one_all + masses_plus_one_all + masses_minus_two_all

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
    #print masses_obs_all

    #print 'the_central = np.array({})'.format(central)
    #plt.clf()
    #f, ax = plt.subplots(figsize=(9, 8))
    #y_range_min = 1e-5#1e-9
    #y_range_max = 1e-1#1e1
    #if not self.do_coupling_scenario:
    #  self.fe = '0.0'
    #  self.fu = '1.0'
    #  self.ft = '0.0'
    #if self.fe == '0p5': fe_label = '1/2'
    #elif self.fe == '0p3': fe_label = '1/3'
    #else: fe_label = self.fe.replace('p', '.')
    #if self.fu == '0p5': fu_label = '1/2'
    #elif self.fu == '0p3': fu_label = '1/3'
    #else: fu_label = self.fu.replace('p', '.')
    #if self.ft == '0p5': ft_label = '1/2'
    #elif self.ft == '0p3': ft_label = '1/3'
    #else: ft_label = self.ft.replace('p', '.')
    ##coupling_scenario = r'(r$_{e}$={fe}, r$_{mu}$={fu}, r$_{tau}$={ft})'.format(e='e', fe=fe_label, mu=r'\mu', fu=fu_label, tau=r'\tau', ft=ft_label)
    ##ax.text(0.1, 0.93, 'CMS', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30, fontweight='bold')
    ##ax.text(0.17, 0.84, 'Preliminary', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=25, fontstyle='italic')
    ##ax.text(0.26, 0.75, coupling_scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
    ##ax.text(0.25, 0.66, 'Lepton universality tests', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='blue', fontsize=18)
    ##ax.text(0.84, 0.93, self.scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='black', fontsize=23, fontweight='bold')
    #plt.axhline(y=1e-2, color='blue', linewidth=3, linestyle='--', zorder=10)

    ##f1 = plt.fill_between(masses_minus_two_sigma, minus_two, boundary_minus_two, color='gold'     , label=r'95% expected')
    ##f2 = plt.fill_between(masses_minus_one_sigma, minus_one, boundary_minus_one, color='forestgreen', label=r'68% expected')
    ##f3 = plt.fill_between(masses_plus_two_sigma, boundary_plus_two, plus_two, color='gold'       , label=r'95% expected')
    ##f4 = plt.fill_between(masses_plus_one_sigma, boundary_plus_one, plus_one, color='forestgreen', label=r'68% expected')
    ##p1, = plt.plot(masses_central, central, color='red', label='Median expected', linewidth=2)
    ##if not self.do_blind:
    ##  p8, = plt.plot(masses_obs, obs, color='black', label='Observed', linewidth=2)

    ##f1 = plt.plot(masses_minus_two_sigma, minus_two, marker='X', linewidth=0,  color='gold', label=r'95% expected')
    ##f2 = plt.plot(masses_minus_one_sigma, minus_one, marker='X', linewidth=0, color='forestgreen', label=r'68% expected')
    ##f3 = plt.plot(masses_plus_two_sigma, plus_two, marker='X', linewidth=0, color='gold'       , label=r'95% expected')
    ##f4 = plt.plot(masses_plus_one_sigma, plus_one, marker='X', linewidth=0, color='forestgreen', label=r'68% expected')

    ## median
    #p1, = plt.plot(masses_central_all, central_all, color='red', label='Median expected', linewidth=2)
    #p1, = plt.plot(masses_central, central, marker='X', color='red', label='Median expected', linewidth=0)
    ###p3, = plt.plot(masses_central_3, central_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    #p2, = plt.plot(masses_central_upper, central_upper, marker='X', color='blue', label='Median expected', linewidth=0)

    ## obs
    ##p1, = plt.plot(masses_obs, obs, marker='X', color='red', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_obs_upper, obs_upper, marker='X', color='blue', label='Median expected', linewidth=0)
    ###p2, = plt.plot(masses_obs_3, obs_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##p3, = plt.plot(masses_obs_missing, obs_missing, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##if not self.do_blind:
    ##  p8, = plt.plot(masses_obs_all, obs_all, color='black', label='Observed', linewidth=2)

    ## minus one
    #p1, = plt.plot(masses_minus_one, minus_one, marker='X', color='red', label='Median expected', linewidth=0)
    #p2, = plt.plot(masses_minus_one_upper, minus_one_upper, marker='X', color='blue', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_minus_one_3, minus_one_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_minus_one_missing, minus_one_missing, marker='X', color='magenta', label='Median expected', linewidth=0)
    #p1, = plt.plot(masses_minus_one_all, minus_one_all, color='forestgreen', label='Median expected', linewidth=2)
    ##f2 = plt.fill_between(masses_minus_one_sigma_all, minus_one_all, central_all, color='forestgreen', label=r'68% expected')

    ## plus one
    #p1, = plt.plot(masses_plus_one, plus_one, marker='X', color='red', label='Median expected', linewidth=0)
    #p2, = plt.plot(masses_plus_one_upper, plus_one_upper, marker='X', color='blue', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_plus_one_3, plus_one_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_plus_one_missing, plus_one_missing, marker='X', color='magenta', label='Median expected', linewidth=0)
    #p1, = plt.plot(masses_plus_one_all, plus_one_all, color='forestgreen', label='Median expected', linewidth=2)

    ## minus two
    #p1, = plt.plot(masses_minus_two, minus_two, marker='X', color='red', label='Median expected', linewidth=0)
    #p2, = plt.plot(masses_minus_two_upper, minus_two_upper, marker='X', color='blue', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_minus_two_3, minus_two_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_minus_two_missing, minus_two_missing, marker='X', color='magenta', label='Median expected', linewidth=0)
    #p1, = plt.plot(masses_minus_two_all, minus_two_all, color='gold', label='Median expected', linewidth=2)

    #p1, = plt.plot(masses_plus_two, plus_two, marker='X', color='red', label='Median expected', linewidth=0)
    #p2, = plt.plot(masses_plus_two_upper, plus_two_upper, marker='X', color='blue', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_plus_two_3, plus_two_3, marker='X', color='magenta', label='Median expected', linewidth=0)
    ##p2, = plt.plot(masses_plus_two_missing, plus_two_missing, marker='X', color='magenta', label='Median expected', linewidth=0)
    #p1, = plt.plot(masses_plus_two_all, plus_two_all, color='gold', label='Median expected', linewidth=2)



    ##p3, = plt.plot(masses_central_3, central_3, marker='X', color='magenta', label='Median expected', linewidth=0)


    ##if not self.do_blind:
    ##  p8, = plt.plot(masses_obs, obs, marker='X', color='black', label='Observed', linewidth=0)

    ##veto_D0 = plt.gca().add_patch(Rectangle((1.74, 1.01e-5), 1.8-1.74, 9e-4-1.01e-5, edgecolor='white', facecolor='white', zorder=3)) 
    ##veto_Jpsi = plt.gca().add_patch(Rectangle((3.05, 1.01e-5), 3.15-3.05, 1e-3-1.01e-5, edgecolor='white', facecolor='white', zorder=3)) 
    ##veto_Psi2S = plt.gca().add_patch(Rectangle((3.65, 1.01e-5), 3.75-3.65, 9e-2-1.01e-5, edgecolor='white', facecolor='white', zorder=3)) 

    ##if not self.do_blind:
    ##  first_legend = plt.legend(handles=[p8, p1, f1, f2], loc='lower right', fontsize=18)
    ##else:
    ##  first_legend = plt.legend(handles=[p1, f2, f1], loc='lower right', fontsize=18)
    ##ax = plt.gca().add_artist(first_legend)

    ##if self.scenario == 'Majorana':
    ##  if not self.do_coupling_scenario:
    ##    #p2, = plt.plot(db.masses_delphidisplaced, db.exp_delphidisplaced, color='darkorange', label='Delphi displaced', linewidth=1.3, linestyle='dashed')
    ##    p2, = plt.plot(db.masses_atlas_lower, db.exp_atlas_lower, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    p2_2, = plt.plot(db.masses_atlas_upper, db.exp_atlas_upper, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    p3, = plt.plot(db.masses_cmsdisplacedmuon, db.exp_cmsdisplacedmuon, color='blueviolet', label='CMS displaced', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    #p3_2, = plt.plot(db.masses_cmsdisplacedmuon_upper, db.exp_cmsdisplacedmuon_upper, color='blueviolet', label='CMS displaced', linewidth=1.3, linestyle='dashed')
    ##    p4, = plt.plot(db.masses_lhcb_peskin, db.exp_lhcb_peskin, color='darkred', label='LHCb', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    p5, = plt.plot(db.masses_belle, db.exp_belle, color='deepskyblue', label='Belle', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    #p6, = plt.plot(db.masses_charm, db.exp_charm, color='magenta', label='CHARM', linewidth=1.3, linestyle='dashed')
    ##    #p7, = plt.plot(db.masses_EXO_22_017_Majorana, db.exp_EXO_22_017_Majorana, color='deepskyblue', label='EXO-22-017', linewidth=1.3, linestyle='dashed', zorder=10)

    ##    second_legend = plt.legend(handles=[p2, p3, p4, p5], loc='lower left', fontsize=18)
    ##    ax = plt.gca().add_artist(second_legend)
    ##  else:
    ##    if self.fe == '0p0' and self.fu == '0p5' and self.ft == '0p5':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Majorana_0p0_0p5_0p5, db.exp_EXO_21_013_Majorana_0p0_0p5_0p5, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    elif self.fe == '0p5' and self.fu == '0p5' and self.ft == '0p0':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Majorana_0p5_0p5_0p0, db.exp_EXO_21_013_Majorana_0p5_0p5_0p0, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    elif self.fe == '0p3' and self.fu == '0p3' and self.ft == '0p3':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Majorana_0p3_0p3_0p3, db.exp_EXO_21_013_Majorana_0p3_0p3_0p3, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)

    ##    second_legend = plt.legend(handles=[p1], loc='lower left', fontsize=18)
    ##    ax = plt.gca().add_artist(second_legend)


    ##elif self.scenario == 'Dirac':
    ##  if not self.do_coupling_scenario:
    ##    p2, = plt.plot(db.masses_atlas_lower_dirac, db.exp_atlas_lower_dirac, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    #p2_2, = plt.plot(db.masses_atlas_upper_dirac, db.exp_atlas_upper_dirac, color='darkorange', label='ATLAS displaced', linewidth=1.3, linestyle='dashed')
    ##    p3, = plt.plot(db.masses_cmsdisplacedmuon_dirac, db.exp_cmsdisplacedmuon_dirac, color='blueviolet', label='CMS displaced', linewidth=1.3, linestyle='dashed', zorder=10)

    ##    second_legend = plt.legend(handles=[p2, p3], loc='lower left', fontsize=18)
    ##  else:
    ##    if self.fe == '0p0' and self.fu == '0p5' and self.ft == '0p5':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Dirac_0p0_0p5_0p5, db.exp_EXO_21_013_Dirac_0p0_0p5_0p5, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    elif self.fe == '0p5' and self.fu == '0p5' and self.ft == '0p0':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Dirac_0p5_0p5_0p0, db.exp_EXO_21_013_Dirac_0p5_0p5_0p0, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)
    ##    elif self.fe == '0p3' and self.fu == '0p3' and self.ft == '0p3':
    ##      p1, = plt.plot(db.masses_EXO_21_013_Dirac_0p3_0p3_0p3, db.exp_EXO_21_013_Dirac_0p3_0p3_0p3, color='darkorange', label='EXO-21-013', linewidth=1.3, linestyle='dashed', zorder=10)

    ##    second_legend = plt.legend(handles=[p1], loc='lower left', fontsize=18)
    ##    ax = plt.gca().add_artist(second_legend)

    #plt.title(lumi + ' (13 TeV)', loc='right', fontsize=23)
    #plt.ylabel(r'$|V|^2$', fontsize=23)
    #plt.yticks(fontsize=17)
    #plt.ylim(y_range_min, y_range_max)
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.xlabel(r'$m_{N}$ (GeV)', fontsize=23)
    ##plt.xlim(min(masses_central), max(masses_central))
    ##plt.xlim(min(masses_tot), max(masses_tot))
    #plt.xlim(2, max(masses_tot))
    #plt.xticks(fontsize=17)
    #plt.yscale('log')
    #plt.xscale('linear')
    #plt.grid(True)
    #if not self.do_coupling_scenario:
    #  name_2d = '2d_hnl_limit_{}'.format(self.scenario) 
    #else:
    #  name_2d = '2d_hnl_limit_scenario_{}_{}_{}_{}'.format(self.scenario, self.fe, self.fu, self.ft) 
    #plt.savefig('{}/{}.pdf'.format(plotDir, name_2d))
    #plt.savefig('{}/{}.png'.format(plotDir, name_2d))
    #print '--> {}/{}.png created'.format(plotDir, name_2d)


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


