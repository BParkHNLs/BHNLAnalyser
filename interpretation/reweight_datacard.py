import os
import sys
from os import path
sys.path.append('../objects')
from coupling_scenarios import CouplingScenario
sys.path.append('../scripts')
from decays_mod import HNLDecays 
import math


class DatacardReweighter(object):
  '''
    Class that aims at reweigthing the signal rate in the datacards according to a 
    given coupling scenario
  '''
  def __init__(self, datacard_name, mass, ctau, fe, fu, ft, path_motherdir, indirlabel, outdirlabel, subdirlabel, flavour_channel):
    self.datacard_name = datacard_name
    self.mass = mass # needed for e- path to workspaces
    self.ctau = ctau # needed for e- path to workspaces
    self.fe = float(fe)
    self.fu = float(fu)
    self.ft = float(ft)
    self.path_motherdir = path_motherdir
    self.indirlabel = indirlabel
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.flavour_channel = flavour_channel
    if flavour_channel not in ['muon', 'electron']:
      raise RuntimeError('Unknown flavour channel. Please choose amongst ["muon", "electron"]')

    self.datacard = '{}/{}/{}/{}'.format(self.path_motherdir, self.flavour_channel, self.indirlabel, self.datacard_name)


  def createOutputDirectory(self):
    '''
      Create output directory
    '''
    fe = str(round(self.fe, 1)).replace('.', 'p')
    fu = str(round(self.fu, 1)).replace('.', 'p')
    ft = str(round(self.ft, 1)).replace('.', 'p')

    output_dirname = './outputs/{}/{}/weighted_datacards/coupling_{}_{}_{}/{}'.format(self.outdirlabel, self.subdirlabel, fe, fu, ft, self.flavour_channel)
    if not path.exists(output_dirname):
      os.system('mkdir -p {}'.format(output_dirname))

    return output_dirname


  def getRateList(self):
    '''
      Get the list of the signal and background rates from the initial datacard
    '''
    f = open(self.datacard)
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
    if self.flavour_channel == 'muon':
      decay_width_ini = HNLDecays(mass=self.mass, fe=0., fu=1., ft=0.).decay_rate['tot']
      decay_width_new = HNLDecays(mass=self.mass, fe=self.fe, fu=self.fu, ft=self.ft).decay_rate['tot']
      coupling_weight = self.fu * self.fu * math.pow(decay_width_ini / decay_width_new, 2) # yields inversely proportional to decay_width
      #coupling_weight = self.fu * self.fu
    elif self.flavour_channel == 'electron':
      decay_width_ini = HNLDecays(mass=self.mass, fe=0., fu=1., ft=0.).decay_rate['tot'] # the approximation that gamma_e = gamma_u is made
      decay_width_new = HNLDecays(mass=self.mass, fe=self.fe, fu=self.fu, ft=self.ft).decay_rate['tot']
      coupling_weight = self.fu * self.fe * math.pow(decay_width_ini / decay_width_new, 2) # yields inversely proportional to decay_width
      coupling_weight = coupling_weight * 0.5 * 0.5 # in the electron datacard, the yields are not normalised to (0.5, 0.5, 0), and this needs to be corrected
      #coupling_weight = self.fu * self.fe

    updated_rate_list = []
    for rate in rate_list:
      if rate != 1.:
        # modify the signal rate according to the given coupling scenario
        updated_rate = rate * coupling_weight
      else:
        # do not modify background rate, keep it to 1
        updated_rate = 1.

      updated_rate_list.append(updated_rate)

    return updated_rate_list


  def updateAlpha(self, alpha):
    '''
      Update the value of alpha
    '''
    if self.flavour_channel == 'muon':
      decay_width_ini = HNLDecays(mass=self.mass, fe=0., fu=1., ft=0.).decay_rate['tot']
      decay_width_new = HNLDecays(mass=self.mass, fe=self.fe, fu=self.fu, ft=self.ft).decay_rate['tot']
      alpha_weight = self.fu * self.fu * math.pow(decay_width_ini / decay_width_new, 2) # yields inversely proportional to decay_width
    elif self.flavour_channel == 'electron':
      decay_width_ini = HNLDecays(mass=self.mass, fe=0., fu=1., ft=0.).decay_rate['tot'] # the approximation that gamma_e = gamma_u is made
      decay_width_new = HNLDecays(mass=self.mass, fe=self.fe, fu=self.fu, ft=self.ft).decay_rate['tot']
      alpha_weight = self.fu * self.fe * math.pow(decay_width_ini / decay_width_new, 2) # yields inversely proportional to decay_width
      alpha_weight = alpha_weight * 0.5 * 0.5 # in the electron datacard, the yields are not normalised to (0.5, 0.5, 0), and this needs to be corrected

    updated_alpha = float(alpha) * alpha_weight
    updated_alpha = str(updated_alpha)

    return updated_alpha


  def getRootfile(self, line):
    '''
      Get the root file name where the worspace is stored
    '''
    idx1 = line.rfind('.root')+5
    idx2 = line.rfind('/')+1
    if idx2 == 0:
      idx2 = line.rfind(' ', 0, idx1-1)+1
    
    return line[idx2:idx1]


  def getPathToRootfile(self, line):
    '''
      Get the rootfile with path
    '''
    idx1 = line.rfind('.root')+5
    idx2 = line.rfind(' ', 0, idx1-1)+1

    return line[idx2:idx1]


  def correctPathToWorkspaces(self, datacard):
    '''
      Make sure that the path to the workspaces is correct
    '''
    f_in = open(datacard)
    f_out = open(datacard+'_tmp', 'w+')

    #if self.flavour_channel == 'muon':
    #  #updated_path = self.path_motherdir + '/muon/' + self.indirlabel + '/workspaces' 
    #  updated_path = self.path_motherdir + '/muon/' + self.indirlabel
    #elif self.flavour_channel == 'electron':
    #  #updated_path = self.path_motherdir + '/electron/' + self.indirlabel + '/workspaces'
    #  if 'multipdf' not in line:
    #    updated_path = self.path_motherdir + '/electron/' + self.indirlabel + '/Mass{mass}/Mass{mass}_ctau{ctau}/SigFits'
    #  else:
    #    updated_path = self.path_motherdir + '/electron/' + self.indirlabel + '/Mass{mass}/ws'

    lines = f_in.readlines()
    for line in lines:
      if '.root' not in line:
        f_out.write(line)
      else:
        if self.flavour_channel == 'muon':
          updated_path = self.path_motherdir + '/muon/' + self.indirlabel
        elif self.flavour_channel == 'electron':
          updated_path = self.path_motherdir + '/electron/' + self.indirlabel
          #if 'multipdf' not in line:
          #  updated_path = self.path_motherdir + '/electron/' + self.indirlabel + '/Mass{mass}/Mass{mass}_ctau{ctau}/SigFits'.format(mass=str(self.mass).replace('.', 'p'), ctau=str(self.ctau).replace('.', 'p'))
          #else:
          #  updated_path = self.path_motherdir + '/electron/' + self.indirlabel + '/Mass{mass}/ws'.format(mass=str(self.mass).replace('.', 'p'))

        rootfile = self.getRootfile(line=line)
        path_to_rootfile = self.getPathToRootfile(line=line)
        line = line.replace(path_to_rootfile, updated_path + '/' + rootfile)
        f_out.write(line)
        
    f_out.close()
    command = 'mv {dc}_tmp {dc}'.format(dc=datacard)
    os.system(command)


  def updateDatacard(self, updated_rate_list):
    '''
      Update the datacard with the updated signal rates
    '''
    updated_datacard_name = self.outdir + '/' + self.datacard_name
    updated_datacard = open(updated_datacard_name, 'w+')
    datacard = open(self.datacard)
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

    # correct the path to the workspaces
    self.correctPathToWorkspaces(datacard=updated_datacard_name)

    print '\n [DatacardReweighter] --> {} created'.format(updated_datacard_name)


  def process(self):
    print '\n [DatacardReweighter] -> create output directory'
    self.outdir = self.createOutputDirectory()

    print '\n [DatacardReweighter] -> get the list of rates'
    rate_list = self.getRateList()

    print '\n [DatacardReweighter] -> update the list of rates'
    updated_rate_list = self.updateRateList(rate_list=rate_list)

    print '\n [DatacardReweighter] -> update the datacard'
    self.updateDatacard(updated_rate_list=updated_rate_list)

    


if __name__ == '__main__':

  coupling_scenario = CouplingScenario(
      fe = 0.5, 
      fu = 0.5, 
      ft = 0.0,
      )

  #datacard_name = 'datacard_combined_m_3p0_v2_4p4em04.txt'
  datacard_name = 'HNL_M2.0_ctau100.0_PF_LxySUnder50_OS.txt'
  path_motherdir = './test'
  #dirlabel = 'study_categorisation_0_50_150_new_features_v2_score0p95'
  dirlabel = '07_10_22_test_pNN_0p99'
  
  reweighter = DatacardReweighter(
      datacard_name = datacard_name,
      coupling_scenario = coupling_scenario,
      path_motherdir = path_motherdir,
      dirlabel = dirlabel,
      flavour_channel = 'electron',
      )

  reweighter.process()
