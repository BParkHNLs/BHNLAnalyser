import os
import sys
from os import path
sys.path.append('../objects')
from coupling_scenarios import CouplingScenario



class FlavourChannelsCombinator(object):
  '''
    Class that performs the datacard combination between the electron and muons channels
  '''
  def __init__(self, datacard_name_muon, datacard_name_electron, fe, fu, ft, outdirlabel, subdirlabel):
    self.fe = float(fe)
    self.fu = float(fu)
    self.ft = float(ft)
    self.datacard_name_muon = datacard_name_muon
    self.datacard_name_electron = datacard_name_electron
    self.datacard_name_combined = self.datacard_name_muon
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel


  def getInputDirectory(self):
    '''
      Get the directory where the weighted datacards are stored
    '''
    fe = str(round(self.fe, 1)).replace('.', 'p')
    fu = str(round(self.fu, 1)).replace('.', 'p')
    ft = str(round(self.ft, 1)).replace('.', 'p')

    input_dirname = './outputs/{}/{}/weighted_datacards/coupling_{}_{}_{}'.format(self.outdirlabel, self.subdirlabel, fe, fu, ft)
    if not path.exists(input_dirname):
      raise RuntimeError('Input directory "{}" not found'.format(input_dirname))

    return input_dirname


  def createOutputDirectory(self):
    '''
      Create output directory
    '''
    output_dirname = '{}/flavour_combined'.format(self.indir)
    if not path.exists(output_dirname):
      os.system('mkdir -p {}'.format(output_dirname))

    return output_dirname


  def correctPathToWorkspaces(self, datacard):
    f_in = open(datacard)
    f_out = open(datacard+'_tmp', 'w+')

    lines = f_in.readlines()
    workspaces = []

    for line in lines:
      if self.indir not in line:
        f_out.write(line)
      else:
        line = line.replace(self.indir + '/muon/', '')
        line = line.replace(self.indir + '/electron/', '')
        f_out.write(line)
        
    f_out.close()
    command = 'mv {dc}_tmp {dc}'.format(dc=datacard)
    os.system(command)


  def combineDatacards(self):
    '''
      Combine the datacards between the flavour channels
    '''
    muon_datacard = '{indir}/muon/{mu}'.format(indir=self.indir, mu=self.datacard_name_muon)
    if not path.exists(muon_datacard):
      raise RuntimeError('Muon datacard "{}" not found'.format(muon_datacard))

    electron_datacard = '{indir}/electron/{el}'.format(indir=self.indir, el=self.datacard_name_electron)
    if not path.exists(electron_datacard):
      raise RuntimeError('Electron datacard "{}" not found'.format(electron_datacard))

    combined_datacard = '{outdir}/{cmb}'.format(outdir=self.outdir, cmb=self.datacard_name_combined)

    command = 'combineCards.py {mu} {el} > {cmb}'.format(
        mu = muon_datacard,
        el = electron_datacard,
        cmb = combined_datacard,
        )

    os.system(command)

    # correct the path to the workspaces
    self.correctPathToWorkspaces(datacard=combined_datacard)

    print '\n [FlavourChannelsCombinator] --> {}/{} created'.format(self.outdir, self.datacard_name_combined)


  def process(self):
    print '\n [FlavourChannelsCombinator] -> get input directory'
    self.indir = self.getInputDirectory()

    print '\n [FlavourChannelsCombinator] -> create output directory'
    self.outdir = self.createOutputDirectory()

    print '\n [FlavourChannelsCombinator] -> will combine the flavour channels into a single datacard'
    self.combineDatacards()



if __name__ == '__main__':

  coupling_scenario = CouplingScenario(
      fe = 0.5, 
      fu = 0.5, 
      ft = 0.0,
      )

  datacard_name_muon = 'datacard_combined_m_3p0_v2_4p4em04.txt'
  datacard_name_electron = 'HNL_M2.0_ctau100.0_PF_LxySUnder50_OS.txt'
  
  combinator = FlavourChannelsCombinator(
      datacard_name_muon = datacard_name_muon,
      datacard_name_electron = datacard_name_electron,
      coupling_scenario = coupling_scenario,
      )

  combinator.process()
