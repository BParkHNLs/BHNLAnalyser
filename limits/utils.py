import re

def getPointGrid(files):
  return re.findall(r'\d+', files.split('/')[len(files.split('/'))-1])
  
#/t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/BHNLAnalyser/datacards/test_scan/datacard_bhnl_m_3_v2_1p2em05.txt
#./datacards_combined/test_categories_selection/datacard_combined_3.0_2.2e-06.txt
def getMassList(files): 
  masses = []

  for limitFile in files:
    #signal_mass = getPointGrid(limitFile) 
    #if '{}.{}'.format(signal_mass[0], signal_mass[1]) not in masses:
    mass = limitFile[limitFile.rfind('_m_')+3:limitFile.rfind('_ctau_', limitFile.rfind('_m_')+3)]
    if mass not in masses:
      #masses.append('{}.{}'.format(signal_mass[0], signal_mass[1])) 
      masses.append(mass) 

  masses.sort()
  return masses


def getCoupling(files):
  signal_coupling = getPointGrid(files) 
  #coupling = '{}.{}Em{}'.format(signal_coupling[2], signal_coupling[3], signal_coupling[4])
  coupling = '{}p{}em{}'.format(signal_coupling[2], signal_coupling[3], signal_coupling[4])
  coupling = '{}p{}e+{}'.format(signal_coupling[2], signal_coupling[3], signal_coupling[4])
  return coupling




