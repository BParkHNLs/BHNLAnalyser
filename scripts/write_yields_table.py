import os
import sys
import ROOT
from collections import OrderedDict

from tools import Tools
sys.path.append('../objects')
from categories import categories

tools = Tools()

output_label = 'V12_08Aug22' 
tag = 'study_categorisation_0_50_150_new_features_v2_score0p95' 
path = '../outputs/{}/datacards/{}'.format(output_label, tag)

masses = [1.0, 1.5, 2.0, 3.0, 4.5]
ctaus = [1.0, 10.0, 100.0, 1000.0]
categories = categories['categories_0_50_150']


def getSignalYields(line):
  yields = float(line[5:line.rfind('1.0')-1])
  return '{:.2e}'.format(yields)


def getBackgroundYields(workspace):
  yields = workspace.data("data_obs").createHistogram("hnl_mass").Integral()
  yields = yields * 0.2 #(2 sigma / 10 sigma)
  return '{:.2e}'.format(yields)


def getCategoryTitle(category):
  title = category.title
  title = title.replace('l_{xy}', 'L$_\mathrm{xy}$')
  title = title.replace('#sigma_{xy}', '$\\sigma$')
  title = title.replace('<=', '$\leq$')
  title = title.replace('<', '$<$')
  title = title.replace('>', '$>$')
  title = title.replace('$$', '')
  return title

table_yields = open('table_yields.txt', 'w+')

for mass in masses:
  coupling_line = '{} & &'.format(mass)
  for ictau, ctau in enumerate(ctaus):
    v2 = tools.getVV(mass, ctau)
    coupling = tools.getCouplingLabel(v2)
    if ictau != len(ctaus)-1:
      coupling_line += ' & {}'.format(coupling)
    else:
      coupling_line += ' & {} \\\ '.format(coupling)
  table_yields.write('\hline')
  table_yields.write('\n' + coupling_line)
  table_yields.write('\n' + '\hline')

  for icat, category in enumerate(categories): 
    if 'incl' in category.label: continue
    signal_yields = OrderedDict()
    data_obs_name = 'workspace_data_obs_bhnl_m_{}_cat_{}.root'.format(str(mass).replace('.', 'p'), category.label)
    try:
      data_obs_file = tools.getRootFile('{}/{}'.format(path, data_obs_name))
      workspace = data_obs_file.Get("workspace")
      background_yields = getBackgroundYields(workspace)
    except:
      background_yields = '-'
    table_entry = ' & {} & {}'.format(getCategoryTitle(category), background_yields)

    for ctau in ctaus:
      v2 = tools.getVV(mass, ctau)
      coupling = tools.getCouplingLabel(v2)

      datacard_name = 'datacard_bhnl_m_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(coupling).replace('.', 'p').replace('-', 'm'), category.label)
      try:
        card = open('{}/{}'.format(path, datacard_name))
        lines = card.readlines()
        for line in lines:
          if 'rate' not in line: continue
          signal_yields[ctau] = getSignalYields(line)
      except:
          signal_yields[ctau] = '-'

    for ictau, ctau in enumerate(ctaus):
      if ictau != len(ctaus)-1: 
        table_entry += ' & {}'.format(signal_yields[ctau])
      else:
        table_entry += ' & {} \\\ '.format(signal_yields[ctau])
    table_yields.write('\n' + table_entry)
    if icat == len(categories)-1:
      table_yields.write('\n' + '\hline')


table_yields.close()
print '-> table_yields.txt created'
os.system('cat table_yields.txt')

    


