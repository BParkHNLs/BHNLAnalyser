import os
import sys
import ROOT
from collections import OrderedDict

from tools import Tools
from compute_yields import ComputeYields
sys.path.append('../objects')
from categories import categories
from samples import data_samples
from baseline_selection import selection
from quantity import Quantity

tools = Tools()


def getSignalYields(line):
  yields = float(line[5:line.rfind('1.0')-1])
  return '{:.2e}'.format(yields)


def getBackgroundYields(workspace, mass, category_label):
  yields = workspace.data("data_obs_bhnl_m_{}_cat_{}".format(str(mass).replace('.', 'p'), category_label)).createHistogram("hnl_mass").Integral()
  yields = yields * 0.2 #(2 sigma / 10 sigma)
  return '{:.2e}'.format(yields)


def getCategoryTitle(category):
  title = category.title
  title = title.replace('L_{xy}', 'L$_\mathrm{xy}$')
  title = title.replace('#sigma', '$\\sigma$')
  title = title.replace('B_{c}', 'B$_{c}$')
  title = title.replace('<=', '$\leq$')
  title = title.replace('<', '$<$')
  title = title.replace('>', '$>$')
  title = title.replace('$$', '')
  return title


def writeYieldsTable():
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
        background_yields = getBackgroundYields(workspace, mass, category.label)
      except:
        background_yields = '-'
      table_entry = ' & {} & {}'.format(getCategoryTitle(category), background_yields)

      for ctau in ctaus:
        v2 = tools.getVV(mass, ctau)
        coupling = tools.getCouplingLabel(v2)

        #datacard_name = 'datacard_bhnl_m_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(coupling).replace('.', 'p').replace('-', 'm'), category.label)
        datacard_name = 'datacard_bhnl_m_{}_ctau_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(ctau).replace('.', 'p'), str(coupling).replace('.', 'p').replace('-', 'm'), category.label)
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


def writeEfficiencyTable():
  table_efficiency = open('table_efficiency.txt', 'w+')

  for mass in masses:
    coupling_line = '{} & &'.format(mass)
    for ictau, ctau in enumerate(ctaus):
      v2 = tools.getVV(mass, ctau)
      coupling = tools.getCouplingLabel(v2)
      if ictau != len(ctaus)-1:
        coupling_line += ' & $|$V$^2|$={}'.format(coupling)
      else:
        coupling_line += ' & $|$V$^2|$={} \\\ '.format(coupling)
    table_efficiency.write('\hline')
    table_efficiency.write('\n' + coupling_line)
    table_efficiency.write('\n' + '\hline')

    for icat, category in enumerate(categories): 
      if 'incl' in category.label: continue

      # get the number of background yields before the pNN cut
      resolution = resolution_p0 + mass * resolution_p1
      hnl_mass = Quantity(name_flat='hnl_mass', label='hnl_mass', nbins=80, bin_min=mass-2*resolution, bin_max=mass+2*resolution)

      f_data = tools.getRootFile(data_sample.filename) 
      tree_data = tools.getTree(f_data, 'signal_tree')
      background_selection = baseline_selection + ' && ' + category.definition_flat + ' && hnl_charge==0'
      hist_data = tools.createHisto(tree=tree_data, quantity=hnl_mass, selection=background_selection)
      
      background_yields_ini = hist_data.Integral() * 41.6 / 5.302

      # get the number of background yields after the pNN cut
      data_obs_name = 'workspace_data_obs_bhnl_m_{}_cat_{}.root'.format(str(mass).replace('.', 'p'), category.label)
      try:
        data_obs_file = tools.getRootFile('{}/{}'.format(path, data_obs_name))
        workspace = data_obs_file.Get("workspace")
        background_yields_sel = float(getBackgroundYields(workspace, mass, category.label))
      except:
        background_yields_sel = '-'
      #print '{} ini: {} sel: {}'.format(category.label, background_yields_ini, background_yields_sel)

      if background_yields_sel != '-':
        efficiency_background = round(background_yields_sel / background_yields_ini * 100, 2)
      else:
        efficiency_background = '-'
      #print efficiency_background

      table_entry = ' & {} & {}\%'.format(getCategoryTitle(category), efficiency_background)

      signal_yields_ini = OrderedDict()
      signal_yields_sel = OrderedDict()
      efficiency_signal = OrderedDict()
      for ctau in ctaus:
        # get the number of signal yields before the pNN cut
        signal_selection = 'ismatched==1 && ' + baseline_selection + ' && ' + category.definition_flat
        lumi = 40.0
        signal_yields_ini[ctau] = ComputeYields(signal_label=signal_label, selection=signal_selection).computeSignalYields(mass=mass, ctau=ctau, lumi=lumi, sigma_B=472.8e9, is_bc=False, strategy='inclusive', add_weight_hlt=True, add_weight_pu=True, add_weight_muid=True, weight_hlt='weight_hlt_fullBpark', weight_pusig='weight_pu_sig_tot', weight_mu0id='weight_mu0_softid', weight_muid='weight_mu_looseid')[0] 

        # get the number of signal yields after the pNN cut
        v2 = tools.getVV(mass, ctau)
        coupling = tools.getCouplingLabel(v2)

        datacard_name = 'datacard_bhnl_m_{}_ctau_{}_v2_{}_cat_{}.txt'.format(str(mass).replace('.', 'p'), str(ctau).replace('.', 'p'), str(coupling).replace('.', 'p').replace('-', 'm'), category.label)
        try:
          card = open('{}/{}'.format(path, datacard_name))
          lines = card.readlines()
          for line in lines:
            if 'rate' not in line: continue
            signal_yields_sel[ctau] = float(getSignalYields(line))
        except:
            signal_yields_sel[ctau] = '-'
        #print '{} {} ini: {} sel: {}'.format(category.label, ctau, signal_yields_ini[ctau], signal_yields_sel[ctau])

        if signal_yields_sel[ctau] != '-':
          efficiency_signal[ctau] = round(signal_yields_sel[ctau] / signal_yields_ini[ctau] * 100, 2)
          efficiency_signal[ctau] = str(efficiency_signal[ctau]) + '\%'
        else:
          efficiency_signal[ctau] = '-'
        print efficiency_signal[ctau]
        #print '{} {} {}'.format(getCategoryTitle(category), ctau, efficiency_signal[ctau])

      for ictau, ctau in enumerate(ctaus):
        if ictau != len(ctaus)-1: 
          table_entry += ' & {}'.format(efficiency_signal[ctau])
        else:
          table_entry += ' & {} \\\ '.format(efficiency_signal[ctau])

      table_efficiency.write('\n' + table_entry)
      if icat == len(categories)-1:
        table_efficiency.write('\n' + '\hline')


  table_efficiency.close()
  print '-> table_efficiency.txt created'
  os.system('cat table_efficiency.txt')
      

if __name__ == '__main__':

  output_label = 'V13_06Feb23' 
  tag = 'unblinding_Bc_fullscan_nobernstein_v2' 
  #path = '../outputs/{}/datacards/{}'.format(output_label, tag)
  path = '/work/anlyon/outputs/{}/datacards/{}'.format(output_label, tag)

  masses = [1.0, 1.5, 2.0, 3.0, 4.5]
  #ctaus = [1.0, 10.0, 100.0, 1000.0]
  ctaus = [0.01, 10.0, 1000.0, 10000.0]
  categories = categories['categories_0_50_150']

  data_label = 'V13_06Feb23'
  data_sample = data_samples[data_label][0]

  baseline_selection = selection['baseline_08Aug22'].flat

  resolution_p0 = 6.98338e-04
  resolution_p1 = 7.78382e-03 

  signal_label = 'V13_06Feb23_training_large'

  writeYieldsTable()
  #writeEfficiencyTable()


