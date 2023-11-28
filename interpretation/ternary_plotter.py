# run within the python-ternary framework (https://github.com/marcharper/python-ternary)
import os
import sys
import ternary
from ternary.ternary_style import TernaryStyle
import matplotlib.pyplot as plt
from os import path
import numpy as np
import math
import glob
import re
sys.path.append('../scripts')
from decays import HNLDecays 
from tools import Tools


# have a look at https://seaborn.pydata.org/tutorial/color_palettes.html for color palettes

from matplotlib.colors import LinearSegmentedColormap
my_map3 = LinearSegmentedColormap.from_list('my_gradient', (
    # Edit this gradient at https://eltos.github.io/gradient/#003072-348DBB-E7DDCA-FFB300-FF7108
    (0.000, (0.000, 0.188, 0.447)),
    (0.250, (0.204, 0.553, 0.733)),
    (0.500, (0.906, 0.867, 0.792)),
    (0.750, (1.000, 0.702, 0.000)),
    (1.000, (1.000, 0.443, 0.031))))

class TernaryPlotter(object):
  def __init__(self, mass, scenario, exclusion, flavour, homedir, outdirlabel, subdirlabel, ternary_style):
    self.mass = mass
    self.scenario = scenario
    if self.scenario not in ['Majorana', 'Dirac']:
      raise RuntimeError('Unrecognised scenario "{}"'.format(self.scenario))
    self.exclusion = exclusion
    if self.exclusion not in ['coupling', 'ctau']:
      raise RuntimeError('Unrecognised exclusion "{}"'.format(self.exclusion))
    self.flavour = flavour
    if self.flavour not in ['combined', 'muon', 'electron']:
      raise RuntimeError('Unrecognised flavour "{}"'.format(self.flavour))
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.plotdir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    self.scale = 10 # monitors the fine graining
    self.cmap = 'RdBu_r' #'Spectral'#'PuOr' #'Blues' #'Oranges'
    #if self.exclusion == 'ctau': self.cmap += '_r' # reverse the map
    if self.exclusion == 'coupling': 
      self.default_value = -99.
      self.z_axis_label = r'Observed $|V|^{2}$'
      ternary_style.exponent = 1e-4
    elif self.exclusion == 'ctau': 
      self.default_value = -99.
      self.z_axis_label = r'Observed $c\tau$ (mm)'
      ternary_style.exponent = 1e3
    self.ternary_style = ternary_style
    self.fontsize = self.ternary_style.fontsize
    self.tmp_dir = './tmp'
    self.tools = Tools()
    

  def get_points(self, boundary=True):
    '''
      Get the list of coupling scenarios
    '''
    points = []
    start = 0
    if not boundary:
        start = 1
    for i in range(start, self.scale + (1 - start)):
      for j in range(start, self.scale + (1 - start) - i):
        k = self.scale - i - j
        points.append([float(i)/self.scale, float(j)/self.scale, float(k)/self.scale])
    return points


  def create_dir(self):
    '''
      Create temporary directory to store the exclusion point
      per coupling scenario
    '''
    if path.exists(self.tmp_dir):
      os.system('rm -r {}'.format(self.tmp_dir))
    os.system('mkdir -p {}'.format(self.tmp_dir)) 


  def remove_dir(self):
    '''
      Delete temporary directory
    '''
    if path.exists(self.tmp_dir):
      os.system('rm -r {}'.format(self.tmp_dir))


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
      if values[icoupling] >= crossing and values[icoupling+1] <= crossing:
        coupling_up = couplings[icoupling+1]
        value_up = values[icoupling+1]
        coupling_down = couplings[icoupling]
        value_down = values[icoupling]
        break

    # do not proceed if crossing is not found
    if coupling_up == 0 or coupling_down == 1e9:
      intersection = self.default_value
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


  def get_coupling_target(self, mass, coupling, fe, fu, ft):
    '''
      Function to convert the value of the coupling
      to the given scenario
    '''
    val_coupling = float(coupling)

    val_fe = float(fe.replace('p', '.'))
    val_fu = float(fu.replace('p', '.'))
    val_ft = float(ft.replace('p', '.'))
    decay_width_ini = HNLDecays(mass=float(mass), fe=0., fu=1., ft=0.).decay_rate['tot']
    decay_width_new = HNLDecays(mass=float(mass), fe=val_fe, fu=val_fu, ft=val_ft).decay_rate['tot']
    corr = decay_width_ini / decay_width_new
    val_coupling = corr * val_coupling 

    if self.scenario == 'Dirac':
      val_coupling = 2.0 * val_coupling

    return val_coupling


  def get_exclusion_value(self, x, fe, fu, ft):
    '''
      Get the value of the exclusion according to the 
      quantity on the z-axis (colorbar)
    '''
    if self.exclusion == 'coupling':
      exclusion_val = x
    elif self.exclusion == 'ctau':
      exclusion_val = x
      if exclusion_val != self.default_value:
        # convert the value of the coupling to ctau
        fe_val = float(fe.replace('p', '.'))
        fu_val = float(fu.replace('p', '.'))
        ft_val = float(ft.replace('p', '.'))
        ismaj = True if self.scenario == 'Majorana' else False
        exclusion_val = self.tools.getCtau(mass=float(self.mass), vv=float(x), fe=fe_val, fu=fu_val, ft=ft_val, ismaj=ismaj)

    return exclusion_val


  def produce_exclusion_files(self, points):
    '''
      Create temporary files reporting the excluded value for
      each coupling scenario
    '''
    count = 1
    for point in points:
      if count%10 == 0: print '   ---> {}% completed'.format(round(float(count)/len(points)*100, 0))
      count += 1
      # define the coupling fractions
      fe = str(point[0]).replace('.', 'p')
      fu = str(point[1]).replace('.', 'p')
      ft = str(point[2]).replace('.', 'p')

      # get the results files 
      if self.flavour == 'combined':
        path_results = '{}/outputs/{}/limits/{}/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft) 
      elif self.flavour == 'muon':
        path_results = '{}/outputs/{}/limits/{}/muon/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft) 
      elif self.flavour == 'electron':
        path_results = '{}/outputs/{}/limits/{}/electron/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, fe, fu, ft) 

      fileName = 'result*{}*.txt'.format(self.scenario)
      files = [f for f in glob.glob(path_results+fileName)]
     
      v2s = []
      obs = []

      for limit_file in files:
        if 'm_{}_'.format(self.mass) not in limit_file: continue

        # get the coupling and ctau from the file name
        ctau = limit_file[limit_file.find('ctau_')+5:limit_file.find('_', limit_file.find('ctau_')+5)]
        coupling_ref = limit_file[limit_file.rfind('v2_')+3:limit_file.find('.txt')] # coupling for the Majorana + (0, 1, 0) scenario
        # correct the value of the coupling for the given scenario
        val_coupling = self.get_coupling_target(mass=self.mass, coupling=coupling_ref, fe=fe, fu=fu, ft=ft)
      
        try:
          the_file = open('{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(path_results, self.scenario, self.mass, ctau, coupling_ref), 'r')
          
          # get the necessary information from the result files
          val_obs = None
          content = the_file.readlines()
          for line in content:
            if 'Observed' in line:
              values = re.findall(r'\d+', line)
              val_obs = values[0] + '.' + values[1]
         
          v2s.append(val_coupling)
          obs.append(float(val_obs))

        except:
          print 'Cannot open {}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(path_results, self.scenario, self.mass, ctau, coupling_ref)

      # sort by coupling
      vector = zip(v2s, obs)
      vector.sort(key = lambda x : float(x[0]))
      v2s = [i[0] for i in vector]
      obs = [i[1] for i in vector]

      # find the intersection (in the coupling)   
      x_obs = self.get_intersection(v2s, obs)

      # translate the exclusion to the desired quantity
      exclusion_val = self.get_exclusion_value(x_obs, fe, fu, ft)

      # create exclusion file 
      exclusion_filename = '{}/exclusion_{}_m_{}_{}_{}_{}.txt'.format(self.tmp_dir, self.scenario, str(self.mass).replace('.', 'p'), fe, fu, ft)
      exclusion_file = open(exclusion_filename, 'w+')
      exclusion_file.write('\n{} {} {} {}'.format(fe.replace('p', '.'), fu.replace('p', '.'), ft.replace('p', '.'), exclusion_val))
      exclusion_file.close()


  def get_exclusion_from_file(self, line):
    '''
      Get exclusion from temporary file
    '''
    idx1 = line.find(' ')+1
    idx2 = line.find(' ', idx1)+1
    idx3 = line.find(' ', idx2)+1

    return float(line[idx3:len(line)])


  def heat_function(self, point):
    '''
      Associate for each coupling scenario the value of the exclusion
    '''
    exclusion_filename = '{}/exclusion_{}_m_{}_{}_{}_{}.txt'.format(self.tmp_dir, self.scenario, str(self.mass).replace('.', 'p'), str(round(point[0], 1)).replace('.', 'p'), str(round(point[1], 1)).replace('.', 'p'), str(round(point[2], 1)).replace('.', 'p'))
    try:
      f = open(exclusion_filename)
      lines = f.readlines()
      for line in lines:
        point_line = '{} {} {}'.format(round(point[0], 1), round(point[1], 1), round(point[2], 1))
        if point_line in line:
          exclusion = self.get_exclusion_from_file(line) 
          if exclusion == self.default_value and '_0p0_' not in exclusion_filename: print '{} missing'.format(exclusion_filename)
      f.close()
    except:
      exclusion = self.default_value
      print '{} missing'.format(exclusion_filename)

    if exclusion == self.default_value: exclusion = None # do not plot missing points

    return exclusion


  def get_boundaries(self):
    '''
      Get the boundaries of the z-axis (colorbar) range
    '''
    # search for min and max values
    value_min = 1e9
    value_max = 1e-9
    for i in np.linspace(1, 0, 10):
      for j in np.linspace(1, 0, 10):
        for k in np.linspace(1, 0, 10):
          if round(i, 1) + round(j, 1) + round(k, 1) != 1.: continue 
          try:
            exclusion_filename = '{}/exclusion_{}_m_{}_{}_{}_{}.txt'.format(self.tmp_dir, self.scenario, str(self.mass).replace('.', 'p'), str(round(i, 1)).replace('.', 'p'), str(round(j, 1)).replace('.', 'p'), str(round(k, 1)).replace('.', 'p'))
            f = open(exclusion_filename)
            lines = f.readlines()
            for line in lines:
              # remove empty lines
              if len(line) == 1: continue
              exclusion_value = self.get_exclusion_from_file(line) 
              if exclusion_value < value_min and exclusion_value != self.default_value: value_min = exclusion_value
              if exclusion_value > value_max and exclusion_value != self.default_value: value_max = exclusion_value
            f.close()
          except:
            continue

    if self.ternary_style.fixed_range:
      range_min = self.ternary_style.range_min
      range_max = self.ternary_style.range_max
      if range_min > value_min: print 'WARNING: min excluded value ({}) is smaller than range_min ({})'.format(value_min, range_min)
      if range_max < value_max: print 'WARNING: max excluded value ({}) is larger than range_max ({})'.format(value_max, range_max)
    else:
      range_min = value_min
      range_max = value_max

    #range_min = range_min - 0.15 * range_min
    #range_max = range_max + 0.15 * range_max

    return range_min, range_max


  def missing_points(self):
    '''
      Unused function
    '''
    grid_points = self.get_points()
    missing_points = grid_points

    f = open(self.exclusion_file)
    lines = f.readlines()

    is_missing = True
    for line in lines:
      for point in grid_points:
        point_line = '{} {} {}'.format(round(point[0], 1), round(point[1], 1), round(point[2], 1))
        if point_line in line:
          is_missing = False
          missing_points.remove(point)
          break
    f.close()

    return missing_points


  def produce_ternary_plot(self):
    '''
      Produce the ternary plot
    '''
    plt.clf()
    fig, ax = plt.subplots(figsize=(6.5/1.5, 5./1.5))
    ax.axis("off")

    # add labels
    ax.text(0.1, 0.96, 'CMS', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=2.0*self.fontsize, fontweight='bold')
    ax.text(0.15, 0.89, 'Preliminary', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=1.5*self.fontsize, fontstyle='italic')
    ax.text(0.85, 0.96, self.scenario, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=1.4*self.fontsize, fontweight='bold')
    ax.text(0.85, 0.89, r'$m_{N}$' + ' = {} GeV'.format(self.mass), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=1.4*self.fontsize)

    # add exponent on colorbar
    if self.ternary_style.scientific:
      exponent = '{:.0e}'.format(self.ternary_style.exponent)
      if self.exclusion == 'coupling':
        exponent = '-' + exponent[exponent.find('e-')+3]
      elif self.exclusion == 'ctau':
        exponent = '+' + exponent[exponent.find('e')+3]
      exponent = 'x$10^{e}$'.format(e='{'+exponent+'}')
      ax.text(1.03, 1.02, r'{}'.format(exponent), transform=ax.transAxes, fontsize=1.3*self.fontsize)

    # define ternary figure
    figure, tax = ternary.figure(ax=ax, scale=self.scale)

    # plot the heatmap
    range_min, range_max = self.get_boundaries()
    tax.heatmapf(self.heat_function, 
                 boundary = True,
                 scale = self.scale,
                 scientific = True,
                 style = "hexagonal", 
                 cmap = plt.cm.get_cmap(self.cmap),
                 cbarlabel = self.z_axis_label,
                 vmin = range_min,
                 vmax = range_max,
                 ternary_style = self.ternary_style,
                 )

    # define axes and grid
    tax.boundary(linewidth=1.0)
    tax.left_axis_label("$r_{\\tau}$", fontsize=1.7*self.fontsize, offset=0.18)
    tax.right_axis_label("$r_{\mu}$", fontsize=1.7*self.fontsize, offset=0.18)
    tax.bottom_axis_label("$r_{e}$", fontsize=1.7*self.fontsize, offset=0.08)
    tax.gridlines(multiple=1, linewidth=1.5, color='black')

    # set and format axes ticks.
    ticks = [i / float(self.scale) for i in range(self.scale+1)]
    tax.ticks(ticks=ticks, axis='rlb', linewidth=1, fontsize=1.4*self.fontsize, clockwise=False, offset=0.023, tick_formats="%0.1f")

    # add crosses on points without exclusion 
    markers = [[0.1, 0.1], [0.5, 0.5], [1, 1], [1, 1.5]]
    x = np.random.rand(10)
    y = np.random.rand(10)
    x = [i+0.5 for i in np.linspace(0, 9, 10)]
    y = np.zeros(10)
    plt.plot(x, y, marker='X', color='crimson', linewidth=0, markersize=10, markeredgecolor='white')

    # format and save
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    plt.tight_layout()
    #tax.show()
    plotname = 'ternary_plot_{}_m_{}'.format(self.scenario, str(self.mass).replace('.', 'p'))
    if self.exclusion == 'coupling': plotname += '_v2'
    elif self.exclusion == 'ctau': plotname += '_ctau'
    if self.flavour == 'combined': plotname += '_combined' 
    elif self.flavour == 'muon': plotname += '_muon' 
    elif self.flavour == 'electron': plotname += '_electron' 
    tax.savefig('{}/{}.png'.format(self.plotdir, plotname))
    tax.savefig('{}/{}.pdf'.format(self.plotdir, plotname))

    print ' -> {}/{}.png created'.format(self.plotdir, plotname)


  def process(self):
    print '\nMass {} GeV vs {} ({}, {})'.format(self.mass, self.exclusion, self.flavour, self.scenario)

    # create temporary directory
    self.create_dir()

    # get the list of coupling scenarios
    points = self.get_points()

    # get the exclusion for each scenario and store results
    # in temporary files
    print '\n -> Produce exclusion files'
    self.produce_exclusion_files(points=points)

    # produce ternary plot
    print '\n -> Produce ternary plot'
    self.produce_ternary_plot()

    # delete temporary directory
    self.remove_dir()



if __name__ == '__main__':

  print '-.-.-.-.-.-.-.-.-'
  print ' Ternary Plotter '
  print '-.-.-.-.-.-.-.-.-'

  # set ternary style
  fontsize = 7.3
  log = True
  scientific = True 
  fixed_range = False 
  range_min = 2.5e-5 #0.5e-4 
  range_max = 0.001 #6.0e-4 

  ternary_style = TernaryStyle(
      fontsize = fontsize,
      log = log,
      scientific = scientific,
      fixed_range = fixed_range,
      range_min = range_min,
      range_max = range_max,
      )

  homedir = '/work/anlyon'
  outdirlabel = 'V13_06Feb23'
  subdirlabel = 'paper-v6_ternary'

  masses = ['1.0', '1.5', '2.0']
  scenarios = ['Majorana', 'Dirac']
  flavours = ['combined', 'muon', 'electron']
  exclusions = ['coupling', 'ctau']

  #masses = ['1.0']
  #scenarios = ['Majorana']
  #flavours = ['combined']
  #exclusions = ['ctau']

  for mass in masses:
    for scenario in scenarios:
      for exclusion in exclusions:
        for flavour in flavours:
          plotter = TernaryPlotter(
            mass = mass,
            scenario = scenario,
            exclusion = exclusion,
            flavour = flavour,
            homedir = homedir,
            outdirlabel = outdirlabel,
            subdirlabel = subdirlabel,
            ternary_style = ternary_style,
            )
          plotter.process()



