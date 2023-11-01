# run within the python-ternary framework (https://github.com/marcharper/python-ternary)
import os
import ternary
from ternary.ternary_style import TernaryStyle
import matplotlib.pyplot as plt
from os import path
import numpy as np


# have a look at https://seaborn.pydata.org/tutorial/color_palettes.html for color palettes


class TernaryPlotter(object):
  def __init__(self, mass, scenario, homedir, outdirlabel, subdirlabel, ternary_style):
    self.mass = mass
    self.scenario = scenario
    if self.scenario not in ['Majorana', 'Dirac']:
      raise RuntimeError('Unrecognised scenario "{}"'.format(self.scenario))
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.plotdir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    self.scale = 10 # monitors the fine graining
    self.cmap = 'Blues' #'Oranges'
    self.ternary_style = ternary_style
    self.fontsize = self.ternary_style.fontsize
    

  def get_points(self, boundary=True):
    points = []
    start = 0
    if not boundary:
        start = 1
    for i in range(start, self.scale + (1 - start)):
      for j in range(start, self.scale + (1 - start) - i):
        k = self.scale - i - j
        points.append([float(i)/self.scale, float(j)/self.scale, float(k)/self.scale])
        #points.append([i, j, k])
    return points


  def get_exclusion(self, line):
    idx1 = line.find(' ')+1
    idx2 = line.find(' ', idx1)+1
    idx3 = line.find(' ', idx2)+1

    return float(line[idx3:len(line)])


  def heat_function(self, point):
    coupling = 1e9 # default value of the coupling used when exclusion is missing
    exclusion_filename = '{}/exclusion_{}_m_{}_{}_{}_{}.txt'.format(self.plotdir, self.scenario, str(self.mass).replace('.', 'p'), str(round(point[0], 1)).replace('.', 'p'), str(round(point[1], 1)).replace('.', 'p'), str(round(point[2], 1)).replace('.', 'p'))
    try:
      f = open(exclusion_filename)
      lines = f.readlines()
      for line in lines:
        point_line = '{} {} {}'.format(round(point[0], 1), round(point[1], 1), round(point[2], 1))
        if point_line in line:
          coupling = self.get_exclusion(line) 
          if coupling == -99 and '_0p0_' not in exclusion_filename: print '{} missing'.format(exclusion_filename)
      f.close()
    except:
      coupling = -99
      print '{} missing'.format(exclusion_filename)

    return coupling


  def get_boundaries(self):
    # search for min and max values
    value_min = 1e9
    value_max = 1e-9
    for i in np.linspace(1, 0, 10):
      for j in np.linspace(1, 0, 10):
        for k in np.linspace(1, 0, 10):
          if round(i, 1) + round(j, 1) + round(k, 1) != 1.: continue 
          try:
            exclusion_filename = '{}/exclusion_{}_m_{}_{}_{}_{}.txt'.format(self.plotdir, self.scenario, str(self.mass).replace('.', 'p'), str(round(i, 1)).replace('.', 'p'), str(round(j, 1)).replace('.', 'p'), str(round(k, 1)).replace('.', 'p'))
            f = open(exclusion_filename)
            lines = f.readlines()
            for line in lines:
              # remove empty lines
              if len(line) == 1: continue
              exclusion_value = self.get_exclusion(line) 
              if exclusion_value < value_min and exclusion_value != -99: value_min = exclusion_value
              if exclusion_value > value_max and exclusion_value != -99: value_max = exclusion_value
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


  def process(self):
    print '\nMass {} GeV ({})'.format(self.mass, self.scenario)

    # define figure
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
      exponent = '-' + exponent[exponent.find('e-')+3]
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
                 cbarlabel = r'Observed $|V|^{2}$',
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
    tax.savefig('{}/ternary_plot_{}_m_{}.png'.format(self.plotdir, self.scenario, str(self.mass).replace('.', 'p')))
    tax.savefig('{}/ternary_plot_{}_m_{}.pdf'.format(self.plotdir, self.scenario, str(self.mass).replace('.', 'p')))

    print ' -> {}/ternary_plot_{}_m_{}.png created'.format(self.plotdir, self.scenario, str(self.mass).replace('.', 'p'))



if __name__ == '__main__':

  print '-.-.-.-.-.-.-.-.-'
  print ' Ternary Plotter '
  print '-.-.-.-.-.-.-.-.-'

  # set ternary style
  fontsize = 7.3
  log = True
  scientific = True 
  fixed_range = True 
  range_min = 2.5e-5 #0.5e-4 
  range_max = 0.001 #6.0e-4 
  exponent = 1e-4

  ternary_style = TernaryStyle(
      fontsize = fontsize,
      log = log,
      scientific = scientific,
      fixed_range = fixed_range,
      range_min = range_min,
      range_max = range_max,
      exponent = exponent,
      )

  homedir = '/work/anlyon'
  outdirlabel = 'V13_06Feb23'
  subdirlabel = 'paper-v5_ternary'

  masses = ['1.0', '1.5', '2.0']
  scenarios = ['Majorana', 'Dirac']

  for scenario in scenarios:
    for mass in masses:
      plotter = TernaryPlotter(
        mass = mass,
        scenario = scenario,
        homedir = homedir,
        outdirlabel = outdirlabel,
        subdirlabel = subdirlabel,
        ternary_style = ternary_style,
        )
      plotter.process()



