# run within the python-ternary framework (https://github.com/marcharper/python-ternary)
import os
import ternary
import matplotlib.pyplot as plt
from os import path
import numpy as np


# have a look at https://seaborn.pydata.org/tutorial/color_palettes.html for color palettes


class TernaryPlotter(object):
  def __init__(self, mass, homedir, outdirlabel, subdirlabel):
    self.mass = mass
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.plotdir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    self.scale = 10 # monitors the fine graining
    self.cmap = 'Blues' #'Oranges'
    

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
    exclusion_filename = '{}/exclusion_m_{}_{}_{}_{}.txt'.format(self.plotdir, str(self.mass).replace('.', 'p'), str(round(point[0], 1)).replace('.', 'p'), str(round(point[1], 1)).replace('.', 'p'), str(round(point[2], 1)).replace('.', 'p'))
    try:
      f = open(exclusion_filename)
      lines = f.readlines()
      for line in lines:
        point_line = '{} {} {}'.format(round(point[0], 1), round(point[1], 1), round(point[2], 1))
        if point_line in line:
          coupling = self.get_exclusion(line) 
          if coupling == -99: print exclusion_filename
      f.close()
    except:
      coupling = -99
      print '{} missing'.format(exclusion_filename)

    return coupling


  def get_boundaries(self):
    coupling_min = 1e9
    coupling_max = 1e-9
    for i in np.linspace(1, 0, 10):
      for j in np.linspace(1, 0, 10):
        for k in np.linspace(1, 0, 10):
          if round(i, 1) + round(j, 1) + round(k, 1) != 1.: continue 
          exclusion_filename = '{}/exclusion_m_{}_{}_{}_{}.txt'.format(self.plotdir, str(self.mass).replace('.', 'p'), str(round(i, 1)).replace('.', 'p'), str(round(j, 1)).replace('.', 'p'), str(round(k, 1)).replace('.', 'p'))
          f = open(exclusion_filename)
          lines = f.readlines()
          for line in lines:
            # remove empty lines
            if len(line) == 1: continue
            coupling = self.get_exclusion(line) 
            if coupling < coupling_min and coupling != -99: coupling_min = coupling
            if coupling > coupling_max and coupling != -99: coupling_max = coupling
          f.close()

    print coupling_min
    print coupling_max
    #coupling_min = coupling_min - 0.15 * coupling_min
    #coupling_max = coupling_max + 0.15 * coupling_max

    return coupling_min, coupling_max


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
    # define figure
    plt.clf()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.axis("off")
    ax.text(0.1, 0.93, 'CMS', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, fontweight='bold')
    ax.text(0.15, 0.86, 'Preliminary', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15, fontstyle='italic')
    ax.text(0.85, 0.89, r'$m_{N}$' + ' = {} GeV'.format(self.mass), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=15)

    figure, tax = ternary.figure(ax=ax, scale=self.scale)

    # plot the heatmap
    tax.heatmapf(self.heat_function, 
                 boundary=True,
                 scale=self.scale,
                 scientific=True,
                 style="hexagonal", 
                 cmap=plt.cm.get_cmap(self.cmap),
                 cbarlabel=r'$|V|^{2}$ (median expected)',
                 ##vmin=1e-4,#self.get_boundaries()[0],
                 ##vmax=5e-4,#self.get_boundaries()[1],
                 vmin=self.get_boundaries()[0],
                 vmax=self.get_boundaries()[1],
                 )

    # define axes and grid
    tax.boundary(linewidth=2.0)
    tax.left_axis_label("$f_{\\tau}$", fontsize=17, offset=0.16)
    tax.right_axis_label("$f_{\mu}$", fontsize=17, offset=0.16)
    tax.bottom_axis_label("$f_{e}$", fontsize=17, offset=0.06)

    tax.gridlines(multiple=1, linewidth=2, color='black')

    # Set and format axes ticks.
    ticks = [i / float(self.scale) for i in range(self.scale+1)]
    tax.ticks(ticks=ticks, axis='rlb', linewidth=1, fontsize=14, clockwise=False, offset=0.02, tick_formats="%0.1f")

    # format and save
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    plt.tight_layout()
    #tax.show()
    tax.savefig('{}/ternary_plot_m_{}.png'.format(self.plotdir, str(self.mass).replace('.', 'p')))
    tax.savefig('{}/ternary_plot_m_{}.pdf'.format(self.plotdir, str(self.mass).replace('.', 'p')))

    print ' -> {}/ternary_plot_m_{}.png created'.format(self.plotdir, str(self.mass).replace('.', 'p'))



if __name__ == '__main__':

  mass = '3.0'
  #homedir = '/t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNLAnalyser'
  homedir = '/work/anlyon'
  outdirlabel = 'V13_06Feb23'
  subdirlabel = 'combination_ANv8'
  plotter = TernaryPlotter(
    mass = mass,
    homedir = homedir,
    outdirlabel = outdirlabel,
    subdirlabel = subdirlabel,
    )
  plotter.process()



