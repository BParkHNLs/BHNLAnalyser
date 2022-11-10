# run within the python-ternary framework (https://github.com/marcharper/python-ternary)
import os
import ternary
import matplotlib.pyplot as plt
from os import path


class TernaryPlotter(object):
  def __init__(self, mass, homedir, outdirlabel, subdirlabel):
    self.mass = mass
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.plotdir = '{}/outputs/{}/limits/{}/plots'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    self.exclusion_file = '{}/exclusion_m_{}.txt'.format(self.plotdir, str(self.mass).replace('.', 'p'))
    # exclusion file should have lines in the format
    # fe fu ft exclusion
    if not path.exists(self.exclusion_file):
      raise RuntimeError('Exclusion file "{}" not found'.format(self.exclusion_file))
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
    idx4 = line.find(' ', idx3)

    return float(line[idx3:idx4])


  def heat_function(self, point):
    coupling = 1e9 # default value of the coupling used when exclusion is missing
    f = open(self.exclusion_file)
    lines = f.readlines()
    for line in lines:
      point_line = '{} {} {}'.format(round(point[0], 1), round(point[1], 1), round(point[2], 1))
      if point_line in line:
        coupling = self.get_exclusion(line) 
    f.close()

    return coupling


  def get_boundaries(self):
    f = open(self.exclusion_file)
    lines = f.readlines()
    coupling_min = 1e9
    coupling_max = 1e-9
    for line in lines:
      # remove empty lines
      if len(line) == 1: continue
      coupling = self.get_exclusion(line) 
      if coupling < coupling_min: coupling_min = coupling
      if coupling > coupling_max: coupling_max = coupling
    f.close()

    coupling_min = coupling_min - 0.15 * coupling_min
    coupling_max = coupling_max + 0.15 * coupling_max

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
    fig, ax = plt.subplots()
    ax.axis("off")
    figure, tax = ternary.figure(ax=ax, scale=self.scale)

    # plot the heatmap
    tax.heatmapf(self.heat_function, 
                 boundary=True,
                 scale=self.scale,
                 scientific=True,
                 style="hexagonal", 
                 cmap=plt.cm.get_cmap(self.cmap),
                 cbarlabel='min(central exclusion coupling)',
                 vmin=self.get_boundaries()[0],
                 vmax=self.get_boundaries()[1],
                 )

    # define axes and grid
    tax.boundary(linewidth=2.0)
    tax.left_axis_label("$f_{\\tau}$", fontsize=15, offset=0.16)
    tax.right_axis_label("$f_{\mu}$", fontsize=15, offset=0.16)
    tax.bottom_axis_label("$f_{e}$", fontsize=15, offset=0.06)

    tax.gridlines(multiple=1, linewidth=2, color='black')

    # Set and format axes ticks.
    ticks = [i / float(self.scale) for i in range(self.scale+1)]
    tax.ticks(ticks=ticks, axis='rlb', linewidth=1, fontsize=13, clockwise=False, offset=0.02, tick_formats="%0.1f")

    # format and save
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    plt.tight_layout()
    #tax.show()
    tax.savefig('{}/ternary_plot_m_{}.png'.format(self.plotdir, self.mass))
    tax.savefig('{}/ternary_plot_m_{}.pdf'.format(self.plotdir, self.mass))

    print ' -> {}/ternary_plot_m_{}.png created'.format(self.plotdir, self.mass)



if __name__ == '__main__':

  mass = '2.0'
  homedir = '/t3home/anlyon/BHNL/BHNLAnalyser/CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNLAnalyser'
  outdirlabel = 'V12_08Aug22'
  subdirlabel = 'test_interpretation_triangle_v1'
  plotter = TernaryPlotter(
    mass = mass,
    homedir = homedir,
    outdirlabel = outdirlabel,
    subdirlabel = subdirlabel,
    )
  plotter.process()



