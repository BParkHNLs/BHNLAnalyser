'''
  List of the ctau ctau_points (mm) to be processed for each mass
'''


class Points(object):
  '''
    Class that assigns to ctau list to a mass list
  '''
  def __init__(self, mass_list, ctau_list):
    self.mass_list = mass_list
    self.ctau_list = ctau_list




ctau_points = {}

ctau_points['generated'] = [
  Points(
      mass_list = [1.0],
      ctau_list = [1000.0, 100.0, 10.0],
      ),
  Points(
      mass_list = [1.5],
      ctau_list = [1000.0, 100.0, 10.0],
      ),
  Points(
      mass_list = [2.0],
      ctau_list = [1000.0, 100.0, 10.0],
      ),
  Points(
      mass_list = [3.0],
      ctau_list = [1000.0, 100.0, 10.0, 1.0],
      ),
  Points(
      mass_list = [4.5],
      ctau_list = [100.0, 10.0, 1.0, 0.1],
      ),
  ]


ctau_points['close_to_exclusion'] = [
  Points(
      mass_list = [1.0],
      ctau_list = [5000.0],
      ),
  Points(
      mass_list = [1.5],
      ctau_list = [1500.0],
      ),
  Points(
      mass_list = [2.0],
      ctau_list = [300.0],
      ),
  Points(
      mass_list = [3.0],
      ctau_list = [15.0],
      ),
  Points(
      mass_list = [4.5],
      ctau_list = [0.02],
      ),
  ]


ctau_points['baseline'] = [
  Points(
      mass_list = [1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4],
      ctau_list = [300.0, 400.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 7000.0, 10000.0, 15000.0, 20000.0],
      ),
  Points(
      mass_list = [1.5, 1.42, 1.44, 1.46, 1.48, 1.53, 1.56, 1.59, 1.62, 1.65, 1.68, 1.71, 1.74, 1.77, 1.8, 1.83, 1.86, 1.89, 1.92, 1.95, 1.98],
      ctau_list = [50.0, 70.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0],
      ),
  Points(
      mass_list = [2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95],
      ctau_list = [10.],
      ),
  Points(
      mass_list = [3.0, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0],
      ctau_list = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7 ,1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0],
      ),
  Points(
      mass_list = [4.5, 4.1, 4.2, 4.3, 4.4, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4],
      ctau_list = [0.0015, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1],
      ),
  Points(
      mass_list = [5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1],
      ctau_list = [0.0002, 0.0003, 0.0004, 0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.015],
      ),
  ]


