import ROOT

class DataSample(object):
  def __init__(self, filename, label):
    self.filename = filename
    self.label = label


class SignalSample(object): 
  def __init__(self, filename='', filename_Bc='', label='', mass='', ctau='', resolution='', filter_efficiency='', filter_efficiency_Bc='', colour=ROOT.kRed):
    self.filename = filename
    self.filename_Bc = filename_Bc
    self.label = label
    self.mass = mass
    self.ctau = ctau
    self.resolution = resolution
    self.filter_efficiency = filter_efficiency
    self.filter_efficiency_Bc = filter_efficiency_Bc
    self.colour = colour


class QCDMCSample(object):
  def __init__(self, filename, label, cross_section, filter_efficiency, colour=ROOT.kBlue):
    self.filename = filename
    self.label = label
    self.cross_section = cross_section
    self.filter_efficiency = filter_efficiency
    self.colour = colour


signal_samples = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      #label = '4.5GeV, 1.2mm (V21)',
      label = '4.5GeV, 1.2mm',
      mass = 4.5,
      ctau = 1.2,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      #label = '3GeV, 184mm (V20_emu)',
      label = '3GeV, 184mm',
      mass = 3,
      ctau = 184.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      filter_efficiency_Bc = 1.45e-01,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      #label = '1GeV, 10^{4}mm (V26)',
      label = '1GeV, 10^{4}mm',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.0104,
      filter_efficiency = 2.24e-04,
      colour = ROOT.kOrange+0
      ),
  ]


signal_samples_tag_and_probe = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_tag_and_probe_v2_tag_fired_DST_DoubleMu1i_sf.root',
      label = 'V15_control, fixed tag trigger',
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_tag_and_probe_v2_sf.root',
      label = 'V15_control, any tag trigger',
      colour = ROOT.kOrange+0,
      ),
]


signal_samples_effscan_m3 = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      label = '3GeV, 184mm (V20_emu)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  ]

signal_samples_effscan_m1 = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 1,
      ctau = 1e4,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 1,
      ctau = 10.0,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 1,
      ctau = 100.0,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 1,
      ctau = 100000.0,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  ]

signal_samples_effscan_m4p5 = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching_ct.root',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  ]

signal_samples_limits_m4p5 = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.1*1e3,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.01*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    label = '4.5GeV, 1.2mm',
  #    mass = 4.5,
  #    ctau = 1.2,
  #    resolution = 0.0373,
  #    filter_efficiency = 3.91e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.001*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.0001*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.01,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.0001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.00001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  ]


signal_samples_limits_m4p5_large = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 7e02,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 5e02,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 3e02,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 1e02,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 7e01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 5e01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 3e01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 1e01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    label = '4.5GeV, 1.2mm',
  #    mass = 4.5,
  #    ctau = 1.2,
  #    resolution = 0.0373,
  #    filter_efficiency = 3.91e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 7e00,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 5e00,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 3e00,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 1e00,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 7e-01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 5e-01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 3e-01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 1e-01,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.01,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.0001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 4.5,
  #    ctau = 0.00001,
  #    resolution = 0.0373,
  #    filter_efficiency = 4.0e-04,
  #    ),
  ]

signal_samples_limits_m3 = [
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    label = '3GeV, 184mm',
    mass = 3.0,
    ctau = 184.0,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03,
    filter_efficiency_Bc = 1.45e-01,
    colour = ROOT.kRed+1,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e00,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.1,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.01,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.001,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
]


signal_samples_limits_m3_large = [
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 7e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 5e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 3e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 7e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 5e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 3e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    label = '3GeV, 184mm',
    mass = 3.0,
    ctau = 184.0,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03,
    filter_efficiency_Bc = 1.45e-01,
    colour = ROOT.kRed+1,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 7e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 5e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 3e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 7e00,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 5e00,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 3e00,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e00,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.1,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.01,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
#SignalSample(
#    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#    mass = 3.0,
#    ctau = 0.001,
#    resolution = 0.0243,
#    filter_efficiency = 5e-03,
#    ),
]


signal_samples_limits_m1 = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e01,
      resolution = 0.0104,
      filter_efficiency = 7.0e-03,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 1,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.1,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.01,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.0001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.00001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  ]


signal_samples_limits_m1_large = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 7e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 5e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 3e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 7e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 5e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 3e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 7e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 5e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 3e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 7e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 5e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 3e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 7e01,
      resolution = 0.0104,
      filter_efficiency = 7.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 5e01,
      resolution = 0.0104,
      filter_efficiency = 7.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 3e01,
      resolution = 0.0104,
      filter_efficiency = 7.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 1.0,
      ctau = 1e01,
      resolution = 0.0104,
      filter_efficiency = 7.0e-03,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 1,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.1,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.01,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.0001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
  #    mass = 1,
  #    ctau = 0.00001,
  #    resolution = 0.0104,
  #    filter_efficiency = 7.0e-03,
  #    ),
  ]


signal_samples_limits = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root',
      label = '3GeV, 184mm (V20_emu)',
      mass = 3,
      ctau = 184,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      mass = 3,
      ctau = 10,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 1,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 0.1,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 0.01,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 0.001,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 1e-3,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 1e-4,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 1e-5,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      mass = 3,
      ctau = 1e-6,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root',
      label = '4.5GeV, 1.2mm',
      mass = 4.5,
      ctau = 1.2,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-1,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-2,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-3,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-5,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-6,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      mass = 4.5,
      ctau = 1e-7,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_updatedgenmatching.root',
      label = '1GeV, 10^{4}mm',
      mass = 1,
      ctau = 1e4,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      mass = 1,
      ctau = 1e3,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e2,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e1,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-1,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-2,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-3,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-4,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-5,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  SignalSample(
      mass = 1,
      ctau = 1e-6,
      resolution = 0.0104,
      filter_efficiency = 2.40e-04,
      ),
  ]

signal_samples_loose = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      label = '4.5GeV, 1.2mm',
      mass = 4.5,
      ctau = 1.2,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      label = '3GeV, 184mm',
      mass = 3,
      ctau = 184,
      filter_efficiency = 4e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root',
      label = '1GeV, 10^{4}mm',
      mass = 1,
      ctau = 1e4,
      colour = ROOT.kOrange+0
      ),
  ]

data_samples_F1 = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/Chunk0_n750/flat/flat_bparknano.root',
      label = 'ParkingBPH1_Run2018A (F1)',
      ),
  ]

data_samples_V02 = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_selCos2d.root',
      label = 'ParkingBPH1_Run2018A (V02)',
      ),
  ]

data_samples_V03 = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'ParkingBPH1_Run2018A (V03)',
      ),
  ]

data_samples_V04 = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V04/ParkingBPH1_Run2018A/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V04/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_updatedmatching.root',
      label = 'ParkingBPH1_Run2018A (V04)',
      ),
  ]

data_samples = [
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH1_Run2018B/merged/flat_bparknano_29Jun21_hlt_weights.root',
  #    label = 'ParkingBPH1_Run2018B (V05_29Jun21)',
  #    ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH1_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH1_Run2018A (V05_29Jun21)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH2_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH2_Run2018A (V05_29Jun21)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH3_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH3_Run2018A (V05_29Jun21)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH4_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH4_Run2018A (V05_29Jun21)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH5_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH5_Run2018A (V05_29Jun21)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH6_Run2018A/merged/flat_bparknano_29Jun21_trigger_sf.root',
      label = 'ParkingBPH6_Run2018A (V05_29Jun21)',
      ),
  ]

data_samples_triggermuon_matching_check = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'Relaxed trigger muon matching',
      ),
  ]

data_samples_small = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH1_Run2018A/Chunk0_n500/flat/flat_bparknano_selected_stdtrgmu_full_partial.root',
      label = 'ParkingBPH1_Run2018A (V03)',
      ),
  ]

data_samples_loose = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH4_Run2018B/Chunk0_n10/bparknano_loosepreselection_v1_nj1.root',
      label = 'background',
      ),
  ]

data_samples_tag_and_probe = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V06_tag_and_probe/ParkingBPH1_Run2018A/merged/flat_bparknano_tag_and_probe_v2_tag_fired_DST_DoubleMu1.root',
      label = 'ParkingBPH1_Run2018A (V06)',
      ),
  ]

qcd_samples = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      label = 'QCD_pt15to20 (V06_29Jun21)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      label = 'QCD_pt20to30 (V06_29Jun21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      label = 'QCD_pt30to50 (V06_29Jun21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      label = 'QCD_pt50to80 (V06_29Jun21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2_extmerged.root',
  #    label = 'QCD_pt80to120 (V06_29Jun21)',
  #    cross_section = 2.75842e+06,
  #    filter_efficiency = 0.03844,
  #    colour = ROOT.kRed-10, 
  #    ),
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2_extmerged.root',
  #    label = 'QCD_pt120to170 (V06_29Jun21)',
  #    cross_section = 469797,
  #    filter_efficiency = 0.05362,
  #    colour = ROOT.kOrange+6, 
  #    ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      label = 'QCD_pt170to300 (V06_29Jun21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples_V05 = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt15to20 (V05)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt20to30 (V05)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt30to50 (V05)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt50to80 (V05)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr_extmerged.root',
      label = 'QCD_pt80to120 (V05)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt80to120_ext (V05)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kOrange-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr_extmerged.root',
      label = 'QCD_pt120to170 (V05)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt120to170_ext (V05)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange-3, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V05/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_newselection_stdtrgmu_v1_corr.root',
      label = 'QCD_pt170to300 (V05)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]


qcd_samples_triggermuon_matching_check = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt15to20 (V02)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt20to30 (V02)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt30to50 (V02)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt50to80 (V02)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study_extmerged.root',
      label = 'QCD_pt80to120 (V02)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt80to120_ext (V02)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kOrange-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt120to170 (V02)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt120to170_ext (V02)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange-3, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'QCD_pt170to300 (V02)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]


qcd_samples_V04 = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt15to20 (V04)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt20to30 (V04)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt30to50 (V04)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt50to80 (V04)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1_extmerged.root',
      label = 'QCD_pt80to120 (V04)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt80to120_ext (V04)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kOrange-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1_extmerged.root',
      label = 'QCD_pt120to170 (V04)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt120to170_ext (V04)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange-3, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V04/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'QCD_pt170to300 (V04)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples_V03 = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt15to20 (V03)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt20to30 (V03)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt30to50 (V03)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt50to80 (V03)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt80to120 (V03)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt80to120_ext (V03)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kOrange-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt120to170 (V03)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt120to170_ext (V03)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange-3, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_selected_selCos2d.root',
      label = 'QCD_pt170to300 (V03)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]


qcd_samples_V02 = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt15to20 (V02)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt20to30 (V02)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt30to50 (V02)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt50to80 (V02)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt80to120 (V02)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano.root',
      label = 'QCD_pt80to120_ext (V02)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kOrange-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt120to170 (V02)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano.root',
      label = 'QCD_pt120to170_ext (V02)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange-3, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V02/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano.root',
      label = 'QCD_pt170to300 (V02)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]



