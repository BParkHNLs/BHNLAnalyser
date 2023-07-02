import ROOT

class DataSample(object):
  def __init__(self, filename, label, lumi=0.074):
    self.filename = filename
    self.label = label
    self.lumi = lumi


class SignalSample(object): 
  def __init__(self, filename=None, filename_Bc=None, label=None, mass=None, ctau=None, resolution=None, filter_efficiency=None, filter_efficiency_Bu=None, filter_efficiency_Bd=None, filter_efficiency_Bs=None, filter_efficiency_Bc=None, n_miniaod_Bu=None, n_miniaod_Bd=None, n_miniaod_Bs=None, muon_rate=None, muon_rate_Bc=None, is_private=None, colour=ROOT.kRed):
    self.filename = filename
    self.filename_Bc = filename_Bc
    self.label = label
    self.mass = mass
    self.ctau = ctau # this is the target ctau, not necessarily the ctau with which the sample was produced
    self.resolution = resolution
    self.filter_efficiency = filter_efficiency
    self.filter_efficiency_Bu = filter_efficiency_Bu
    self.filter_efficiency_Bd = filter_efficiency_Bd
    self.filter_efficiency_Bs = filter_efficiency_Bs
    self.filter_efficiency_Bc = filter_efficiency_Bc
    self.n_miniaod_Bu = n_miniaod_Bu
    self.n_miniaod_Bd = n_miniaod_Bd
    self.n_miniaod_Bs = n_miniaod_Bs
    self.muon_rate = muon_rate
    self.muon_rate_Bc = muon_rate_Bc
    self.is_private = is_private
    self.colour = colour


class QCDMCSample(object):
  def __init__(self, filename, label, cross_section, filter_efficiency, colour=ROOT.kBlue):
    self.filename = filename
    self.label = label
    self.cross_section = cross_section
    self.filter_efficiency = filter_efficiency
    self.colour = colour

signal_samples = {}
signal_samples['private'] = [
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


signal_samples['private_29Sep21'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      label = '4.5GeV, 1.2mm (V21_29Sep21)',
      mass = 4.5,
      ctau = 1.2,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 184mm (V20_emu_29Sep21)',
      mass = 3,
      ctau = 184.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      filter_efficiency_Bc = 1.45e-01,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      label = '1GeV, 10^{4}mm (V26_29Sep21)',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.0104,
      filter_efficiency = 2.24e-04,
      colour = ROOT.kOrange+0
      ),
  ]


signal_samples['central_V11_24Apr22_benchmark'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  ]

signal_samples['central_V11_24Apr22_benchmark_loose'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30files.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n324/bparknano_nj1.root',
      label = '1 GeV, 1000 mm',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_5files.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/bparknano_nj1.root',
      label = '3 GeV, 100 mm',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_2files.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n15/bparknano_nj1.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  ]

signal_samples['V44_study_normalisation_Bu'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V44_Bu/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.00888,
      filter_efficiency = 4.22e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V44_study_normalisation_Bd'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V44_Bd/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.00888,
      filter_efficiency = 3.99e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V44_study_normalisation_Bs'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V44_Bs/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.00888,
      filter_efficiency = 1.09e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_study_normalisation'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.00888,
      filter_efficiency = 6.67e-03, # from central,  #6.82e-03, from V44_muon
      filter_efficiency_Bu = 4.22e-03,
      filter_efficiency_Bd = 3.99e-03,
      filter_efficiency_Bs = 1.09e-03,
      n_miniaod_Bu = 38789.,
      n_miniaod_Bd = 37730.,
      n_miniaod_Bs = 8113.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_study_normalisation_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      filter_efficiency = 7.75e-03, # taken from central 8.20e-03, # taken from V44_muon, to check
      filter_efficiency_Bu = 5.03e-03,
      filter_efficiency_Bd = 5.05e-03,
      filter_efficiency_Bs = 1.30e-03,
      n_miniaod_Bu = 12713.,
      n_miniaod_Bd = 12884.,
      n_miniaod_Bs = 2782.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_study_normalisation_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      filter_efficiency = 7.01e-03, # taken from central, 7.24e-03, # taken from V44_muon, to check
      filter_efficiency_Bu = 4.57e-03,
      filter_efficiency_Bd = 4.55e-03,
      filter_efficiency_Bs = 1.20e-03,
      n_miniaod_Bu = 44812.,
      n_miniaod_Bd = 44858.,
      n_miniaod_Bs = 9531.,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_study_normalisation_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02, # taken from central 1.55e-02, # taken from V44_muon, to check
      filter_efficiency_Bu = 1.42e-02,
      filter_efficiency_Bd = 6.16e-03,
      filter_efficiency_Bs = 1.63e-03,
      n_miniaod_Bu = 51125.,
      n_miniaod_Bd = 24599.,
      n_miniaod_Bs = 5185.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_study_normalisation_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      filter_efficiency = 2.55e-02, #2.38e-02, # take from central 2.55e-02, # taken from V44_muon, to check
      filter_efficiency_Bu = 2.08e-02,
      filter_efficiency_Bd = 1.25e-02,
      filter_efficiency_Bs = 2.03e-03,
      n_miniaod_Bu = 152046.,
      n_miniaod_Bd = 99224.,
      n_miniaod_Bs = 14190.,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_preselection'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      filter_efficiency = 7.75e-03, # taken from central 8.20e-03, # taken from V44_muon, to check
      filter_efficiency_Bu = 5.03e-03,
      filter_efficiency_Bd = 5.05e-03,
      filter_efficiency_Bs = 1.30e-03,
      n_miniaod_Bu = 12713.,
      n_miniaod_Bd = 12884.,
      n_miniaod_Bs = 2782.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      filter_efficiency = 7.01e-03, # taken from central, 7.24e-03, # taken from V44_muon, to check
      filter_efficiency_Bu = 4.57e-03,
      filter_efficiency_Bd = 4.55e-03,
      filter_efficiency_Bs = 1.20e-03,
      n_miniaod_Bu = 44812.,
      n_miniaod_Bd = 44858.,
      n_miniaod_Bs = 9531.,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.00888,
      filter_efficiency = 6.67e-03, # from central,  #6.82e-03, from V44_muon
      filter_efficiency_Bu = 4.22e-03,
      filter_efficiency_Bd = 3.99e-03,
      filter_efficiency_Bs = 1.09e-03,
      n_miniaod_Bu = 38789.,
      n_miniaod_Bd = 37730.,
      n_miniaod_Bs = 8113.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kBlue
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02, # taken from central 1.55e-02, # taken from V44_muon, to check
      filter_efficiency_Bu = 1.42e-02,
      filter_efficiency_Bd = 6.16e-03,
      filter_efficiency_Bs = 1.63e-03,
      n_miniaod_Bu = 51125.,
      n_miniaod_Bd = 24599.,
      n_miniaod_Bs = 5185.,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      filter_efficiency = 2.55e-02, #2.38e-02, # take from central 2.55e-02, # taken from V44_muon, to check
      filter_efficiency_Bu = 2.08e-02,
      filter_efficiency_Bd = 1.25e-02,
      filter_efficiency_Bs = 2.03e-03,
      n_miniaod_Bu = 152046.,
      n_miniaod_Bd = 99224.,
      n_miniaod_Bs = 14190.,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+8
      ),
  ]


signal_samples['V13_06Feb23_preselection'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_norm.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_norm.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample( # reweighted to 0.01
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.5 GeV, 0.01 mm',
      mass = 4.5,
      ctau = 0.01,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  ]

signal_samples['V13_06Feb23_trackid'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kGreen+8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kBlue+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample( # reweighted to 0.01
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 0.01 mm',
      mass = 4.5,
      ctau = 0.01,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  ]

signal_samples['V13_06Feb23_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V13_06Feb23_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+3
      ),
  ]

signal_samples['V13_06Feb23_m2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+3
      ),
  ]

signal_samples['V13_06Feb23_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root', #FIXME take version with gen-matching
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      filter_efficiency_Bc = 8.72e-02,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.82e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_06Feb23.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  ]

signal_samples['V13_06Feb23_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.30e-1,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V13_06Feb23_m5p5'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '5.5 GeV, 10 mm',
      mass = 5.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '5.5 GeV, 1 mm',
      mass = 5.5,
      ctau = 1.0,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '5.5 GeV, 0.1 mm',
      mass = 5.5,
      ctau = 0.1,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p01mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '5.5 GeV, 0.01 mm',
      mass = 5.5,
      ctau = 0.01,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  ]


signal_samples['V13_06Feb23_training_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root', #FIXME take version with gen-matching
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      filter_efficiency_Bc = 8.72e-02,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.82e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.30e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V12_08Aug22_benchmark'] = [
  #SignalSample(
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst_v2.root',
  #    label = '2 GeV, 1000 mm',
  #    mass = 2.0,
  #    ctau = 1000.0,
  #    resolution = 0.0167,
  #    filter_efficiency = 1.66e-03,
  #    muon_rate = 0.44,
  #    is_private = False,
  #    colour = ROOT.kMagenta-10
  #    ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst_v2.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_deltaR0p1_deltaPt3.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst_v2.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_deltaR0p1_deltaPt3.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst_v2.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_deltaR0p1_deltaPt3.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  ]


signal_samples['V12_08Aug22_sensitivity'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V12_08Aug22_luminorm'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_08Aug22.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  ]

signal_samples['V12_08Aug22_training'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n191/flat/flat_bparknano_08Aug22_sr_nj1.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/flat/flat_bparknano_08Aug22_nj1.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  ]

signal_samples['V12_08Aug22_training_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root', #FIXME take version with gen-matching
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      filter_efficiency_Bc = 8.72e-02,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.82e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.30e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V12_08Aug22_genmatching'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root', #FIXME take version with gen-matching
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      filter_efficiency_Bc = 8.72e-02,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.82e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.30e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_v2.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V12_08Aug22_m1'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V12_08Aug22_m1p5'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kGreen-8
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen-3
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kGreen+3
      ),
  ]

signal_samples['V12_08Aug22_m2'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kMagenta-10
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+0
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kMagenta+3
      ),
  ]

signal_samples['V12_08Aug22_m3'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root', #FIXME take version with gen-matching
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      filter_efficiency_Bc = 8.72e-02,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.82e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.99e-01,
      muon_rate = 0.47, #TODO check
      is_private = False,
      colour = ROOT.kRed-4
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      filter_efficiency_Bc = 1.86e-01,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  ]

signal_samples['V12_08Aug22_m4p5'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.30e-1,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_syst.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['V12_08Aug22_m5p5'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '5.5 GeV, 10 mm',
      mass = 5.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '5.5 GeV, 1 mm',
      mass = 5.5,
      ctau = 1.0,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '5.5 GeV, 0.1 mm',
      mass = 5.5,
      ctau = 0.1,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.39e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p01mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_sr.root',
      label = '5.5 GeV, 0.01 mm',
      mass = 5.5,
      ctau = 0.01,
      resolution = 0.0391,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  ]

signal_samples['V42_06Feb23_m1p0_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.0 GeV, 10.0 mm',
     mass = 1.0,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00509474610858,
     filter_efficiency_Bd = 0.00505787918375,
     filter_efficiency_Bs = 0.00131789931582,
     n_miniaod_Bu = 12713.0,
     n_miniaod_Bd = 12884.0,
     n_miniaod_Bs = 2782.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p02_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p02_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.02 GeV, 10.0 mm',
     mass = 1.02,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00503064774422,
     filter_efficiency_Bd = 0.00506791966914,
     filter_efficiency_Bs = 0.00131093058305,
     n_miniaod_Bu = 39039.0,
     n_miniaod_Bd = 39815.0,
     n_miniaod_Bs = 8397.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p04_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p04_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.04 GeV, 10.0 mm',
     mass = 1.04,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00497473668818,
     filter_efficiency_Bd = 0.00497376345838,
     filter_efficiency_Bs = 0.00130765369091,
     n_miniaod_Bu = 39299.0,
     n_miniaod_Bd = 40508.0,
     n_miniaod_Bs = 8603.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p06_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p06_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.06 GeV, 10.0 mm',
     mass = 1.06,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00496663523959,
     filter_efficiency_Bd = 0.00494376336892,
     filter_efficiency_Bs = 0.00130398351908,
     n_miniaod_Bu = 39525.0,
     n_miniaod_Bd = 40301.0,
     n_miniaod_Bs = 8556.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p08_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p08_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.08 GeV, 10.0 mm',
     mass = 1.08,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00496104730855,
     filter_efficiency_Bd = 0.00491215859269,
     filter_efficiency_Bs = 0.00129988556157,
     n_miniaod_Bu = 43448.0,
     n_miniaod_Bd = 44038.0,
     n_miniaod_Bs = 9389.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p1_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.1 GeV, 10.0 mm',
     mass = 1.1,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00489328097605,
     filter_efficiency_Bd = 0.00492073824034,
     filter_efficiency_Bs = 0.00130637184596,
     n_miniaod_Bu = 45303.0,
     n_miniaod_Bd = 46349.0,
     n_miniaod_Bs = 9818.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p12_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p12_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.12 GeV, 10.0 mm',
     mass = 1.12,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00497150527898,
     filter_efficiency_Bd = 0.00492388160179,
     filter_efficiency_Bs = 0.00130170494182,
     n_miniaod_Bu = 42692.0,
     n_miniaod_Bd = 43544.0,
     n_miniaod_Bs = 9199.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p14_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p14_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.14 GeV, 10.0 mm',
     mass = 1.14,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00487669172932,
     filter_efficiency_Bd = 0.00487266307699,
     filter_efficiency_Bs = 0.00129788707957,
     n_miniaod_Bu = 63419.0,
     n_miniaod_Bd = 64911.0,
     n_miniaod_Bs = 13846.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p16_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p16_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.16 GeV, 10.0 mm',
     mass = 1.16,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00489325377664,
     filter_efficiency_Bd = 0.00485449528867,
     filter_efficiency_Bs = 0.00128512889534,
     n_miniaod_Bu = 35302.0,
     n_miniaod_Bd = 36196.0,
     n_miniaod_Bs = 7672.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p18_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p18_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.18 GeV, 10.0 mm',
     mass = 1.18,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00486403226575,
     filter_efficiency_Bd = 0.00485147080133,
     filter_efficiency_Bs = 0.00129342095237,
     n_miniaod_Bu = 47811.0,
     n_miniaod_Bd = 48890.0,
     n_miniaod_Bs = 10331.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p2_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.2 GeV, 10.0 mm',
     mass = 1.2,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00487101018949,
     filter_efficiency_Bd = 0.00483910883797,
     filter_efficiency_Bs = 0.0012855087506,
     n_miniaod_Bu = 45241.0,
     n_miniaod_Bd = 45908.0,
     n_miniaod_Bs = 9708.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p22_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p22_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.22 GeV, 10.0 mm',
     mass = 1.22,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00483813781701,
     filter_efficiency_Bd = 0.00482652901426,
     filter_efficiency_Bs = 0.00128016335403,
     n_miniaod_Bu = 42962.0,
     n_miniaod_Bd = 42901.0,
     n_miniaod_Bs = 9172.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p24_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p24_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.24 GeV, 10.0 mm',
     mass = 1.24,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00481950077112,
     filter_efficiency_Bd = 0.0048593407215,
     filter_efficiency_Bs = 0.00127626002776,
     n_miniaod_Bu = 41498.0,
     n_miniaod_Bd = 41804.0,
     n_miniaod_Bs = 8930.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p26_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p26_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.26 GeV, 10.0 mm',
     mass = 1.26,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00474695977549,
     filter_efficiency_Bd = 0.00478819188192,
     filter_efficiency_Bs = 0.0012683006733,
     n_miniaod_Bu = 32262.0,
     n_miniaod_Bd = 32854.0,
     n_miniaod_Bs = 7103.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p28_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p28_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.28 GeV, 10.0 mm',
     mass = 1.28,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00474720718646,
     filter_efficiency_Bd = 0.00475448190896,
     filter_efficiency_Bs = 0.00125915580583,
     n_miniaod_Bu = 35107.0,
     n_miniaod_Bd = 35139.0,
     n_miniaod_Bs = 7580.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p3_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.3 GeV, 10.0 mm',
     mass = 1.3,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00477730107438,
     filter_efficiency_Bd = 0.00478081439869,
     filter_efficiency_Bs = 0.00126792126792,
     n_miniaod_Bu = 45362.0,
     n_miniaod_Bd = 45775.0,
     n_miniaod_Bs = 9723.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p32_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p32_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.32 GeV, 10.0 mm',
     mass = 1.32,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00474419226229,
     filter_efficiency_Bd = 0.0047281207565,
     filter_efficiency_Bs = 0.00124868822035,
     n_miniaod_Bu = 42053.0,
     n_miniaod_Bd = 42961.0,
     n_miniaod_Bs = 9148.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p34_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p34_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.34 GeV, 10.0 mm',
     mass = 1.34,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00467727728765,
     filter_efficiency_Bd = 0.00473455463304,
     filter_efficiency_Bs = 0.00124408999679,
     n_miniaod_Bu = 38201.0,
     n_miniaod_Bd = 38316.0,
     n_miniaod_Bs = 8260.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p36_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p36_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.36 GeV, 10.0 mm',
     mass = 1.36,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00464313869077,
     filter_efficiency_Bd = 0.00469492807606,
     filter_efficiency_Bs = 0.00124441686102,
     n_miniaod_Bu = 34269.0,
     n_miniaod_Bd = 34330.0,
     n_miniaod_Bs = 7342.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p38_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p38_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.38 GeV, 10.0 mm',
     mass = 1.38,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00465482419903,
     filter_efficiency_Bd = 0.00463492176895,
     filter_efficiency_Bs = 0.0012315530123,
     n_miniaod_Bu = 36772.0,
     n_miniaod_Bd = 37159.0,
     n_miniaod_Bs = 7965.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p4_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.4 GeV, 10.0 mm',
     mass = 1.4,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00467405867406,
     filter_efficiency_Bd = 0.00456491696008,
     filter_efficiency_Bs = 0.00123610224106,
     n_miniaod_Bu = 44870.0,
     n_miniaod_Bd = 45541.0,
     n_miniaod_Bs = 9728.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p42_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p42_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.42 GeV, 10.0 mm',
     mass = 1.42,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00467338622345,
     filter_efficiency_Bd = 0.00459370362398,
     filter_efficiency_Bs = 0.00122140940897,
     n_miniaod_Bu = 42425.0,
     n_miniaod_Bd = 42312.0,
     n_miniaod_Bs = 9273.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p44_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p44_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.44 GeV, 10.0 mm',
     mass = 1.44,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00460797411126,
     filter_efficiency_Bd = 0.00454100854101,
     filter_efficiency_Bs = 0.0012187713285,
     n_miniaod_Bu = 31874.0,
     n_miniaod_Bd = 31807.0,
     n_miniaod_Bs = 6830.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p46_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p46_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.46 GeV, 10.0 mm',
     mass = 1.46,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0045368906059,
     filter_efficiency_Bd = 0.00458469877217,
     filter_efficiency_Bs = 0.00119427293004,
     n_miniaod_Bu = 38419.0,
     n_miniaod_Bd = 38760.0,
     n_miniaod_Bs = 8212.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p48_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p48_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.48 GeV, 10.0 mm',
     mass = 1.48,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00454621056717,
     filter_efficiency_Bd = 0.00450559826092,
     filter_efficiency_Bs = 0.00119275169644,
     n_miniaod_Bu = 37345.0,
     n_miniaod_Bd = 37698.0,
     n_miniaod_Bs = 8090.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p5_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.5 GeV, 10.0 mm',
     mass = 1.5,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00452781780607,
     filter_efficiency_Bd = 0.00448656954069,
     filter_efficiency_Bs = 0.0011658951363,
     n_miniaod_Bu = 58515.0,
     n_miniaod_Bd = 58314.0,
     n_miniaod_Bs = 12316.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p53_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p53_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.53 GeV, 10.0 mm',
     mass = 1.53,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00446668848156,
     filter_efficiency_Bd = 0.00451702475324,
     filter_efficiency_Bs = 0.00116064930841,
     n_miniaod_Bu = 43054.0,
     n_miniaod_Bd = 43139.0,
     n_miniaod_Bs = 9213.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p56_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p56_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.56 GeV, 10.0 mm',
     mass = 1.56,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0044525619068,
     filter_efficiency_Bd = 0.00454871237095,
     filter_efficiency_Bs = 0.00115548548402,
     n_miniaod_Bu = 43052.0,
     n_miniaod_Bd = 43346.0,
     n_miniaod_Bs = 9376.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p59_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p59_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.59 GeV, 10.0 mm',
     mass = 1.59,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00442377687951,
     filter_efficiency_Bd = 0.00449831619899,
     filter_efficiency_Bs = 0.00114489563965,
     n_miniaod_Bu = 33698.0,
     n_miniaod_Bd = 33795.0,
     n_miniaod_Bs = 7444.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p62_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p62_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.62 GeV, 10.0 mm',
     mass = 1.62,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00438951137808,
     filter_efficiency_Bd = 0.00434858520694,
     filter_efficiency_Bs = 0.0011339436347,
     n_miniaod_Bu = 34015.0,
     n_miniaod_Bd = 34135.0,
     n_miniaod_Bs = 7417.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p65_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.65 GeV, 10.0 mm',
     mass = 1.65,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00435525067503,
     filter_efficiency_Bd = 0.0043864293154,
     filter_efficiency_Bs = 0.00111130755866,
     n_miniaod_Bu = 28489.0,
     n_miniaod_Bd = 28710.0,
     n_miniaod_Bs = 6178.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p68_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p68_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.68 GeV, 10.0 mm',
     mass = 1.68,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00431371593656,
     filter_efficiency_Bd = 0.00428054376739,
     filter_efficiency_Bs = 0.00111399325335,
     n_miniaod_Bu = 36111.0,
     n_miniaod_Bd = 36422.0,
     n_miniaod_Bs = 7815.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p71_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p71_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.71 GeV, 10.0 mm',
     mass = 1.71,
     ctau = 10.0,
     filter_efficiency_Bu = 0.004318667588,
     filter_efficiency_Bd = 0.00426780324991,
     filter_efficiency_Bs = 0.00111814290331,
     n_miniaod_Bu = 28064.0,
     n_miniaod_Bd = 27706.0,
     n_miniaod_Bs = 5987.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p74_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p74_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.74 GeV, 10.0 mm',
     mass = 1.74,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00432341784992,
     filter_efficiency_Bd = 0.00421893096565,
     filter_efficiency_Bs = 0.0011104034509,
     n_miniaod_Bu = 36977.0,
     n_miniaod_Bd = 36920.0,
     n_miniaod_Bs = 7809.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p77_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p77_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.77 GeV, 10.0 mm',
     mass = 1.77,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00428765264586,
     filter_efficiency_Bd = 0.00421204908595,
     filter_efficiency_Bs = 0.00109929812752,
     n_miniaod_Bu = 37748.0,
     n_miniaod_Bd = 37679.0,
     n_miniaod_Bs = 8230.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p8_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.8 GeV, 10.0 mm',
     mass = 1.8,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00423519281372,
     filter_efficiency_Bd = 0.00424472968378,
     filter_efficiency_Bs = 0.001097398713,
     n_miniaod_Bu = 41709.0,
     n_miniaod_Bd = 42010.0,
     n_miniaod_Bs = 8842.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p83_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p83_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.83 GeV, 10.0 mm',
     mass = 1.83,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00426373441423,
     filter_efficiency_Bd = 0.00420024728873,
     filter_efficiency_Bs = 0.00109027478344,
     n_miniaod_Bu = 35142.0,
     n_miniaod_Bd = 35647.0,
     n_miniaod_Bs = 7530.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p86_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p86_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.86 GeV, 10.0 mm',
     mass = 1.86,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00424279608156,
     filter_efficiency_Bd = 0.00415812939532,
     filter_efficiency_Bs = 0.0010705977424,
     n_miniaod_Bu = 37301.0,
     n_miniaod_Bd = 37089.0,
     n_miniaod_Bs = 7912.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p89_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p89_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.89 GeV, 10.0 mm',
     mass = 1.89,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00419456002391,
     filter_efficiency_Bd = 0.00409016093339,
     filter_efficiency_Bs = 0.00107245829655,
     n_miniaod_Bu = 29212.0,
     n_miniaod_Bd = 28767.0,
     n_miniaod_Bs = 6138.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p92_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p92_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.92 GeV, 10.0 mm',
     mass = 1.92,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0042312418826,
     filter_efficiency_Bd = 0.00407497863727,
     filter_efficiency_Bs = 0.00106447262649,
     n_miniaod_Bu = 59727.0,
     n_miniaod_Bd = 59294.0,
     n_miniaod_Bs = 12894.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p95_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.95 GeV, 10.0 mm',
     mass = 1.95,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00416420902132,
     filter_efficiency_Bd = 0.00410193270368,
     filter_efficiency_Bs = 0.00106006871021,
     n_miniaod_Bu = 36548.0,
     n_miniaod_Bd = 36222.0,
     n_miniaod_Bs = 7834.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m1p98_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p98_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '1.98 GeV, 10.0 mm',
     mass = 1.98,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00423035922688,
     filter_efficiency_Bd = 0.00406869558737,
     filter_efficiency_Bs = 0.0010639727502,
     n_miniaod_Bu = 35610.0,
     n_miniaod_Bd = 35102.0,
     n_miniaod_Bs = 7349.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p0_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.0 GeV, 10.0 mm',
     mass = 2.0,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00421887614228,
     filter_efficiency_Bd = 0.0040456931854,
     filter_efficiency_Bs = 0.00104751025391,
     n_miniaod_Bu = 56785.0,
     n_miniaod_Bd = 55220.0,
     n_miniaod_Bs = 11874.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p05_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.05 GeV, 10.0 mm',
     mass = 2.05,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00425872017684,
     filter_efficiency_Bd = 0.0041043342972,
     filter_efficiency_Bs = 0.00104556939967,
     n_miniaod_Bu = 29212.0,
     n_miniaod_Bd = 28754.0,
     n_miniaod_Bs = 6076.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p1_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.1 GeV, 10.0 mm',
     mass = 2.1,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00428024745382,
     filter_efficiency_Bd = 0.00413516521868,
     filter_efficiency_Bs = 0.00107198400355,
     n_miniaod_Bu = 40285.0,
     n_miniaod_Bd = 38847.0,
     n_miniaod_Bs = 8378.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p15_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.15 GeV, 10.0 mm',
     mass = 2.15,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00436772410742,
     filter_efficiency_Bd = 0.00423450826027,
     filter_efficiency_Bs = 0.00107076441316,
     n_miniaod_Bu = 39674.0,
     n_miniaod_Bd = 38699.0,
     n_miniaod_Bs = 8085.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p2_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.2 GeV, 10.0 mm',
     mass = 2.2,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00440639669216,
     filter_efficiency_Bd = 0.00419186816935,
     filter_efficiency_Bs = 0.00107536722647,
     n_miniaod_Bu = 37255.0,
     n_miniaod_Bd = 35590.0,
     n_miniaod_Bs = 7661.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p25_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.25 GeV, 10.0 mm',
     mass = 2.25,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0045108464992,
     filter_efficiency_Bd = 0.00418046759546,
     filter_efficiency_Bs = 0.00109038945813,
     n_miniaod_Bu = 32317.0,
     n_miniaod_Bd = 30814.0,
     n_miniaod_Bs = 6462.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p3_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.3 GeV, 10.0 mm',
     mass = 2.3,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00458287783319,
     filter_efficiency_Bd = 0.00426012772726,
     filter_efficiency_Bs = 0.00108311069391,
     n_miniaod_Bu = 42821.0,
     n_miniaod_Bd = 40272.0,
     n_miniaod_Bs = 8566.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p35_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.35 GeV, 10.0 mm',
     mass = 2.35,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00468138388476,
     filter_efficiency_Bd = 0.00431074291049,
     filter_efficiency_Bs = 0.00110304610549,
     n_miniaod_Bu = 33028.0,
     n_miniaod_Bd = 30180.0,
     n_miniaod_Bs = 6512.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p4_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.4 GeV, 10.0 mm',
     mass = 2.4,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00488761790477,
     filter_efficiency_Bd = 0.00440275062976,
     filter_efficiency_Bs = 0.00111904566472,
     n_miniaod_Bu = 33143.0,
     n_miniaod_Bd = 30138.0,
     n_miniaod_Bs = 6539.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p45_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.45 GeV, 10.0 mm',
     mass = 2.45,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00505945537142,
     filter_efficiency_Bd = 0.00445445098139,
     filter_efficiency_Bs = 0.00113978237202,
     n_miniaod_Bu = 42072.0,
     n_miniaod_Bd = 37337.0,
     n_miniaod_Bs = 8048.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p5_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.5 GeV, 10.0 mm',
     mass = 2.5,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00534595734532,
     filter_efficiency_Bd = 0.00458240421606,
     filter_efficiency_Bs = 0.00115163598091,
     n_miniaod_Bu = 45556.0,
     n_miniaod_Bd = 39872.0,
     n_miniaod_Bs = 8564.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p55_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.55 GeV, 10.0 mm',
     mass = 2.55,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00554257744241,
     filter_efficiency_Bd = 0.00469250911583,
     filter_efficiency_Bs = 0.00118177607215,
     n_miniaod_Bu = 40255.0,
     n_miniaod_Bd = 34335.0,
     n_miniaod_Bs = 7246.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p6_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.6 GeV, 10.0 mm',
     mass = 2.6,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00581773044709,
     filter_efficiency_Bd = 0.00472393779266,
     filter_efficiency_Bs = 0.00119962460002,
     n_miniaod_Bu = 48817.0,
     n_miniaod_Bd = 40470.0,
     n_miniaod_Bs = 8567.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p65_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.65 GeV, 10.0 mm',
     mass = 2.65,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00621758061463,
     filter_efficiency_Bd = 0.00479392716523,
     filter_efficiency_Bs = 0.00122806029205,
     n_miniaod_Bu = 47953.0,
     n_miniaod_Bd = 37794.0,
     n_miniaod_Bs = 8019.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p7_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.7 GeV, 10.0 mm',
     mass = 2.7,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00667708481318,
     filter_efficiency_Bd = 0.00486746544259,
     filter_efficiency_Bs = 0.00125570762671,
     n_miniaod_Bu = 52768.0,
     n_miniaod_Bd = 39698.0,
     n_miniaod_Bs = 8409.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p75_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.75 GeV, 10.0 mm',
     mass = 2.75,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00741341477924,
     filter_efficiency_Bd = 0.00497503840246,
     filter_efficiency_Bs = 0.00128220596322,
     n_miniaod_Bu = 44104.0,
     n_miniaod_Bd = 31428.0,
     n_miniaod_Bs = 6654.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p8_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.8 GeV, 10.0 mm',
     mass = 2.8,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00821676504996,
     filter_efficiency_Bd = 0.00508682643659,
     filter_efficiency_Bs = 0.00129320020691,
     n_miniaod_Bu = 52216.0,
     n_miniaod_Bd = 34824.0,
     n_miniaod_Bs = 7266.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p85_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.85 GeV, 10.0 mm',
     mass = 2.85,
     ctau = 10.0,
     filter_efficiency_Bu = 0.00910408359179,
     filter_efficiency_Bd = 0.00523293299147,
     filter_efficiency_Bs = 0.00134463947702,
     n_miniaod_Bu = 47929.0,
     n_miniaod_Bd = 29050.0,
     n_miniaod_Bs = 6087.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p9_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.9 GeV, 10.0 mm',
     mass = 2.9,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0106178308509,
     filter_efficiency_Bd = 0.00554579960523,
     filter_efficiency_Bs = 0.0014208522654,
     n_miniaod_Bu = 177690.0,
     n_miniaod_Bd = 98760.0,
     n_miniaod_Bs = 20853.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m2p95_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '2.95 GeV, 10.0 mm',
     mass = 2.95,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0123944809849,
     filter_efficiency_Bd = 0.00581035831256,
     filter_efficiency_Bs = 0.00153500701049,
     n_miniaod_Bu = 47607.0,
     n_miniaod_Bd = 23827.0,
     n_miniaod_Bs = 5103.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p0_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.0 GeV, 10.0 mm',
     mass = 3.0,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0141719975756,
     filter_efficiency_Bd = 0.00624397658164,
     filter_efficiency_Bs = 0.00166144439192,
     n_miniaod_Bu = 51125.0,
     n_miniaod_Bd = 24599.0,
     n_miniaod_Bs = 5185.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.0 GeV, 1.0 mm',
     mass = 3.0,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0142412589035,
     filter_efficiency_Bd = 0.0062782647922,
     filter_efficiency_Bs = 0.00166766042776,
     filter_efficiency_Bc = 0.187197580645,
     n_miniaod_Bu = 55643.0,
     n_miniaod_Bd = 26006.0,
     n_miniaod_Bs = 5615.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p05_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.05 GeV, 1.0 mm',
     mass = 3.05,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0159491285916,
     filter_efficiency_Bd = 0.00692491035194,
     filter_efficiency_Bs = 0.00187123062859,
     filter_efficiency_Bc = 0.189521772426,
     n_miniaod_Bu = 54096.0,
     n_miniaod_Bd = 24633.0,
     n_miniaod_Bs = 5361.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.05 GeV, 10.0 mm',
     mass = 3.05,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0160755712231,
     filter_efficiency_Bd = 0.00694755267112,
     filter_efficiency_Bs = 0.0018611018611,
     n_miniaod_Bu = 44970.0,
     n_miniaod_Bd = 20314.0,
     n_miniaod_Bs = 4357.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p1_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.1 GeV, 1.0 mm',
     mass = 3.1,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0177642756058,
     filter_efficiency_Bd = 0.00788671903104,
     filter_efficiency_Bs = 0.00206567913941,
     filter_efficiency_Bc = 0.18789696542,
     n_miniaod_Bu = 62086.0,
     n_miniaod_Bd = 28734.0,
     n_miniaod_Bs = 6261.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.1 GeV, 10.0 mm',
     mass = 3.1,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0177642756058,
     filter_efficiency_Bd = 0.0078640824671,
     filter_efficiency_Bs = 0.00208354733778,
     n_miniaod_Bu = 48314.0,
     n_miniaod_Bd = 22225.0,
     n_miniaod_Bs = 4769.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p15_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.15 GeV, 10.0 mm',
     mass = 3.15,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0189858746005,
     filter_efficiency_Bd = 0.00896869763289,
     filter_efficiency_Bs = 0.00228411047268,
     n_miniaod_Bu = 49198.0,
     n_miniaod_Bd = 23505.0,
     n_miniaod_Bs = 5041.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.15 GeV, 1.0 mm',
     mass = 3.15,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0190004762877,
     filter_efficiency_Bd = 0.00890733097694,
     filter_efficiency_Bs = 0.00227963027963,
     filter_efficiency_Bc = 0.186738304285,
     n_miniaod_Bu = 211513.0,
     n_miniaod_Bd = 101169.0,
     n_miniaod_Bs = 21515.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p2_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.2 GeV, 1.0 mm',
     mass = 3.2,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0195597598973,
     filter_efficiency_Bd = 0.00967682723967,
     filter_efficiency_Bs = 0.00246601680235,
     filter_efficiency_Bc = 0.183691926885,
     n_miniaod_Bu = 43193.0,
     n_miniaod_Bd = 21527.0,
     n_miniaod_Bs = 4378.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.2 GeV, 10.0 mm',
     mass = 3.2,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0195133963751,
     filter_efficiency_Bd = 0.00968274508054,
     filter_efficiency_Bs = 0.0024672067108,
     n_miniaod_Bu = 49133.0,
     n_miniaod_Bd = 24521.0,
     n_miniaod_Bs = 5000.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p25_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.25 GeV, 10.0 mm',
     mass = 3.25,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0199739338032,
     filter_efficiency_Bd = 0.0100663371718,
     filter_efficiency_Bs = 0.00257393269881,
     n_miniaod_Bu = 61444.0,
     n_miniaod_Bd = 31959.0,
     n_miniaod_Bs = 6827.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.25 GeV, 1.0 mm',
     mass = 3.25,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0200413711584,
     filter_efficiency_Bd = 0.0100894540254,
     filter_efficiency_Bs = 0.00261927347804,
     filter_efficiency_Bc = 0.180910673485,
     n_miniaod_Bu = 40710.0,
     n_miniaod_Bd = 21328.0,
     n_miniaod_Bs = 4414.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p3_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.3 GeV, 1.0 mm',
     mass = 3.3,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0199165004695,
     filter_efficiency_Bd = 0.0103399860636,
     filter_efficiency_Bs = 0.00269713781016,
     filter_efficiency_Bc = 0.180257860951,
     n_miniaod_Bu = 48703.0,
     n_miniaod_Bd = 26410.0,
     n_miniaod_Bs = 5548.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.3 GeV, 10.0 mm',
     mass = 3.3,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0198138798139,
     filter_efficiency_Bd = 0.0102746075563,
     filter_efficiency_Bs = 0.00267900986563,
     n_miniaod_Bu = 43187.0,
     n_miniaod_Bd = 23408.0,
     n_miniaod_Bs = 4788.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p35_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.35 GeV, 1.0 mm',
     mass = 3.35,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0195469478188,
     filter_efficiency_Bd = 0.01043789958,
     filter_efficiency_Bs = 0.00269878443181,
     filter_efficiency_Bc = 0.177778355609,
     n_miniaod_Bu = 176935.0,
     n_miniaod_Bd = 97513.0,
     n_miniaod_Bs = 20067.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.35 GeV, 10.0 mm',
     mass = 3.35,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0194258208699,
     filter_efficiency_Bd = 0.0104928403432,
     filter_efficiency_Bs = 0.00270045022011,
     n_miniaod_Bu = 54148.0,
     n_miniaod_Bd = 29715.0,
     n_miniaod_Bs = 6133.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p4_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.4 GeV, 10.0 mm',
     mass = 3.4,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0193971692167,
     filter_efficiency_Bd = 0.010729388702,
     filter_efficiency_Bs = 0.00272362495301,
     n_miniaod_Bu = 54642.0,
     n_miniaod_Bd = 30659.0,
     n_miniaod_Bs = 6404.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.4 GeV, 1.0 mm',
     mass = 3.4,
     ctau = 1.0,
     filter_efficiency_Bu = 0.019723047856,
     filter_efficiency_Bd = 0.0106557776746,
     filter_efficiency_Bs = 0.00271805580626,
     filter_efficiency_Bc = 0.177884615385,
     n_miniaod_Bu = 45684.0,
     n_miniaod_Bd = 25793.0,
     n_miniaod_Bs = 5288.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p45_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.45 GeV, 1.0 mm',
     mass = 3.45,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0194041692361,
     filter_efficiency_Bd = 0.0109283127107,
     filter_efficiency_Bs = 0.00271634960111,
     filter_efficiency_Bc = 0.170892185753,
     n_miniaod_Bu = 33089.0,
     n_miniaod_Bd = 19328.0,
     n_miniaod_Bs = 3957.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.45 GeV, 10.0 mm',
     mass = 3.45,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0194851220076,
     filter_efficiency_Bd = 0.0109468309468,
     filter_efficiency_Bs = 0.0027269987947,
     n_miniaod_Bu = 73286.0,
     n_miniaod_Bd = 42514.0,
     n_miniaod_Bs = 8651.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p5_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.5 GeV, 1.0 mm',
     mass = 3.5,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0192479537107,
     filter_efficiency_Bd = 0.0109944703646,
     filter_efficiency_Bs = 0.00273849792371,
     filter_efficiency_Bc = 0.16882110869,
     n_miniaod_Bu = 37476.0,
     n_miniaod_Bd = 21964.0,
     n_miniaod_Bs = 4515.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.5 GeV, 10.0 mm',
     mass = 3.5,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0194441791615,
     filter_efficiency_Bd = 0.0109390607193,
     filter_efficiency_Bs = 0.0027654892301,
     n_miniaod_Bu = 40246.0,
     n_miniaod_Bd = 23646.0,
     n_miniaod_Bs = 4829.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p55_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.55 GeV, 1.0 mm',
     mass = 3.55,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0193233358424,
     filter_efficiency_Bd = 0.0112114592214,
     filter_efficiency_Bs = 0.00278110950757,
     filter_efficiency_Bc = 0.163474140538,
     n_miniaod_Bu = 72755.0,
     n_miniaod_Bd = 43606.0,
     n_miniaod_Bs = 8993.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.55 GeV, 10.0 mm',
     mass = 3.55,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0192608891639,
     filter_efficiency_Bd = 0.0112400704311,
     filter_efficiency_Bs = 0.0027924474528,
     n_miniaod_Bu = 45891.0,
     n_miniaod_Bd = 27684.0,
     n_miniaod_Bs = 5590.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p6_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.6 GeV, 1.0 mm',
     mass = 3.6,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0191522578987,
     filter_efficiency_Bd = 0.0114710286182,
     filter_efficiency_Bs = 0.00284412323617,
     filter_efficiency_Bc = 0.165737619912,
     n_miniaod_Bu = 156439.0,
     n_miniaod_Bd = 96392.0,
     n_miniaod_Bs = 19653.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.6 GeV, 10.0 mm',
     mass = 3.6,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0189676709952,
     filter_efficiency_Bd = 0.0114319328038,
     filter_efficiency_Bs = 0.0028441769853,
     n_miniaod_Bu = 58117.0,
     n_miniaod_Bd = 35821.0,
     n_miniaod_Bs = 7298.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p65_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.65 GeV, 1.0 mm',
     mass = 3.65,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0189421689422,
     filter_efficiency_Bd = 0.0116580251812,
     filter_efficiency_Bs = 0.00288852325631,
     filter_efficiency_Bc = 0.164650537634,
     n_miniaod_Bu = 49189.0,
     n_miniaod_Bd = 31486.0,
     n_miniaod_Bs = 6260.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.65 GeV, 10.0 mm',
     mass = 3.65,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0190904017857,
     filter_efficiency_Bd = 0.0115989686819,
     filter_efficiency_Bs = 0.00288412985651,
     n_miniaod_Bu = 48926.0,
     n_miniaod_Bd = 31057.0,
     n_miniaod_Bs = 6241.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p7_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.7 GeV, 1.0 mm',
     mass = 3.7,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0188423547557,
     filter_efficiency_Bd = 0.0118532818533,
     filter_efficiency_Bs = 0.00285555318264,
     filter_efficiency_Bc = 0.161351655,
     n_miniaod_Bu = 35650.0,
     n_miniaod_Bd = 22629.0,
     n_miniaod_Bs = 4624.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.7 GeV, 10.0 mm',
     mass = 3.7,
     ctau = 10.0,
     filter_efficiency_Bu = 0.018798349616,
     filter_efficiency_Bd = 0.0118409937074,
     filter_efficiency_Bs = 0.00286289787761,
     n_miniaod_Bu = 44340.0,
     n_miniaod_Bd = 28299.0,
     n_miniaod_Bs = 5797.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p75_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.75 GeV, 1.0 mm',
     mass = 3.75,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0189051875621,
     filter_efficiency_Bd = 0.0119911553923,
     filter_efficiency_Bs = 0.00284818916546,
     filter_efficiency_Bc = 0.155869175627,
     n_miniaod_Bu = 137263.0,
     n_miniaod_Bd = 90065.0,
     n_miniaod_Bs = 18051.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.75 GeV, 10.0 mm',
     mass = 3.75,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0187871981175,
     filter_efficiency_Bd = 0.0120029418275,
     filter_efficiency_Bs = 0.00284706915399,
     n_miniaod_Bu = 47774.0,
     n_miniaod_Bd = 31165.0,
     n_miniaod_Bs = 6222.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p8_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.8 GeV, 1.0 mm',
     mass = 3.8,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0187352739986,
     filter_efficiency_Bd = 0.0122707853303,
     filter_efficiency_Bs = 0.00292524468648,
     filter_efficiency_Bc = 0.157507176436,
     n_miniaod_Bu = 68335.0,
     n_miniaod_Bd = 45573.0,
     n_miniaod_Bs = 9171.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.8 GeV, 10.0 mm',
     mass = 3.8,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0186646854954,
     filter_efficiency_Bd = 0.0122522548469,
     filter_efficiency_Bs = 0.00293427655406,
     n_miniaod_Bu = 46551.0,
     n_miniaod_Bd = 31295.0,
     n_miniaod_Bs = 6185.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p85_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.85 GeV, 1.0 mm',
     mass = 3.85,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0187473185041,
     filter_efficiency_Bd = 0.0123973398451,
     filter_efficiency_Bs = 0.00294007604255,
     filter_efficiency_Bc = 0.156446311936,
     n_miniaod_Bu = 88351.0,
     n_miniaod_Bd = 59419.0,
     n_miniaod_Bs = 12003.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.85 GeV, 10.0 mm',
     mass = 3.85,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0187686410597,
     filter_efficiency_Bd = 0.012383573243,
     filter_efficiency_Bs = 0.00291990363054,
     n_miniaod_Bu = 81556.0,
     n_miniaod_Bd = 54833.0,
     n_miniaod_Bs = 10846.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p9_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.9 GeV, 10.0 mm',
     mass = 3.9,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0188068241001,
     filter_efficiency_Bd = 0.0124350932632,
     filter_efficiency_Bs = 0.00298247164458,
     n_miniaod_Bu = 106513.0,
     n_miniaod_Bd = 73289.0,
     n_miniaod_Bs = 14639.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.9 GeV, 1.0 mm',
     mass = 3.9,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0188783482143,
     filter_efficiency_Bd = 0.0124350932632,
     filter_efficiency_Bs = 0.00298247348775,
     filter_efficiency_Bc = 0.152636530621,
     n_miniaod_Bu = 209597.0,
     n_miniaod_Bd = 143174.0,
     n_miniaod_Bs = 28581.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m3p95_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.95 GeV, 1.0 mm',
     mass = 3.95,
     ctau = 1.0,
     filter_efficiency_Bu = 0.0186972498554,
     filter_efficiency_Bd = 0.0126752692542,
     filter_efficiency_Bs = 0.00299297580776,
     filter_efficiency_Bc = 0.147621070518,
     n_miniaod_Bu = 82571.0,
     n_miniaod_Bd = 57337.0,
     n_miniaod_Bs = 11172.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '3.95 GeV, 10.0 mm',
     mass = 3.95,
     ctau = 10.0,
     filter_efficiency_Bu = 0.0186700388151,
     filter_efficiency_Bd = 0.0126876476415,
     filter_efficiency_Bs = 0.00299359565817,
     n_miniaod_Bu = 122070.0,
     n_miniaod_Bd = 84607.0,
     n_miniaod_Bs = 16535.0,
     muon_rate = 1.,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p0_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.0 GeV, 0.1 mm',
     mass = 4.0,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0188104217821,
     filter_efficiency_Bd = 0.0129517487822,
     filter_efficiency_Bs = 0.00295230158316,
     filter_efficiency_Bc = 0.14425627423,
     n_miniaod_Bu = 136200.0,
     n_miniaod_Bd = 94574.0,
     n_miniaod_Bs = 18242.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p1_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.1 GeV, 0.1 mm',
     mass = 4.1,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0190526215261,
     filter_efficiency_Bd = 0.0133916856858,
     filter_efficiency_Bs = 0.00294565088525,
     filter_efficiency_Bc = 0.139813981398,
     n_miniaod_Bu = 223438.0,
     n_miniaod_Bd = 158480.0,
     n_miniaod_Bs = 29820.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p2_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.2 GeV, 0.1 mm',
     mass = 4.2,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0196614196614,
     filter_efficiency_Bd = 0.0133351915853,
     filter_efficiency_Bs = 0.00294980694981,
     filter_efficiency_Bc = 0.13772753401,
     n_miniaod_Bu = 214510.0,
     n_miniaod_Bd = 153011.0,
     n_miniaod_Bs = 28057.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p3_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.3 GeV, 0.1 mm',
     mass = 4.3,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0202104880414,
     filter_efficiency_Bd = 0.013700761848,
     filter_efficiency_Bs = 0.00290366749968,
     filter_efficiency_Bc = 0.133676767677,
     n_miniaod_Bu = 316460.0,
     n_miniaod_Bd = 222978.0,
     n_miniaod_Bs = 39525.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p4_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.4 GeV, 0.1 mm',
     mass = 4.4,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0201683072372,
     filter_efficiency_Bd = 0.013130124722,
     filter_efficiency_Bs = 0.0025990836764,
     filter_efficiency_Bc = 0.133914222093,
     n_miniaod_Bu = 212441.0,
     n_miniaod_Bd = 146541.0,
     n_miniaod_Bs = 24048.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p5_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.5 GeV, 0.1 mm',
     mass = 4.5,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0208806569981,
     filter_efficiency_Bd = 0.0126426264263,
     filter_efficiency_Bs = 0.00210758721628,
     filter_efficiency_Bc = 0.13245508982,
     n_miniaod_Bu = 222515.0,
     n_miniaod_Bd = 145221.0,
     n_miniaod_Bs = 20792.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p6_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.6 GeV, 0.1 mm',
     mass = 4.6,
     ctau = 0.1,
     filter_efficiency_Bu = 0.0205202769041,
     filter_efficiency_Bd = 0.0117000117,
     filter_efficiency_Bs = 0.0014095504338,
     filter_efficiency_Bc = 0.13116617413,
     n_miniaod_Bu = 146480.0,
     n_miniaod_Bd = 88107.0,
     n_miniaod_Bs = 9668.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p7_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.7 GeV, 0.1 mm',
     mass = 4.7,
     ctau = 0.1,
     filter_efficiency_Bu = 0.019966345678,
     filter_efficiency_Bd = 0.0099225610291,
     filter_efficiency_Bs = 0.000866038490439,
     filter_efficiency_Bc = 0.132296966271,
     n_miniaod_Bu = 266641.0,
     n_miniaod_Bd = 142966.0,
     n_miniaod_Bs = 9968.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p8_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.8 GeV, 0.1 mm',
     mass = 4.8,
     ctau = 0.1,
     filter_efficiency_Bu = -99,
     filter_efficiency_Bd = -99,
     filter_efficiency_Bs = -99,
     filter_efficiency_Bc = 0.132886032886,
     n_miniaod_Bu = 212869.0,
     n_miniaod_Bd = 94490.0,
     n_miniaod_Bs = 0.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m4p9_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '4.9 GeV, 0.1 mm',
     mass = 4.9,
     ctau = 0.1,
     filter_efficiency_Bu = -99,
     filter_efficiency_Bd = -99,
     filter_efficiency_Bs = -99,
     filter_efficiency_Bc = 0.133171960099,
     n_miniaod_Bu = 253967.0,
     n_miniaod_Bd = 87127.0,
     n_miniaod_Bs = 0.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p0_norm'] = [
 SignalSample(
     filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass5p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.0 GeV, 0.1 mm',
     mass = 5.0,
     ctau = 0.1,
     filter_efficiency_Bu = -99,
     filter_efficiency_Bd = -99,
     filter_efficiency_Bs = -99,
     filter_efficiency_Bc = 0.137866799461,
     n_miniaod_Bu = 304452.0,
     n_miniaod_Bd = 50739.0,
     n_miniaod_Bs = 0.0,
     muon_rate = 1.,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p1_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.1 GeV, 0.1 mm',
     mass = 5.1,
     ctau = 0.1,
     filter_efficiency_Bc = 0.1397320904,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p2_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.2 GeV, 0.1 mm',
     mass = 5.2,
     ctau = 0.1,
     filter_efficiency_Bc = 0.141518235793,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p3_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.3 GeV, 0.1 mm',
     mass = 5.3,
     ctau = 0.1,
     filter_efficiency_Bc = 0.143748568812,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p4_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.4 GeV, 0.1 mm',
     mass = 5.4,
     ctau = 0.1,
     filter_efficiency_Bc = 0.143049503469,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p5_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.5 GeV, 0.1 mm',
     mass = 5.5,
     ctau = 0.1,
     filter_efficiency_Bc = 0.143539647705,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p6_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.6 GeV, 0.1 mm',
     mass = 5.6,
     ctau = 0.1,
     filter_efficiency_Bc = 0.142474758409,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p7_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.7 GeV, 0.1 mm',
     mass = 5.7,
     ctau = 0.1,
     filter_efficiency_Bc = 0.133590971769,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p8_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.8 GeV, 0.1 mm',
     mass = 5.8,
     ctau = 0.1,
     filter_efficiency_Bc = 0.1245474179,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m5p9_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '5.9 GeV, 0.1 mm',
     mass = 5.9,
     ctau = 0.1,
     filter_efficiency_Bc = 0.103264414759,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m6p0_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '6.0 GeV, 0.1 mm',
     mass = 6.0,
     ctau = 0.1,
     filter_efficiency_Bc = 0.0704549817623,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]

signal_samples['V42_06Feb23_m6p1_norm'] = [
 SignalSample(
     filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_15Jun23.root',
     label = '6.1 GeV, 0.1 mm',
     mass = 6.1,
     ctau = 0.1,
     filter_efficiency_Bc = 0.0442446136992,
     muon_rate_Bc = 0.47,
     colour = ROOT.kOrange+0,
     ),
  ]




signal_samples['V42_06Feb23_m1p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.0 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 8.20e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p02'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p02_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.02 GeV, 10 mm',
      mass = 1.02,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 8.02e-03, #7.75e-03, #5.45e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p04'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p04_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.04 GeV, 10 mm',
      mass = 1.04,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.97e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p06'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p06_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.06 GeV, 10 mm',
      mass = 1.06,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.81e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p08'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p08_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.08 GeV, 10 mm',
      mass = 1.08,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.95e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.1 GeV, 10 mm',
      mass = 1.1,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.88e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p12'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p12_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.12 GeV, 10 mm',
      mass = 1.12,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.80e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p14'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p14_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.14 GeV, 10 mm',
      mass = 1.14,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.90e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p16'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p16_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.16 GeV, 10 mm',
      mass = 1.16,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.81e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p18'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p18_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.18 GeV, 10 mm',
      mass = 1.18,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.83e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.2 GeV, 10 mm',
      mass = 1.2,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.84e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p22'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p22_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.22 GeV, 10 mm',
      mass = 1.22,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.78e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p24'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p24_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.24 GeV, 10 mm',
      mass = 1.24,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.79e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p26'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p26_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.26 GeV, 10 mm',
      mass = 1.26,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.83e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p28'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p28_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.28 GeV, 10 mm',
      mass = 1.28,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.67e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.3 GeV, 10 mm',
      mass = 1.3,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.54e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p32'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p32_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.32 GeV, 10 mm',
      mass = 1.32,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.53e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p34'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p34_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.34 GeV, 10 mm',
      mass = 1.34,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.54e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p36'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p36_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.36 GeV, 10 mm',
      mass = 1.36,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.56e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p38'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p38_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.38 GeV, 10 mm',
      mass = 1.38,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.56e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.4 GeV, 10 mm',
      mass = 1.4,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.52e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p42'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p42_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.42 GeV, 10 mm',
      mass = 1.42,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.41e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p44'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p44_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.44 GeV, 10 mm',
      mass = 1.44,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.40e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p46'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p46_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.46 GeV, 10 mm',
      mass = 1.46,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.29e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p48'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p48_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.48 GeV, 10 mm',
      mass = 1.48,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.27e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.24e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p53'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p53_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.53 GeV, 10 mm',
      mass = 1.53,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.24e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p56'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p56_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.56 GeV, 10 mm',
      mass = 1.56,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.26e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p59'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p59_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.59 GeV, 10 mm',
      mass = 1.59,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.16e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p62'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p62_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.62 GeV, 10 mm',
      mass = 1.62,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.17e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.65 GeV, 10 mm',
      mass = 1.65,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.10e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p68'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p68_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.68 GeV, 10 mm',
      mass = 1.68,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.08e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p71'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p71_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.71 GeV, 10 mm',
      mass = 1.71,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.94e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p74'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p74_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.74 GeV, 10 mm',
      mass = 1.74,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.90e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p77'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p77_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.77 GeV, 10 mm',
      mass = 1.77,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.76e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.8 GeV, 10 mm',
      mass = 1.8,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.86e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p83'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p83_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.83 GeV, 10 mm',
      mass = 1.83,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.82e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p86'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p86_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.86 GeV, 10 mm',
      mass = 1.86,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.84e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p89'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p89_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.89 GeV, 10 mm',
      mass = 1.89,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.80e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p92'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p92_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.92 GeV, 10 mm',
      mass = 1.92,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.76e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.95 GeV, 10 mm',
      mass = 1.95,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.76e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m1p98'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p98_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '1.98 GeV, 10 mm',
      mass = 1.98,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.79e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.0 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.82e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p05'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.05 GeV, 10 mm',
      mass = 2.05,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.80e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.1 GeV, 10 mm',
      mass = 2.1,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.96e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p15'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.15 GeV, 10 mm',
      mass = 2.15,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.83e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.2 GeV, 10 mm',
      mass = 2.2,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.89e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p25'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.25 GeV, 10 mm',
      mass = 2.25,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 6.99e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.3 GeV, 10 mm',
      mass = 2.3,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.16e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p35'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.35 GeV, 10 mm',
      mass = 2.35,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.46e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.4 GeV, 10 mm',
      mass = 2.4,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.57e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p45'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.45 GeV, 10 mm',
      mass = 2.45,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.74e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.5 GeV, 10 mm',
      mass = 2.5,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 7.94e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p55'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.55 GeV, 10 mm',
      mass = 2.55,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 8.31e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.6 GeV, 10 mm',
      mass = 2.6,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 8.60e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.65 GeV, 10 mm',
      mass = 2.65,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 8.94e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.7 GeV, 10 mm',
      mass = 2.7,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 9.43e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p75'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.75 GeV, 10 mm',
      mass = 2.75,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 9.92e-03,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.8 GeV, 10 mm',
      mass = 2.8,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.08e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p85'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.85 GeV, 10 mm',
      mass = 2.85,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.16e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p9'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.9 GeV, 10 mm',
      mass = 2.9,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m2p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '2.95 GeV, 10 mm',
      mass = 2.95,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.39e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.0 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.55e-02, #1.48e-02,#1.06e-02,
      filter_efficiency_Bc = 1.87e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.0 GeV, 10 mm',
      mass = 3.0,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.55e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p05'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.05 GeV, 1 mm',
      mass = 3.05,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.77e-02,
      filter_efficiency_Bc = 1.4e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.05 GeV, 10 mm',
      mass = 3.05,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.77e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.1 GeV, 1 mm',
      mass = 3.1,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.96e-02,
      filter_efficiency_Bc = 1.88e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.1 GeV, 10 mm',
      mass = 3.1,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 1.96e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p15'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.15 GeV, 1 mm',
      mass = 3.15,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.06e-02,
      filter_efficiency_Bc = 1.87e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.15 GeV, 10 mm',
      mass = 3.15,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.06e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.2 GeV, 1 mm',
      mass = 3.2,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.13e-02,
      filter_efficiency_Bc = 1.84e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.2 GeV, 10 mm',
      mass = 3.2,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.13e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p25'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.25 GeV, 1 mm',
      mass = 3.25,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.21e-02,
      filter_efficiency_Bc = 1.81e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.25 GeV, 10 mm',
      mass = 3.25,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.21e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.3 GeV, 1 mm',
      mass = 3.3,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      filter_efficiency_Bc = 1.80e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.3 GeV, 10 mm',
      mass = 3.3,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p35'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.35 GeV, 1 mm',
      mass = 3.35,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.22e-02,
      filter_efficiency_Bc = 1.78e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.35 GeV, 10 mm',
      mass = 3.35,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.22e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.4 GeV, 1 mm',
      mass = 3.4,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.24e-02,
      filter_efficiency_Bc = 1.78e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.4 GeV, 10 mm',
      mass = 3.4,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.24e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p45'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.45 GeV, 1 mm',
      mass = 3.45,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      filter_efficiency_Bc = 1.71e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.45 GeV, 10 mm',
      mass = 3.45,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.5 GeV, 1 mm',
      mass = 3.5,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      filter_efficiency_Bc = 1.69e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.5 GeV, 10 mm',
      mass = 3.5,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p55'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.55 GeV, 1 mm',
      mass = 3.55,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      filter_efficiency_Bc = 1.63e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.55 GeV, 10 mm',
      mass = 3.55,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.6 GeV, 1 mm',
      mass = 3.6,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      filter_efficiency_Bc = 1.66e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.6 GeV, 10 mm',
      mass = 3.6,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.65 GeV, 1 mm',
      mass = 3.65,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      filter_efficiency_Bc = 1.65e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.65 GeV, 10 mm',
      mass = 3.65,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.7 GeV, 1 mm',
      mass = 3.7,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.28e-02,
      filter_efficiency_Bc = 1.61e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.7 GeV, 10 mm',
      mass = 3.7,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.28e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p75'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.75 GeV, 1 mm',
      mass = 3.75,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.56e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.75 GeV, 10 mm',
      mass = 3.75,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.8 GeV, 1 mm',
      mass = 3.8,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.27e-02,
      filter_efficiency_Bc = 1.58e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.8 GeV, 10 mm',
      mass = 3.8,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.27e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p85'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.85 GeV, 1 mm',
      mass = 3.85,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.31e-02,
      filter_efficiency_Bc = 1.56e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.85 GeV, 10 mm',
      mass = 3.85,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.31e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p9'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.9 GeV, 1 mm',
      mass = 3.9,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.34e-02,
      filter_efficiency_Bc = 1.53e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.9 GeV, 10 mm',
      mass = 3.9,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.34e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m3p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.95 GeV, 1 mm',
      mass = 3.95,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.35e-02,
      filter_efficiency_Bc = 1.48e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '3.95 GeV, 10 mm',
      mass = 3.95,
      ctau = 10.0,
      resolution = None,
      filter_efficiency = 2.35e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.0 GeV, 0.1 mm',
      mass = 4.0,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.37e-02,
      filter_efficiency_Bc = 1.44e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.1 GeV, 0.1 mm',
      mass = 4.1,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.42e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.2 GeV, 0.1 mm',
      mass = 4.2,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.47e-02,
      filter_efficiency_Bc = 1.38e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.3 GeV, 0.1 mm',
      mass = 4.3,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.49e-02,
      filter_efficiency_Bc = 1.34e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.4 GeV, 0.1 mm',
      mass = 4.4,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.55e-02,
      filter_efficiency_Bc = 1.34e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.55e-02, #2.38e-02,#1.75e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.6 GeV, 0.1 mm',
      mass = 4.6,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.47e-02,
      filter_efficiency_Bc = 1.31e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.7 GeV, 0.1 mm',
      mass = 4.7,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.28e-02,
      filter_efficiency_Bc = 1.32e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m4p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.8 GeV, 0.1 mm',
      mass = 4.8,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 1.46e-02,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]


signal_samples['V42_06Feb23_m4p9'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '4.9 GeV, 0.1 mm',
      mass = 4.9,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.33e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p0'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.0 GeV, 0.1 mm',
      mass = 5.0,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.38e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p1'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.1 GeV, 0.1 mm',
      mass = 5.1,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p2'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.2 GeV, 0.1 mm',
      mass = 5.2,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.42e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p3'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.3 GeV, 0.1 mm',
      mass = 5.3,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.44e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p4'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.4 GeV, 0.1 mm',
      mass = 5.4,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.43e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p5'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.5 GeV, 0.1 mm',
      mass = 5.5,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.44e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p6'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.6 GeV, 0.1 mm',
      mass = 5.6,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.42e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p7'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.7 GeV, 0.1 mm',
      mass = 5.7,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.34e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p8'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.8 GeV, 0.1 mm',
      mass = 5.8,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.25e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m5p9'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '5.9 GeV, 0.1 mm',
      mass = 5.9,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 1.03e-01,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m6p0'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '6.0 GeV, 0.1 mm',
      mass = 6.0,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 7.05e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_06Feb23_m6p1'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_norm.root',
      label = '6.1 GeV, 0.1 mm',
      mass = 6.1,
      ctau = 0.1,
      resolution = None,
      filter_efficiency_Bc = 4.42e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]


signal_samples['central_V11_24Apr22_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 1000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['central_V11_24Apr22_m1p5'] = [
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
  #    label = '1.5 GeV, 1000 mm',
  #    mass = 1.5,
  #    ctau = 1000.0,
  #    resolution = 0.01285,
  #    filter_efficiency = 1.66e-03,
  #    muon_rate = 0.45,
  #    is_private = False,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.01265,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m2'] = [
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
  #    label = '2 GeV, 1000 mm',
  #    mass = 2.0,
  #    ctau = 1000.0,
  #    resolution = 0.0167,
  #    filter_efficiency = 1.66e-03,
  #    muon_rate = 0.44,
  #    is_private = False,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kRed-9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+1
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kRed+2
      ),
  ]


signal_samples['V10_30Dec21_benchmark'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 1mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 10^{3}mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1GeV, 1000mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1GeV, 100mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1GeV, 10mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V10_30Dec21_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1.5GeV, 1000mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1.5GeV, 100mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '1.5GeV, 10mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '2GeV, 1000mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04, # is this one correct?
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '2GeV, 100mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '2GeV, 10mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '3GeV, 1000mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed-9,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '3GeV, 10mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '3GeV, 1mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed+2,
      ),
  #SignalSample( # reweighted
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    label = '3GeV, 1mm (central, V10_30Dec21)',
  #    mass = 3,
  #    ctau = 0.1,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private= False,
  #    colour = ROOT.kRed+3,
  #    ),
  ]

signal_samples['V10_30Dec21_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '4.5GeV, 100mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '4.5GeV, 10mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed-5,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '4.5GeV, 1mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21_sr.root',
      label = '4.5GeV, 1mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kRed-2,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    label = '4.5GeV, 0.1mm (central, V10_30Dec21)',
  #    mass = 4.5,
  #    ctau = 0.1,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.62e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  ]

signal_samples['central_benchmark_looseselection'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0373, # to be checked
      filter_efficiency = 5.53e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3 GeV, 100 mm',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0243, # to be checked
      filter_efficiency = 5.87e-03,
      filter_efficiency_Bc = 1.54e-01, 
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n38/merged/bparknano.root',
      label = '1 GeV, 1000 mm',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.0104, # to be checked
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V00_looseselection'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n38/merged/bparknano.root',
      label = '1GeV, 1000mm (central, loose)',
      mass = 1,
      ctau = 1000.0,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '1GeV, 100mm (central, loose)',
      mass = 1,
      ctau = 100.0,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '1GeV, 10mm (central, loose)',
      mass = 1,
      ctau = 10.0,
      colour = ROOT.kOrange+9
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '1.5GeV, 1000mm (central, loose)',
      mass = 1.5,
      ctau = 1000.0,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '1.5GeV, 100mm (central, loose)',
      mass = 1.5,
      ctau = 100.0,
      colour = ROOT.kYellow-2
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '1.5GeV, 10mm (central, loose)',
      mass = 1.5,
      ctau = 10.0,
      colour = ROOT.kYellow+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '2GeV, 1000mm (central, loose)',
      mass = 2,
      ctau = 1000.0,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '2GeV, 100mm (central, loose)',
      mass = 2,
      ctau = 100.0,
      colour = ROOT.kBlue-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '2GeV, 10mm (central, loose)',
      mass = 2,
      ctau = 10.0,
      colour = ROOT.kBlue+3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '3GeV, 1000mm (central, loose)',
      mass = 3,
      ctau = 1000.0,
      colour = ROOT.kRed-9,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '3GeV, 100mm (central, loose)',
      mass = 3,
      ctau = 100.0,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '3GeV, 10mm (central, loose)',
      mass = 3,
      ctau = 10.0,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '3GeV, 1mm (central, loose)',
      mass = 3,
      ctau = 1.0,
      colour = ROOT.kRed+2,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '4.5GeV, 100mm (central, loose)',
      mass = 4.5,
      ctau = 100.0,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '4.5GeV, 10mm (central, loose)',
      mass = 4.5,
      ctau = 10.0,
      colour = ROOT.kRed-5,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '4.5GeV, 1mm (central, loose)',
      mass = 4.5,
      ctau = 1.0,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      label = '4.5GeV, 0.1mm (central, loose)',
      mass = 4.5,
      ctau = 0.1,
      colour = ROOT.kRed-2,
      ),
  ]

signal_samples['V00_looseselection_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '1GeV, 1000mm (central, loose)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '1GeV, 100mm (central, loose)',
      mass = 1,
      ctau = 100.0,
      filter_efficiency = 5.17e-03,
      resolution = 0.00853,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '1GeV, 10mm (central, loose)',
      mass = 1,
      ctau = 10.0,
      filter_efficiency = 5.85e-03,
      resolution = 0.00853,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V00_looseselection_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '3GeV, 1000mm (central, loose)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed-9,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 70.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 50.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 40.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 30.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 20.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 15.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '3GeV, 100mm (central, loose)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 7.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 5.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 4.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 3.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 2.0,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 1.5,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '3GeV, 10mm (central, loose)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.7,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.5,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.4,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.3,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.2,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      mass = 3,
      ctau = 0.15,
      resolution = 0.0251,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '3GeV, 1mm (central, loose)',
      mass = 3,
      ctau = 0.1,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+2,
      ),
  ]

signal_samples['V00_looseselection_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '4.5GeV, 100mm (central, loose)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 4.58e-04,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '4.5GeV, 10mm (central, loose)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 5.53e-04,
      colour = ROOT.kRed-5,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '4.5GeV, 1mm (central, loose)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 5.53e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano.root',
      label = '4.5GeV, 0.1mm (central, loose)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 5.53e-04,
      colour = ROOT.kRed-2,
      ),
  ]

signal_samples['tag_and_probe'] = [
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


signal_samples['tag_and_probe_BToJPsiKstar'] = [
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V0/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_signal_sf_study15Sep21_A1_v1.root',
  #    label = 'BToJPsiKstar (V0), B->#mu#mu#pi',
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V0/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_sf_study15Sep21_A1_B1_v1.root',
  #    label = 'BToJPsiKstar (V0), JPsi->#mu#mu',
  #    colour = ROOT.kOrange+0,
  #    ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V12_08Aug22/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V12_08Aug22/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_smearing.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V12_08Aug22/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_mcweights.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V12_08Aug22/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_scale.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V12_08Aug22/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_08Aug22_linearscale.root',
      label = 'BToJPsiKstar (V12_08Aug22)',
      colour = ROOT.kOrange+0,
      ),
]

signal_samples['V15_control_06Feb23'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_06Feb23.root',
      colour = ROOT.kRed+1,
      label = 'control',
      ),
]


data_samples = {}
data_samples['V10_30Dec21'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  ]

data_samples['V10_30Dec21_small'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018D (V10_30Dec21)',
      lumi = 5.302, #0.075, # approximate, 5.302 / 70
      ),
  ]

data_samples['V10_30Dec21_fullBPark'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018A/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH2_Run2018A/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH2_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH3_Run2018A/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH3_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH4_Run2018A/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH4_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH5_Run2018A/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH5_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH6_Run2018A/merged/flat_bparknano_30Dec21_partial.root',
      label = 'ParkingBPH6_Run2018A (V10_30Dec21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018B (V10_30Dec21)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH2_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH2_Run2018B (V10_30Dec21)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH3_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH3_Run2018B (V10_30Dec21)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH4_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH4_Run2018B (V10_30Dec21)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH5_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH5_Run2018B (V10_30Dec21)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH6_Run2018B/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH6_Run2018B (V10_30Dec21)',
      lumi = 0.377,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018C/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018C (V10_30Dec21)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH2_Run2018C/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH2_Run2018C (V10_30Dec21)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH3_Run2018C/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH3_Run2018C (V10_30Dec21)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH4_Run2018C/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH4_Run2018C (V10_30Dec21)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH5_Run2018C/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH5_Run2018C (V10_30Dec21)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH2_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH2_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH3_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH3_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH4_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH4_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH5_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH5_Run2018D (V10_30Dec21)',
      lumi = 5.302,
      ),
  ]

data_samples['V11_24Apr22'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/merged/flat_bparknano_24Apr22.root',
      label = 'ParkingBPH1_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  ]

data_samples['V11_24Apr22_small'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_24Apr22.root',
      label = 'ParkingBPH1_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  ]

data_samples['V11_24Apr22_fullBPark'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH1_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH2_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH2_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH3_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH3_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH4_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH4_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH5_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH5_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH6_Run2018A/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH6_Run2018A (V11_24Apr22)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH1_Run2018B (V11_24Apr22)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH2_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH2_Run2018B (V11_24Apr22)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH3_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH3_Run2018B (V11_24Apr22)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH4_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH4_Run2018B (V11_24Apr22)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH5_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH5_Run2018B (V11_24Apr22)',
      lumi = 0.991,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH6_Run2018B/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH6_Run2018B (V11_24Apr22)',
      lumi = 0.377,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018C/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH1_Run2018C (V11_24Apr22)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH2_Run2018C/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH2_Run2018C (V11_24Apr22)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH3_Run2018C/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH3_Run2018C (V11_24Apr22)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH4_Run2018C/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH4_Run2018C (V11_24Apr22)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH5_Run2018C/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH5_Run2018C (V11_24Apr22)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/merged/flat_bparknano_24Apr22.root',
      label = 'ParkingBPH1_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH2_Run2018D/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH2_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH3_Run2018D/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH3_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH4_Run2018D/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH4_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH5_Run2018D/merged/flat_bparknano_24Apr22_partial.root',
      label = 'ParkingBPH5_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  ]

data_samples['V12_08Aug22'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_sr.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302,
      ),
  ]

data_samples['V12_08Aug22_CR'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302,
      ),
  ]

data_samples['V12_08Aug22_small'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_sr.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_nj1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk1_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk2_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk3_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk4_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk5_n500/flat/flat_bparknano_08Aug22.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      lumi = 5.302/75.,
      ),
  ]

data_samples['V13_06Feb23'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH1_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  ]

data_samples['V13_06Feb23_small'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH1_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  ]

data_samples['V13_06Feb23_fullBPark'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH1_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH2_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH3_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH4_Run2018D (V13_06Feb23)',
      lumi = 5.302 * 0.9894,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018D/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH5_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018C/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH1_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018C/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH2_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018C/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH3_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018C/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH4_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018C/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH5_Run2018C (V13_06Feb23)',
      lumi = 1.103 * 0.9997,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH1_Run2018B (V13_06Feb23)',
      lumi = 0.911 * 0.9997,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH2_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH3_Run2018B (V13_06Feb23)',
      lumi = 0.911, # * 0.0122,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH4_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH5_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH6_Run2018B/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH6_Run2018B (V13_06Feb23)',
      lumi = 0.377,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH1_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH2_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH3_Run2018A (V13_06Feb23)',
      lumi = 0.774,# * 0.3259,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH4_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH5_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH6_Run2018A/merged/flat_bparknano_06Feb23_15Jun23.root',
      label = 'ParkingBPH6_Run2018A (V13_06Feb23)',
      lumi = 0.774 *  0.9968, # * 0.9843,
      ),
  ]

data_samples['V13_06Feb23_control'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_control.root',
      label = 'ParkingBPH1_Run2018A',
      lumi = 0.774,
      ),
  ]

data_samples['track_id'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_06Feb23_norm.root',
      label = 'ParkingBPH1_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  ]

data_samples['triggermuon_matching_check'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/merged/flat_bparknano_for_triggermuon_matching_study.root',
      label = 'Relaxed trigger muon matching',
      ),
  ]

data_samples['small'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH1_Run2018A/Chunk0_n500/flat/flat_bparknano_selected_stdtrgmu_full_partial.root',
      label = 'ParkingBPH1_Run2018A (V03)',
      ),
  ]

data_samples['loose'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH4_Run2018B/Chunk0_n10/bparknano_loosepreselection_v2_nj1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/bparknano_data_1file_looseselection.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/merged/bparknano_loosepreselection_1file.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_data_loosedr_1file.root',
      label = 'background',
      ),
  ]

data_samples['tag_and_probe'] = [
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V06_tag_and_probe/ParkingBPH1_Run2018A/merged/flat_bparknano_tag_and_probe_v2_tag_fired_DST_DoubleMu1.root',
  #    label = 'ParkingBPH1_Run2018A (V06)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V06_tag_and_probe/ParkingBPH1_Run2018B/merged/flat_bparknano_tag_and_probe_v2.root',
  #    label = 'ParkingBPH1_Run2018B (V06)',
  #    ),
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/flat_bparknano_08Aug22_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/flat/flat_bparknano_08Aug22_smearing.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/splots_Chunk0_smearing.root',
      label = 'ParkingBPH1_Run2018D (V12_08Aug22)',
      ),
  ]

data_samples['control_Bc'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH1_Run2018D/merged/flat_bparknano_24Apr22_control_Bc.root',
      label = 'ParkingBPH1_Run2018D (V11_24Apr22)',
      lumi = 5.302,
      ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH2_Run2018D/merged/flat_bparknano_24Apr22_control_Bc.root',
  #    label = 'ParkingBPH2_Run2018D (V11_24Apr22)',
  #    lumi = 5.302,
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH3_Run2018D/merged/flat_bparknano_24Apr22_control_Bc.root',
  #    label = 'ParkingBPH3_Run2018D (V11_24Apr22)',
  #    lumi = 5.302,
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH4_Run2018D/merged/flat_bparknano_24Apr22_control_Bc.root',
  #    label = 'ParkingBPH4_Run2018D (V11_24Apr22)',
  #    lumi = 5.302,
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V11_24Apr22/ParkingBPH5_Run2018D/merged/flat_bparknano_24Apr22_control_Bc_partial.root',
  #    label = 'ParkingBPH5_Run2018D (V11_24Apr22)',
  #    lumi = 5.302,
  #    ),
  ]


qcd_samples = {}
qcd_samples['V13_06Feb23'] = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Feb23.root',
      label = 'QCD_pt15to20 (V13_06Feb23)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Feb23.root',
      label = 'QCD_pt20to30 (V13_06Feb23)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Feb23.root',
      label = 'QCD_pt30to50 (V13_06Feb23)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Feb23.root',
      label = 'QCD_pt50to80 (V13_06Feb23)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_06Feb23_extmerged.root',
      label = 'QCD_pt80to120 (V13_06Feb23)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_06Feb23_extmerged.root',
      label = 'QCD_pt120to170 (V13_06Feb23)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Feb23.root',
      label = 'QCD_pt170to300 (V13_06Feb23)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V12_08Aug22'] = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_sf_v1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt15to20 (V12_08Aug22)',
      #label = 'QCD_pt15to20',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt20to30 (V12_08Aug22)',
      #label = 'QCD_pt20to30',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk1_n9/flat/flat_bparknano_08Aug22_sf_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_v4.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt30to50 (V12_08Aug22)',
      #label = 'QCD_pt30to50',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt50to80 (V12_08Aug22)',
      #label = 'QCD_pt50to80',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt80to120 (V12_08Aug22)',
      #label = 'QCD_pt80to120',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt120to170 (V12_08Aug22)',
      #label = 'QCD_pt120to170',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_mccorrections_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V12_08Aug22/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_08Aug22_withhlt.root',
      label = 'QCD_pt170to300 (V12_08Aug22)',
      #label = 'QCD_pt170to300',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V11_24Apr22'] = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt15to20 (V11_24Apr22)',
      #label = 'QCD_pt15to20',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt20to30 (V11_24Apr22)',
      #label = 'QCD_pt20to30',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt30to50 (V11_24Apr22)',
      #label = 'QCD_pt30to50',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt50to80 (V11_24Apr22)',
      #label = 'QCD_pt50to80',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt80to120 (V11_24Apr22)',
      #label = 'QCD_pt80to120',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = 'QCD_pt120to170 (V11_24Apr22)',
      #label = 'QCD_pt120to170',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sf_v3.root',
  #    label = 'QCD_pt170to300 (V11_24Apr22)',
  #    #label = 'QCD_pt170to300',
  #    cross_section = 117989,
  #    filter_efficiency = 0.07335,
  #    colour = ROOT.kOrange+8, 
  #    ),
  ]

qcd_samples['V11_24Apr22_sources'] = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt15to20 (V11_24Apr22)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt20to30 (V11_24Apr22)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sources.root',
      #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano_pt30to50.root',
      label = 'QCD_pt30to50 (V11_24Apr22)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt50to80 (V11_24Apr22)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt80to120 (V11_24Apr22)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt120to170 (V11_24Apr22)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V11_24Apr22/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_24Apr22_sources.root',
      label = 'QCD_pt170to300 (V11_24Apr22)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V10_30Dec21'] = [
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
  #    label = 'QCD_pt15to20 (V10_30Dec21)',
  #    cross_section = 1.27319e+09,
  #    filter_efficiency = 0.00300,
  #    colour = ROOT.kBlue+2, 
  #    ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_sf_v2.root',
      label = 'QCD_pt20to30 (V10_30Dec21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_sf_v2.root',
      label = 'QCD_pt30to50 (V10_30Dec21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_sf_v2.root',
      label = 'QCD_pt50to80 (V10_30Dec21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_v2_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_sf_v2_extmerged.root',
      label = 'QCD_pt80to120 (V10_30Dec21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_v2_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_sf_v2_extmerged.root',
      label = 'QCD_pt120to170 (V10_30Dec21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  #QCDMCSample(
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_v2.root',
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21_sf_v2.root',
  #    label = 'QCD_pt170to300 (V10_30Dec21)',
  #    cross_section = 117989,
  #    filter_efficiency = 0.07335,
  #    colour = ROOT.kOrange+8, 
  #    ),
  ]

qcd_samples['V10_30Dec21_sources'] = [
  QCDMCSample(
      filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano_pt30to50.root',
      #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano.root',
      label = 'QCD_pt30to50 (V10_30Dec21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/plugins/dumper/flat_bparknano_pt170to300.root',
      label = 'QCD_pt170to300 (V10_30Dec21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]


