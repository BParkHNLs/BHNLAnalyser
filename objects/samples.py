import ROOT

class DataSample(object):
  def __init__(self, filename, label, lumi=0.074):
    self.filename = filename
    self.label = label
    self.lumi = lumi


class SignalSample(object): 
  def __init__(self, filename=None, filename_Bc=None, label=None, mass=None, ctau=None, resolution=None, filter_efficiency=None, filter_efficiency_Bc=None, muon_rate=None, is_private=None, colour=ROOT.kRed):
    self.filename = filename
    self.filename_Bc = filename_Bc
    self.label = label
    self.mass = mass
    self.ctau = ctau # this is the target ctau, not necessarily the ctau with which the sample was produced
    self.resolution = resolution
    self.filter_efficiency = filter_efficiency
    self.filter_efficiency_Bc = filter_efficiency_Bc
    self.muon_rate = muon_rate
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


signal_samples['V13_06Feb23_preselection'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root', #FIXME take version with gen-matching
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN5p5_ctau0p01mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root', #FIXME take version with gen-matching
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0m_TuneCP5_13TeV_pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V13_06Feb23/BcToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Feb23_partial.root',
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


signal_samples['V42_06Feb23_m1p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p02_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p04_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p06_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p08_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p12_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m1p14'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p14_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '1.14 GeV, 10 mm',
#      mass = 1.14,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.90e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_06Feb23_m1p16'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p16_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p18_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p22_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p24_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p26_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m1p28'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p28_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '1.28 GeV, 10 mm',
#      mass = 1.28,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.67e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_06Feb23_m1p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p32_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p34_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p36_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p38_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p42_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p44_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p46_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p48_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p53_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p56_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p59_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m1p62'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p62_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '1.62 GeV, 10 mm',
#      mass = 1.62,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.17e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_06Feb23_m1p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m1p68'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p68_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '1.68 GeV, 10 mm',
#      mass = 1.68,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.08e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_06Feb23_m1p71'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p71_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p74_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p77_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p83_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p86_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p89_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p92_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p98_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m4p3'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '4.3 GeV, 0.1 mm',
#      mass = 4.3,
#      ctau = 0.1,
#      resolution = None,
#      filter_efficiency = 2.49e-02,
#      filter_efficiency_Bc = 1.34e-01,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_06Feb23_m4p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_06Feb23_m4p8'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
#      label = '4.8 GeV, 0.1 mm',
#      mass = 4.8,
#      ctau = 0.1,
#      resolution = None,
#      filter_efficiency = 1.46e-02,
#      filter_efficiency_Bc = 1.33e-01,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]


signal_samples['V42_06Feb23_m4p9'] = [
  SignalSample(
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass4p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p2_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p3_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p4_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p5_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p6_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p7_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p8_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass5p9_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p0_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42_Bc/mass6p1_ctau0p1/nanoFiles/merged/flat_bparknano_06Feb23_v2.root',
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

#signal_samples['V42_08Aug22_m0p5'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass0p5_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '0.5 GeV, 10 mm',
#      mass = 0.5,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = None,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

#signal_samples['V42_08Aug22_m0p6'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass0p6_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '0.6 GeV, 10 mm',
#      mass = 0.6,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 8.55e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]
#
#signal_samples['V42_08Aug22_m0p7'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass0p7_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '0.7 GeV, 10 mm',
#      mass = 0.7,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 8.36e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]
#
#signal_samples['V42_08Aug22_m0p8'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass0p8_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '0.8 GeV, 10 mm',
#      mass = 0.8,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 8.44e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]
#
#signal_samples['V42_08Aug22_m0p9'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass0p9_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '0.9 GeV, 10 mm',
#      mass = 0.9,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 8.28e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m1p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p0_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p02'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p02_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p04'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p04_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p06'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p06_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p08'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p08_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p1_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p12'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p12_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

#signal_samples['V42_08Aug22_m1p14'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p14_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '1.14 GeV, 10 mm',
#      mass = 1.14,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.90e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m1p16'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p16_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p18'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p18_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p2_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p22'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p22_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p24'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p24_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p26'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p26_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

#signal_samples['V42_08Aug22_m1p28'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p28_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '1.28 GeV, 10 mm',
#      mass = 1.28,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.67e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m1p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p3_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p32'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p32_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p34'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p34_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p36'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p36_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p38'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p38_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p4_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p42'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p42_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p44'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p44_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p46'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p46_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p48'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p48_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p5_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p53'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p53_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p56'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p56_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p59'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p59_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

#signal_samples['V42_08Aug22_m1p62'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p62_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '1.62 GeV, 10 mm',
#      mass = 1.62,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.17e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m1p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p65_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

#signal_samples['V42_08Aug22_m1p68'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p68_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '1.68 GeV, 10 mm',
#      mass = 1.68,
#      ctau = 10.0,
#      resolution = None,
#      filter_efficiency = 7.08e-03,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m1p71'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p71_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p74'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p74_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p77'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p77_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p8_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p83'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p83_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p86'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p86_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p89'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p89_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p92'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p92_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p95_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m1p98'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass1p98_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p05'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p05_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p1_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p15'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p15_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p2_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p25'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p25_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p3_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p35'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p35_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p4_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p45'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p45_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p5_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p55'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p55_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p6_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p65_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p7_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p75'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p75_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p8_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p85'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p85_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p9'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p9_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m2p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p95_ctau10p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
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

signal_samples['V42_08Aug22_m3p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p0_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.0 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.55e-02, #1.48e-02,#1.06e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p05'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p05_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.05 GeV, 1 mm',
      mass = 3.05,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.77e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p1_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.1 GeV, 1 mm',
      mass = 3.1,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 1.96e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p15'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p15_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.15 GeV, 1 mm',
      mass = 3.15,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.06e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p2_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.2 GeV, 1 mm',
      mass = 3.2,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.13e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p25'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p25_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.25 GeV, 1 mm',
      mass = 3.25,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.21e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p3_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.3 GeV, 1 mm',
      mass = 3.3,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p35'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p35_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.35 GeV, 1 mm',
      mass = 3.35,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.22e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p4_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.4 GeV, 1 mm',
      mass = 3.4,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.24e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p45'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p45_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.45 GeV, 1 mm',
      mass = 3.45,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p5_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.5 GeV, 1 mm',
      mass = 3.5,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p55'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p55_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.55 GeV, 1 mm',
      mass = 3.55,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p6_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.6 GeV, 1 mm',
      mass = 3.6,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.23e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p65'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p65_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.65 GeV, 1 mm',
      mass = 3.65,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.25e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p7_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.7 GeV, 1 mm',
      mass = 3.7,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.28e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p75'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p75_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.75 GeV, 1 mm',
      mass = 3.75,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p8'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p8_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.8 GeV, 1 mm',
      mass = 3.8,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.27e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p85'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p85_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.85 GeV, 1 mm',
      mass = 3.85,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.31e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p9'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p9_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.9 GeV, 1 mm',
      mass = 3.9,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.34e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m3p95'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass3p95_ctau1p0/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '3.95 GeV, 1 mm',
      mass = 3.95,
      ctau = 1.0,
      resolution = None,
      filter_efficiency = 2.35e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p0'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p0_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.0 GeV, 0.1 mm',
      mass = 4.0,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.37e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p1_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.1 GeV, 0.1 mm',
      mass = 4.1,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.42e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p2_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.2 GeV, 0.1 mm',
      mass = 4.2,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.47e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

#signal_samples['V42_08Aug22_m4p3'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p3_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '4.3 GeV, 0.1 mm',
#      mass = 4.3,
#      ctau = 0.1,
#      resolution = None,
#      filter_efficiency = 2.49e-02,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

signal_samples['V42_08Aug22_m4p4'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p4_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.4 GeV, 0.1 mm',
      mass = 4.4,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.55e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p5_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.55e-02, #2.38e-02,#1.75e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p6'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p6_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.6 GeV, 0.1 mm',
      mass = 4.6,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.47e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V42_08Aug22_m4p7'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p7_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
      label = '4.7 GeV, 0.1 mm',
      mass = 4.7,
      ctau = 0.1,
      resolution = None,
      filter_efficiency = 2.28e-02,
      muon_rate = 0.47, #FIXME
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

#signal_samples['V42_08Aug22_m4p8'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p8_ctau0p1/nanoFiles/merged/flat_bparknano_08Aug22.root',
#      label = '4.8 GeV, 0.1 mm',
#      mass = 4.8,
#      ctau = 0.1,
#      resolution = None,
#      filter_efficiency = 1.46e-02,
#      muon_rate = 0.47, #FIXME
#      is_private = True,
#      colour = ROOT.kOrange+0
#      ),
#  ]

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

signal_samples['central_V11_24Apr22_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-8
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-5
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed+4
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22_sf_v3.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kRed-1
      ),
  ]

signal_samples['central_V11_24Apr22_m1_large_v2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 100000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 50000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 10000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1.0,
      ctau = 5000.0,
      resolution = 0.00888,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 500.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 100 mm',
      mass = 1.0,
      ctau = 100.0,
      resolution = 0.00873,
      filter_efficiency = 6.24e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 50.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 10 mm',
      mass = 1.0,
      ctau = 10.0,
      resolution = 0.00861,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m1p5_large_v2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 50000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 10000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 5000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.01285,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.45,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 50.0,
      resolution = 0.0129,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
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

signal_samples['central_V11_24Apr22_m2_large_v2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 10000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 5000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 500.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 50.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
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

signal_samples['central_V11_24Apr22_m3_large_v2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 10000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 5000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 1000.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1000 mm',
      mass = 3.0,
      ctau = 500.0,
      resolution = 0.0255,
      filter_efficiency = 4.97e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 100.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 50.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 10.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3.0,
      ctau = 5.0,
      resolution = 0.0248,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 1.0,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1 mm',
      mass = 3.0,
      ctau = 0.5,
      resolution = 0.0244,
      filter_efficiency = 1.49e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m4p5_large_v2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 1000.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 500.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 50.0,
      resolution = 0.0392,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 5.0,
      resolution = 0.0391,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 0.5,
      resolution = 0.0396,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.05,
      resolution = 0.0382,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['crab_V11_24Apr22_m2_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 10000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 5000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 1000.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 1000 mm',
      mass = 2.0,
      ctau = 500.0,
      resolution = 0.0167,
      filter_efficiency = 1.66e-03,
      muon_rate = 0.44,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 100.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 100 mm',
      mass = 2.0,
      ctau = 50.0,
      resolution = 0.0169,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_partial_v3.root',
      label = '2 GeV, 10 mm',
      mass = 2.0,
      ctau = 10.0,
      resolution = 0.0166,
      filter_efficiency = 6.67e-03,
      muon_rate = 0.46,
      is_private = True,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m1_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 100000 mm',
      mass = 1,
      ctau = 100000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 50000 mm',
      mass = 1,
      ctau = 50000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 10000 mm',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 5000 mm',
      mass = 1,
      ctau = 5000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 500 mm',
      mass = 1,
      ctau = 500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 100 mm',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 50 mm',
      mass = 1,
      ctau = 50.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1 GeV, 1000 mm',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 7.75e-03,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m1p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 10000 mm',
      mass = 1.5,
      ctau = 10000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 5000 mm',
      mass = 1.5,
      ctau = 5000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 1000 mm',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 500 mm',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 100 mm',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 500 mm',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '1.5 GeV, 10 mm',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 7.01e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m2_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 10000 mm',
      mass = 2.,
      ctau = 10000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 5000 mm',
      mass = 2.,
      ctau = 5000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 1000 mm',
      mass = 2.,
      ctau = 1000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 500 mm',
      mass = 2.,
      ctau = 500.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 100 mm',
      mass = 2.,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 50 mm',
      mass = 2.,
      ctau = 50.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '2 GeV, 10 mm',
      mass = 2.,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 5.86e-03,
      muon_rate = 0.46,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m3_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 10000 mm',
      mass = 3.,
      ctau = 10000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 5000 mm',
      mass = 3.,
      ctau = 5000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 1000 mm',
      mass = 3.,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 500 mm',
      mass = 3.,
      ctau = 500.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3.,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 50 mm',
      mass = 3.,
      ctau = 50.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 10 mm',
      mass = 3.,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      muon_rate = 0.47,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V11_24Apr22_m4p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 5000 mm',
      mass = 4.5,
      ctau = 5000.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 1000 mm',
      mass = 4.5,
      ctau = 1000.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 500 mm',
      mass = 4.5,
      ctau = 500.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 100 mm',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 50 mm',
      mass = 4.5,
      ctau = 50.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 10 mm',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 5 mm',
      mass = 4.5,
      ctau = 5.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 0.5 mm',
      mass = 4.5,
      ctau = 0.5,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 0.1 mm',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V11_24Apr22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_24Apr22.root',
      label = '4.5 GeV, 0.05 mm',
      mass = 4.5,
      ctau = 0.05,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      muon_rate = 0.49,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V11_24Apr22_benchmark'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1 GeV, 1000 mm',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3 GeV, 100 mm',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '4.5 GeV, 1 mm',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      colour = ROOT.kRed+4,
      ),
  ]

signal_samples['V11_24Apr22_sliding_window'] = [
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
  #    label = '1GeV, 1000mm (private, V11_24Apr22)',
  #    mass = 1,
  #    ctau = 1000.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 1.59e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 10mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 7.75e-03,
      colour = ROOT.kOrange+9
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
  #    label = '1.5GeV, 100mm (private, V11_24Apr22)',
  #    mass = 1.5,
  #    ctau = 100.0,
  #    resolution = 0.0126,
  #    filter_efficiency = 5.99e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 10mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 7.01e-03,
      colour = ROOT.kOrange+0
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
  #    label = '2GeV, 10mm (private, V11_24Apr22)',
  #    mass = 2,
  #    ctau = 10.0,
  #    resolution = 0.0158,
  #    filter_efficiency = 6.67e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  ]

signal_samples['V11_24Apr22_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '1GeV, 10mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 7.75e-03,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V11_24Apr22_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '1.5GeV, 10mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 7.01e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V11_24Apr22_m2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '2GeV, 10mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 6.67e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V11_24Apr22_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1.,
      is_private = True,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 10mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.91e-01,
      muon_rate = 1.,
      is_private = True,
      colour = ROOT.kRed-4,
      ),
  ]

signal_samples['V11_24Apr22_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1.,
      colour = ROOT.kRed-8,
      ),
  SignalSample( # temporarilly reweighted
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1.,
      colour = ROOT.kRed-8,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 1.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    colour = ROOT.kRed+4,
  #    ),
  ]

signal_samples['V11_24Apr22_m1_large'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 100000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 50000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 5000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau1000.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 1000mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 50.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1GeV, 10mm (private, V11_24Apr22)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 7.75e-03,
      colour = ROOT.kOrange+9,
      muon_rate = 1,
      is_private = True,
      ),
  ]

signal_samples['V11_24Apr22_m1p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 10000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 5000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      colour = ROOT.kOrange+0,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 100mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.0126,
      filter_efficiency = 5.99e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass1.5_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '1.5GeV, 10mm (private, V11_24Apr22)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 7.01e-03,
      muon_rate = 1,
      is_private = True,
      colour = ROOT.kOrange+0,
      ),
  ]

signal_samples['V11_24Apr22_m2_large'] = [
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 10000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 5000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 500.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      colour = ROOT.kOrange+0,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 100mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 50.0,
      resolution = 0.0164,
      filter_efficiency = 5.86e-03,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass2.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '2GeV, 10mm (private, V11_24Apr22)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 6.67e-03,
      muon_rate = 1,
      is_private = True,
      colour = ROOT.kOrange+0,
      ),
  ]

signal_samples['V11_24Apr22_m3_large'] = [
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 10000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 5000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 500.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      colour = ROOT.kRed+1,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( 
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 100mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 50.0,
      resolution = 0.0251,
      filter_efficiency = 1.40e-02,
      filter_efficiency_Bc = 1.93e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 10mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.91e-01,
      colour = ROOT.kRed-4,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 10mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 5.0,
      resolution = 0.0251,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.91e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass3.0_ctau10.0/nanoFiles/merged/flat_bparknano_24Apr22.root',
      label = '3GeV, 10mm (private, V11_24Apr22)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 1.48e-02,
      filter_efficiency_Bc = 1.91e-01,
      muon_rate = 1,
      is_private = True,
      ),
  ]

signal_samples['V11_24Apr22_m4p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 5000.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 1000.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 500.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      colour = ROOT.kRed-8,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 50.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 100mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 5.0,
      resolution = 0.0368,
      filter_efficiency = 2.30e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      colour = ROOT.kRed+4,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 0.5,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 0.05,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 0.01,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau100.0/nanoFiles/merged/flat_bparknano_24Apr22_v2.root',
      filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau100.0/nanoFiles/Chunk0_n219/flat/flat_bparknano_24Apr22.root',
      label = '4.5GeV, 1mm (private, V11_24Apr22)',
      mass = 4.5,
      ctau = 0.005,
      resolution = 0.0368,
      filter_efficiency = 2.38e-02,
      filter_efficiency_Bc = 1.40e-01,
      muon_rate = 1,
      is_private = True,
      ),
  #SignalSample( # generated
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 1.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 0.5,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 0.1,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 0.05,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass4.5_ctau1.0/nanoFiles/merged/flat_bparknano_24Apr22_merged.root',
  #    label = '4.5GeV, 1mm (private, V11_24Apr22)',
  #    mass = 4.5,
  #    ctau = 0.01,
  #    resolution = 0.0368,
  #    filter_efficiency = 2.38e-02,
  #    ),
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

#signal_samples['V10_30Dec21_m1_large'] = [
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 100000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 70000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 50000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 40000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 30000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 20000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 15000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 10000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 7000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 5000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 4000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 3000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 2000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 1500.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '1GeV, 1000mm (central, V10_30Dec21)',
#      mass = 1,
#      ctau = 1000.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 700.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 500.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 400.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 300.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 200.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 150.0,
#      resolution = 0.00853,
#      filter_efficiency = 1.59e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '1GeV, 100mm (central, V10_30Dec21)',
#      mass = 1,
#      ctau = 100.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+7
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 70.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 50.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 40.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 30.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 20.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 1,
#      ctau = 15.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.17e-03,
#      colour = ROOT.kOrange+0
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '1GeV, 10mm (central, V10_30Dec21)',
#      mass = 1,
#      ctau = 10.0,
#      resolution = 0.00853,
#      filter_efficiency = 5.85e-03,
#      colour = ROOT.kOrange+9
#      ),
#  ]

signal_samples['V10_30Dec21_m1_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 100000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 70000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 50000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 40000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 30000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 20000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 15000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 7000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 5000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 4000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 3000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 2000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 1500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 1000mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 700.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 400.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 300.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 200.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 150.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 100mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 70.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 50.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 40.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 30.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 20.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1,
      ctau = 15.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 10mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V10_30Dec21_m1p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 10000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 7000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 5000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 4000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 3000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 2000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 1500.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kYellow-3
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 700.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 500.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 400.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 300.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 200.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 150.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 70.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 50.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 40.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 30.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 20.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 1.5,
      ctau = 15.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      muon_rate = 0.5,
      is_private= False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
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

signal_samples['V10_30Dec21_m2_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 10000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 7000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 5000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 4000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 3000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 2000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 1500.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 1000mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04, # is this one correct?
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 700.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 500.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 400.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 300.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 200.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 150.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 100mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 70.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 50.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 40.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 30.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 20.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 2,
      ctau = 15.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 10mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m3_large'] = [
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1000mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-9,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 700.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 500.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 400.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    label = '3GeV, 1000mm (central, V10_30Dec21)',
  #    mass = 3,
  #    ctau = 300.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    label = '3GeV, 1000mm (central, V10_30Dec21)',
  #    mass = 3,
  #    ctau = 200.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    label = '3GeV, 1000mm (central, V10_30Dec21)',
  #    mass = 3,
  #    ctau = 150.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 2.40e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+1,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 70.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 50.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 40.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 30.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 20.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 15.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 10mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-4,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 7.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 5.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 4.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 3.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
  #    mass = 3,
  #    ctau = 2.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    muon_rate = 0.5,
  #    is_private = False,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 3,
      ctau = 1.5,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+1,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+2,
      ),
  ]

signal_samples['V10_30Dec21_m3_large_generated'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1000mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-9,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 100mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 10mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+2,
      ),
  ]

#signal_samples['V10_30Dec21_m3_large'] = [
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 1000mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 1000.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed-9,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 700.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 500.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 400.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 1000mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 300.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 1000mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 200.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 1000mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 150.0,
#      resolution = 0.0251,
#      filter_efficiency = 2.40e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
#      label = '3GeV, 100mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 100.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      filter_efficiency_Bc = 1.54e-01, # to modify
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 70.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 50.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 40.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 30.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 20.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 15.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.58e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 10mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 10.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed-4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 7.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 5.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 4.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 3.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 2.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 3,
#      ctau = 1.5,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+1,
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '3GeV, 1mm (central, V10_30Dec21)',
#      mass = 3,
#      ctau = 1.0,
#      resolution = 0.0251,
#      filter_efficiency = 5.87e-03,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+2,
#      ),
#  ]


#signal_samples['V10_30Dec21_m4p5_large'] = [
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '4.5GeV, 100mm (central, V10_30Dec21)',
#      mass = 4.5,
#      ctau = 100.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed-8,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 70.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 50.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 40.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 30.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 20.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 15.0,
#      resolution = 0.0368,
#      filter_efficiency = 6.10e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '4.5GeV, 10mm (central, V10_30Dec21)',
#      mass = 4.5,
#      ctau = 10.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed-5,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 7.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 5.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 4.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 3.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 2.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 1.5,
#      resolution = 0.0368,
#      filter_efficiency = 8.00e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample( # generated
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      label = '4.5GeV, 1mm (central, V10_30Dec21)',
#      mass = 4.5,
#      ctau = 1.0,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.7,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.5,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.4,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.3,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.2,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.15,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.1,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.07,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.05,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.04,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  SignalSample(
#      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau*mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
#      mass = 4.5,
#      ctau = 0.03,
#      resolution = 0.0368,
#      filter_efficiency = 8.35e-04,
#      muon_rate = 0.5,
#      is_private = False,
#      colour = ROOT.kRed+4,
#      ),
#  ]

signal_samples['V10_30Dec21_m4p5_large'] = [
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 100mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 70.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 50.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 40.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 30.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 20.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 15.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 10mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed-5,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 7.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 5.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 4.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 3.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 2.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 1.5,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 1mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.7,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.5,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.4,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.3,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.2,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.15,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.07,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.05,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.04,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      mass = 4.5,
      ctau = 0.03,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      muon_rate = 0.5,
      is_private = False,
      colour = ROOT.kRed+4,
      ),
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH1_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH2_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH3_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH4_Run2018D (V13_06Feb23)',
      lumi = 5.302 * 0.9854,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018D/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH5_Run2018D (V13_06Feb23)',
      lumi = 5.302,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018C/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH1_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018C/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH2_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018C/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH3_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018C/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH4_Run2018C (V13_06Feb23)',
      lumi = 1.103,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018C/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH5_Run2018C (V13_06Feb23)',
      lumi = 1.103 * 0.9997,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH1_Run2018B (V13_06Feb23)',
      lumi = 0.911 * 0.9997,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH2_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH3_Run2018B (V13_06Feb23)',
      lumi = 0.911 * 0.0122,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH4_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH5_Run2018B (V13_06Feb23)',
      lumi = 0.911,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH6_Run2018B/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH6_Run2018B (V13_06Feb23)',
      lumi = 0.377,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH1_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH1_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH2_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH2_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH3_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH3_Run2018A (V13_06Feb23)',
      lumi = 0.774 * 0.3259,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH4_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH4_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH5_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH5_Run2018A (V13_06Feb23)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V13_06Feb23/ParkingBPH6_Run2018A/merged/flat_bparknano_06Feb23_partial.root',
      label = 'ParkingBPH6_Run2018A (V13_06Feb23)',
      lumi = 0.774 * 0.9843,
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


