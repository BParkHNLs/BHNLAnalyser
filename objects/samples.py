import ROOT

class DataSample(object):
  def __init__(self, filename, label, lumi=0.074):
    self.filename = filename
    self.label = label
    self.lumi = lumi


class SignalSample(object): 
  def __init__(self, filename='', filename_Bc='', label='', mass='', ctau='', resolution='', filter_efficiency='', filter_efficiency_Bc='', colour=ROOT.kRed):
    self.filename = filename
    self.filename_Bc = filename_Bc
    self.label = label
    self.mass = mass
    self.ctau = ctau # this is the target ctau, not necessarily the ctau with which the sample was produced
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 1000mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 100mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1GeV, 10mm (central, V10_30Dec21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['V10_30Dec21_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1.5GeV, 1000mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1.5GeV, 100mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '1.5GeV, 10mm (central, V10_30Dec21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 1000mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04, # is this one correct?
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 100mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '2GeV, 10mm (central, V10_30Dec21)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['V10_30Dec21_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1000mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed-9,
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 10mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+2,
      ),
  SignalSample( # reweighted
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '3GeV, 1mm (central, V10_30Dec21)',
      mass = 3,
      ctau = 0.1,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+3,
      ),
  ]

signal_samples['V10_30Dec21_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 100mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 10mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed-5,
      ),
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_30Dec21.root',
      label = '4.5GeV, 1mm (central, V10_30Dec21)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
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
      label = '4.5GeV, 1mm (central, loose)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0373, # to be checked
      filter_efficiency = 5.53e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V00_looseselection/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_3files.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, loose)',
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
      label = '1GeV, 10^{3}mm (central, loose)',
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

signal_samples['central_V09_06Nov21_benchmark'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 10^{3}mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_significance'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 10mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m1'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 10mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      colour = ROOT.kOrange+9
      ),
  ]

signal_samples['central_V09_06Nov21_m1p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 1000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 100mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m2'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04, # is this one correct?
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m3'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed-9,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+2,
      ),
  SignalSample( # reweighted
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 0.1,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+3,
      ),
  ]

signal_samples['central_V09_06Nov21_m4p5'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed-8,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed-5,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed-2,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 0.1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 0.1,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.62e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  ]

signal_samples['central_V09_06Nov21_m1_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 100000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 70000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 50000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 40000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 30000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 20000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 15000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 7000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 7000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 5000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 5000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 4000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 4000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 3000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 2000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 1500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 1000.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 700.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 500.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 400.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 300.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 200.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 1000mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 150.0,
      resolution = 0.00853,
      filter_efficiency = 1.59e-03,
      colour = ROOT.kOrange+0
      ),
  #
  #
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '1GeV, 100mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 400.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.17e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 100mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 300.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.17e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 200mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 200.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.17e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 200mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 150.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.17e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 100.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 70.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 50.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 40.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 30.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 20.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 100mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 15.0,
      resolution = 0.00853,
      filter_efficiency = 5.17e-03,
      colour = ROOT.kOrange+0
      ),
  #
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '1GeV, 10mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 40.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.85e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 10mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 30.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.85e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 10mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 20.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.85e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '1GeV, 10mm (central, V09_06Nov21)',
  #    mass = 1,
  #    ctau = 15.0,
  #    resolution = 0.00853,
  #    filter_efficiency = 5.85e-03,
  #    colour = ROOT.kOrange+0
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1GeV, 10mm (central, V09_06Nov21)',
      mass = 1,
      ctau = 10.0,
      resolution = 0.00853,
      filter_efficiency = 5.85e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m1p5_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 1000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 1000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 7000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 5000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 4000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 3000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 2000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 1500mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 1000mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 1000.0,
      resolution = 0.0126,
      filter_efficiency = 1.03e-03,
      colour = ROOT.kYellow-3
      ),
  #
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 700mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 500mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 400mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 300mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 200mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 100mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 150.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 100mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 100.0,
      resolution = 0.0126,
      filter_efficiency = 3.35e-03,
      colour = ROOT.kOrange+0
      ),
  #
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 70mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 70.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 50.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 40.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 30.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 20mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 15mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '1.5GeV, 10mm (central, V09_06Nov21)',
      mass = 1.5,
      ctau = 10.0,
      resolution = 0.0122,
      filter_efficiency = 3.72e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m2_large'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 10000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 7000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 5000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 4000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 3000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 2000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 1500.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04,
      colour = ROOT.kAzure+7
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 1000mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 1000.0,
      resolution = 0.0170,
      filter_efficiency = 8.95e-04, # is this one correct?
      colour = ROOT.kAzure+7
      ),
  #
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 700.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 500.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 400.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 300.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 200.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 150.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 100mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 100.0,
      resolution = 0.0164,
      filter_efficiency = 1.57e-03,
      colour = ROOT.kOrange+0
      ),
    #
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 70.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 50.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 40.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 30.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 20.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 15.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '2GeV, 10mm (central, V09_06Nov21)',
      mass = 2,
      ctau = 10.0,
      resolution = 0.0158,
      filter_efficiency = 1.78e-03,
      colour = ROOT.kOrange+0
      ),
  ]

signal_samples['central_V09_06Nov21_m3_large'] = [
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 1000.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 700.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 500.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 400.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 300.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 200.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1000mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 150.0,
      resolution = 0.0251,
      filter_efficiency = 2.40e-03,
      colour = ROOT.kRed+1,
      ),
  #
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '3GeV, 100mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 400.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 100mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 300.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 100mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 200.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 100mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 150.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.58e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 100.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      filter_efficiency_Bc = 1.54e-01, # to modify
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 70.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 50.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 40.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 30.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 20.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 100mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 15.0,
      resolution = 0.0251,
      filter_efficiency = 5.58e-03,
      colour = ROOT.kRed+1,
      ),
  #
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '3GeV, 10mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 40.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 10mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 30.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 10mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 20.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample( 
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 10mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 15.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 10.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 7.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 5.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 4.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 3.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 2.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 10mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 1.5,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  #
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '3GeV, 1mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 4.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 1mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 3.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 1mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 2.0,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 1mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 1.5,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '3GeV, 1mm (central, V09_06Nov21)',
      mass = 3,
      ctau = 1.0,
      resolution = 0.0251,
      filter_efficiency = 5.87e-03,
      colour = ROOT.kRed+1,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '3GeV, 1mm (central, V09_06Nov21)',
  #    mass = 3,
  #    ctau = 0.1,
  #    resolution = 0.0251,
  #    filter_efficiency = 5.87e-03,
  #    colour = ROOT.kRed+1,
  #    ),
  ]

signal_samples['central_V09_06Nov21_m4p5_large'] = [
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 100.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 70.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 50.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 40.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 30.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 20.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 100mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 15.0,
      resolution = 0.0368,
      filter_efficiency = 6.10e-04,
      colour = ROOT.kRed+4,
      ),
  #
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 10mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 40.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.00e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 10mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 30.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.00e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 10mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 20.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.00e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 10mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 15.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.00e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 10.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 7.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 5.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 4.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 3.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 2.0,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 10mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 1.5,
      resolution = 0.0368,
      filter_efficiency = 8.00e-04,
      colour = ROOT.kRed+4,
      ),
  #
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 4.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 3.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 2.0,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 1.5,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  SignalSample( # generated
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 1.0,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.7,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.5,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.4,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.3,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.2,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.15,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21_ctaumerged.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.1,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.07,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.05,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.04,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
      label = '4.5GeV, 1mm (central, V09_06Nov21)',
      mass = 4.5,
      ctau = 0.03,
      resolution = 0.0368,
      filter_efficiency = 8.35e-04,
      colour = ROOT.kRed+4,
      ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 0.02,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 0.015,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  #SignalSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V09_06Nov21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_06Nov21.root',
  #    label = '4.5GeV, 1mm (central, V09_06Nov21)',
  #    mass = 4.5,
  #    ctau = 0.01,
  #    resolution = 0.0368,
  #    filter_efficiency = 8.35e-04,
  #    colour = ROOT.kRed+4,
  #    ),
  ]

signal_samples['loose_dsaonly'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      label = '4.5GeV, 1.2mm (V21)',
      mass = 4.5,
      ctau = 1.2,
      resolution = 0.0373,
      filter_efficiency = 3.91e-04,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      label = '3GeV, 184mm (V20_emu)',
      mass = 3,
      ctau = 184.0,
      resolution = 0.0243,
      filter_efficiency = 3.93e-03,
      filter_efficiency_Bc = 1.45e-01,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/bparknano_looseselection_matched_dsaonly.root',
      label = '1GeV, 10^{4}mm (V26)',
      mass = 1,
      ctau = 10000.0,
      resolution = 0.0104,
      filter_efficiency = 2.24e-04,
      colour = ROOT.kOrange+0
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
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V0/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_sf_study15Sep21_A1_B1_v1.root',
      label = 'BToJPsiKstar (V0), JPsi->#mu#mu',
      colour = ROOT.kOrange+0,
      ),
]

signal_samples['limits_m4p5_29Sep21'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.1*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04, #3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.01*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.001*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.0001*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  ]

signal_samples['limits_m4p5'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_18Aug21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.1*1e3,
      resolution = 0.0373,
      filter_efficiency = 3.53e-04,
      filter_efficiency_Bc = 3.0e-02,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_18Aug21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
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
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_18Aug21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      mass = 4.5,
      ctau = 0.001*1e3,
      resolution = 0.0373,
      filter_efficiency = 4.0e-04,
      filter_efficiency_Bc = 3.43e-02,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_18Aug21.root',
      #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21_Bc/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_29Jun21.root',
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


signal_samples['limits_m4p5_large'] = [
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

signal_samples['limits_m3_29Sep21'] = [
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e03,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03, #1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    label = '3GeV, 184mm',
    mass = 3.0,
    ctau = 184.0,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03,
    filter_efficiency_Bc = 1.45e-01,
    colour = ROOT.kRed+1,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e02,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03, #5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e01,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03, #5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e00,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03, #5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
]

signal_samples['limits_m3'] = [
SignalSample(
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e03,
    resolution = 0.0243,
    filter_efficiency = 1.93e-03,
    filter_efficiency_Bc = 4.8e-02,
    ),
SignalSample(
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    label = '3GeV, 184mm',
    mass = 3.0,
    ctau = 184.0,
    resolution = 0.0243,
    filter_efficiency = 3.93e-03,
    filter_efficiency_Bc = 1.45e-01,
    colour = ROOT.kRed+1,
    ),
SignalSample(
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e02,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    mass = 3.0,
    ctau = 1e01,
    resolution = 0.0243,
    filter_efficiency = 5e-03,
    filter_efficiency_Bc = 2.5e-01,
    ),
SignalSample(
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
    #filename_Bc = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_Bc/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
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


signal_samples['limits_m3_large'] = [
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

signal_samples['limits_m1_29Sep21'] = [
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1e05,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4, #2.34e-05,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1e03,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4, #1.80e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1e02,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4, #6.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1e01,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4, #7.0e-03,
      ),
  SignalSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Sep21_v2.root',
      mass = 1.0,
      ctau = 1,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4, #7.0e-03,
      ),
  ]


signal_samples['limits_m1'] = [
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
      mass = 1.0,
      ctau = 1e05,
      resolution = 0.0104,
      filter_efficiency = 2.34e-05,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
      mass = 1.0,
      ctau = 1e04,
      resolution = 0.0104,
      filter_efficiency = 2.24e-4,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
      mass = 1.0,
      ctau = 1e03,
      resolution = 0.0104,
      filter_efficiency = 1.80e-03,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
      mass = 1.0,
      ctau = 1e02,
      resolution = 0.0104,
      filter_efficiency = 6.0e-03,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_29Jun21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_18Aug21.root',
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


signal_samples['limits_m1_large'] = [
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


signal_samples['loose'] = [
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

data_samples = {}
data_samples['F1'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano_largestPt_SS.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/Chunk0_n750/flat/flat_bparknano.root',
      label = 'ParkingBPH1_Run2018A (F1)',
      ),
  ]

data_samples['V02'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_selCos2d_nonNullChargeCR.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_selCos2d.root',
      label = 'ParkingBPH1_Run2018A (V02)',
      ),
  ]

data_samples['V03'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V03/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      label = 'ParkingBPH1_Run2018A (V03)',
      ),
  ]

data_samples['V04'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V04/ParkingBPH1_Run2018A/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V04/ParkingBPH1_Run2018A/merged/flat_bparknano_selected_updatedmatching.root',
      label = 'ParkingBPH1_Run2018A (V04)',
      ),
  ]

data_samples['V05_29Jun21'] = [
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

data_samples['V07_18Aug21'] = [
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH1_Run2018A/merged/flat_bparknano_18Aug21_nodsa.root',
  #    label = 'ParkingBPH1_Run2018A (V07_18Aug21)',
  #    ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH1_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
      label = 'ParkingBPH1_Run2018A (V07_18Aug21)',
      ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH2_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH2_Run2018A (V07_18Aug21)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH3_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH3_Run2018A (V07_18Aug21)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH4_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH4_Run2018A (V07_18Aug21)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH5_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH5_Run2018A (V07_18Aug21)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH6_Run2018A/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH6_Run2018A (V07_18Aug21)',
  #    ),
  #DataSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH1_Run2018B/merged/flat_bparknano_18Aug21_v2.root',
  #    label = 'ParkingBPH1_Run2018B (V07_18Aug21)',
  #    ),
  ]

data_samples['V08_29Sep21'] = [
  DataSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V08/ParkingBPH1_Run2018A/merged/flat_bparknano_29Sep21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V08/ParkingBPH1_Run2018A/Chunk0_n500/flat/flat_bparknano_29Sep21_v3.root', # 10.8% of A1
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V08/ParkingBPH1_Run2018A/merged/flat_bparknano_29Sep21_v3.root',
      label = 'ParkingBPH1_Run2018A (V08_29Sep21)',
      ),
  ]

data_samples['V09_06Nov21'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH1_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH2_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH2_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH3_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH3_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH4_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH4_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH5_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH5_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH6_Run2018A/merged/flat_bparknano_06Nov21.root',
      label = 'ParkingBPH6_Run2018A (V09_06Nov21)',
      lumi = 0.774,
      ),
  ]

data_samples['V10_30Dec21'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21.root',
      label = 'ParkingBPH1_Run2018D (V10_30Dec21)',
      lumi = 5.302,
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/bparknano_data_1file_looseselection.root',
      label = 'background',
      ),
  ]

data_samples['tag_and_probe'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V06_tag_and_probe/ParkingBPH1_Run2018A/merged/flat_bparknano_tag_and_probe_v2_tag_fired_DST_DoubleMu1.root',
      label = 'ParkingBPH1_Run2018A (V06)',
      ),
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V06_tag_and_probe/ParkingBPH1_Run2018B/merged/flat_bparknano_tag_and_probe_v2.root',
      label = 'ParkingBPH1_Run2018B (V06)',
      ),
  ]

data_samples['loose_dsaonly'] = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V07_18Aug21/ParkingBPH1_Run2018A/merged/bparknano_looseselection_dsaonly.root',
      label = 'ParkingBPH1_Run2018A (1file)',
      ),
  ]


qcd_samples = {}
qcd_samples['V10_30Dec21'] = [
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
  #    label = 'QCD_pt15to20 (V10_30Dec21)',
  #    cross_section = 1.27319e+09,
  #    filter_efficiency = 0.00300,
  #    colour = ROOT.kBlue+2, 
  #    ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
      label = 'QCD_pt20to30 (V10_30Dec21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
      label = 'QCD_pt30to50 (V10_30Dec21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
      label = 'QCD_pt50to80 (V10_30Dec21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_extmerged.root',
      label = 'QCD_pt80to120 (V10_30Dec21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_30Dec21_extmerged.root',
      label = 'QCD_pt120to170 (V10_30Dec21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_30Dec21.root',
      label = 'QCD_pt170to300 (V10_30Dec21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V09_06Nov21'] = [
  #QCDMCSample(
  #    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Nov21_v3.root',
  #    label = 'QCD_pt15to20 (V09_06Nov21)',
  #    cross_section = 1.27319e+09,
  #    filter_efficiency = 0.00300,
  #    colour = ROOT.kBlue+2, 
  #    ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Nov21_v4.root',
      label = 'QCD_pt20to30 (V09_06Nov21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Nov21_v4.root',
      label = 'QCD_pt30to50 (V09_06Nov21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Nov21_v4.root',
      label = 'QCD_pt50to80 (V09_06Nov21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_06Nov21_v4_extmerged.root',
      label = 'QCD_pt80to120 (V09_06Nov21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_06Nov21_v4_extmerged.root',
      label = 'QCD_pt120to170 (V09_06Nov21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V09_06Nov21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_06Nov21_v4.root',
      label = 'QCD_pt170to300 (V09_06Nov21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V08_29Sep21'] = [
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Sep21_v5.root',
      label = 'QCD_pt15to20 (V08_29Sep21)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Sep21_v5.root',
      label = 'QCD_pt20to30 (V08_29Sep21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Sep21_v5.root',
      label = 'QCD_pt30to50 (V08_29Sep21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Sep21_v5.root',
      label = 'QCD_pt50to80 (V08_29Sep21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Sep21_v5_extmerged.root',
      label = 'QCD_pt80to120 (V08_29Sep21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Sep21_v5_extmerged.root',
      label = 'QCD_pt120to170 (V08_29Sep21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V08/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Sep21_v5.root',
      label = 'QCD_pt170to300 (V08_29Sep21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V07_18Aug21'] = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1.root',
      label = 'QCD_pt15to20 (V07_18Aug21)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1.root',
      label = 'QCD_pt20to30 (V07_18Aug21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1.root',
      label = 'QCD_pt30to50 (V07_18Aug21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1.root',
      label = 'QCD_pt50to80 (V07_18Aug21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_extmerged.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1_extmerged.root',
      label = 'QCD_pt80to120 (V07_18Aug21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_extmerged.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1_extmerged.root',
      label = 'QCD_pt120to170 (V07_18Aug21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_v4.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_18Aug21_sf_study15Sep21_A1_6_B1_v1.root',
      label = 'QCD_pt170to300 (V07_18Aug21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V06_29Jun21'] = [
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study31Jul21.root',
      label = 'QCD_pt15to20 (V06_29Jun21)',
      cross_section = 1.27319e+09,
      filter_efficiency = 0.00300,
      colour = ROOT.kBlue+2, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study31Jul21.root',
      label = 'QCD_pt20to30 (V06_29Jun21)',
      cross_section = 5.58528e+08,
      filter_efficiency = 0.00530,
      colour = ROOT.kBlue-4, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study31Jul21.root',
      label = 'QCD_pt30to50 (V06_29Jun21)',
      cross_section = 1.39803e+08, 
      filter_efficiency = 0.01182,
      colour = ROOT.kBlue-9, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study31Jul21.root',
      label = 'QCD_pt50to80 (V06_29Jun21)',
      cross_section = 1.92225e+07,
      filter_efficiency = 0.02276,
      colour = ROOT.kBlue-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1_extmerged.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_sf_study31Jul21_extmerged.root',
      label = 'QCD_pt80to120 (V06_29Jun21)',
      cross_section = 2.75842e+06,
      filter_efficiency = 0.03844,
      colour = ROOT.kRed-10, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2_extmerged.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1_extmerged.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_ext/merged/flat_bparknano_29Jun21_sf_study31Jul21_extmerged.root',
      label = 'QCD_pt120to170 (V06_29Jun21)',
      cross_section = 469797,
      filter_efficiency = 0.05362,
      colour = ROOT.kOrange+6, 
      ),
  QCDMCSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_trigger_scale_factor_A1_6_v2.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study15Sep21_A1_6_v1.root',
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/merged/flat_bparknano_29Jun21_sf_study31Jul21.root',
      label = 'QCD_pt170to300 (V06_29Jun21)',
      cross_section = 117989,
      filter_efficiency = 0.07335,
      colour = ROOT.kOrange+8, 
      ),
  ]

qcd_samples['V05'] = [
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


qcd_samples['triggermuon_matching_check'] = [
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


qcd_samples['V04'] = [
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

qcd_samples['V03'] = [
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


qcd_samples['V02'] = [
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



