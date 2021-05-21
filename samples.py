import ROOT

class DataSample(object):
  def __init__(self, filename, label):
    self.filename = filename
    self.label = label


class SignalSample(object): 
  def __init__(self, filename, label, mass, ctau, filter_efficiency='', colour=ROOT.kRed):
    self.filename = filename
    self.label = label
    self.mass = mass
    self.ctau = ctau
    self.filter_efficiency = filter_efficiency
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
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V21/mass4.5_ctau1.2/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #label = '4.5GeV, 1.2mm (V21)',
      label = '4.5GeV, 1.2mm',
      mass = 4.5,
      ctau = 1.2,
      colour = ROOT.kRed+4,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #label = '3GeV, 184mm (V20_emu)',
      label = '3GeV, 184mm',
      mass = 3,
      ctau = 184,
      filter_efficiency = 4e-03,
      colour = ROOT.kRed+1,
      ),
  SignalSample(
      #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_selected_stdtrgmu_full_v1.root',
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V26/mass1.0_ctau10000.0/nanoFiles/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      #label = '1GeV, 10^{4}mm (V26)',
      label = '1GeV, 10^{4}mm',
      mass = 1,
      ctau = 1e4,
      colour = ROOT.kOrange+0
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

data_samples = [
  DataSample(
      filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V04/ParkingBPH1_Run2018A/merged/flat_bparknano_newselection_stdtrgmu_v1.root',
      label = 'ParkingBPH1_Run2018A (V04)',
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

qcd_samples = [
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



