'''
  This class defines extra selection to be applied at analysis level
  Both the 'flat' and 'nano' syntaxes can be reported
'''

class Selection(object):
  def __init__(self, flat='', nano=''):
    self.flat = flat
    self.nano = nano


selection = {}
selection['none'] = Selection(
    flat = 'hnl_pt > 0',
    #flat = 'mu0_triggermatching_dr<0.05',
    nano = 'BToMuMuPi_hnl_pt > 0',
)

selection['control_Bc'] = Selection(
    #flat = 'b_pt > 0',
    #flat = 'abs(dimu_mass-3.097)<0.07 && b_cos2d>0.997 && sv_prob>0.05 && abs(l1_dxy)<0.01 && abs(l2_dxy)<0.015 && mindr>0.3 && sv_lxy>0.02 && sv_lxy<0.075 && b_mass>5.5 && b_mass<7.5',
    #flat = 'abs(dimu_mass-3.097)<0.07 && b_cos2d>0.998 && sv_prob>0.1 && abs(l1_dxy)<0.015 && abs(l2_dxy)<0.015 && mindr>0.3 && maxdr<1.5 && sv_lxy<0.075 && b_mass>5.5 && b_mass<7.5'
    flat = 'abs(dimu_mass-3.097)<0.07 && l2_pt>3 && k_pt>1.5 && b_cos2d>0.999 && l1_mediumid==1 && l2_mediumid==1 && sv_prob>0.1 && abs(l1_dxy)<0.05 && abs(l2_dxy)<0.05 && mindr>0.1 && maxdr<1.5 && sv_lxysig>3 && b_mass>5.5 && b_mass<7.5'
)

selection['nodsa'] = Selection(
    flat = 'mu_isdsa!=1',
)

selection['standard'] = Selection(
    flat = ' && '.join([
      'mu_isdsa !=1',
      'trgmu_softid == 1',
      'mu_looseid == 1',
      'mu_intimemuon == 1',
      'mu_trackerhighpurityflag == 1'
    ])
)

selection['preselection'] = Selection(
    nano = ' && '.join([
      'Muon_pt[BToMuMuPi_trg_mu_idx]>1.5',
      'abs(Muon_eta[BToMuMuPi_trg_mu_idx])<2',
      'BToMuMuPi_pi_pt>0.7',
      'abs(BToMuMuPi_pi_eta)<2',
      'abs(BToMuMuPi_pi_dz)>0.005',
      'abs(BToMuMuPi_pi_dxy)>0.005',
      'abs(BToMuMuPi_pi_dzS)>1.5',
      'abs(BToMuMuPi_pi_dxyS)>3',
      'abs(BToMuMuPi_pi_DCASig)>5',
      'Muon_pt[BToMuMuPi_sel_mu_idx]>1.5',
      'abs(Muon_eta[BToMuMuPi_sel_mu_idx])<2',
      'abs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.0015',
      'abs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.001',
      'abs(Muon_dzS[BToMuMuPi_sel_mu_idx])>1.',
      'abs(Muon_dxyS[BToMuMuPi_sel_mu_idx])>1.5',
      'BToMuMuPi_sv_prob>0.001',
      'BToMuMuPi_hnl_cos2D>0.995',
      'BToMuMuPi_sv_lxy_sig>15',
      'BToMuMuPi_mass<8',
      'Muon_softId[BToMuMuPi_trg_mu_idx] == 1',
      'Muon_looseId[BToMuMuPi_sel_mu_idx] == 1',
      #'ProbeTracks_highPurityFlag[BToMuMuPi_pi_idx] == 1',
    ]),
    flat = ' && '.join([
      'mu0_pt>1.5',
      'abs(mu0_eta)<2',
      'pi_pt>0.7',
      'abs(pi_eta)<2',
      'pi_dz>0.005',
      'pi_dxy>0.005',
      'pi_dzsig>1.5',
      'pi_dxysig>3',
      'pi_dcasig>5',
      'mu_pt>1.5',
      'abs(mu_eta)<2',
      'mu_dz>0.0015',
      'mu_dxy>0.001',
      'mu_dzsig>1.',
      'mu_dxysig>1.5',
      'sv_prob>0.001',
      'hnl_cos2d>0.995',
      'sv_lxysig>15',
      'b_mass<8',
      'mu0_softid == 1',
      'mu_looseid == 1',
      'pi_packedcandhashighpurity == 1',
      '((mu0_charge!=mu_charge && abs(mu0_mu_mass-3.097)>0.15 && abs(mu0_mu_mass-3.686)>0.08 && abs(mu0_mu_mass-1.019)>0.01) || (mu0_charge==mu_charge))',
      '((mu0_charge==mu_charge && abs(mu0_pi_mass-3.097)>0.05 && abs(mu0_pi_mass-1.76)>0.05) || (mu0_charge!=mu_charge))',
      'sv_lxy < 100',
      'sv_lxy/sv_lxysig>0',
      #'hnl_mass>4.4247326 && hnl_mass<4.5752674',
    ])
)

selection['study_Nov21'] = Selection(
    flat = ' && '.join([
      'mu_isdsa !=1',
      'trgmu_softid == 1',
      'mu_looseid == 1',
      #'mu_intimemuon==1',
      #'mu_trackerhighpurityflag==1',
      #'((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))',
      #'trgmu_pi_mass < 4.5',
      #'trgmu_mu_mass < 4.5',
      '(trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)',
      #'deltar_trgmu_pi < 1.5',
      #'deltar_trgmu_mu < 1',
    ])
)

selection['baseline_30Dec21'] = Selection(
    flat = ' && '.join([
      'mu_isdsa !=1',
      'trgmu_softid == 1',
      'mu_looseid == 1',
      #'mu_customisedid == 1',
      'pi_packedcandhashighpurity == 1',
      #'trgmu_pi_mass < 4.5',
      #'trgmu_mu_mass < 4.5',
      '((trgmu_charge!=mu_charge && abs(trgmu_mu_mass-3.097)>0.15 && abs(trgmu_mu_mass-3.686)>0.08 && abs(trgmu_mu_mass-1.019)>0.01) || (trgmu_charge==mu_charge))',
      '((trgmu_charge==mu_charge && abs(trgmu_pi_mass-3.097)>0.05 && abs(trgmu_pi_mass-1.76)>0.05) || (trgmu_charge!=mu_charge))',
      #'deltar_trgmu_pi < 1.5',
      #'deltar_trgmu_mu < 1',
      'sv_lxy < 100',
    ])
)

selection['baseline_08Aug22'] = Selection(
    flat = ' && '.join([
      'hnl_charge==0',
      'mu0_softid==1',
      'mu_looseid==1',
      'pi_packedcandhashighpurity==1',
      '((mu0_charge!=mu_charge && abs(mu0_mu_mass-3.097)>0.15 && abs(mu0_mu_mass-3.686)>0.08 && abs(mu0_mu_mass-1.019)>0.01) || (mu0_charge==mu_charge))',
      '((mu0_charge==mu_charge && abs(mu0_pi_mass-3.097)>0.05 && abs(mu0_pi_mass-1.76)>0.05) || (mu0_charge!=mu_charge))',
      'sv_lxy<100',
      'sv_lxy/sv_lxysig>0',
      '((mu_istriggering==1 && mu_pt>7 && abs(mu_eta)<1.5) || (mu_istriggering==0))',
      '((mu0_istriggering==1 && mu0_pt>7 && abs(mu0_eta)<1.5) || (mu0_istriggering==0))',
    ])
)

selection['baseline_08Aug22_CR'] = Selection(
    flat = ' && '.join([
      'mu0_softid==1',
      'mu_looseid==1',
      'pi_packedcandhashighpurity==1',
      '((mu0_charge!=mu_charge && abs(mu0_mu_mass-3.097)>0.15 && abs(mu0_mu_mass-3.686)>0.08 && abs(mu0_mu_mass-1.019)>0.01) || (mu0_charge==mu_charge))',
      '((mu0_charge==mu_charge && abs(mu0_pi_mass-3.097)>0.05 && abs(mu0_pi_mass-1.76)>0.05) || (mu0_charge!=mu_charge))',
      'sv_lxy<100',
      'sv_lxy/sv_lxysig>0',
      '((mu_istriggering==1 && mu_pt>7 && abs(mu_eta)<1.5) || (mu_istriggering==0))',
      '((mu0_istriggering==1 && mu0_pt>7 && abs(mu0_eta)<1.5) || (mu0_istriggering==0))',
      #'((mu0_istriggering==1 && mu0_pt>7. && abs(mu0_eta)<1.5 && abs(mu0_dxysig_bs)>3.) || (mu_istriggering==1 && mu_pt>7. && abs(mu_eta)<1.5 && abs(mu_dxysig_bs)>3.))',
      #'((mu0_istriggering==1 && mu0_pt>8. && abs(mu0_eta)<1.5 && abs(mu0_dxysig_bs)>3.) || (mu_istriggering==1 && mu_pt>8. && abs(mu_eta)<1.5 && abs(mu_dxysig_bs)>3.))',
      #'mu_istriggering==0 && mu0_istriggering==1 && mu0_pt>8. && abs(mu0_eta)<1.5 && abs(mu0_dxysig_bs)>3.',
    ])
)

selection['tag_and_probe'] = Selection(
    flat = ' && '.join([
      #'mass > 3.05 && mass < 3.15',
      'pt > 12',
      'probe_pt > 2.5',
      'lxy_sig > 10',
      #'abs(probe_dxy_bs) > 0.015',
      'probe_isloose == 1',
    ])
)

selection['baseline_sources'] = Selection(
    flat = ' && '.join([
      'hnl_charge == 0',
      'sv_lxy < 100',
    ])
)

selection['mass_charge'] = Selection(
    flat = ' && '.join([
      'b_mass < 6.4',
      'hnl_charge == 0',
    ])
)

selection['mass_charge_ID_nodsa'] = Selection(
    flat = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && trgmu_softid==1 && mu_looseid==1'
)

selection['mass_charge_ID_pipt_nodsa'] = Selection(
    flat = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && pi_pt>0.8 && trgmu_softid==1 && mu_looseid==1'
)

selection['sensitivity_study_July'] = Selection(
    #flat = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'
    flat = 'hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))'
)

selection['sensitivity_study_July_bmass_hnlcharge'] = Selection(
    flat = 'hnl_cos2d>0.993 && sv_prob>0.05 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))'
)

selection['sensitivity_study_July_TF'] = Selection(
    flat = 'mu_isdsa!=1 && b_mass<6.4 && hnl_cos2d>0.993 && sv_prob>0.05 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))'
)

selection['13Oct21'] = Selection(
    flat = 'hnl_charge==0 && b_mass<6.4 && ((mu_isdsa!=1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))) || (mu_isdsa==1 && mu_passdsaid==1))'
    #flat = '((mu_isdsa!=1 && trgmu_softid==1 && mu_looseid==1) || (mu_isdsa==1 && mu_passdsaid==1))'
)





