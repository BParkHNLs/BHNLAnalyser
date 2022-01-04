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
    nano = 'BToMuMuPi_hnl_pt > 0',
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





