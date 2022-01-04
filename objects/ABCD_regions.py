'''
  This class defines the different ABCD categorisations
'''

class ABCD(object):
  def __init__(self, label='',  SR_selection='', CR_B_selection='', CR_C_selection='', CR_D_selection=''):
    self.label = label
    self.SR_selection = SR_selection
    self.CR_B_selection = CR_B_selection
    self.CR_C_selection = CR_C_selection
    self.CR_D_selection = CR_D_selection



ABCD_regions = {}
ABCD_regions['cos2d_svprob'] = ABCD(
    label = 'cos2d_svprob',
    SR_selection = 'hnl_cos2d > 0.993 && sv_prob > 0.05',
    CR_B_selection = 'hnl_cos2d > 0.993 && sv_prob < 0.05',
    CR_C_selection = 'hnl_cos2d < 0.993 && sv_prob > 0.05',
    CR_D_selection = 'hnl_cos2d < 0.993 && sv_prob < 0.05'
    )

ABCD_regions['cos2d_svprob_0p996'] = ABCD(
    label = 'cos2d_svprob_0p996',
    SR_selection = 'hnl_cos2d > 0.996 && sv_prob > 0.05',
    CR_B_selection = 'hnl_cos2d > 0.996 && sv_prob < 0.05',
    CR_C_selection = 'hnl_cos2d < 0.996 && sv_prob > 0.05',
    CR_D_selection = 'hnl_cos2d < 0.996 && sv_prob < 0.05'
    )

ABCD_regions['bmass_hnlcharge'] = ABCD(
    label = 'bmass_hnlcharge',
    SR_selection = 'b_mass < 6.27 && hnl_charge == 0',
    CR_B_selection = 'b_mass < 6.27 && hnl_charge != 0',
    CR_C_selection = 'b_mass > 6.27 && hnl_charge == 0',
    CR_D_selection = 'b_mass > 6.27 && hnl_charge != 0'
    )

ABCD_regions['cos2d_hnlcharge'] = ABCD(
    label = 'cos2d_hnlcharge',
    SR_selection = 'hnl_cos2d > 0.998 && hnl_charge == 0',
    CR_B_selection = 'hnl_cos2d > 0.998 && hnl_charge != 0',
    CR_C_selection = 'hnl_cos2d < 0.998 && hnl_charge == 0',
    CR_D_selection = 'hnl_cos2d < 0.998 && hnl_charge != 0'
    )
