'''
  This class defines the final state categorisation
  As for the definition of the category, both the 'flat' and 'nano' syntaxes can be used
'''

class Category(object):
  def __init__(self, label='', title='', definition_flat='', definition_nano='', cutbased_selection=''):
    self.label = label
    self.title = title
    self.definition_flat = definition_flat
    self.definition_nano = definition_nano
    self.cutbased_selection = cutbased_selection


categories = {}
categories['inclusive'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
          ),
]

categories['standard'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
]
