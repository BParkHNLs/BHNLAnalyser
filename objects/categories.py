import os


class Category(object):
  def __init__(self, label='', title='', definition_flat='', definition_nano=''):
    self.label = label
    self.title = title
    self.definition_flat = definition_flat
    self.definition_nano = definition_nano




categories = {}
categories['standard'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
          ),
]
