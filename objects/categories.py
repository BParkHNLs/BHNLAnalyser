'''
  This class defines the final state categorisation
  As for the definition of the category, both the 'flat' and 'nano' syntaxes can be used
'''

class Category(object):
  def __init__(self, label='', title='', definition_flat='', definition_nano='', cutbased_selection='', cutbased_selection_lowmass='', cutbased_selection_mediummass='', cutbased_selection_highmass=''):
    self.label = label
    self.title = title
    self.definition_flat = definition_flat
    self.definition_nano = definition_nano
    self.cutbased_selection = cutbased_selection
    self.cutbased_selection_lowmass = cutbased_selection_lowmass
    self.cutbased_selection_mediummass = cutbased_selection_mediummass
    self.cutbased_selection_highmass = cutbased_selection_highmass


categories = {}
categories['control_Bc'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'b_pt > 0',
           cutbased_selection = 'b_pt > 0',
          ),
]

categories['inclusive'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0',
          ),
]

categories['tag_and_probe'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'mass > 0',
           cutbased_selection = 'mass>0',
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
           cutbased_selection = 'mu_isdsa!=1 && b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_isdsa!=1 && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_isdsa!=1 && b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_isdsa!=1 && b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
]

categories['standard_sources'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30',
           cutbased_selection = 'hnl_charge==0',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>75',
           cutbased_selection = 'hnl_charge==0',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>150',
           cutbased_selection = 'hnl_charge==0',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30',
           cutbased_selection = 'hnl_charge==0',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>75',
           cutbased_selection = 'hnl_charge==0',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>150',
           cutbased_selection = 'hnl_charge==0',
          ),
]

categories['all_categories'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'sv_lxysig>25 && pi_dcasig>9'
          ),
  Category(label = 'lxy1to3_OS',
           title = '(1<l_{xy}<=3)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_dxy>0.01 && sv_lxysig>40 && pi_dcasig>20'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxy1to10_OS',
           title = '(1<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy1to20_OS',
           title = '(1<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy3to5_OS',
           title = '(3<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to10_OS',
           title = '(3<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to20_OS',
           title = '(3<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxy5to10_OS',
           title = '(5<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>5 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy5to20_OS',
           title = '(5<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>5 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_pt > 0' # 'mu_isdsa!=1 && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_OS',
           title = 'l_{xy}>10cm, OS',
           definition_flat = 'sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<3 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<3 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && pi_pt>1.3'
          ),
  Category(label = 'lxygt20_OS',
           title = 'l_{xy}>20cm, OS',
           definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<2 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<2 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && mu_trkisoid==1 && pi_pt>1.8'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'sv_lxysig>25 && pi_dcasig>9'
          ),
  Category(label = 'lxy1to3_SS',
           title = '(1<l_{xy}<=3)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_dxy>0.01 && sv_lxysig>40 && pi_dcasig>20'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxy1to10_SS',
           title = '(1<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy1to20_SS',
           title = '(1<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy3to5_SS',
           title = '(3<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to10_SS',
           title = '(3<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to20_SS',
           title = '(3<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxy5to10_SS',
           title = '(5<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>5 && sv_lxy<=10 && trgmu_charge==mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy5to20_SS',
           title = '(5<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>5 && sv_lxy<=20 && trgmu_charge==mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_SS',
           title = 'l_{xy}>10cm, SS',
           definition_flat = 'sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<3 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<3 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && pi_pt>1.3'
          ),
  Category(label = 'lxygt20_SS',
           title = 'l_{xy}>20cm, SS',
           definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<2 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<2 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && mu_trkisoid==1 && pi_pt>1.8'
          ),
]

categories['study_Nov21'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'sv_lxysig>25 && pi_dcasig>9'
          ),
  Category(label = 'lxy1to3_OS',
           title = '(1<l_{xy}<=3)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_dxy>0.01 && sv_lxysig>40 && pi_dcasig>20'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxy1to10_OS',
           title = '(1<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy1to20_OS',
           title = '(1<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy3to5_OS',
           title = '(3<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to10_OS',
           title = '(3<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to20_OS',
           title = '(3<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_pt > 0' # 'mu_isdsa!=1 && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_OS',
           title = 'l_{xy}>10cm, OS',
           definition_flat = 'sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<3 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<3 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && pi_pt>1.3'
          ),
  Category(label = 'lxygt20_OS',
           title = 'l_{xy}>20cm, OS',
           definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<2 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<2 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && mu_trkisoid==1 && pi_pt>1.8'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'sv_lxysig>25 && pi_dcasig>9'
          ),
  Category(label = 'lxy1to3_SS',
           title = '(1<l_{xy}<=3)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_dxy>0.01 && sv_lxysig>40 && pi_dcasig>20'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxy1to10_SS',
           title = '(1<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy1to20_SS',
           title = '(1<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'abs(cos_theta_star_pion)<0.87 && mu_dxy>0.01 && pi_pt>0.9 && sv_lxysig>60 && pi_dcasig>25'
          ),
  Category(label = 'lxy3to5_SS',
           title = '(3<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to10_SS',
           title = '(3<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && mu_numberofvalidpixelhits<5 && abs(cos_theta_star_pion)<0.85 && mu_dxy>0.03 && pi_pt>0.9 && sv_lxysig>50'
          ),
  Category(label = 'lxy3to20_SS',
           title = '(3<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<4 && abs(cos_theta_star_pion)<0.87 && mu_dxy>0.03 && pi_pt>0.8 && sv_lxysig>50'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_pt > 0' #'mu_isdsa!=1 && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_SS',
           title = 'l_{xy}>10cm, SS',
           definition_flat = 'sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<3 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<3 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && pi_pt>1.3'
          ),
  Category(label = 'lxygt20_SS',
           title = 'l_{xy}>20cm, SS',
           definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'mu_numberofpixellayers<2 && mu_numberoftrackerlayers<13 && mu_numberofvalidpixelhits<2 && abs(cos_theta_star_pion)<0.8 && mu_dxy>0.03 && mu_trkisoid==1 && pi_pt>1.8'
          ),
]

categories['3cat_0_1_5_benchmark'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9997 && pi_pt>1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9997 && pi_pt>1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
]

categories['3cat_0_1_10_benchmark'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to10_OS',
           title = '(1<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9998 && pi_pt>1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.04 && sv_lxysig>70'
          ),
  Category(label = 'lxygt10_OS',
           title = 'l_{xy}>10cm, OS',
           definition_flat = 'sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to10_SS',
           title = '(1<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9998 && pi_pt>1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.04 && sv_lxysig>70'
          ),
  Category(label = 'lxygt10_SS',
           title = 'l_{xy}>10cm, SS',
           definition_flat = 'sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
]

categories['3cat_0_1_20_benchmark'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to20_OS',
           title = '(1<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9999 && pi_pt>1.1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.04 && sv_lxysig>120'
          ),
  Category(label = 'lxygt20_OS',
           title = 'l_{xy}>20cm, OS',
           definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99998 && pi_pt>1.8 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to20_SS',
           title = '(1<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9999 && pi_pt>1.1 && abs(mu_dxy)>0.03 && abs(pi_dxy)>0.04 && sv_lxysig>120'
          ),
  Category(label = 'lxygt20_SS',
           title = 'l_{xy}>20cm, SS',
           definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99998 && pi_pt>1.8 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
]

categories['4cat_0_1_3_10_benchmark'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to3_OS',
           title = '(1<l_{xy}<=3)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9994 && pi_pt>0.9 && abs(mu_dxy)>0.02 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxy3to10_OS',
           title = '(3<l_{xy}<=10)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99991 && pi_pt>0.9 && abs(mu_dxy)>0.06 && abs(pi_dxy)>0.17 && sv_lxysig>70'
          ),
  Category(label = 'lxygt10_OS',
           title = 'l_{xy}>10cm, OS',
           definition_flat = 'sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to3_SS',
           title = '(1<l_{xy}<=3)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9994 && pi_pt>0.9 && abs(mu_dxy)>0.02 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxy3to10_SS',
           title = '(3<l_{xy}<=10)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=10 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99991 && pi_pt>0.9 && abs(mu_dxy)>0.06 && abs(pi_dxy)>0.17 && sv_lxysig>70'
          ),
  Category(label = 'lxygt10_SS',
           title = 'l_{xy}>10cm, SS',
           definition_flat = 'sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
          ),
]

categories['4cat_0_1_3_20_benchmark'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to3_OS',
           title = '(1<l_{xy}<=3)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9994 && pi_pt>0.9 && abs(mu_dxy)>0.02 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxy3to20_OS',
           title = '(3<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99994 && pi_pt>1.0 && abs(mu_dxy)>0.07 && abs(pi_dxy)>0.2 && sv_lxysig>70'
          ),
  Category(label = 'lxygt20_OS',
           title = 'l_{xy}>20cm, OS',
           definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99998 && pi_pt>1.8 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to3_SS',
           title = '(1<l_{xy}<=3)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=3 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9994 && pi_pt>0.9 && abs(mu_dxy)>0.02 && abs(pi_dxy)>0.03 && sv_lxysig>70'
          ),
  Category(label = 'lxy3to20_SS',
           title = '(3<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>3 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99994 && pi_pt>1.0 && abs(mu_dxy)>0.07 && abs(pi_dxy)>0.2 && sv_lxysig>70'
          ),
  Category(label = 'lxygt20_SS',
           title = 'l_{xy}>20cm, SS',
           definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.99998 && pi_pt>1.8 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
]

categories['3cat_0_1_5_significance'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.8 && sv_lxysig>25'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10'
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.9997 && pi_pt>0.9 && abs(mu_dxy)>0.04 && abs(pi_dxy)>0.08 && sv_lxysig>90'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25'
           #cutbased_selection = 'hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>11 && abs(pi_dxysig)>24' # used for gt5forGeV
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20'
           #cutbased_selection = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>150 && abs(mu_dxysig)>20'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.8 && sv_lxysig>25'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10'
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.9997 && pi_pt>0.9 && abs(mu_dxy)>0.04 && abs(pi_dxy)>0.08 && sv_lxysig>90'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25'
           #cutbased_selection = 'hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>11 && abs(pi_dxysig)>24'
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           #cutbased_selection = 'hnl_cos2d>0.99997 && pi_pt>1.3 && abs(mu_dxy)>0.1 && abs(pi_dxy)>0.2'
           cutbased_selection = 'hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20'
           #cutbased_selection = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>150 && abs(mu_dxysig)>20'
          ),
]

categories['3cat_0_1_5_significance_permass'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<1e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1 && sv_lxysig>90 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<1e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>6 && sv_lxysig>30',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>20 && (1-hnl_cos2d)<1e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>270 && min(abs(mu_dxysig), abs(pi_dxysig))>25 && (1-hnl_cos2d)<1e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4.5 && sv_lxysig>75',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>230 && min(abs(mu_dxysig), abs(pi_dxysig))>60 && (1-hnl_cos2d)<1e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4 && sv_lxysig>150',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<1e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1 && sv_lxysig>90 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<1e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>6 && sv_lxysig>30',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>20 && (1-hnl_cos2d)<1e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>25 && (1-hnl_cos2d)<1e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4 && sv_lxysig>75',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>230 && min(abs(mu_dxysig), abs(pi_dxysig))>60 && (1-hnl_cos2d)<1e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>20 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4 && sv_lxysig>150',
          ),
]

categories['V12_08Aug22_permass'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && mu0_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.15 && pi_pt>0.8 && sv_lxysig>90 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.15 && pi_pt>2. && sv_lxysig>30 && (1-hnl_cos2d)<2e-3',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.2 && pi_pt>1.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.2 && pi_pt>3. && sv_lxysig>50 && (1-hnl_cos2d)<2e-4',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && mu0_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.35 && pi_pt>2.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.35 && pi_pt>3.5 && sv_lxysig>50 && (1-hnl_cos2d)<2e-5',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && mu0_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.15 && pi_pt>0.8 && sv_lxysig>100 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.15 && pi_pt>2. && sv_lxysig>30 && (1-hnl_cos2d)<2e-3',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && mu0_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.2 && pi_pt>1.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.2 && pi_pt>3. && sv_lxysig>50 && (1-hnl_cos2d)<2e-4',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && mu0_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && abs(b_mass-5.3)<0.35 && pi_pt>2.5 && sv_lxysig>150 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && abs(b_mass-5.3)<0.35 && pi_pt>3.5 && sv_lxysig>50 && (1-hnl_cos2d)<2e-5',
          ),
]


categories['V11_24Apr22_permass'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>90 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>100',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>60 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>200',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>100 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>70',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<3e-4',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>300',
          ),
]

categories['V11_24Apr22_permass_Bc'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>90 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>100 && b_mass<5.9',
          ),
  Category(label = 'lxy1to5_OS',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<2e-4 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100 && b_mass<5.9',
          ),
  Category(label = 'lxygt5_OS',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>60 && (1-hnl_cos2d)<2e-5 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>200 && b_mass<5.9',
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>100 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>70 && b_mass<5.9',
          ),
  Category(label = 'lxy1to5_SS',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<3e-4 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100 && b_mass<5.9',
          ),
  Category(label = 'lxygt5_SS',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5 && b_mass<5.9',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-5 && b_mass<5.9',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>300 && b_mass<5.9',
          ),
  #
  Category(label = 'lxy0to1_OS_Bc',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>30 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>90 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>100 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && sv_lxysig>80 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>50 && fabs(b_mass-6.27)<0.5',
          ),
  Category(label = 'lxy1to5_OS_Bc',
           title = '(1<l_{xy}<=5)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<2e-4 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>200 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && sv_lxysig>30 && fabs(b_mass-6.27)<0.5',
          ),
  Category(label = 'lxygt5_OS_Bc',
           title = 'l_{xy}>5cm, OS',
           definition_flat = 'sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>60 && (1-hnl_cos2d)<2e-5 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>200 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>200 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>3 && sv_lxysig>30 && fabs(b_mass-6.27)<0.5',
          ),
  Category(label = 'lxy0to1_SS_Bc',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>50 && min(abs(mu_dxysig), abs(pi_dxysig))>0 && (1-hnl_cos2d)<2e-3 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>0.8 && sv_lxysig>100 && min(abs(mu_dxysig), abs(pi_dxysig))>5 && (1-hnl_cos2d)<2e-3 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>70 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && sv_lxysig>80 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>50 && fabs(b_mass-6.27)<0.5',
          ),
  Category(label = 'lxy1to5_SS_Bc',
           title = '(1<l_{xy}<=5)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>0.7 && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>10 && (1-hnl_cos2d)<3e-4 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>300 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-4 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>4. && sv_lxysig>100 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>1.5 && sv_lxysig>200 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && sv_lxysig>30 && fabs(b_mass-6.27)<0.5',
          ),
  Category(label = 'lxygt5_SS_Bc',
           title = 'l_{xy}>5cm, SS',
           definition_flat = 'sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection_lowmass = 'hnl_charge==0 && pi_pt>2. && sv_lxysig>150 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<3e-5 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>3.5 && sv_lxysig>200 && min(abs(mu_dxysig), abs(pi_dxysig))>40 && (1-hnl_cos2d)<2e-5 && fabs(b_mass-6.27)<0.5',
           #cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>5. && sv_lxysig>300 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_mediummass = 'hnl_charge==0 && pi_pt>2 && sv_lxysig>200 && fabs(b_mass-6.27)<0.5',
           cutbased_selection_highmass = 'hnl_charge==0 && pi_pt>3 && sv_lxysig>30 && fabs(b_mass-6.27)<0.5',
          ),
]


categories['3cat_0_1_20_significance'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS',
           title = 'l_{xy}<=1cm, OS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to20_OS',
           title = '(1<l_{xy}<=20)cm, OS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9999 && pi_pt>1.1 && abs(mu_dxy)>0.05 && abs(pi_dxy)>0.12 && sv_lxysig>120'
          ),
  Category(label = 'lxygt20_OS',
           title = 'l_{xy}>20cm, OS',
           definition_flat = 'sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_cos2d>0.999993 && pi_pt>2 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
  Category(label = 'lxy0to1_SS',
           title = 'l_{xy}<=1cm, SS',
           definition_flat = 'sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.998 && pi_pt>0.9 && sv_lxysig>25'
          ),
  Category(label = 'lxy1to20_SS',
           title = '(1<l_{xy}<=20)cm, SS',
           definition_flat = 'sv_lxy>1 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.9999 && pi_pt>1.1 && abs(mu_dxy)>0.05 && abs(pi_dxy)>0.12 && sv_lxysig>120'
          ),
  Category(label = 'lxygt20_SS',
           title = 'l_{xy}>20cm, SS',
           definition_flat = 'sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_cos2d>0.999993 && pi_pt>2 && abs(mu_dxy)>0.25 && abs(pi_dxy)>0.5'
          ),
]


categories['13Oct21_slimmed'] = [
  Category(label = 'incl_slimmed',
           title = 'inclusive, slimmed',
           definition_flat = 'hnl_pt > 0 && mu_isdsa!=1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS_slimmed',
           title = 'l_{xy}<=1cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>0.75 && b_mass>2.0 && pi_dcasig>6 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxy1to5_OS_slimmed',
           title = '(1<l_{xy}<=5)cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>0.9 && b_mass>1.7 && pi_dcasig>8 && sv_lxysig>30 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxygt5_OS_slimmed',
           title = 'l_{xy}>5cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.1 && pi_dcasig>8 && sv_lxysig>40 && deltar_mu_pi>0.1'
          ),
  #Category(label = 'lxygt20_OS_slimmed',
  #         title = 'l_{xy}>20cm, OS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'pi_pt>2 && pi_dcasig>8 && sv_lxysig>130 && deltar_mu_pi>0.1'
  #        ),
  Category(label = 'lxy0to1_SS_slimmed',
           title = 'l_{xy}<=1cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.8 && b_mass>2.8 && pi_dcasig>6 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxy1to5_SS_slimmed',
           title = '(1<l_{xy}<=5)cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.9 && b_mass>1.7 && pi_dcasig>9 && sv_lxysig>30 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxygt5_SS_slimmed',
           title = 'l_{xy}>5cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  #Category(label = 'lxygt20_SS_slimmed',
  #         title = 'l_{xy}>20cm, SS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'pi_pt>2 && pi_dcasig>8 && sv_lxysig>130'
  #        ),
]

categories['combined_dsa'] = [
  #Category(label = 'incl',
  #         title = 'inclusive',
  #         definition_flat = 'hnl_pt > 0',
  #         definition_nano = 'BToMuMuPi_hnl_pt > 0',
  #         cutbased_selection = 'hnl_pt > 0'
  #        ),
  #Category(label = 'lxy0to1_OS_dsa',
  #         title = 'l_{xy}<=1cm, OS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'hnl_pt>0'
  #        ),
  #Category(label = 'lxy1to5_OS_dsa',
  #         title = '(1<l_{xy}<=5)cm, OS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'deltaeta_mu_pi>0.05 && deltaeta_trgmu_mu>0.05 && deltar_mu_pi>0.15 && deltar_trgmu_mu>0.05 && pi_pt>1 && sv_prob>0.03'
  #        ),
  #Category(label = 'lxygt5_OS_dsa',
  #         title = 'l_{xy}>5cm, OS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
  #        ),
  #Category(label = 'lxy0to1_SS_dsa',
  #         title = 'l_{xy}<=1cm, SS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'hnl_pt>0'
  #        ),
  #Category(label = 'lxy1to5_SS_dsa',
  #         title = '(1<l_{xy}<=5)cm, SS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'deltaeta_mu_pi>0.05 && deltaeta_trgmu_mu>0.05 && deltar_mu_pi>0.15 && deltar_trgmu_mu>0.05 && pi_pt>1 && sv_prob>0.03'
  #        ),
  #Category(label = 'lxygt5_SS_dsa',
  #         title = 'l_{xy}>5cm, SS, dsa',
  #         definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
  #        ),
  Category(label = 'lxy0to1_OS_nodsa',
           title = 'l_{xy}<=1cm, OS, no dsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_OS_nodsa',
           title = '(1<l_{xy}<=5)cm, OS, no dsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_OS_nodsa',
           title = 'l_{xy}>5cm, OS, no dsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxy0to1_SS_nodsa',
           title = 'l_{xy}<=1cm, SS, no dsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_SS_nodsa',
           title = '(1<l_{xy}<=5)cm, SS, no dsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_SS_nodsa',
           title = 'l_{xy}>5cm, SS, nodsa',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
]

categories['combined_dsa_add'] = [
  Category(label = 'lxy5to20_OS_slimmed',
           title = '(5<l_{xy}<=20)cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
  Category(label = 'lxy5to20_SS_slimmed',
           title = '(5<l_{xy}<=20)cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
  Category(label = 'lxy5to20_OS_dsa',
           title = '(5<l_{xy}<=20)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
  Category(label = 'lxy5to20_SS_dsa',
           title = '(5<l_{xy}<=20)cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
]


categories['category_study_combined_dsa_slimmed'] = [
  Category(label = 'incl_slimmed',
           title = 'inclusive, slimmed',
           definition_flat = 'hnl_pt > 0 && mu_isdsa!=1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  #Category(label = 'lxy0to1_OS_slimmed',
  #         title = 'l_{xy}<=1cm, OS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'pi_pt>0.75 && b_mass>2.0 && pi_dcasig>6 && deltar_mu_pi>0.1'
  #        ),
  #Category(label = 'lxy1to5_OS_slimmed',
  #         title = '(1<l_{xy}<=5)cm, OS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
  #         cutbased_selection = 'pi_pt>0.9 && b_mass>1.7 && pi_dcasig>8 && sv_lxysig>30 && deltar_mu_pi>0.1'
  #        ),
  Category(label = 'lxygt5_OS_slimmed',
           title = 'l_{xy}>5cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.1 && pi_dcasig>8 && sv_lxysig>40 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxygt10_OS_slimmed',
           title = 'l_{xy}>10cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxygt15_OS_slimmed',
           title = 'l_{xy}>15cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>15 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxygt20_OS_slimmed',
           title = 'l_{xy}>20cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>2 && pi_dcasig>8 && sv_lxysig>130 && deltar_mu_pi>0.1'
          ),
  Category(label = 'lxygt25_OS_slimmed',
           title = 'l_{xy}>25cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>25 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt > 0'
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt30_OS_slimmed',
           title = 'l_{xy}>30cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt > 0'
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  #Category(label = 'lxy0to1_SS_slimmed',
  #         title = 'l_{xy}<=1cm, SS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'pi_pt>0.8 && b_mass>2.8 && pi_dcasig>6 && deltar_mu_pi>0.1'
  #        ),
  #Category(label = 'lxy1to5_SS_slimmed',
  #         title = '(1<l_{xy}<=5)cm, SS, slimmed',
  #         definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
  #         cutbased_selection = 'pi_pt>0.9 && b_mass>1.7 && pi_dcasig>9 && sv_lxysig>30 && deltar_mu_pi>0.1'
  #        ),
  Category(label = 'lxygt5_SS_slimmed',
           title = 'l_{xy}>5cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_SS_slimmed',
           title = 'l_{xy}>10cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxygt15_SS_slimmed',
           title = 'l_{xy}>15cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>15 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxygt20_SS_slimmed',
           title = 'l_{xy}>20cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>2 && pi_dcasig>8 && sv_lxysig>130'
          ),
  Category(label = 'lxygt25_SS_slimmed',
           title = 'l_{xy}>25cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>25 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxygt30_SS_slimmed',
           title = 'l_{xy}>30cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt > 0'
          ),
]

categories['category_study_combined_dsa_dsa'] = [
  Category(label = 'incl_dsa',
           title = 'inclusive, dsa',
           definition_flat = 'hnl_pt > 0 && mu_isdsa==1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS_dsa',
           title = 'l_{xy}<=1cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>2.7 && b_mass>4.75 && pi_dcasig>8 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy1to5_OS_dsa',
           title = '(1<l_{xy}<=5)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.5 && b_mass>2.3 && pi_dcasig>3 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt5_OS_dsa',
           title = 'l_{xy}>5cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>0.9 && mu_dxy>0.2 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt10_OS_dsa',
           title = 'l_{xy}>10cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>10 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.1 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt15_OS_dsa',
           title = 'l_{xy}>15cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>15 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.1 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt20_OS_dsa',
           title = 'l_{xy}>20cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.2 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt25_OS_dsa',
           title = 'l_{xy}>25cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>25 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.2 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt30_OS_dsa',
           title = 'l_{xy}>30cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>30 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.2 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy0to1_SS_dsa',
           title = 'l_{xy}<=1cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>2.7 && b_mass>4.75 && pi_dcasig>8 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy1to5_SS_dsa',
           title = '(1<l_{xy}<=5)cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1.5 && b_mass>2.3 && pi_dcasig>3 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt5_SS_dsa',
           title = 'l_{xy}>5cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.9 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt10_SS_dsa',
           title = 'l_{xy}>10cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>10 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.9 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt15_SS_dsa',
           title = 'l_{xy}>15cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>15 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.9 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt20_SS_dsa',
           title = 'l_{xy}>20cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt25_SS_dsa',
           title = 'l_{xy}>25cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>25 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1 && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt30_SS_dsa',
           title = 'l_{xy}>30cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>30 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1 && mu_ismatchedtoslimmedmuon==0'
          ),
]

categories['category_study_combined_dsa'] = [
  Category(label = 'incl',
           title = 'inclusive',
           definition_flat = 'hnl_pt > 0',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'incl_slimmed',
           title = 'inclusive, slimmed',
           definition_flat = 'hnl_pt > 0 && mu_isdsa!=1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           #cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS_slimmed',
           title = 'l_{xy}<=1cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_OS_slimmed',
           title = '(1<l_{xy}<=5)cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_OS_slimmed',
           title = 'l_{xy}>5cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_OS_slimmed',
           title = 'l_{xy}>10cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>10 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt15_OS_slimmed',
           title = 'l_{xy}>15cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>15 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt20_OS_slimmed',
           title = 'l_{xy}>20cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt25_OS_slimmed',
           title = 'l_{xy}>25cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>25 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt30_OS_slimmed',
           title = 'l_{xy}>30cm, OS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxy0to1_SS_slimmed',
           title = 'l_{xy}<=1cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
          ),
  Category(label = 'lxy1to5_SS_slimmed',
           title = '(1<l_{xy}<=5)cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
          ),
  Category(label = 'lxygt5_SS_slimmed',
           title = 'l_{xy}>5cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt10_SS_slimmed',
           title = 'l_{xy}>10cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>10 && trgmu_charge==mu_charge',
           #cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt15_SS_slimmed',
           title = 'l_{xy}>15cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>15 && trgmu_charge==mu_charge',
           #cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt20_SS_slimmed',
           title = 'l_{xy}>20cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>20 && trgmu_charge==mu_charge',
           #cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt25_SS_slimmed',
           title = 'l_{xy}>25cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>25 && trgmu_charge==mu_charge',
           #cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'lxygt30_SS_slimmed',
           title = 'l_{xy}>30cm, SS, slimmed',
           definition_flat = 'mu_isdsa!=1 && sv_lxy>30 && trgmu_charge==mu_charge',
           #cutbased_selection = 'b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
          ),
  Category(label = 'incl_dsa',
           title = 'inclusive, dsa',
           definition_flat = 'hnl_pt > 0 && mu_isdsa==1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS_dsa',
           title = 'l_{xy}<=1cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
  Category(label = 'lxy1to5_OS_dsa',
           title = '(1<l_{xy}<=5)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'deltaeta_mu_pi>0.05 && deltaeta_trgmu_mu>0.05 && deltar_mu_pi>0.15 && deltar_trgmu_mu>0.05 && pi_pt>1 && sv_prob>0.03'
          ),
  Category(label = 'lxygt5_OS_dsa',
           title = 'l_{xy}>5cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt10_OS_dsa',
           title = 'l_{xy}>10cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>10 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt15_OS_dsa',
           title = 'l_{xy}>15cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>15 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt20_OS_dsa',
           title = 'l_{xy}>20cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt25_OS_dsa',
           title = 'l_{xy}>25cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>25 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt30_OS_dsa',
           title = 'l_{xy}>30cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>30 && trgmu_charge!=mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxy0to1_SS_dsa',
           title = 'l_{xy}<=1cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'hnl_pt>0'
          ),
  Category(label = 'lxy1to5_SS_dsa',
           title = '(1<l_{xy}<=5)cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'deltaeta_mu_pi>0.05 && deltaeta_trgmu_mu>0.05 && deltar_mu_pi>0.15 && deltar_trgmu_mu>0.05 && pi_pt>1 && sv_prob>0.03'
          ),
  Category(label = 'lxygt5_SS_dsa',
           title = 'l_{xy}>5cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && trgmu_charge==mu_charge',
           cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt10_SS_dsa',
           title = 'l_{xy}>10cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>10 && trgmu_charge==mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt15_SS_dsa',
           title = 'l_{xy}>15cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>15 && trgmu_charge==mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt20_SS_dsa',
           title = 'l_{xy}>20cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge==mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt25_SS_dsa',
           title = 'l_{xy}>25cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>25 && trgmu_charge==mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
  Category(label = 'lxygt30_SS_dsa',
           title = 'l_{xy}>30cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>30 && trgmu_charge==mu_charge',
           #cutbased_selection = 'deltaeta_trgmu_mu>0.06 && deltaeta_mu_pi>0.05 && deltaeta_trgmu_pi>0.03 && deltaphi_trgmu_pi>0.01 && deltar_trgmu_mu>0.1 && mu_dxy>0.2 && pi_pt>0.9 && sv_prob>0.05'
          ),
]

categories['sensitivity_study_dsa_only'] = [
  Category(label = 'incl_dsa',
           title = 'inclusive, dsa',
           definition_flat = 'hnl_pt > 0 && mu_isdsa==1',
           definition_nano = 'BToMuMuPi_hnl_pt > 0',
           cutbased_selection = 'hnl_pt > 0'
          ),
  Category(label = 'lxy0to1_OS_dsa',
           title = 'l_{xy}<=1cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>2.7 && b_mass>4.75 && pi_dcasig>8'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy1to5_OS_dsa',
           title = '(1<l_{xy}<=5)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.5 && b_mass>2.3 && pi_dcasig>3'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy5to20_OS_dsa',
           title = '(5<l_{xy}<=20)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>0.8 && b_mass>1.8 && pi_dcasig>3'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt20_OS_dsa',
           title = 'l_{xy}>20cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge!=mu_charge',
           cutbased_selection = 'pi_pt>1.2'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy0to1_SS_dsa',
           title = 'l_{xy}<=1cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy<=1 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>2.7 && b_mass>4.75 && pi_dcasig>8'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy1to5_SS_dsa',
           title = '(1<l_{xy}<=5)cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1.5 && b_mass>2.3 && pi_dcasig>3'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxy5to20_SS_dsa',
           title = '(5<l_{xy}<=20)cm, OS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>5 && sv_lxy<=20 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>0.8 && b_mass>1.8 && pi_dcasig>3'# && mu_ismatchedtoslimmedmuon==0'
          ),
  Category(label = 'lxygt20_SS_dsa',
           title = 'l_{xy}>20cm, SS, dsa',
           definition_flat = 'mu_isdsa==1 && sv_lxy>20 && trgmu_charge==mu_charge',
           cutbased_selection = 'pi_pt>1'# && mu_ismatchedtoslimmedmuon==0'
          ),
]


