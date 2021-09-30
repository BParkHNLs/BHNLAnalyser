import os
from os import path

import ROOT
import math

from computeYields import ComputeYields
from tools import Tools
from samples import data_samples, signal_samples, signal_samples_limits, qcd_samples

#import sys
#sys.path.append('../../../../BHNLGen/CMSSW_10_2_15/src/HNLsGen/python/.')
#from my_common import getVV


class DatacardsMaker(Tools):
  def __init__(self, data_file='', signal_files='', do_shape_analysis=True, do_categories=True, categories=None, add_Bc=False, outdirlabel=''):
    self.tools = Tools()
    self.data_file = data_file
    self.signal_files = signal_files 
    self.do_shape_analysis = do_shape_analysis
    self.do_categories = do_categories
    self.categories = categories
    if do_categories and categories == None:
      raise RuntimeError('Please fill the dictionnary with the different categories')
    self.add_Bc = add_Bc
    self.outputdir = './datacards/{}'.format(outdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))


  def getCouplingLabel(self, v2):
    coupling = "{:e}".format(v2)
    part1 = coupling[:coupling.find('e')]
    part1 = str(round(float(part1), 1))
    part2 = coupling[coupling.find('e'):]
    return (part1+part2) #.replace('.', 'p').replace('-', 'm')


  def getBackgroundYields(self, data_file, signal_file, baseline_selection='', category=''):
    background_selection = baseline_selection
    if self.do_categories:
      category_selection = self.categories[category]
      background_selection += ' && {}'.format(category_selection)

    # ABCD method
    background_yields = ComputeYields(data_file=data_file, signal_file=signal_file, selection=background_selection).computeBkgYieldsFromABCDData()[0] 

    # TF method
    #from samples import qcd_samples
    #white_list_15to300 = ['QCD_pt15to20 (V06_29Jun21)', 'QCD_pt20to30 (V06_29Jun21)', 'QCD_pt30to50 (V06_29Jun21)', 'QCD_pt50to80 (V06_29Jun21)', 'QCD_pt80to120 (V06_29Jun21)', 'QCD_pt120to170 (V06_29Jun21)', 'QCD_pt170to300 (V06_29Jun21)']
    #background_yields = ComputeYields(data_file=data_file, qcd_files=qcd_samples, signal_file=signal_file, selection=background_selection, white_list=white_list_15to300).computeBkgYieldsFromMC()[0]

    if background_yields == 0.: background_yields = 1e-9

    if background_yields != 1e-9: background_yields = background_yields * 41.6/(0.774 * 0.1) # added 10percent since we only use Chunk2

    return background_yields


  def getSignalYields(self, signal_file, baseline_selection='', category=''):
    signal_selection = 'ismatched==1 && {}'.format(baseline_selection)
    if self.do_categories:
      category_selection = self.categories[category]
      signal_selection += ' && {}'.format(category_selection)

    signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=False)[0] 
    #signal_yields = ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields()[0]
    if self.add_Bc and signal_file.filename_Bc != '':
      signal_yields += ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=True)[0]

      #print 'yields {} + {} = {}'.format(ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=False)[0], ComputeYields(signal_file=signal_file, selection=signal_selection).computeSignalYields(lumi=41.6, isBc=True), signal_yields)[0]

    return signal_yields


  def getSignalMassCoupling(self, signal_point):
    signal_mass = signal_point.mass
    signal_ctau = signal_point.ctau
    signal_v2 = self.tools.getVV(mass=signal_mass, ctau=signal_ctau, ismaj=True)
    signal_coupling = self.getCouplingLabel(signal_v2)

    return signal_mass, signal_coupling


  def getLabel(self, signal_mass='', signal_coupling='', category=''):
    if not self.do_categories:
      label = 'bhnl_m_{}_v2_{}'.format(signal_mass, signal_coupling)
    else:
      label = 'bhnl_cat_{}_m_{}_v2_{}'.format(category, signal_mass, signal_coupling)
    return label


  def writeCard(self, label, signal_yields, background_yields):
    datacard_name = 'datacard_{}.txt'.format(label)
    datacard = open('{}/{}'.format(self.outputdir, datacard_name), 'w+')
    #--------------------------------------------------------------------------------------------------------------------------------------------
    #shapes *          {lbl}   shape_{lbl}.root   $PROCESS $PROCESS_$SYSTEMATIC
    datacard.write(
'''\
imax 1 number of bins
jmax 1 number of backgrounds
kmax * number of nuisance parameters
--------------------------------------------------------------------------------------------------------------------------------------------
bin               {lbl}
observation       {obs}
--------------------------------------------------------------------------------------------------------------------------------------------
bin                                                      {lbl}                          {lbl}             
process                                                  sig                            qcd              
process                                                  -1                             1                   
rate                                                     {sig_yields}                   {bkg_yields}    
--------------------------------------------------------------------------------------------------------------------------------------------
lumi                                       lnN           1.025                          -   
syst_sig_{lbl}                             lnN           1.3                            -    
syst_bkg_{lbl}                             lnN           -                              1.3   
--------------------------------------------------------------------------------------------------------------------------------------------
{lbl} autoMCStats 0 0 1
'''.format(
            lbl = label,
            obs =  -1, # for the moment, we only look at blinded data
            sig_yields = signal_yields,
            bkg_yields = background_yields,
        )
      )

    datacard.close()
    print '--> {}/{} created'.format(self.outputdir, datacard_name)


  def produceDatacard(self, signal_file='', category='', baseline_selection=''):
    # get the signal mass/coupling
    signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)

    # get the process label
    label = self.getLabel(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)

    # get the background yields
    background_yields = self.getBackgroundYields(data_file=self.data_file, signal_file=signal_file, baseline_selection=baseline_selection, category=category)

    # get the signal yields
    signal_yields = self.getSignalYields(signal_file=signal_file, baseline_selection=baseline_selection, category=category)

    #print '{} {}'.format(signal_mass, background_yields)

    # create the datacard
    self.writeCard(label=label, signal_yields=signal_yields, background_yields=background_yields)


  def process(self):
    #baseline_selection = '(trgmu_mu_mass<3.03 || trgmu_mu_mass>3.15) && hnl_charge==0 && b_mass<5.5'
    #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'
    #baseline_selection = 'trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))'
    # for bkg TF method
    #baseline_selection = 'mu_isdsa!=1 && b_mass<6.4 && hnl_cos2d>0.993 && sv_prob>0.05 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'

    #baseline_selection = 'hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8 && trgmu_looseid==1 && trgmu_softid==1 && mu_looseid==1 && mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))' # used for sensitivity
    #baseline_selection = 'mu_isdsa!=1 && hnl_charge==0 && b_mass<6.4 && deltar_mu_pi>0.1 && deltar_mu_pi<1.7 && deltar_trgmu_mu<1 && deltar_trgmu_pi<1.5 && pi_pt>0.8'
    baseline_selection = 'hnl_charge==0 && b_mass<6.4'

    if not self.do_categories:
      # loop on the signal points
      for signal_file in signal_files:
        self.produceDatacard(signal_file=signal_file, baseline_selection=baseline_selection)


        if self.do_shape_analysis:
          signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)
          label = self.getLabel(signal_mass=signal_mass, signal_coupling=signal_coupling)

          rootfile_name = 'shape_{}.root'.format(label)
          root_file = ROOT.TFile.Open('{}/{}'.format(self.outputdir, rootfile_name), 'RECREATE')  
          root_file.cd()

          from quantity import Quantity
          sigma = signal_file.resolution
          quantity = Quantity(name_flat='hnl_mass', nbins=80, bin_min=signal_mass-2*sigma, bin_max=signal_mass+2*sigma)

          # data
          data_hist = ROOT.TH1D('data_obs', 'data_obs', quantity.nbins, quantity.bin_min, quantity.bin_max)
          root_file.cd()
          data_hist.Write()

          # signal shape
          f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+signal_file.filename, 'READ')
          
          filename = signal_file.filename
          original_ctau = filename[filename.find('ctau')+4:filename.find('/', filename.find('ctau')+1)]
          target_ctau = signal_file.ctau
          ctau_weight = -99
          if float(original_ctau) != float(target_ctau):
            ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

          signal_selection = 'ismatched==1 && {}'.format(baseline_selection)
          hist_sig = self.tools.createHisto(f, 'signal_tree', quantity, hist_name='sig', branchname='flat', selection=signal_selection, weight=ctau_weight)
          
          signal_yields = self.getSignalYields(signal_file=signal_file, baseline_selection=baseline_selection)
          hist_sig.Scale(signal_yields/hist_sig.Integral())
          root_file.cd() 
          hist_sig.Write()


          # background shape
          int_mc_tot = 0.
          hist_mc_tot = ROOT.TH1D('qcd', 'qcd', quantity.nbins, quantity.bin_min, quantity.bin_max)
          hist_mc_tot.Sumw2()

          background_selection = baseline_selection
          for ifile, qcd_file in enumerate(qcd_samples):
            f_mc = ROOT.TFile.Open(qcd_file.filename, 'READ')
            
            weight_mc = self.tools.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
            hist_mc_name = 'hist'
            hist_mc = self.tools.createHisto(f_mc, 'signal_tree', quantity, hist_name=hist_mc_name, branchname='flat', selection=background_selection, weight=weight_mc) 
            hist_mc.Sumw2()

            int_mc_tot += hist_mc.Integral()
          
            hist_mc_tot.Add(hist_mc)

          background_yields = self.getBackgroundYields(data_file=self.data_file, signal_file=signal_file, baseline_selection=baseline_selection)
          hist_mc_tot.Scale(background_yields/int_mc_tot)
          root_file.cd()
          hist_mc_tot.Write()

          #root_file.Write()
          root_file.Close()
    
          print '--> {}/{} created'.format(self.outputdir, rootfile_name)
          

    else:
      for category in categories:
        # loop on the signal points
        for signal_file in signal_files:
          self.produceDatacard(signal_file=signal_file, category=category, baseline_selection=baseline_selection)

          #label = self.getLabel(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)
          #if self.do_shape_analysis:
          #  root_file = ROOT.TFile('shape_{}.root'.format(label), 'RECREATE')  
          #  root_file.cd()

            # signal shape

          if self.do_shape_analysis:
            signal_mass, signal_coupling = self.getSignalMassCoupling(signal_file)
            label = self.getLabel(signal_mass=signal_mass, signal_coupling=signal_coupling, category=category)

            rootfile_name = 'shape_{}.root'.format(label)
            root_file = ROOT.TFile.Open('{}/{}'.format(self.outputdir, rootfile_name), 'RECREATE')  
            root_file.cd()

            from quantity import Quantity
            sigma = signal_file.resolution
            quantity = Quantity(name_flat='hnl_mass', nbins=80, bin_min=signal_mass-2*sigma, bin_max=signal_mass+2*sigma)

            # data
            data_hist = ROOT.TH1D('data_obs', 'data_obs', quantity.nbins, quantity.bin_min, quantity.bin_max)
            root_file.cd()
            data_hist.Write()

            # signal shape
            f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+signal_file.filename, 'READ')
            
            filename = signal_file.filename
            original_ctau = filename[filename.find('ctau')+4:filename.find('/', filename.find('ctau')+1)]
            target_ctau = signal_file.ctau
            ctau_weight = -99
            if float(original_ctau) != float(target_ctau):
              ctau_weight = '({ctau0} / {ctau1} * exp((1./{ctau0} - 1./{ctau1}) * gen_hnl_ct))'.format(ctau0=original_ctau, ctau1=target_ctau)

            signal_selection = 'ismatched==1 && {}'.format(baseline_selection)
            category_selection = self.categories[category]
            signal_selection += ' && {}'.format(category_selection)
            hist_sig = self.tools.createHisto(f, 'signal_tree', quantity, hist_name='sig', branchname='flat', selection=signal_selection, weight=ctau_weight)
            
            signal_yields = self.getSignalYields(signal_file=signal_file, baseline_selection=baseline_selection, category=category)
            hist_sig.Scale(signal_yields/hist_sig.Integral())
            root_file.cd() 
            hist_sig.Write()


            # background shape
            int_mc_tot = 0.
            hist_mc_tot = ROOT.TH1D('qcd', 'qcd', quantity.nbins, quantity.bin_min, quantity.bin_max)
            hist_mc_tot.Sumw2()

            background_selection = baseline_selection
            background_selection += ' && {}'.format(category_selection)
            for ifile, qcd_file in enumerate(qcd_samples):
              f_mc = ROOT.TFile.Open(qcd_file.filename, 'READ')
              
              weight_mc = self.tools.computeQCDMCWeight(f_mc, qcd_file.cross_section, qcd_file.filter_efficiency)
              hist_mc_name = 'hist'
              hist_mc = self.tools.createHisto(f_mc, 'signal_tree', quantity, hist_name=hist_mc_name, branchname='flat', selection=background_selection, weight=weight_mc) 
              hist_mc.Sumw2()

              int_mc_tot += hist_mc.Integral()
            
              hist_mc_tot.Add(hist_mc)

            background_yields = self.getBackgroundYields(data_file=self.data_file, signal_file=signal_file, baseline_selection=baseline_selection, category=category)
            #if background_yields == 1e-09 and int_mc_tot==0.:
            #  print 'ok'
            #elif background_yields == 1e-09 and int_mc_tot!=0.:
            #  print 'not ok'

            if background_yields == 1e-09: 
              hist_mc_tot.Fill(signal_file.mass)
            if int_mc_tot != 0.: hist_mc_tot.Scale(background_yields/int_mc_tot)
            if int_mc_tot==0. and background_yields == 1e-09: hist_mc_tot.Scale(1e-09)
            root_file.cd()
            hist_mc_tot.Write()

            #root_file.Write()
            root_file.Close()
      
            print '--> {}/{} created'.format(self.outputdir, rootfile_name)


if __name__ == '__main__':

  data_file = data_samples[0]

  categories = {}
  #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.7'
  #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_pi_mass<4.5 '
  #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_pi_mass<4.5'
  #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.7'
  #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && pi_pt>1 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'
  #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>1.5 && deltar_mu_pi<1.5 && deltar_trgmu_pi<1 && pi_pt>2 && trgmu_mu_mass>0.5 && trgmu_pi_mass<4.5'

  # used for sensitivity
  #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2 '
  #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
  #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
  #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
  #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
  #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'

  categories['lxy0to1_OS_dsa'] = 'mu_isdsa==1 && sv_lxy<1 && trgmu_charge!=mu_charge '
  categories['lxy1to5_OS_dsa'] = 'mu_isdsa==1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge'
  categories['lxygt5_OS_dsa'] = 'mu_isdsa==1 && sv_lxy>=5 && trgmu_charge!=mu_charge'
  categories['lxy0to1_SS_dsa'] = 'mu_isdsa==1 && sv_lxy<1 && trgmu_charge==mu_charge'
  categories['lxy1to5_SS_dsa'] = 'mu_isdsa==1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge'
  categories['lxygt5_SS_dsa'] = 'mu_isdsa==1 && sv_lxy>=5 && trgmu_charge==mu_charge'
  categories['lxy0to1_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.0 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2 '
  categories['lxy1to5_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && b_pt>10 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
  categories['lxygt5_OS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.2 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'
  categories['lxy0to1_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.75 && deltaphi_trgmu_hnl>0.015 && deltaphi_trgmu_hnl<0.45 && deltaphi_trgmu_pi<2'
  categories['lxy1to5_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && b_pt>10 && deltaphi_trgmu_hnl>0.015 && pi_dcasig>8 && sv_lxysig>30 && pi_pt>1'
  categories['lxygt5_SS_nodsa'] = 'mu_isdsa!=1 && sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>0.8 && pi_dcasig>8 && sv_lxysig>40 && pi_pt>1.3'

  #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge && b_mass>2.7'
  #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge && b_mass>1.7 && pi_pt>1'
  #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge && b_mass>1.5 && pi_pt>2'
  #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge && b_mass>2.7'
  #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge && b_mass>1.7 && pi_pt>1'
  #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge && b_mass>1.5 && pi_pt>2'
  #categories['lxy0to1_OS'] = 'sv_lxy<1 && trgmu_charge!=mu_charge'
  #categories['lxy1to5_OS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge!=mu_charge '
  #categories['lxygt5_OS'] = 'sv_lxy>=5 && trgmu_charge!=mu_charge'
  #categories['lxy0to1_SS'] = 'sv_lxy<1 && trgmu_charge==mu_charge'
  #categories['lxy1to5_SS'] = 'sv_lxy>=1 && sv_lxy<5 && trgmu_charge==mu_charge'
  #categories['lxygt5_SS'] = 'sv_lxy>=5 && trgmu_charge==mu_charge'


  from samples import signal_samples_limits_m1, signal_samples_limits_m1_large
  from samples import signal_samples_limits_m3, signal_samples_limits_m3_large
  from samples import signal_samples_limits_m4p5, signal_samples_limits_m4p5_large

  #dirlabel = 'test_categories_selection_29Jun21_fulllumi_muonid'
  dirlabel = 'categories_selection_19Aug21_fulllumi_combineddsa'

  signal_files = signal_samples_limits_m1
  datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, do_categories=True, do_shape_analysis=False, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 

  #signal_files = signal_samples_limits_m3
  #datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, do_categories=True, do_shape_analysis=False, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 

  #signal_files = signal_samples_limits_m4p5
  #datacards = DatacardsMaker(data_file=data_file, signal_files=signal_files, do_categories=True, do_shape_analysis=False, categories=categories, add_Bc=False, outdirlabel=dirlabel).process() 
