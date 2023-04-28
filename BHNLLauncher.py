import os
from os import path
import importlib
import subprocess
import sys
sys.path.append('./objects')
from categories import categories
from ctau_points import ctau_points
from samples import signal_samples
from coupling_scenarios import coupling_scenarios
sys.path.append('./cfgs')
sys.path.append('./scripts')
sys.path.append('./data')
from tools import Tools


#"----------------User's decision board-----------------"

output_label = 'V13_06Feb23'
tag = 'unblinding_Bc_fullscan_nobernstein_v2'
cfg_filename = 'V13_06Feb23_cfg.py'
submit_batch = True
do_plotter = False
do_datacards = True
do_limits = True
do_interpretation = False

### in do_limits
do_combine_datacards = True
do_produce_limits = True
do_plot_limits = True
do_plot_fits = True

### in do_interpretation
do_interpretation_limits = True
do_interpretation_plotter = True

### other
do_check_config = False
do_standard_queue = True
memory = 20000 #5500
scenarios = ['Majorana']

#'------------------------------------------------------'



class BHNLLauncher(object):
  def __init__(self, cfg_name, outlabel, tag, submit_batch, do_plotter, do_datacards, do_limits, do_combine_datacards, do_produce_limits, do_plot_limits, do_plot_fits, do_interpretation, do_interpretation_limits, do_interpretation_plotter, do_check_config, do_standard_queue, memory, scenarios):
    self.cfg_name = cfg_name
    self.outlabel = outlabel
    self.tag = tag
    self.submit_batch = submit_batch
    self.do_plotter = do_plotter
    self.do_datacards = do_datacards
    self.do_limits = do_limits
    self.do_combine_datacards = do_combine_datacards
    self.do_produce_limits = do_produce_limits
    self.do_plot_limits = do_plot_limits
    self.do_plot_fits = do_plot_fits
    self.do_interpretation = do_interpretation
    self.do_interpretation_limits = do_interpretation_limits
    self.do_interpretation_plotter = do_interpretation_plotter
    self.do_check_config = do_check_config
    self.do_standard_queue = do_standard_queue
    self.memory = memory
    self.scenarios = scenarios
    # get current directory
    self.homedir = '/work/anlyon/'
    #self.homedir = os.getcwd()

    self.logdir_name = '/work/anlyon/outputs/{}/logs/{}'.format(self.outlabel, self.tag)
    if not path.exists(self.logdir_name):
      #print 'creating directory'
      os.system('mkdir -p {}'.format(self.logdir_name))


  def printHeader(self):
    print '            #########################   '
    print '                  BHNL Launcher         '    
    print '            #########################   '

    print '\n'
    print ' -> Will process config file:  ./cfgs/{}.py'.format(self.getConfigName())
    print ' -> Output directory:          ./outputs/{}'.format(self.outlabel)
    print ' -> Tag subdirectory:          ./outputs/{}'.format(self.tag)
    print '\n-.-.-.-.-.-.-.\n' 


  def getConfigName(self):
    if self.cfg_name[len(self.cfg_name)-3:] == '.py': self.cfg_name = self.cfg_name[:len(self.cfg_name)-3]
    return self.cfg_name


  def getConfig(self):
    # check and adapt syntax if necessary
    self.cfg_name = self.getConfigName()

    # check if the config exists
    if not path.exists('./cfgs/{}.py'.format(self.cfg_name)):
      raise RuntimeError('Config file of name "{}" not found in ./cfgs'.format(self.cfg_name))

    # import the config
    cfg_file = importlib.import_module(self.cfg_name)
    cfg = cfg_file.cfg
    return cfg


  def copyConfig(self):
    outputdir = './outputs/{}/cfgs'.format(self.outlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    command_cp = 'cp ./cfgs/{name}.py {out}/{name}_{tag}.py'.format(name=self.cfg_name, out=outputdir, tag=self.tag)
    os.system(command_cp)

  
  def getMassList(self, signal_label=''):
    masses = []
    for signal_sample in signal_samples[signal_label]:
      if signal_sample.mass not in masses: masses.append(signal_sample.mass)
    return masses


  def getCtauList(self, mass):
    from ctau_points import ctau_points
    for ctau_points_list in ctau_points[self.cfg.ctau_points_label]:
      if mass not in ctau_points_list.mass_list: continue
      ctau_points = ctau_points_list.ctau_list
    return ctau_points


  def getJobId(self, job):
    return int(job[job.find('job')+4:])


  def getJobIdsList(self, jobIds):
    listIds = ''
    for jobId in jobIds:
      listIds += '{}:'.format(jobId)
    return listIds[:len(listIds)-1]


  def writeSubmitter(self, command, label, mass=None):
    content = '\n'.join([
        '#!/bin/bash',
        'homedir="$PWD"',
        'workdir="/scratch/{}/{}/job_nj{}"'.format(os.environ["USER"], label, '${SLURM_JOB_ID}' if self.submit_batch else '0'),
        'echo ""',
        'echo " --> creating workdir "$workdir',
        'mkdir -p $workdir',
        'echo " --> copying scripts"',
        'cp -r ./scripts/*py $workdir',
        'cp -r ./limits/*py $workdir',
        'cp -r ./interpretation/*py $workdir',
        'cp -r ./objects/*py $workdir',
        'cp -r ./data $workdir/..',
        '{}'.format('cp -r ./flashgg_plugin $workdir' if self.do_datacards and self.cfg.use_discrete_profiling else ''),
        'cp -r ./scripts/mva/outputs/{} $workdir'.format(self.cfg.training_label),
        '{}'.format('cp -r {}/outputs/{}/datacards/{}/*m_{}* $workdir'.format(self.homedir, self.outlabel, self.tag, str(mass).replace('.', 'p')) if self.do_limits and (self.do_combine_datacards or self.do_produce_limits) else ''),
        '{}'.format('cp -r {}/outputs/{}/datacards_combined/{}/ $workdir'.format(self.homedir, self.outlabel, self.tag) if self.do_limits and self.do_produce_limits else ''),
        '{}'.format('cp -r {}/outputs/{}/limits/{}/results* $workdir'.format(self.homedir, self.outlabel, self.tag) if (self.do_limits and (self.do_plot_limits or self.do_plot_fits)) or (self.do_interpretation and self.do_interpretation_plotter) else ''),
        'cd $workdir',
        'DATE_START=`date +%s`',
        command,
        #'echo " --> content of the wordir"',
        #'ls -l',
        'echo " --> coyping the files to the home directory $homedir"',
        #FIXME we probably want to remove the copying and save directly in the good location
        #'cp -r ./outputs $homedir',
        'DATE_END=`date +%s`',
        'runtime=$((DATE_END-DATE_START))',
        'echo " --> Wallclock running time: $runtime s"',
        'cd $homedir',
        'echo " --> removing workdir"',
        'rm -r $workdir/../data',
        'rm -r $workdir',
        'echo "Done"'
        ])
    
    submitter = open('submitter_{}.sh'.format(label), 'w+')
    submitter.write(content)
    submitter.close()


  def launchPlotter(self, category='', quantity_label=''):
    ## get current directory
    #homedir = os.getcwd()

    # get the command to run the plotting script
    command_plotter = ' '.join([
        'python plotter.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--homedir {}'.format(self.homedir),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(self.cfg.signal_labels[0]), #FIXME
        '--quantities_label {}'.format(quantity_label),
        '--selection_label {}'.format(self.cfg.baseline_selection_label),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--sample_type {}'.format(self.cfg.sample_type),
        '--tree_name {}'.format(self.cfg.tree_name),
        '--qcd_white_list {}'.format(self.cfg.qcd_white_list),
        '--CMStag {}'.format(self.cfg.CMStag),
        '--weight_hlt {}'.format(self.cfg.branch_weight_hlt),
        '--weight_puqcd {}'.format(self.cfg.branch_weight_puqcd),
        '--weight_pusig {}'.format(self.cfg.branch_weight_pusig),
        '--weight_muid {}'.format(self.cfg.branch_weight_muid),
        '--weight_mu0id {}'.format(self.cfg.branch_weight_mu0id),
        '{}'.format('--add_weight_hlt' if self.cfg.add_weight_hlt else ''),
        '{}'.format('--add_weight_pu' if self.cfg.add_weight_pu else ''),
        '{}'.format('--add_weight_muid' if self.cfg.add_weight_muid else ''),
        '{}'.format('--plot_CR' if self.cfg.plot_CR else ''),
        '{}'.format('--plot_SR' if self.cfg.plot_SR else ''),
        '{}'.format('--plot_dataSig' if self.cfg.plot_dataSig else ''),
        '{}'.format('--plot_ratio' if self.cfg.plot_ratio else ''),
        '{}'.format('--do_shape' if self.cfg.do_shape else ''),
        '{}'.format('--do_luminorm' if self.cfg.do_luminorm else ''),
        '{}'.format('--do_stack' if self.cfg.do_stack else ''),
        '{}'.format('--do_log' if self.cfg.do_log else ''),
        '{}'.format('--add_overflow' if self.cfg.add_overflow else ''),
        '{}'.format('--add_CMSlabel' if self.cfg.add_CMSlabel else ''),
        '{}'.format('--do_tdrstyle' if self.cfg.do_tdrstyle else ''),
        ])

    # write the submitter
    label = 'plotter_' + self.outlabel + '_' + self.tag + '_' + category.label + '_' + quantity_label
    self.writeSubmitter(command_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      command_submit = 'sbatch -p {que} --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnlplt_{lbl} submitter_{lbl}.sh'.format(
          que = 'standard' if self.do_standard_queue else 'short',
          ld = self.logdir_name,
          lbl=label,
          ) 

      os.system(command_submit)

      print '--> Plotter job has been submitted'
      print '---> The plots will be stored in ./outputs/{}/plots/{}'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    

  def launchDatacards(self, signal_label='', category=''):
    # get the command to run the datacards script
    command_datacards = ' '.join([
        'python create_datacards.py',
        '--homedir {}'.format(self.homedir),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(signal_label),
        '--ctau_points_label {}'.format(self.cfg.ctau_points_label),
        '--selection_label {}'.format(self.cfg.baseline_selection_label),
        '--training_label {}'.format(self.cfg.training_label),
        '--cut_score {}'.format(self.cfg.cut_score),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--reweighting_strategy {}'.format(self.cfg.reweighting_strategy),
        '--ABCD_label {}'.format(self.cfg.ABCD_label),
        '--signal_model_label {}'.format(self.cfg.signal_model_label),
        '--background_model_label {}'.format(self.cfg.background_model_label),
        '--fit_window_size {}'.format(self.cfg.fit_window_size),
        '--mass_window_size {}'.format(self.cfg.mass_window_size),
        '--nbins {}'.format(self.cfg.nbins),
        '--lumi_target {}'.format(self.cfg.lumi_target),
        '--sigma_B {}'.format(self.cfg.sigma_B),
        '--lhe_efficiency {}'.format(self.cfg.lhe_efficiency),
        '--sigma_mult {}'.format(self.cfg.sigma_mult_window),
        '--resolution_p0 {}'.format(self.cfg.resolution_p0),
        '--resolution_p1 {}'.format(self.cfg.resolution_p1),
        '--weight_hlt {}'.format(self.cfg.branch_weight_hlt),
        '--weight_pusig {}'.format(self.cfg.branch_weight_pusig),
        '--weight_muid {}'.format(self.cfg.branch_weight_muid),
        '--weight_mu0id {}'.format(self.cfg.branch_weight_mu0id),
        '--CMStag {}'.format(self.cfg.CMStag),
        '{}'.format('--do_cutbased' if self.cfg.do_cutbased else ''),
        '{}'.format('--do_mva' if self.cfg.do_mva else ''),
        '{}'.format('--do_parametric' if self.cfg.do_parametric else ''),
        '{}'.format('--add_weight_hlt' if self.cfg.add_weight_hlt else ''),
        '{}'.format('--add_weight_pu' if self.cfg.add_weight_pu else ''),
        '{}'.format('--add_weight_muid' if self.cfg.add_weight_muid else ''),
        '{}'.format('--do_ABCD' if self.cfg.do_ABCD else ''),
        '{}'.format('--do_ABCDHybrid' if self.cfg.do_ABCDHybrid else ''),
        '{}'.format('--do_TF' if self.cfg.do_TF else ''),
        '{}'.format('--do_realData' if self.cfg.do_realData else ''),
        '{}'.format('--do_counting' if self.cfg.do_counting else ''),
        '{}'.format('--do_shape_analysis' if self.cfg.do_shape_analysis else ''),
        '{}'.format('--do_shape_TH1' if self.cfg.do_shape_TH1 else ''),
        '{}'.format('--use_discrete_profiling' if self.cfg.use_discrete_profiling else ''),
        '{}'.format('--do_binned_fit' if self.cfg.do_binned_fit else ''),
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        '{}'.format('--plot_pulls' if self.cfg.plot_pulls else ''),
        '{}'.format('--do_categories' if self.cfg.do_categories else ''),
        '{}'.format('--add_Bc' if self.cfg.add_Bc else ''),
        '{}'.format('--plot_prefit' if self.cfg.plot_prefit else ''),
        '{}'.format('--add_CMSlabel' if self.cfg.add_CMSlabel else ''),
        '{}'.format('--add_lumilabel' if self.cfg.add_lumilabel else ''),
        '{}'.format('--do_tdrstyle' if self.cfg.do_tdrstyle else ''),
        ])

    # write the submitter
    label = 'datacards_' + self.outlabel + '_' + self.tag + '_' + signal_label + '_' + category.label
    self.writeSubmitter(command_datacards, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      command_submit = 'sbatch -p {que} --mem {mem} --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnldcs_{lbl} submitter_{lbl}.sh'.format(
          que = 'standard' if self.do_standard_queue else 'short',
          mem = self.memory,
          ld = self.logdir_name,
          lbl=label,
          ) 

      job = subprocess.check_output(command_submit, shell=True)
      job_id = self.getJobId(job)

      print '--> Datacards job has been submitted'
      print '---> The datacards will be stored in ./outputs/{}/datacards/{}'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    
    return job_id


  def launchCombineDatacards(self, mass='', scenario='', do_dependency=False, job_id=''):
    # get the command to run the datacards combination script
    command_datacards_combine = ' '.join([
        'python combine_datacards.py',
        '--homedir {}'.format(self.homedir),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--mass {}'.format(mass),
        '--scenario {}'.format(scenario),
        #'--wildcard {}'.format(self.cfg.datacards_wildcard), #TODO add that in config
        #'--mass_whitelist {}'.format(mass), #TODO adapt
        #'--mass_blacklist {}'.format(self.getParserString(self.cfg.mass_black_list)), #TODO adapt
        #'--coupling_whitelist {}'.format(self.getParserString(self.cfg.coupling_white_list)), #TODO adapt
        #'--coupling_blacklist {}'.format(self.getParserString(self.cfg.coupling_black_list)), #TODO adapt
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        ])

    # write the submitter
    label = 'combine_datacards_' + self.outlabel + '_' + self.tag + '_' + str(mass).replace('.', 'p')
    self.writeSubmitter(command_datacards_combine, label, mass)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      command_submit = 'sbatch -p short --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnldcscmb_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = self.logdir_name,
          lbl = label,
          dpd = '--dependency=afterany:{}'.format(job_id) if do_dependency else ''
          ) 

      #os.system(command_submit)
      job = subprocess.check_output(command_submit, shell=True)
      job_id = self.getJobId(job)

      print '--> Datacards Combination job has been submitted'
      print '---> The combined datacards will be stored in ./outputs/{}/datacards_combined/{}'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)

    return job_id


  def launchLimitsProducer(self, mass='', ctau='', scenario='', do_dependency=False, job_id=''):
    # get the command to run the limits producer script
    command_limits_producer = ' '.join([
        'python produce_limits.py',
        '--homedir {}'.format(self.homedir),
        '--indirlabel {}'.format(self.tag),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--mass {}'.format(mass),
        '--ctau {}'.format(ctau), 
        '--scenario {}'.format(scenario),
        '{}'.format('--use_discrete_profiling' if self.cfg.use_discrete_profiling else ''),
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        ])

    # write the submitter
    label = 'limits_producer_' + self.outlabel + '_' + self.tag + '_' + str(mass).replace('.', 'p') + '_' + str(ctau).replace('.', 'p')
    self.writeSubmitter(command_limits_producer, label, mass)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      command_submit = 'sbatch -p {que} --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnllimits_{lbl} {dpd} submitter_{lbl}.sh'.format(
          #que = 'short',
          que = 'standard' if self.do_standard_queue else 'short',
          ld = self.logdir_name,
          lbl = label,
          dpd = '--dependency=afterany:{}'.format(job_id) if do_dependency else ''
          ) 

      job = subprocess.check_output(command_submit, shell=True)
      job_id = self.getJobId(job)

      print '--> Limits production job has been submitted'
      print '---> The limits results will be stored in ./outputs/{}/limits/{}/results'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    
    return job_id


  def launchLimitsPlotter(self, scenario='', do_dependency=False, job_id='', fe=None, fu=None, ft=None):
    # get the command to run the limits plotter script
    command_limits_plotter = ' '.join([
        #'python pvalue.py',
        'python limit_plotter.py',
        '--scenario {}'.format(scenario),
        '--homedir {}'.format(self.homedir),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        #'--mass_whitelist {}'.format(self.getParserString(self.cfg.mass_white_list)), #FIXME 
        #'--mass_whitelist 2.0',
        #'--mass_blacklist 2.55,3.15,3.3,3.35,3.85,3.95,5.0',
        #'--coupling_whitelist {}'.format(self.getParserString(self.cfg.coupling_white_list)), #FIXME 
        #'--coupling_blacklist {}'.format(self.getParserString(self.cfg.coupling_black_list)), #FIXME  
        #'--coupling_blacklist 0.00054,0.00017',
        '--fe {}'.format(fe),
        '--fu {}'.format(fu),
        '--ft {}'.format(ft),
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        ])

    # write the submitter
    label = 'limits_plotter_' + self.outlabel + '_' + self.tag
    self.writeSubmitter(command_limits_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      command_submit = 'sbatch -p short --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=limitplotter_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = self.logdir_name,
          lbl = label,
          dpd = '--dependency=afterany:{}'.format(job_id) if do_dependency else ''
          ) 

      os.system(command_submit)

      print '--> Limits plotter job has been submitted'
      print '---> The limits plots will be stored in ./outputs/{}/limits/{}/plots'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)


  def launchFitsPlotter(self, mass='', scenario='', do_dependency=False, job_id='', fe=None, fu=None, ft=None):
    # get the command to run the limits plotter script
    command_fits_plotter = ' '.join([
        'python fit_plotter.py',
        '--scenario {}'.format(scenario),
        '--mass {}'.format(mass),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--homedir {}'.format(self.homedir),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--fe {}'.format(fe),
        '--fu {}'.format(fu),
        '--ft {}'.format(ft),
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        ])

    # write the submitter
    label = 'fits_plotter_' + self.outlabel + '_' + self.tag
    self.writeSubmitter(command_fits_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      command_submit = 'sbatch -p short --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=fitplotter_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = self.logdir_name,
          lbl = label,
          dpd = '--dependency=afterany:{}'.format(job_id) if do_dependency else ''
          ) 

      os.system(command_submit)

      print '--> Fits plotter job has been submitted'
      print '---> The fits plots will be stored in ./outputs/{}/limits/{}/plots'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)


  def launchInterpretationLimits(self, mass='', ctau='', fe='', fu='', ft=''):
    # get the command to run the limits producer script
    command_interpretation_limits_producer = ' '.join([
        'python interpretation_launcher.py',
        '--homedir {}'.format(self.homedir),
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--mass {}'.format(mass),
        '--ctau {}'.format(ctau), 
        '--fe {}'.format(fe),
        '--fu {}'.format(fu),
        '--ft {}'.format(ft),
        '--muon_label {}'.format(self.cfg.muon_eoslabel),
        '--electron_label {}'.format(self.cfg.electron_eoslabel),
        '{}'.format('--use_discrete_profiling' if self.cfg.use_discrete_profiling else ''),
        '{}'.format('--do_blind' if self.cfg.do_blind else ''),
        ])

    # write the submitter
    label = 'interpretation_limits_producer_' + self.outlabel + '_' + self.tag + '_' + str(mass).replace('.', 'p') + '_' + str(ctau).replace('.', 'p')
    self.writeSubmitter(command_interpretation_limits_producer, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      command_submit = 'sbatch -p {que} --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnllimits_{lbl} submitter_{lbl}.sh'.format(
          #que = 'standard' if self.do_standard_queue else 'short',
          que = 'short' if self.do_standard_queue else 'short',
          ld = self.logdir_name,
          lbl = label,
          ) 

      job = subprocess.check_output(command_submit, shell=True)
      job_id = self.getJobId(job)

      print '--> Interpretation limits production job has been submitted'
      print '---> The interpretation limits results will be stored in ./outputs/{}/limits/{}/results'.format(self.outlabel, self.tag) #add coupling

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    
    return job_id


  def process(self):
    self.printHeader()

    print ' -> Getting the config'
    self.cfg = self.getConfig()

    if self.do_check_config:
      print '\n -> Checking the config'
      self.cfg.checkConfig()

    print '\n -> Copying the config to the output directory'
    self.copyConfig()

    sys.path.append('./objects')
    from categories import categories
    categories = categories[self.cfg.categories_label]

    if self.do_plotter: 
      for category in categories:
        for quantity_label in self.cfg.quantities_label:
          print '\n -> Launching the plotter for the category "{}" and quantities "{}"'.format(category.title, quantity_label)
          self.launchPlotter(category=category, quantity_label=quantity_label)

    if self.do_datacards: 
      dependency_datacards = {}
      for signal_label in self.cfg.signal_labels:
        dependency_datacards[signal_label] = []
        for category in categories:
          print '\n -> Launching the datacards producer for the signal samples "{}" in the category "{}"'.format(signal_label, category.title)
          job_id = self.launchDatacards(signal_label=signal_label, category=category)
          dependency_datacards[signal_label].append(job_id)

    if self.do_limits:
      if self.do_combine_datacards:
        dependency_combined_datacards = {}
        for signal_label in self.cfg.signal_labels:
          dependency_combined_datacards[signal_label] = []
          print '\n -> Launching the combination of the datacards'
          for mass in self.getMassList(signal_label):
            for scenario in self.scenarios:
              if self.submit_batch and self.do_datacards:
                job_id = self.launchCombineDatacards(mass=str(mass).replace('.', 'p'), scenario=scenario, do_dependency=True, job_id=self.getJobIdsList(dependency_datacards[signal_label]))
              else:
                job_id = self.launchCombineDatacards(mass=str(mass).replace('.', 'p'), scenario=scenario, do_dependency=False)
              dependency_combined_datacards[signal_label].append(job_id)

      if self.do_produce_limits:
        dependency_limits = []
        for signal_label in self.cfg.signal_labels:
          for mass in self.getMassList(signal_label):
            for scenario in self.scenarios:
              for ctau in self.getCtauList(mass=mass):
                print '\n -> Launching the limits production for mass {} and ctau {} ({})'.format(mass, ctau, scenario)
                if self.submit_batch and self.do_combine_datacards:
                  job_id = self.launchLimitsProducer(mass=mass, ctau=ctau, scenario=scenario, do_dependency=True, job_id=self.getJobIdsList(dependency_combined_datacards[signal_label]))
                else:
                  job_id = self.launchLimitsProducer(mass=mass, ctau=ctau, scenario=scenario, do_dependency=False)
                dependency_limits.append(job_id)

      if self.do_plot_limits:
        for scenario in self.scenarios:
          print '\n -> Launching the limits plotter ({})'.format(scenario)
          if self.submit_batch and self.do_produce_limits:
            job_id = self.launchLimitsPlotter(scenario=scenario, do_dependency=True, job_id=self.getJobIdsList(dependency_limits))
          else:
            job_id = self.launchLimitsPlotter(scenario=scenario, do_dependency=False)

      if self.do_plot_fits:
        for signal_label in self.cfg.signal_labels:
          for mass in self.getMassList(signal_label):
            for scenario in self.scenarios:
              print '\n -> Launching the fits plotter ({})'.format(scenario)
              if self.submit_batch and self.do_produce_limits:
                job_id = self.launchFitsPlotter(mass=mass, scenario=scenario, do_dependency=True, job_id=self.getJobIdsList(dependency_limits))
              else:
                job_id = self.launchFitsPlotter(mass=mass, scenario=scenario, do_dependency=False)

    if self.do_interpretation:
      if self.do_interpretation_limits:
        dependency_interpretation_limits = []
        for coupling_scenario in coupling_scenarios[self.cfg.coupling_scenarios_label]:
          fe = coupling_scenario.fe
          fu = coupling_scenario.fu
          ft = coupling_scenario.ft
          for signal_label in self.cfg.signal_labels:
            for mass in self.getMassList(signal_label):
              for ctau in self.getCtauList(mass=mass):
                # add loop on coupling scenario
                signal_v2 = Tools().getVV(mass=mass, ctau=ctau, ismaj=True)
                signal_coupling = Tools().getCouplingLabel(signal_v2)
                print '\n -> Launching the interpretation limits production for (fe={}, fu={}, ft={}) and mass {} and ctau {}'.format(fe, fu, ft, mass, ctau)
                job_id = self.launchInterpretationLimits(mass=mass, ctau=ctau, fe=fe, fu=fu, ft=ft)
                dependency_interpretation_limits.append(job_id)
          os.system('sleep 30s')

      if self.do_interpretation_plotter:
        for coupling_scenario in coupling_scenarios[self.cfg.coupling_scenarios_label]:
          fe = coupling_scenario.fe
          fu = coupling_scenario.fu
          ft = coupling_scenario.ft
          print '\n -> Launching the limits plotter for (fe={}, fu={}, ft={})'.format(fe, fu, ft)
          if self.submit_batch and self.do_interpretation_limits:
            job_id = self.launchLimitsPlotter(scenario='Majorana', do_dependency=True, job_id=self.getJobIdsList(dependency_interpretation_limits), fe=fe, fu=fu, ft=ft)
          else:
            job_id = self.launchLimitsPlotter(scenario='Majorana', do_dependency=False, fe=fe, fu=fu, ft=ft)

    print '\nDone'
  

if __name__ == '__main__':

  launcher = BHNLLauncher(cfg_name=cfg_filename, 
                          outlabel=output_label, 
                          tag=tag, 
                          submit_batch=submit_batch, 
                          do_plotter=do_plotter, 
                          do_datacards=do_datacards, 
                          do_limits=do_limits, 
                          do_combine_datacards=do_combine_datacards, 
                          do_produce_limits=do_produce_limits, 
                          do_plot_limits=do_plot_limits,
                          do_plot_fits=do_plot_fits,
                          do_interpretation=do_interpretation,
                          do_interpretation_limits=do_interpretation_limits,
                          do_interpretation_plotter=do_interpretation_plotter,
                          do_check_config=do_check_config,
                          do_standard_queue=do_standard_queue,
                          memory=memory,
                          scenarios=scenarios,
                          )
  launcher.process()





