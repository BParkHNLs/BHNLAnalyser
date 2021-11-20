import os
from os import path
import importlib
import subprocess
import sys
sys.path.append('./objects')
from categories import categories
from samples import signal_samples
sys.path.append('./cfgs')


#"----------------User's decision board-----------------"

output_label = 'V09_06Nov21'
#tag = 'ABCDHybrid_bmasshnlcharge_3cat0120significance_fullA'
tag = 'test'
#cfg_filename = 'example_cfg.py'
cfg_filename = '06Nov21_cfg.py'
submit_batch = True
do_plotter = False
do_datacards = True
do_limits = True

###
do_combine_datacards = True
do_produce_limits = True
do_plot_limits = True


#'------------------------------------------------------'



class BHNLLauncher(object):
  def __init__(self, cfg_name, outlabel, tag, submit_batch, do_plotter, do_datacards, do_limits, do_combine_datacards, do_produce_limits, do_plot_limits):
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


  def getCtauList(self, signal_label=''):
    ctaus = []
    for signal_sample in signal_samples[signal_label]:
      if signal_sample.ctau not in ctaus: ctaus.append(signal_sample.ctau)
    return ctaus


  def getJobId(self, job):
    return int(job[job.find('job')+4:])


  def getJobIdsList(self, jobIds):
    listIds = ''
    for jobId in jobIds:
      listIds += '{}:'.format(jobId)
    return listIds[:len(listIds)-1]


  def writeSubmitter(self, command, label):
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
        'cp -r ./objects/*py $workdir',
        'cp -r ./data $workdir',
        '{}'.format('cp -r ./outputs/{}/datacards/{}/ $workdir'.format(self.outlabel, self.tag) if self.do_limits and self.do_combine_datacards else ''),
        '{}'.format('cp -r ./outputs/{}/datacards_combined/{}/ $workdir'.format(self.outlabel, self.tag) if self.do_limits and self.do_produce_limits else ''),
        '{}'.format('cp -r ./outputs/{}/limits/{}/results $workdir'.format(self.outlabel, self.tag) if self.do_limits and self.do_plot_limits else ''),
        'cd $workdir',
        'DATE_START=`date +%s`',
        command,
        'echo " --> content of the wordir"',
        'ls -l',
        'echo " --> coyping the files to the home directory $homedir"',
        'cp -r ./outputs $homedir',
        'DATE_END=`date +%s`',
        'runtime=$((DATE_END-DATE_START))',
        'echo " --> Wallclock running time: $runtime s"',
        'cd $homedir',
        'echo " --> removing workdir"',
        #'rm -r $workdir',
        'echo "Done"'
        ])
    
    submitter = open('submitter_{}.sh'.format(label), 'w+')
    submitter.write(content)
    submitter.close()


  def launchPlotter(self, category='', quantity_label=''):
    # get the command to run the plotting script
    command_plotter = ' '.join([
        'python plotter.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(self.cfg.signal_labels[0]), #FIXME
        '--quantities_label {}'.format(quantity_label),
        '--selection_label {}'.format(self.cfg.selection_label),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--sample_type {}'.format(self.cfg.sample_type),
        '--tree_name {}'.format(self.cfg.tree_name),
        '--qcd_white_list {}'.format(self.cfg.qcd_white_list),
        '--CMStag {}'.format(self.cfg.CMStag),
        '--weight_hlt {}'.format(self.cfg.branch_weight_hlt),
        '--weight_pu {}'.format(self.cfg.branch_weight_pu),
        '{}'.format('--add_weight_hlt' if self.cfg.add_weight_hlt else ''),
        '{}'.format('--add_weight_pu' if self.cfg.add_weight_pu else ''),
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
        ])

    # write the submitter
    label = 'plotter_' + self.outlabel + '_' + self.tag + '_' + category.label + '_' + quantity_label
    self.writeSubmitter(command_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      # create logdir
      logdir_name = './outputs/{}/logs/{}'.format(self.outlabel, self.tag)
      if not path.exists(logdir_name):
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnlplt_{lbl} submitter_{lbl}.sh'.format(
          ld = logdir_name,
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
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(signal_label),
        '--selection_label {}'.format(self.cfg.selection_label),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--ABCD_label {}'.format(self.cfg.ABCD_label),
        '--lumi_target {}'.format(self.cfg.lumi_target),
        '--sigma_B {}'.format(self.cfg.sigma_B),
        '--weight_hlt {}'.format(self.cfg.branch_weight_hlt),
        #'--weight_pu {}'.format(self.cfg.branch_weight_pu),
        '{}'.format('--add_weight_hlt' if self.cfg.add_weight_hlt else ''),
        '{}'.format('--do_ABCD' if self.cfg.do_ABCD else ''),
        '{}'.format('--do_ABCDHybrid' if self.cfg.do_ABCDHybrid else ''),
        '{}'.format('--do_TF' if self.cfg.do_TF else ''),
        '{}'.format('--do_categories' if self.cfg.do_categories else ''),
        '{}'.format('--add_Bc' if self.cfg.add_Bc else ''),
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
      # create logdir
      logdir_name = './outputs/{}/logs/{}'.format(self.outlabel, self.tag)
      if not path.exists(logdir_name):
        print 'creating directory'
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnldcs_{lbl} submitter_{lbl}.sh'.format(
          ld = logdir_name,
          lbl=label,
          ) 

      job = subprocess.check_output(command_submit, shell=True)
      job_id = self.getJobId(job)

      print '--> Datacards job has been submitted'
      print '---> The datacards will be stored in ./outputs/{}/datacards/{}'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    
    return job_id


  def launchCombineDatacards(self, mass='', do_dependency=False, job_id=''):
    # get the command to run the datacards combination script
    command_datacards_combine = ' '.join([
        'python combine_datacards.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--categories_label {}'.format(self.cfg.categories_label),
        #'--wildcard {}'.format(self.cfg.datacards_wildcard), #TODO add that in config
        '--mass_whitelist {}'.format(mass), #TODO adapt
        #'--mass_blacklist {}'.format(self.getParserString(self.cfg.mass_black_list)), #TODO adapt
        #'--coupling_whitelist {}'.format(self.getParserString(self.cfg.coupling_white_list)), #TODO adapt
        #'--coupling_blacklist {}'.format(self.getParserString(self.cfg.coupling_black_list)), #TODO adapt
        '{}'.format('--run_blind' if self.cfg.run_blind else ''),
        ])

    # write the submitter
    label = 'combine_datacards_' + self.outlabel + '_' + self.tag + '_' + str(mass).replace('.', 'p')
    self.writeSubmitter(command_datacards_combine, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      # create logdir
      logdir_name = './outputs/{}/logs/{}'.format(self.outlabel, self.tag)
      if not path.exists(logdir_name):
        print 'creating directory'
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnldcscmb_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = logdir_name,
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


  def launchLimitsProducer(self, mass='', ctau='', do_dependency=False, job_id=''):
    # get the command to run the limits producer script
    command_limits_producer = ' '.join([
        'python produce_limits.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--mass {}'.format(mass),
        '--ctau {}'.format(ctau), 
        '{}'.format('--run_blind' if self.cfg.run_blind else ''),
        ])

    # write the submitter
    label = 'limits_producer_' + self.outlabel + '_' + self.tag + '_' + str(mass).replace('.', 'p') + '_' + str(ctau).replace('.', 'p')
    self.writeSubmitter(command_limits_producer, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
      job_id = 0
    else:
      # create logdir
      logdir_name = './outputs/{}/logs/{}'.format(self.outlabel, self.tag)
      if not path.exists(logdir_name):
        print 'creating directory'
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnllimits_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = logdir_name,
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


  def launchLimitsPlotter(self, do_dependency=False, job_id=''):
    # get the command to run the limits plotter script
    command_limits_plotter = ' '.join([
        'python limit_plotter.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        #'--mass_whitelist {}'.format(self.getParserString(self.cfg.mass_white_list)), #FIXME 
        #'--mass_blacklist {}'.format(self.getParserString(self.cfg.mass_black_list)), #FIXME  
        #'--coupling_whitelist {}'.format(self.getParserString(self.cfg.coupling_white_list)), #FIXME 
        #'--coupling_blacklist {}'.format(self.getParserString(self.cfg.coupling_black_list)), #FIXME  
        '{}'.format('--run_blind' if self.cfg.run_blind else ''),
        ])

    # write the submitter
    label = 'limits_plotter_' + self.outlabel + '_' + self.tag
    self.writeSubmitter(command_limits_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      # create logdir
      logdir_name = './outputs/{}/logs/{}'.format(self.outlabel, self.tag)
      if not path.exists(logdir_name):
        print 'creating directory'
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=limitplotter_{lbl} {dpd} submitter_{lbl}.sh'.format(
          ld = logdir_name,
          lbl = label,
          dpd = '--dependency=afterany:{}'.format(job_id) if do_dependency else ''
          ) 

      os.system(command_submit)

      print '--> Limits plotter job has been submitted'
      print '---> The limits plots will be stored in ./outputs/{}/limits/{}/plots'.format(self.outlabel, self.tag)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)


  def process(self):
    self.printHeader()

    print ' -> Getting the config'
    self.cfg = self.getConfig()

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
            if self.submit_batch and self.do_datacards:
              job_id = self.launchCombineDatacards(mass=mass, do_dependency=True, job_id=self.getJobIdsList(dependency_datacards[signal_label]))
            else:
              job_id = self.launchCombineDatacards(mass=mass, do_dependency=False)
          dependency_combined_datacards[signal_label].append(job_id)

      if self.do_produce_limits:
        for signal_label in self.cfg.signal_labels:
          dependency_limits = []
          for mass in self.getMassList(signal_label):
            for ctau in self.getCtauList(signal_label):
              print '\n -> Launching the limits production for mass {} and ctau {}'.format(mass, ctau)
              if self.submit_batch and self.do_combine_datacards:
                job_id = self.launchLimitsProducer(mass=mass, ctau=ctau, do_dependency=True, job_id=self.getJobIdsList(dependency_combined_datacards[signal_label]))
              else:
                job_id = self.launchLimitsProducer(mass=mass, ctau=ctau, do_dependency=False)
              dependency_limits.append(job_id)

      if self.do_plot_limits:
        print '\n -> Launching the limits plotter'
        if self.submit_batch and self.do_produce_limits:
          job_id = self.launchLimitsPlotter(do_dependency=True, job_id=self.getJobIdsList(dependency_limits))
        else:
          job_id = self.launchLimitsPlotter(do_dependency=False)

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
                          do_plot_limits=do_plot_limits
                          )
  launcher.process()





