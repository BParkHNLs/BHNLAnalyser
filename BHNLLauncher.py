import os
import importlib
from os import path
import sys
sys.path.append('./objects')
from categories import categories
sys.path.append('./cfgs')


#"----------------User's decision board-----------------"

output_label = 'V09_06Nov21'
tag = 'with_weight_pu_qcd_ntrueint'
#cfg_filename = 'example_cfg.py'
cfg_filename = '06Nov21_cfg.py'
submit_batch = False
do_plotter = True
do_datacards = False

#'------------------------------------------------------'



class BHNLLauncher(object):
  def __init__(self, cfg_name, outlabel, tag, submit_batch, do_plotter, do_datacards):
    self.cfg_name = cfg_name
    self.outlabel = outlabel
    self.tag = tag
    self.submit_batch = submit_batch
    self.do_plotter = do_plotter
    self.do_datacards = do_datacards


  def printHeader(self):
    print '            #########################   '
    print '                  BHNL Launcher         '    
    print '            #########################   '

    print '\n'
    print ' -> Will process config file:  ./cfgs/{}.py'.format(self.getConfigName())
    print ' -> Output directory:          ./outputs/{}'.format(self.outlabel)
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
    outputdir = './outputs/{}'.format(self.outlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    command_cp = 'cp ./cfgs/{}.py {}'.format(self.cfg_name, outputdir)
    os.system(command_cp)


  def writeSubmitter(self, command, label):
    #TODO create workdir on scratch?
    # and add time stamp
    content = '\n'.join([
        '#!/bin/bash',
        #'homedir="$PWD/outputs/{}'
        #'workdir="/scratch/{}/{}/job_nj{}_{}"'.format(os.environ["USER"], label, '${SLURM_JOB_ID}' if self.submit_batch else '0', '${SLURM_ARRAY_TASK_ID}' if self.submit_batch else '0'),
        #'echo "creating workdir "$workdir',
        #'mkdir -p $workdir',
        #'echo "copying scripts"',
        #'cp -r ./scripts/*py $workdir',
        #'cp -r ./objects/*py $workdir',
        #'cd $workdir',
        'cd ./scripts',
        command,
        'cd ..',
        'echo "Done"'
        ])
    
    submitter = open('submitter_{}.sh'.format(label), 'w+')
    submitter.write(content)
    submitter.close()


  def launchPlotter(self, category=''):
    # get the command to run the plotting script
    command_plotter = ' '.join([
        'python plotter.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(self.cfg.signal_label),
        '--quantities_label {}'.format(self.cfg.quantities_label),
        '--selection_label {}'.format(self.cfg.selection_label),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--sample_type {}'.format(self.cfg.sample_type),
        '--tree_name {}'.format(self.cfg.tree_name),
        '--qcd_white_list {}'.format(self.cfg.qcd_white_list),
        '--CMStag {}'.format(self.cfg.CMStag),
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
    label = 'plotter_' + self.outlabel + '_' + category.label
    self.writeSubmitter(command_plotter, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      # create logdir
      logdir_name = './outputs/{}/logs'.format(self.outlabel)
      if not path.exists(logdir_name):
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnlplt submitter_{lbl}.sh'.format(
          ld = logdir_name,
          lbl=label,
          ) 

      os.system(command_submit)

      print '--> Plotter job has been submitted'
      print '---> The plots will be stored in ./outputs/{}/plots'.format(self.outlabel)

    command_clean = 'rm submitter_{}.sh'.format(label)
    os.system(command_clean)
    

  def launchDatacards(self, category=''):
    # get the command to run the datacards script
    command_datacards = ' '.join([
        'python create_datacards.py',
        '--outdirlabel {}'.format(self.outlabel),
        '--subdirlabel {}'.format(self.tag),
        '--data_label {}'.format(self.cfg.data_label),
        '--qcd_label {}'.format(self.cfg.qcd_label),
        '--signal_label {}'.format(self.cfg.signal_label),
        '--selection_label {}'.format(self.cfg.selection_label),
        '--categories_label {}'.format(self.cfg.categories_label),
        '--category_label {}'.format(category.label),
        '--ABCD_label {}'.format(self.cfg.ABCD_label),
        '{}'.format('--do_categories' if self.cfg.do_categories else ''),
        '{}'.format('--add_weight_hlt' if self.cfg.add_weight_hlt else ''),
        '{}'.format('--add_Bc' if self.cfg.add_Bc else ''),
        ])

    # write the submitter
    label = 'datacards_' + self.outlabel + '_' + category.label
    self.writeSubmitter(command_datacards, label)

    # launch submitter
    if not self.submit_batch:
      command_submit = 'sh submitter_{}.sh'.format(label)
      os.system(command_submit)
    else:
      # create logdir
      logdir_name = './outputs/{}/logs'.format(self.outlabel)
      if not path.exists(logdir_name):
        print 'creating directory'
        os.system('mkdir -p {}'.format(logdir_name))

      command_submit = 'sbatch -p standard --account t3 -o {ld}/{lbl}.txt -e {ld}/{lbl}.txt --job-name=bhnldcs submitter_{lbl}.sh'.format(
          ld = logdir_name,
          lbl=label,
          ) 

      os.system(command_submit)

      print '--> Datacards job has been submitted'
      print '---> The datacards will be stored in ./outputs/{}/datacards'.format(self.outlabel)

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
        print '\n -> Launching the plotter for the category "{}"'.format(category.title)
        self.launchPlotter(category=category)

    if self.do_datacards: 
      for category in categories:
        print '\n -> Launching the datacards producer for the category "{}"'.format(category.title)
        self.launchDatacards(category=category)

    print '\nDone'
  

if __name__ == '__main__':

  launcher = BHNLLauncher(cfg_name=cfg_filename, outlabel=output_label, tag=tag, submit_batch=submit_batch, do_plotter=do_plotter, do_datacards=do_datacards)
  launcher.process()





