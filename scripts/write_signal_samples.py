import os
from os import path
from glob import glob
import ROOT


class SignalSampleWriter(object):
  def __init__(self, nano_pl, nano_pl_Bc, ntuple_name, ntuple_name_Bc, nano_general_filename, gen_pl, gen_pl_Bc):
    self.nano_pl = nano_pl
    self.nano_pl_Bc = nano_pl_Bc
    self.ntuple_name = ntuple_name
    self.ntuple_name_Bc = ntuple_name_Bc
    self.nano_general_filename = nano_general_filename
    self.gen_pl = gen_pl
    self.gen_pl_Bc = gen_pl_Bc

    self.sample_collection_templatename = 'V42_06Feb23_m{mass}_norm'
    self.force_overwrite = False

    self.path_SE = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen'
    self.ntuple_templatepath = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/{pl}/mass{mass}_ctau{ctau}/nanoFiles/merged'
    self.nano_templatepath = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass{mass}_ctau{ctau}/nanoFiles/Chunk*'
    self.gen_path = '/work/anlyon/request_common/CMSSW_10_2_28_patch1/src/McRequest/slurm/'


  def getMass(self, point):
    idx1 = point.rfind('mass')+4
    idx2 = point.rfind('_', idx1+1)
    mass = point[idx1:idx2]

    return mass


  def getCtau(self, point):
    idx1 = point.rfind('ctau')+4
    idx2 = len(point)
    ctau = point[idx1:idx2]

    return ctau


  def getPoint(self, directory):
    idx1 = directory.rfind('/')+1
    idx2 = len(directory)
    point = directory[idx1:idx2]

    return point


  def getMiniaAODEventsPerBSpecies(self, mass, ctau):
    point = 'mass{}_ctau{}'.format(mass, ctau)

    fetch_n_miniaod_fromfile = False

    for line in self.lines_miniaod:
      if point in line:
        fetch_n_miniaod_fromfile = True
        idx1 = line.find(' ', 1)+1
        idx2 = line.find(' ', idx1+1)
        idx3 = line.find(' ', idx2+1)
        idx4 = len(line)
        count_Bu_muon = line[idx1:idx2]
        count_Bd_muon = line[idx2+1:idx3]
        count_Bs_muon = line[idx3+1:idx4]
    
    if not fetch_n_miniaod_fromfile:
      inputfiles_name = self.nano_templatepath.format(mass=mass, ctau=ctau) + '/' + self.nano_general_filename
      inputfiles = [f for f in glob(inputfiles_name)]

      tree = ROOT.TChain('Events')
      #for i, inputfile in enumerate(inputfiles):
      for inputfile in inputfiles:
        #if i > 3: continue
        tree.Add(inputfile)

      n_events = tree.GetEntries()
      print '   ---> number of events: {}'.format(n_events)

      count_muon = 0
      count_electron = 0
      count_Bu_muon = 0
      count_Bd_muon = 0
      count_Bs_muon = 0
      count_Bu_electron = 0
      count_Bd_electron = 0
      count_Bs_electron = 0
      count_tot = 0

      for entry in tree:
        if count_tot%10000 == 0: print '   ---> processed {}% of the files'.format(round(float(count_tot)/n_events*100, 1))
        count_tot += 1

        # searching for hnl
        has_hnl = False
        hnl_idx = -1
        b_idx = -1
        for icand in range(0, entry.nGenPart):
          if abs(entry.GenPart_pdgId[icand]) == 9900015:
            has_hnl = True
            hnl_idx = icand
            b_idx = entry.GenPart_genPartIdxMother[icand]
            break
        if not has_hnl: print 'WARNING: no hnl found'

        if has_hnl:
          # searching for hnl daughter
          daughter_ismuon = False
          daughter_iselectron = False
          for icand in range(0, entry.nGenPart):
            if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == hnl_idx:
              if abs(entry.GenPart_pdgId[icand]) == 13:
                daughter_ismuon = True
              elif abs(entry.GenPart_pdgId[icand]) == 11: 
                daughter_iselectron = True

          # searching for hnl sister
          sister_ismuon = False
          sister_iselectron = False
          for icand in range(0, entry.nGenPart):
            if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == b_idx:
              if abs(entry.GenPart_pdgId[icand]) == 13:
                sister_ismuon = True
              elif abs(entry.GenPart_pdgId[icand]) == 11: 
                sister_iselectron = True

          if daughter_ismuon and sister_ismuon: 
            count_muon = count_muon + 1
            if abs(entry.GenPart_pdgId[b_idx]) == 521: count_Bu_muon += 1
            elif abs(entry.GenPart_pdgId[b_idx]) == 511: count_Bd_muon += 1
            elif abs(entry.GenPart_pdgId[b_idx]) == 531: count_Bs_muon += 1
          elif (daughter_ismuon and sister_iselectron) or (daughter_iselectron and sister_ismuon):
            count_electron = count_electron + 1
            if abs(entry.GenPart_pdgId[b_idx]) == 521: count_Bu_electron += 1
            elif abs(entry.GenPart_pdgId[b_idx]) == 511: count_Bd_electron += 1
            elif abs(entry.GenPart_pdgId[b_idx]) == 531: count_Bs_electron += 1
      
      #print 'filename: {}'.format(inputfile)
      #print 'muon rate: {}'.format(round(float(count_muon) / float(count_tot), 2))
      #print 'electron rate: {}'.format(round(float(count_electron) / float(count_tot), 2))
      #print 'number of Bu (muon): {}'.format(count_Bu_muon)
      #print 'number of Bd (muon): {}'.format(count_Bd_muon)
      #print 'number of Bs (muon): {}'.format(count_Bs_muon)
      #print 'number of Bu (electron): {}'.format(count_Bu_electron)
      #print 'number of Bd (electron): {}'.format(count_Bd_electron)
      #print 'number of Bs (electron): {}'.format(count_Bs_electron)
      #if count_Bu_muon + count_Bd_muon + count_Bs_muon + count_Bu_electron + count_Bd_electron + count_Bs_electron != count_tot: print '--> WARNING: numbers do not close'

      self.miniaod_file.write('{} {} {} {}\n'.format(point, float(count_Bu_muon), float(count_Bd_muon), float(count_Bs_muon)))

    return float(count_Bu_muon), float(count_Bd_muon), float(count_Bs_muon)


  def getFilterEfficiency(self, logs):
    list_n_gen = []
    list_n_tot = []
    for log in logs:
      with open(log, 'r') as f:
        for line in f:
          if ('Filter efficiency' in line) and ('TO BE USED IN MC' in line):
            n_gen = float(line.split('= (')[1].split(') / (')[0])
            n_tot = float(line.split('= (')[1].split(') / (')[1].split(') = ')[0])
            list_n_gen.append(n_gen)
            list_n_tot.append(n_tot)

    filter_efficiency = sum(list_n_gen) / sum(list_n_tot) if sum(list_n_tot) != 0 else -99

    return filter_efficiency


  def getFilterEfficienciesPerBSpecies(self, mass, ctau, is_bc):
    if not is_bc:
      path_log_Bu = self.gen_path + '/' + self.gen_pl + '_Bu_muon/logs'
      logs_Bu = [f for f in glob('{}/prod_mass{}_ctau{}_*.log'.format(path_log_Bu, mass.replace('p', '.'), ctau.replace('p', '.')))]
      filter_efficiency_Bu = self.getFilterEfficiency(logs_Bu)

      path_log_Bd = self.gen_path + '/' + self.gen_pl + '_Bd_muon/logs'
      logs_Bd = [f for f in glob('{}/prod_mass{}_ctau{}_*.log'.format(path_log_Bd, mass.replace('p', '.'), ctau.replace('p', '.')))]
      filter_efficiency_Bd = self.getFilterEfficiency(logs_Bd)

      path_log_Bs = self.gen_path + '/' + self.gen_pl + '_Bs_muon/logs'
      #path_log_Bs = self.gen_path + '/test_Bs_muon/logs'
      logs_Bs = [f for f in glob('{}/prod_mass{}_ctau{}_*.log'.format(path_log_Bs, mass.replace('p', '.'), ctau.replace('p', '.')))]
      filter_efficiency_Bs = self.getFilterEfficiency(logs_Bs)

      return filter_efficiency_Bu, filter_efficiency_Bd, filter_efficiency_Bs

    else:
      path_log_Bc = self.gen_path + '/' + self.gen_pl_Bc + '_Bc_muon/logs'
      logs_Bc = [f for f in glob('{}/prod_mass{}_ctau{}_*.log'.format(path_log_Bc, mass.replace('p', '.'), ctau.replace('p', '.')))]
      filter_efficiency_Bc = self.getFilterEfficiency(logs_Bc)

      return filter_efficiency_Bc


  def process(self):
    # get signal points
    point_directories = [f for f in glob('{}/{}/*'.format(self.path_SE, self.nano_pl))]  
    point_directories_Bc = [f for f in glob('{}/{}/*'.format(self.path_SE, self.nano_pl_Bc))]  

    points = [self.getPoint(f) for f in point_directories]
    points_Bc = [self.getPoint(f) for f in point_directories_Bc]

    # get mass points
    masses = []
    points_full = [] 
    for point in points:
      if 'None' in point or 'mass0p' in point: 
        continue
      mass = self.getMass(point)
      if mass not in masses: masses.append(mass)
      if point not in points_full: points_full.append(point)

    for point in points_Bc:
      if 'None' in point or 'mass0p' in point: 
        continue
      mass = self.getMass(point)
      if mass not in masses: masses.append(mass)
      if point not in points_full: points_full.append(point)

    masses.sort()

    # create samples outfile
    outfile = open('signal_samples.txt', 'w+')

    # create miniaod file
    miniaod_filename = 'miniaod.txt'
    if not path.exists(miniaod_filename) or self.force_overwrite:
      self.miniaod_file = open(miniaod_filename, 'w+')
    else:
      self.miniaod_file = open(miniaod_filename, 'r')
    self.lines_miniaod = self.miniaod_file.readlines()

    # create signal samples template
    signal_samples = []
    mass_count = 1
    for mass in masses:
      #if mass != '1p0' and mass != '3p0' and mass != '6p0': continue

      print '\n-> Processing mass {} ({} / {})'.format(mass, mass_count, len(masses))
      mass_count += 1
      sample_name = self.sample_collection_templatename.format(mass=mass)
      template = []
      template.append("signal_samples['{}'] = [".format(sample_name))

      for point in points_full:
        if 'mass{}_ctau'.format(mass) not in point: continue

        is_point_maingrid = False
        is_point_Bcgrid = False

        print '\n --> getting sample information' 
        ctau = self.getCtau(point)

        directory = '{}/{}/mass{}_ctau{}'.format(self.path_SE, self.nano_pl, mass, ctau)
        directory_Bc = '{}/{}/mass{}_ctau{}'.format(self.path_SE, self.nano_pl_Bc, mass, ctau)

        if directory in point_directories:
          is_point_maingrid = True

        if directory_Bc in point_directories_Bc:
          is_point_Bcgrid = True

        if is_point_maingrid: sample_filename = self.ntuple_templatepath.format(pl=self.nano_pl, mass=mass, ctau=ctau) + '/' + self.ntuple_name
        if is_point_Bcgrid: sample_filename_Bc = self.ntuple_templatepath.format(pl=self.nano_pl_Bc, mass=mass, ctau=ctau) + '/' + self.ntuple_name_Bc

        sample_label = '{mass} GeV, {ctau} mm'.format(mass=mass.replace('p', '.'), ctau=ctau.replace('p', '.'))
        sample_mass = mass.replace('p', '.')
        sample_ctau = ctau.replace('p', '.')

        if is_point_maingrid:
          print '\n --> getting number of miniAOD events per B species' 
          n_miniaod_Bu, n_miniaod_Bd, n_miniaod_Bs = self.getMiniaAODEventsPerBSpecies(mass=mass, ctau=ctau)

        print '\n --> getting filter efficiencies per B species' 
        if is_point_maingrid:
          filter_efficiency_Bu, filter_efficiency_Bd, filter_efficiency_Bs = self.getFilterEfficienciesPerBSpecies(mass=mass, ctau=ctau, is_bc=False)

        if is_point_Bcgrid:
          filter_efficiency_Bc = self.getFilterEfficienciesPerBSpecies(mass=mass, ctau=ctau, is_bc=True)

        print '\n --> writing template' 
        template.append(" SignalSample(")
        if is_point_maingrid:
          template.append("     filename = '{}',".format(sample_filename))
        if is_point_Bcgrid:
          template.append("     filename_Bc = '{}',".format(sample_filename_Bc))
        template.append("     label = '{}',".format(sample_label))
        template.append("     mass = {},".format(sample_mass))
        template.append("     ctau = {},".format(sample_ctau))
        if is_point_maingrid:
          template.append("     filter_efficiency_Bu = {},".format(filter_efficiency_Bu))
          template.append("     filter_efficiency_Bd = {},".format(filter_efficiency_Bd))
          template.append("     filter_efficiency_Bs = {},".format(filter_efficiency_Bs))
        if is_point_Bcgrid:
          template.append("     filter_efficiency_Bc = {},".format(filter_efficiency_Bc))
        if is_point_maingrid:
          template.append("     n_miniaod_Bu = {},".format(n_miniaod_Bu))
          template.append("     n_miniaod_Bd = {},".format(n_miniaod_Bd))
          template.append("     n_miniaod_Bs = {},".format(n_miniaod_Bs))
        if is_point_maingrid:
          template.append("     muon_rate = 1.,")
        if is_point_Bcgrid:
          template.append("     muon_rate_Bc = 0.47,") #FIXME
        template.append("     colour = ROOT.kOrange+0,")
        template.append("     ),")

      template.append("  ]")
      template.append("")

      signal_sample = '\n'.join(template)
      outfile.write(signal_sample)
      signal_samples.append(signal_sample)

    for signal_sample in signal_samples:
      print signal_sample

    outfile.close()
    self.miniaod_file.close()
    print '\n-> signal_samples.txt created'

    print '\nDone'


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  nano_pl = 'V42'
  nano_pl_Bc = 'V42_Bc'

  ntuple_name = 'flat_bparknano_06Feb23_norm.root'
  ntuple_name_Bc = ntuple_name

  nano_general_filename = 'bparknano_generalstep_nj*.root'

  gen_pl = 'V44'
  gen_pl_Bc = 'V41'

  writer = SignalSampleWriter(
    nano_pl = nano_pl,
    nano_pl_Bc = nano_pl_Bc,
    ntuple_name = ntuple_name,
    ntuple_name_Bc = ntuple_name_Bc, 
    nano_general_filename = nano_general_filename,
    gen_pl = gen_pl,
    gen_pl_Bc = gen_pl_Bc,
    )

  writer.process()

