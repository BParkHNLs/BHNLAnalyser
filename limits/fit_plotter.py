import os
import sys
import glob
import re
from os import path
import ROOT

sys.path.append('../scripts')
sys.path.append('../objects')
from tools import Tools
from categories import categories


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Run limit on a single mass/coupling point', add_help=True)
  parser.add_argument('--scenario'          , type=str, dest='scenario'          , help='signal under consideration'                   , default='Majorana', choices=['Majorana', 'Dirac'])
  parser.add_argument('--categories_label ' , type=str, dest='categories_label'  , help='label of the list of categories'              , default='standard')
  parser.add_argument('--mass'              , type=str, dest='mass'              , help='mass'                                         , default='1.0')
  parser.add_argument('--homedir'           , type=str, dest='homedir'           , help='name of the homedir'                          , default=None)
  parser.add_argument('--outdirlabel'       , type=str, dest='outdirlabel'       , help='name of the outdir'                           , default=None)
  parser.add_argument('--subdirlabel'       , type=str, dest='subdirlabel'       , help='name of the subdir'                           , default=None)
  parser.add_argument('--fe'                , type=str, dest='fe'                , help='electron coupling fraction'                   , default='1.0')
  parser.add_argument('--fu'                , type=str, dest='fu'                , help='muon coupling fraction'                       , default='1.0')
  parser.add_argument('--ft'                , type=str, dest='ft'                , help='tau coupling fraction'                        , default='1.0')
  return parser.parse_args()

#TODO adapt to Dirac case


class FitPlotter(object):
  def __init__(self, mass, categories_label, homedir, outdirlabel, subdirlabel, scenario, fe, fu, ft):
    self.tools = Tools()
    self.mass = float(mass)
    self.categories_label = categories_label
    self.categories = categories[self.categories_label]
    self.homedir = homedir
    self.outdirlabel = outdirlabel
    self.subdirlabel = subdirlabel
    self.scenario = scenario
    self.fe = fe
    self.fu = fu
    self.ft = ft
    if self.fe != 'None' and self.fu != 'None' and self.ft != 'None':
      self.do_coupling_scenario = True
    else:
      self.do_coupling_scenario = False
    if self.do_coupling_scenario:
      self.fe = str(round(float(fe), 1)).replace('.', 'p')
      self.fu = str(round(float(fu), 1)).replace('.', 'p')
      self.ft = str(round(float(ft), 1)).replace('.', 'p')
      
    self.inputdir = '{}/outputs/{}/limits/{}/results/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    self.outputdir = '{}/outputs/{}/limits/{}/plots/fits/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))

    # hardcoded options
    self.print_tag = True
    self.CMS_tag = '' #'Preliminary'
    self.lumi = 41.6
    self.flavour_channel = 'dimuon channel'
    # if true, only plot for the ctau the closest to the exclusion
    self.plot_exclusion = True


  def getEventLabel(self, axis):
    label = axis.GetTitle()
    idx1 = label.find('(') + 2
    idx2 = label.find(')', idx1) -1
    bin_size = float(label[idx1:idx2])
    bin_size = round(bin_size, 4)
    event_label = 'Events / {} GeV'.format(bin_size)
    return event_label


  def setStyle(self):
    #ROOT.gStyle.SetPadTopMargin(0.05) 
    #ROOT.gStyle.SetPadBottomMargin(0.13) 
    ROOT.gStyle.SetPadRightMargin(0.04) 
    ROOT.gStyle.SetPadLeftMargin(0.11) 


  def getCoupling(self, couplings, values, crossing=1):
    '''
      Function that returns the coupling at which the limit intersects with 1 (crossing) in the log-log plane
    '''
    # first, search for couplings whose associated limit is the closest to up and down 1
    coupling_up = 0
    coupling_down = 1e9
    value_up = 0
    value_down = 0
    for icoupling, coupling in enumerate(couplings):
      if icoupling+1 >= len(couplings): continue
      if values[icoupling] >= crossing and values[icoupling+1] <= crossing:
        coupling_up = couplings[icoupling+1]
        value_up = values[icoupling+1]
        coupling_down = couplings[icoupling]
        value_down = values[icoupling]
        break

    if abs(value_down -1) < abs(value_up -1):
      coupling = coupling_down
    else:
      coupling = coupling_up

    return coupling


  def getExclusionPoint(self):
    '''
     Search for the coupling point the closest to the exclusion and return 
     the associated fitDiagnostics file
    '''

    # get the files 
    if not self.do_coupling_scenario:
      pathToResults = '{}/outputs/{}/limits/{}/results/'.format(self.homedir, self.outdirlabel, self.subdirlabel) 
    else:
      pathToResults = '{}/outputs/{}/limits/{}/results_{}_{}_{}/'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.fe, self.fu, self.ft) 

    fileName = 'result*{}*m_{}*.txt'.format(self.scenario, self.mass)

    files = [f for f in glob.glob(pathToResults+fileName)]

    v2s       = []
    obs       = []
    for limitFile in files:
      if 'm_{}_'.format(self.mass) not in limitFile: continue
   
      # for each mass, get the list of the couplings and ctaus from the file name
      ctau = limitFile[limitFile.find('ctau_')+5:limitFile.find('_', limitFile.find('ctau_')+5)]
      coupling = limitFile[limitFile.rfind('v2_')+3:limitFile.find('.txt')]
      val_coupling = float(coupling)
    
      try:
        thefile = open('{}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, self.scenario, self.mass, ctau, coupling), 'r')

        # get the necessary information from the result files
        val_obs       = None
        content = thefile.readlines()
        for line in content:
          if 'Observed' in line:
            values = re.findall(r'\d+', line)
            val_obs = values[0] + '.' + values[1]
        obs.append(float(val_obs))
        v2s.append(val_coupling)

      except:
        print 'Cannot open {}/result_{}_m_{}_ctau_{}_v2_{}.txt'.format(pathToResults, self.scenario, self.mass, ctau, coupling)

    graph = zip(v2s, obs)
    graph.sort(key = lambda x : float(x[0])) # sort by coupling
  
    v2s = [jj[0] for jj in graph]
    obs = [jj[1] for jj in graph]
    
    # find the coupling the closest to the exclusion
    coupling_obs = self.getCoupling(v2s, obs)
    coupling_label = self.tools.getCouplingLabel(coupling_obs)
    coupling_label = coupling_label.replace('.', 'p').replace('-', 'm')

    # get the corresponding file
    filename = glob.glob('{}/outputs/{}/limits/{}/results/fitDiagnostics_{}_m_{}_ctau_*_v2_{}.root'.format(self.homedir, self.outdirlabel, self.subdirlabel, self.scenario, str(self.mass).replace('.', 'p'), coupling_label))[0]

    return filename


  def plot(self, filename, category, fit):
    if fit not in ['prefit', 'fit_s', 'fit_b']:
      raise RuntimeError('Unknown fit type "{}". Choose amongst ["prefit", "fit_s", "fit_b"]'.format(fit))

    # setting style
    self.setStyle()

    signal_mass = str(self.mass).replace('.', 'p')
    signal_ctau = filename[filename.rfind('ctau_')+5:filename.rfind('_v2')]
    ctau = float(signal_ctau.replace('p', '.'))
    signal_v2 = self.tools.getVV(mass=self.mass, ctau=ctau, ismaj=True)
    signal_coupling = self.tools.getCouplingLabel(signal_v2)
    signal_coupling = signal_coupling.replace('.', 'p').replace('-', 'm')

    f = ROOT.TFile.Open(filename)
    rooplot_name = '{}_hnl_mass_muon_channel_m_{}_{}'.format(category.label, str(self.mass).replace('.', 'p'), fit)

    try:
      rooplot = f.Get(rooplot_name)

      canv_name = 'canv_{}_{}_{}'.format(category.label, fit, signal_ctau)
      canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
      pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
      pad.Draw()
      pad.cd()

      if fit == 'prefit':
        # do not plot the signal+background component in the prefit plot
        rooplot.remove(rooplot.getObject(4).GetName())
        rooplot.remove(rooplot.getObject(1).GetName())
        idx_sig = 1
        idx_bkg = 2
      else:
        idx_sig = 2
        idx_bkg = 3
        idx_sig_plus_bkg = 4

      rooplot.SetTitle('')
      rooplot.GetXaxis().SetTitle('m(#mu^{#pm}#pi^{#mp}) (GeV)')
      rooplot.GetXaxis().SetLabelSize(0.04)
      rooplot.GetXaxis().SetTitleSize(0.047)
      rooplot.GetXaxis().SetTitleOffset(1.0)
      #rooplot.GetYaxis().SetTitle(self.getEventLabel(rooplot.GetYaxis()))
      rooplot.GetYaxis().SetTitle('Events / Bin')
      rooplot.GetYaxis().SetLabelSize(0.04)
      rooplot.GetYaxis().SetTitleSize(0.047)
      rooplot.GetYaxis().SetTitleOffset(1.1)
      rooplot.GetYaxis().SetRangeUser(0, rooplot.GetMaximum() + 0.7*rooplot.GetMaximum())
      rooplot.Draw()

      self.tools.printLatexBox(0.15, 0.7, category.title, size=0.04, pos='left', font=42)
      if 'Bc' in category.label:
        b_mass_label = '#mu_{P}#mu#pi mass > 5.7 GeV'
      else:
        b_mass_label = '#mu_{P}#mu#pi mass #leq 5.7 GeV'
      self.tools.printLatexBox(0.15, 0.65, b_mass_label, size=0.04, pos='left', font=42)
      self.tools.printLatexBox(0.15, 0.6, self.flavour_channel, size=0.04, pos='left', font=42)

      # print the CMS tag
      self.tools.printInnerCMSTag(pad, self.CMS_tag, self.print_tag, x_pos=0.15, y_pos=0.83, size=0.55)
    
      # print lumi tag
      self.tools.printLumiTag(pad, self.lumi, size=0.5, offset=0.57)
          
      # print the legend
      if fit == 'fit_s':
        leg = self.tools.getRootTLegend(xmin=0.55, ymin=0.63, xmax=0.87, ymax=0.88, size=0.04)
        leg.AddEntry(rooplot.findObject(rooplot.getObject(0).GetName()), 'data')
        leg.AddEntry(rooplot.findObject(rooplot.getObject(idx_sig).GetName()), 'signal fit')
        leg.AddEntry(rooplot.findObject(rooplot.getObject(idx_bkg).GetName()), 'background fit')
        leg.AddEntry(rooplot.findObject(rooplot.getObject(4).GetName()), 'signal+background fit')
      else:
        leg = self.tools.getRootTLegend(xmin=0.47, ymin=0.63, xmax=0.83, ymax=0.85, size=0.04)
        leg.AddEntry(rooplot.findObject(rooplot.getObject(0).GetName()), 'data')
        leg.AddEntry(rooplot.findObject(rooplot.getObject(idx_bkg).GetName()), 'background prediction')
        leg.AddEntry(rooplot.findObject(rooplot.getObject(idx_sig).GetName()), 'signal - {} GeV, {} mm'.format(self.mass, ctau))
      leg.Draw()

      canv.cd()
      plot_name = 'hnl_mass_m_{}_ctau_{}_cat_{}'.format(signal_mass, signal_ctau, category.label)
      if fit == 'fit_s':
        plot_name += '_postfit'
      else:
        plot_name += '_prefit'
      canv.SaveAs('{}/{}.png'.format(self.outputdir, plot_name))
      canv.SaveAs('{}/{}.pdf'.format(self.outputdir, plot_name))
      canv.SaveAs('{}/{}.C'.format(self.outputdir, plot_name))

    except:
      print 'Rooplot "{}" not found'.format(rooplot_name)


  def process(self):
    if self.plot_exclusion:
      filenames = [self.getExclusionPoint()]
    else:
      filenames = glob.glob('{}/fitDiagnostics_{}_m_{}_*.root'.format(self.inputdir, self.scenario, str(self.mass).replace('.', 'p')))

    if len(filenames) == 0:
      print 'no files found'

    for filename in filenames:
      for category in self.categories:
        if 'incl' in category.label: continue
        if self.mass < 3 and 'Bc' in category.label: continue
        if 'gt150_OS' not in category.label: continue

        self.plot(filename=filename, category=category, fit='prefit')
        self.plot(filename=filename, category=category, fit='fit_s')



if __name__ == '__main__': 
  ROOT.gROOT.SetBatch(True)

  opt = getOptions()

  plotter = FitPlotter(
      mass = opt.mass,
      categories_label = opt.categories_label,
      homedir = opt.homedir,
      outdirlabel= opt.outdirlabel,
      subdirlabel = opt.subdirlabel,
      scenario = opt.scenario,
      fe = opt.fe,
      fu = opt.fu, 
      ft = opt.ft,
      )
  plotter.process()

