from common import Tools
import ROOT

class MatchingChecks(Tools):
  def __init__(self):
    self.tools = Tools()

  def plotMatchingDiagram(self, filename, treename, outdirlabel, particle):
    f = ROOT.TFile.Open(filename, 'READ')
    tree = self.tools.getTree(f, treename)

    outputdir = self.tools.getOutDir('./myPlots/matching_checks', outdirlabel)

    for ientry, entry in enumerate(tree):
      if ientry > (100 if particle=='muon' else 200): continue
      # plot the reco particles (for events where a candidate exists)
      pt_box_reco = []
      deltaR_limits = []
      canv = self.tools.createTCanvas('canv_{}'.format(ientry))
      graph_reco = ROOT.TGraph() 
      flag_reco = False
      if entry.nBToMuMuPi < 1: continue
      flag_reco = True
      #for imuon in range(0, entry.nMuon):
      for ipart in range(0, entry.nBToMuMuPi):
        #print 'reco pt {} eta {} phi {}'.format(entry.Muon_pt[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]])
        #point = graph_reco.GetN()
        if particle == 'muon':
          point = graph_reco.GetN()
          graph_reco.SetPoint(point, entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]])
          pt_info = ROOT.TLatex(entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]]+0.05, entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]]+0.15, '{}'.format(round(entry.Muon_pt[entry.BToMuMuPi_sel_mu_idx[ipart]], 3)))
          deltaR_circle = ROOT.TEllipse(entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]], 0.3, 0.3)
          pt_box_reco.append(pt_info)
          deltaR_limits.append(deltaR_circle)
        else:
          point = graph_reco.GetN()
          graph_reco.SetPoint(point, entry.ProbeTracks_eta[entry.BToMuMuPi_pi_idx[ipart]], entry.ProbeTracks_phi[entry.BToMuMuPi_pi_idx[ipart]])
          pt_info = ROOT.TLatex(entry.ProbeTracks_eta[entry.BToMuMuPi_pi_idx[ipart]]+0.05, entry.ProbeTracks_phi[entry.BToMuMuPi_pi_idx[ipart]]+0.15, '{}'.format(round(entry.ProbeTracks_pt[entry.BToMuMuPi_pi_idx[ipart]], 3)))
          deltaR_circle = ROOT.TEllipse(entry.ProbeTracks_eta[entry.BToMuMuPi_pi_idx[ipart]], entry.ProbeTracks_phi[entry.BToMuMuPi_pi_idx[ipart]], 0.1, 0.1)
          pt_box_reco.append(pt_info)
          deltaR_limits.append(deltaR_circle)
        
      # trigger muon
      if particle == 'muon':
        graph_reco_trg = ROOT.TGraph() 
        for imuon in range(0, entry.nBToMuMuPi):
          point = graph_reco_trg.GetN()
          graph_reco_trg.SetPoint(point, entry.Muon_eta[entry.BToMuMuPi_trg_mu_idx[imuon]], entry.Muon_phi[entry.BToMuMuPi_trg_mu_idx[imuon]])
          pt_info = ROOT.TLatex(entry.Muon_eta[entry.BToMuMuPi_trg_mu_idx[imuon]]+0.05, entry.Muon_phi[entry.BToMuMuPi_trg_mu_idx[imuon]]+0.15, '{}'.format(round(entry.Muon_pt[entry.BToMuMuPi_trg_mu_idx[imuon]], 3)))
          pt_box_reco.append(pt_info)

      # plot the matched reco particles
      pt_box_matched = []
      graph_matched = ROOT.TGraph() 
      flag_matched = False
      for ipart in range(0, entry.nBToMuMuPi):
        if particle == 'muon':
          if entry.BToMuMuPi_sel_mu_isMatched[ipart] != 1: continue
          flag_matched = True
          #if entry.BToMuMuPi_sel_mu_isMatched[entry.BToMuMuPi_sel_mu_idx] != 1: continue
          point = graph_matched.GetN()
          graph_matched.SetPoint(point, entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]])
          pt_info = ROOT.TLatex(entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]]+0.05, entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]]+0.15, '{}'.format(round(entry.Muon_pt[entry.BToMuMuPi_sel_mu_idx[ipart]], 3)))
        else:
          if entry.BToMuMuPi_pi_isMatched[ipart] != 1: continue
          flag_matched = True
          point = graph_matched.GetN()
          graph_matched.SetPoint(point, entry.ProbeTracks_eta[entry.BToMuMuPi_pi_idx[ipart]], entry.ProbeTracks_phi[entry.BToMuMuPi_pi_idx[ipart]])
          pt_info = ROOT.TLatex(entry.ProbeTracks_eta[entry.BToMuMuPi_pi_idx[ipart]]+0.05, entry.ProbeTracks_phi[entry.BToMuMuPi_pi_idx[ipart]]+0.15, '{}'.format(round(entry.ProbeTracks_pt[entry.BToMuMuPi_pi_idx[ipart]], 3)))
        pt_box_matched.append(pt_info)

      # trigger muon
      if particle == 'muon':
        graph_matched_trg = ROOT.TGraph() 
        flag_matched_trg = False
        for imuon in range(0, entry.nBToMuMuPi):
          if entry.BToMuMuPi_trg_mu_isMatched[imuon] != 1: continue
          flag_matched_trg = True
          point = graph_matched_trg.GetN()
          graph_matched_trg.SetPoint(point, entry.Muon_eta[entry.BToMuMuPi_trg_mu_idx[imuon]], entry.Muon_phi[entry.BToMuMuPi_trg_mu_idx[imuon]])
          pt_info = ROOT.TLatex(entry.Muon_eta[entry.BToMuMuPi_trg_mu_idx[imuon]]+0.05, entry.Muon_phi[entry.BToMuMuPi_trg_mu_idx[imuon]]+0.15, '{}'.format(round(entry.Muon_pt[entry.BToMuMuPi_trg_mu_idx[imuon]], 3)))
          pt_box_matched.append(pt_info)

      # plot the gen particles
      pt_box_gen = []
      flag_gen = False
      graph_gen = ROOT.TGraph() 
      for igen in range(0, entry.nGenPart):
        pdgId = 13 if particle == 'muon' else 211
        if abs(entry.GenPart_pdgId[igen]) != pdgId: continue
        #if abs(entry.GenPart_pdgId[entry.GenPart_genPartIdxMother[igen]]) not in [9900015, 521, 511]: continue
        flag_gen = True
        #print 'gen pt {} eta {} phi {}'.format(entry.GenPart_pt[igen], entry.GenPart_eta[igen], entry.GenPart_phi[igen])
        point = graph_gen.GetN()
        graph_gen.SetPoint(point, entry.GenPart_eta[igen], entry.GenPart_phi[igen])
        pt_info = ROOT.TLatex(entry.GenPart_eta[igen]-0.05, entry.GenPart_phi[igen]-0.25, '{}, {}'.format(round(entry.GenPart_pt[igen], 3), entry.GenPart_pdgId[entry.GenPart_genPartIdxMother[igen]]))
        pt_box_gen.append(pt_info)

      if not flag_reco or not flag_gen: continue
      graph_reco.Draw('AP')  
      if particle == 'muon': graph_reco_trg.Draw('P same')  
      if flag_matched: graph_matched.Draw('P same')  
      if particle == 'muon' and flag_matched_trg: graph_matched_trg.Draw('P same')  
      graph_gen.Draw('P same')  

      for pt_info in pt_box_reco:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kBlue+2)
        pt_info.SetTextAlign(21)

      for pt_info in pt_box_matched:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kOrange+2)
        pt_info.SetTextAlign(21)

      for pt_info in pt_box_gen:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kGreen+2)
        pt_info.SetTextAlign(21)

      for deltaR_circle in deltaR_limits:
        #deltaR_circle.Draw()
        deltaR_circle.SetLineColor(2)
        deltaR_circle.SetFillColorAlpha(0, 0)

      graph_reco.SetMarkerStyle(41)
      graph_reco.SetMarkerSize(4)
      graph_reco.SetMarkerColor(ROOT.kBlue+2)

      if particle == 'muon':
        graph_reco_trg.SetMarkerStyle(33)
        graph_reco_trg.SetMarkerSize(4)
        graph_reco_trg.SetMarkerColor(ROOT.kBlue+2)

      graph_matched.SetMarkerStyle(41)
      graph_matched.SetMarkerSize(4)
      graph_matched.SetMarkerColor(ROOT.kOrange+2)

      if particle == 'muon':
        graph_matched_trg.SetMarkerStyle(33)
        graph_matched_trg.SetMarkerSize(4)
        graph_matched_trg.SetMarkerColor(ROOT.kOrange+2)

      graph_gen.SetMarkerStyle(43)
      graph_gen.SetMarkerSize(4)
      graph_gen.SetMarkerColor(ROOT.kGreen+2)

      graph_reco.SetTitle('Event {}'.format(ientry))
      graph_reco.GetXaxis().SetTitle('#eta')
      graph_reco.GetXaxis().SetLabelSize(0.037)
      graph_reco.GetXaxis().SetTitleSize(0.042)
      graph_reco.GetXaxis().SetTitleOffset(1.1)
      #graph_reco.GetXaxis().SetRangeUser(-5, 5)
      graph_reco.GetXaxis().SetLimits(-5.,5.)
      graph_reco.GetYaxis().SetTitle('#phi')
      graph_reco.GetYaxis().SetLabelSize(0.037)
      graph_reco.GetYaxis().SetTitleSize(0.042)
      graph_reco.GetYaxis().SetTitleOffset(1.1)
      graph_reco.GetYaxis().SetRangeUser(-3.2, 3.2)

      canv.SaveAs('{}/event_{}.png'.format(outputdir, ientry))


  def plotTriggerMatchingDiagram(self, filename, treename, outdirlabel):
    f = ROOT.TFile.Open(filename, 'READ')
    tree = self.tools.getTree(f, treename)

    outputdir = self.tools.getOutDir('./myPlots/trigger_matching_checks', outdirlabel)

    for ientry, entry in enumerate(tree):
      if ientry > 100: continue
      # plot the unmatched muon
      pt_box_reco = []
      #deltaR_limits = []
      canv = self.tools.createTCanvas('canv_{}'.format(ientry))
      graph_reco = ROOT.TGraph() 
      flag_reco = False
      #if entry.nBToMuMuPi < 1: continue
      for ipart in range(0, entry.nMuon):
        #print 'reco pt {} eta {} phi {}'.format(entry.Muon_pt[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]])
        #point = graph_reco.GetN()
        if entry.Muon_isTriggering[ipart] != 1 and entry.Muon_isDSAMuon[ipart] != 1:
          flag_reco = True
          point = graph_reco.GetN()
          graph_reco.SetPoint(point, entry.Muon_eta[ipart], entry.Muon_phi[ipart])
          pt_info = ROOT.TLatex(entry.Muon_eta[ipart]+0.05, entry.Muon_phi[ipart]+0.15, '{}'.format(round(entry.Muon_pt[ipart], 3)))
          #deltaR_circle = ROOT.TEllipse(entry.Muon_eta[entry.BToMuMuPi_sel_mu_idx[ipart]], entry.Muon_phi[entry.BToMuMuPi_sel_mu_idx[ipart]], 0.3, 0.3)
          pt_box_reco.append(pt_info)
          #deltaR_limits.append(deltaR_circle)
      if not flag_reco:
        graph_reco.SetPoint(0, 0, 0)
        
      # plot the matched muon
      pt_box_matched = []
      graph_matched = ROOT.TGraph() 
      flag_matched = False
      for ipart in range(0, entry.nMuon):
        if entry.Muon_isTriggering[ipart] == 1 and entry.Muon_isDSAMuon[ipart] != 1:
          flag_matched = True
          point = graph_matched.GetN()
          graph_matched.SetPoint(point, entry.Muon_eta[ipart], entry.Muon_phi[ipart])
          pt_info = ROOT.TLatex(entry.Muon_eta[ipart]+0.05, entry.Muon_phi[ipart]+0.15, '{}'.format(round(entry.Muon_pt[ipart], 3)))
          pt_box_matched.append(pt_info)

      # plot the trigger objects
      pt_box_gen = []
      flag_gen = False
      graph_trgobj = ROOT.TGraph() 
      for itrg in range(0, entry.nTrigObj):
        if entry.TrigObj_filterBits[itrg] != 1: continue
        flag_gen = True
        point = graph_trgobj.GetN()
        graph_trgobj.SetPoint(point, entry.TrigObj_eta[itrg], entry.TrigObj_phi[itrg])
        pt_info = ROOT.TLatex(entry.TrigObj_eta[itrg]-0.05, entry.TrigObj_phi[itrg]-0.25, '{}'.format(round(entry.TrigObj_pt[itrg], 3)))
        pt_box_gen.append(pt_info)

      #if not flag_reco or not flag_gen: continue
      #if not flag_reco: continue
      graph_reco.Draw('AP')  
      if flag_matched: graph_matched.Draw('P same')  
      if flag_gen: graph_trgobj.Draw('P same')  

      for pt_info in pt_box_reco:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kBlue+2)
        pt_info.SetTextAlign(21)

      for pt_info in pt_box_matched:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kOrange+2)
        pt_info.SetTextAlign(21)

      for pt_info in pt_box_gen:
        pt_info.Draw()
        pt_info.SetTextSize(0.02)
        pt_info.SetTextFont(62)
        pt_info.SetTextColor(ROOT.kGreen+2)
        pt_info.SetTextAlign(21)

      #for deltaR_circle in deltaR_limits:
      #  #deltaR_circle.Draw()
      #  deltaR_circle.SetLineColor(2)
      #  deltaR_circle.SetFillColorAlpha(0, 0)

      graph_reco.SetMarkerStyle(41)
      graph_reco.SetMarkerSize(4)
      graph_reco.SetMarkerColor(ROOT.kBlue+2)

      graph_matched.SetMarkerStyle(41)
      graph_matched.SetMarkerSize(4)
      graph_matched.SetMarkerColor(ROOT.kOrange+2)

      graph_trgobj.SetMarkerStyle(43)
      graph_trgobj.SetMarkerSize(4)
      graph_trgobj.SetMarkerColor(ROOT.kGreen+2)

      graph_reco.SetTitle('Event {}'.format(ientry))
      graph_reco.GetXaxis().SetTitle('#eta')
      graph_reco.GetXaxis().SetLabelSize(0.037)
      graph_reco.GetXaxis().SetTitleSize(0.042)
      graph_reco.GetXaxis().SetTitleOffset(1.1)
      #graph_reco.GetXaxis().SetRangeUser(-5, 5)
      graph_reco.GetXaxis().SetLimits(-3.1,3.1)
      graph_reco.GetYaxis().SetTitle('#phi')
      graph_reco.GetYaxis().SetLabelSize(0.037)
      graph_reco.GetYaxis().SetTitleSize(0.042)
      graph_reco.GetYaxis().SetTitleOffset(1.1)
      graph_reco.GetYaxis().SetRangeUser(-3.2, 3.2)

      canv.SaveAs('{}/event_{}.png'.format(outputdir, ientry))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  plot_matchingDiagram = False
  plot_triggerMatchingDiagram = True

  if plot_matchingDiagram:
    plotter = MatchingChecks()
    #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/bparknano.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_loosepreselection_stdtrgmu_v1.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_mu_0p15_0p25.root'
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_mu_0p15_0p25_pi_0p15_0p5.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_updatedgenmatching_mu_0p1_0p25_pi_0p15_0p5_massreldiff_0p1.root'
    #filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged//bparknano_looseselection_stdmatching_0p5.root'
    #plotter.plotMatchingDiagram(filename=filename, treename='Events', outdirlabel='displaced_dsamuon_0p5', particle='muon')
    plotter.plotMatchingDiagram(filename=filename, treename='Events', outdirlabel='displaced_pion_0p15_0p5', particle='pion')


  if plot_triggerMatchingDiagram:
    plotter = MatchingChecks()
    #filename = '/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/bparknano.root'
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_allevts_trgmumatchingchecks_false.root'
    plotter.plotTriggerMatchingDiagram(filename=filename, treename='Events', outdirlabel='mc_matchcond_false')
