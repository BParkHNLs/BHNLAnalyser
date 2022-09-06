#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "tdrstyle.C"
#include "CMS_lumi.C"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = false;
bool runFtestCheckWithToys=false;

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
/* Construct the pdf model using the PdfModelBuilder */
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Polynomial") return pdfsModel.getPolynomial(Form("%s_poly%d",ext,order),order); 
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, TString fit_range){
/* Basic fitting routine, fit is not extended */

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  //params_test->Print("v");
  int stat=1;
  double minnll=10e8;
  while (stat!=0){
    if (ntries>=MaxTries) break;
    RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Range(fit_range),RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize")); 
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name, int nBinsForFit, TString fit_range){
/* Get probability of the f-test. Currently toys are not used and only simple TMath::Prob is used. */
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
  
  // fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitNullData = pdfNull->fitTo(*data, RooFit::Range(fit_range), RooFit::Save(1), RooFit::Strategy(1), RooFit::Minimizer("Minuit2","minimize"), RooFit::PrintLevel(-1)); //FIXME
  RooFitResult *fitTestData = pdfTest->fitTo(*data, RooFit::Range(fit_range), RooFit::Save(1), RooFit::Strategy(1), RooFit::Minimizer("Minuit2","minimize"), RooFit::PrintLevel(-1)); 

  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys = 5000;
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<toyhist.GetNbinsX();b++){
    double x = toyhist.GetBinCenter(b+1);
    if (x>0){
      gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
      ipoint++;
    }
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForFit);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);
    RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

    int stat_n=1;
    int stat_t=1;
    int ntries = 0;
    double nllNull,nllTest;
    // Iterate on the fit 
    int MaxTries = 2;
    while (stat_n!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy, RooFit::Range(fit_range), RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
      //,RooFit::Optimize(0));

      nllNull = fitNull->minNll();
      stat_n = fitNull->status();
      if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
      ntries++; 
    }
  
    ntries = 0;
    while (stat_t!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy, RooFit::Range(fit_range), RooFit::Save(1),RooFit::Strategy(1)
                                             ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
      nllTest = fitTest->minNll();
      stat_t = fitTest->status();
      if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
      ntries++; 
    }
       
    toyhistStatN.Fill(stat_n);
    toyhistStatT.Fill(stat_t);

    if (stat_t !=0 || stat_n !=0) continue;
    nsuccesst++;
    double chi2_t = 2*(nllNull-nllTest);
    if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);

  } // end loop over toys

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  std::cout << "debug get probability f test: " << name.c_str() << std::endl; 
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, int nBinsForFit){
/* Get goodness of fit, based on chi-square, using binned dataset and fitted pdf 
   use toys or chi-square distributions depending on avg number of events in bin */

  double prob;
  int ntoys = 500;
  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForFit),Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForFit < 5 ){

    // for the moment, do not compute p-value when statistics is too low, 
    // and remove criterium on the p-value to enter the envelope for those cases
    prob = -1.;

    /*
    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
    //  std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0)); 

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForFit-np));
      //TCanvas *can = new TCanvas();
      //plot_t->Draw();
      //can->SaveAs("test.png");
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForFit-np),toyhist.GetMaximum(),chi2*(nBinsForFit-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit   
    params->assignValueOnly(preParams);
    */
  } else {
    prob = TMath::Prob(chi2*(nBinsForFit-np),nBinsForFit-np);
  }
  std::cout << "[INFO] GOF Chi2 in Observed =  " << chi2*(nBinsForFit-np) << std::endl;
  std::cout << "[INFO] GOF p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> category_label, float fit_window_min, float fit_window_max, float mass_window_min, float mass_window_max, int nBinsForFit, int nBinsForPlot, int status, double *prob){
/* Plot single pdf vs data, with pulls */
    
  // Chi2 taken from full range fit
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForFit));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name,nBinsForFit);
  RooPlot *plot = mass->frame();
  mass->setRange("sideband_left", fit_window_min, mass_window_min);
  mass->setRange("sideband_right", mass_window_max, fit_window_max);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_left"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_right"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForPlot));

  TCanvas *canv = new TCanvas();
  canv->Divide(1,2);
  canv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.10);
  gPad->SetPad(0.01,0.2,0.99,0.99);
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitle("m_{#mu#pi} (GeV)");
  plot->GetYaxis()->SetTitleSize(0.04);;
  //plot->GetYaxis()->SetTitle("Entries");
  plot->GetYaxis()->SetTitleOffset(1.35);

  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.34,0.85,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  plot->getAttText()->SetTextSize(0.025);
  plot->SetMaximum(plot->GetMaximum()*1.4);
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.034);
  lat->DrawLatex(0.15,0.94,Form("#chi^{2}/ndof = %.3f, Prob = %.2f, Fit Status = %d ",chi2,*prob,status));

  TLatex *lab = new TLatex();
  lab->SetNDC();
  lab->SetTextFont(42);
  lab->SetTextSize(0.034);
  std::string delimiter = "/";
  size_t pos = 0;
  std::string s = name;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
  }
  //delimiter = ".";
  //token = s.substr(0, s.find(delimiter));
  //std::cout << s << std::endl;
  lab->DrawLatex(0.55,0.25,s.c_str());

  // pulls
  canv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetPad(0.01,0.01,0.99,0.2);
  gPad->SetGridy();
  RooHist* hpull = plot->pullHist();
  RooPlot *plot2 = mass->frame();
  plot2->GetYaxis()->SetNdivisions(504);
  plot2->GetYaxis()->SetLabelSize(0.17);
  plot2->GetYaxis()->SetTitleSize(0.17);
  plot2->GetYaxis()->SetTitleOffset(0.24);
  //plot2->GetYaxis()->SetRangeUser(-3.0,3.0);
  plot2->SetMinimum(-3.);
  plot2->SetMaximum(3.);
  plot2->GetYaxis()->SetTitle("Pulls") ;
  plot2->GetXaxis()->SetTitle("");
  plot2->GetXaxis()->SetLabelOffset(5);
  plot2->addPlotable(hpull,"P"); 
  plot2->Draw();

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  //plot_chi2->Draw();
  //canv->SaveAs((name+"debug").c_str());

  delete canv;
  delete lat;
}


void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> category_label, int cat, float fit_window_min, float fit_window_max, float mass_window_min, float mass_window_max, int nBinsForPlot, TString fit_range, int bestFitPdf=-1){
/* Plot MultiPdf vs data */
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TLegend *leg = new TLegend(0.6,0.65,0.95,0.90);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();

  mass->setRange("sideband_left", fit_window_min, mass_window_min);
  mass->setRange("sideband_right", mass_window_max, fit_window_max);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_left"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_right"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForPlot)); 
  TCanvas *canv = new TCanvas();
  //TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
  //pad1->SetBottomMargin(0.18);
  //pad1->Draw();
  //pad1->cd();

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - %s",category_label[cat].c_str()),"LEP");
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  int bestcol= -1;
  int style=1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Range(fit_range),RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"));  
    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
      ext=" (Best Fit Pdf) ";
      pdf= pdfs->getCurrentPdf();
      nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
      bestcol = col;
    }
    string pdfName = pdfs->getCurrentPdf()->GetName(); 
    std::string token;
    std::string delimiter = "_";
    size_t pos = 0;
    while ((pos = pdfName.find(delimiter)) != std::string::npos) {
      token = pdfName.substr(0, pos);
      pdfName.erase(0, pos + delimiter.length());
    }
    leg->AddEntry(pdfLeg,Form("%s%s",pdfName.c_str(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %s",category_label[cat].c_str()));
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->GetXaxis()->SetTitle("m_{#mu#pi} (GeV)");
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}


void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> category_label, int cat, float fit_window_min, float fit_window_max, float mass_window_min, float mass_window_max, int nBinsForPlot, int bestFitPdf=-1){
/* Plot several Pdfs vs data, without ratio plot, (used for the "truth") */
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("sideband_left", fit_window_min, mass_window_min);
  mass->setRange("sideband_right", mass_window_max, fit_window_max);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_left"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("sideband_right"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForPlot));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  if(category_label.size() >0){
  leg->AddEntry(datLeg,Form("Data - %s",category_label[cat].c_str()),"LEP");
  } else {
  leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
  }
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->SetTitle(Form(" %s",category_label[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, TString fit_range, bool silent=false){
/* Get index of the best fit pdf (minimum NLL, including correction) among functions in the multipdf.
   All fits are performed again. */

  double global_minNll = 1E10;
  int best_index = 0;
  int number_of_indeces = cat->numTypes();
    
  RooArgSet snap,clean;
  RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
  params->remove(*cat);  // pdf_index is RooCategory, removed from parameters
  params->snapshot(snap);
  params->snapshot(clean);
  if (!silent) {
    //params->Print("V");
  }
 
  // Uncomment to try to make converge a failed fit
  //bkg->setDirtyInhibit(1);
  //RooAbsReal *nllm = bkg->createNLL(*data);
  //RooMinimizer minim(*nllm);
  //minim.setStrategy(1);
  
  for (int id=0;id<number_of_indeces;id++){    
    params->assignValueOnly(clean);
    cat->setIndex(id);

    //RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

    if (!silent) {
      /*
      std::cout << "BEFORE  MAKING FIT" << std::endl;
      params->Print("V");
      std::cout << "-----------------------" << std::endl;    
      */
    }
    
    //minim.minimize("Minuit2","minimize");
    double minNll=0; //(nllm->getVal())+bkg->getCorrection();
    int fitStatus=1;    
    runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/7, fit_range);
    // Add the penalty

    minNll=minNll+bkg->getCorrection();

    if (!silent) {
      /*
      std::cout << "After Minimization ------------------  " <<std::endl;
      std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
      bkg->Print("v");
      bkg->getCurrentPdf()->getParameters(*data)->Print("V");
      std::cout << " ------------------------------------  " << std::endl;
  
      params->Print("V");
      */
      std::cout << "[INFO] AFTER FITTING" << std::endl;
      std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
      std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
      std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
      std::cout << "-----------------------" << std::endl;
    }
      
    if (minNll < global_minNll){
      global_minNll = minNll;
      snap.assignValueOnly(*params);
      best_index=id;
    }
  } // end loop over pdf_index 
  cat->setIndex(best_index);
  params->assignValueOnly(snap);

  std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
  std::cout << "[INFO] Best fit parameters " << std::endl;
  params->Print("V");
  
  return best_index;
}

int main(int argc, char* argv[]){
 
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  //TODO add lumi_target?

  string fileName;
  int ncats = 1;
  int catOffset = 0;
  //string datfile;
  string outDir;
  string outfilename;
  bool verbose=false;
  bool saveMultiPdf=false;
  string category_label_str;
  vector<string> category_label;
  float mN;
  string mN_label;
  float resolution;
  int mass_window_size;
  int fit_window_size;
  int nbins;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    //("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),                                           "Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys",                                                                   "When running the F-test, use toys to calculate pvals (and make plots) ")
    ("blind",                                                                                   "blind plots")
    ("category_label", po::value<string>(&category_label_str),                                  "Category label")
    ("mN", po::value<float>(&mN),                                                               "Mass of the peak, for center of window")
    ("mN_label", po::value<string>(&mN_label),                                                  "Mass label")
    ("resolution", po::value<float>(&resolution),                                               "Resolution (=sigma) of the signal hypothesis")
    ("mass_window_size", po::value<int>(&mass_window_size),                                     "Sigma multiplier for the mass window")
    ("fit_window_size", po::value<int>(&fit_window_size),                                       "Sigma multiplier for the fit window")
    ("nbins", po::value<int>(&nbins),                                                           "Number of bins for the fit and plotting")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("blind")) BLIND=true;
  saveMultiPdf = vm.count("saveMultiPdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  std::cout << "DEBUG mN=" << mN << std::endl; 

  float fit_window_min = mN - fit_window_size * resolution;
  float fit_window_max = mN + fit_window_size * resolution;
  float mass_window_min = mN - mass_window_size * resolution;
  float mass_window_max = mN + mass_window_size * resolution;

  // apply same binning for fit and plotting
  int nBinsForFit = nbins;
  int nBinsForPlot = nbins;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
  split(category_label,category_label_str,boost::is_any_of(","));
  
  int startingCategory=0;

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  outputws = new RooWorkspace(); outputws->SetName("workspace");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("fTest_workspace");
  if (verbose) std::cout << "[INFO]  inWS open " << inWS << std::endl;
  if (saveMultiPdf){
    transferMacros(inFile,outputfile);

    RooRealVar *intL; 
    RooRealVar *sqrts;

    intL  = intLumi_;
    sqrts = (RooRealVar*)inWS->var("SqrtS");
    if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;

    outputws->import(*intL);
    outputws->import(*sqrts);
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
  }

  // Set up which families of functions you want to test
  vector<string> functionClasses;
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  functionClasses.push_back("PowerLaw");
  functionClasses.push_back("Laurent");
  //functionClasses.push_back("Chebychev");
  //functionClasses.push_back("Polynomial");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Bernstein","pol"));
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("PowerLaw","pow"));
  namingMap.insert(pair<string,string>("Laurent","lau"));
  //namingMap.insert(pair<string,string>("Chebychev","che"));
  //namingMap.insert(pair<string,string>("Polynomial","pol"));

  FILE *resFile ;
  resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("hnl_mass");
  std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.1; // upper threshold on prob_ftest to include function in envelope (looser than truth function)
  double minGofThreshold = 0.01;  // minimal goodness of fit to include function in envelope

  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
  fprintf(resFile,"\\hline\n");

  std::string ext = "13TeV";

  // define fit range
  mass->setRange("full", fit_window_min, fit_window_max);
  mass->setRange("sideband_left", fit_window_min, mass_window_min);
  mass->setRange("sideband_right", mass_window_max, fit_window_max);
  TString fit_range;
  if(BLIND){
    fit_range = "sideband_left,sideband_right";
  }
  else{
    fit_range = "full";
  }

  std::cout << "[INFO] Number of categories to process: " << ncats << std::endl;
  for (int cat=startingCategory; cat<ncats; cat++){

    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    string catname = Form("%s",category_label[cat].c_str()); //TODO keep?

    // Option 1: Use as input an unbinned RooDataSet and bin it
    /*
    RooDataSet *dataFull;
    RooDataSet *dataFull0;
    if (isData_) {
    dataFull = (RooDataSet*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("Data_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }
    else 
    {dataFull = (RooDataSet*)inWS->data(Form("data_mass_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("data_mass_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }

    mass->setBins(nBinsForFit);
    RooDataSet *data;
    string thisdataBinned_name;

    if ( isFlashgg_){
      thisdataBinned_name =Form("CAT_roohist_data_mass_%s",category_label[cat].c_str());
    } else {
      thisdataBinned_name= Form("CAT_roohist_data_mass_cat%d",cat);
    }
    RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
    data = (RooDataSet*)&thisdataBinned; 
    // both "data" and "thisdataBinned" are binned (number of bins is given by "mass" bins) 
    // "dataFull" is unbinned
    */

    // Option 2 (equivalent): Use as input a binned RooDataHist as is
    //string data_name = Form("CAT_roohist_data_mass_%s",category_label[cat].c_str()); 
    //RooDataHist *thisdataBinned = (RooDataHist*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
    RooDataHist *thisdataBinned = (RooDataHist*)inWS->data("hnl_mass_rdh");
    RooDataSet *data = (RooDataSet*)thisdataBinned;

    RooArgList storedPdfs("store");

    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");

    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;

    // Standard F-Test to find the truth functions
    for (vector<string>::iterator funcType=functionClasses.begin(); 
        funcType!=functionClasses.end(); funcType++){

      std::cout << "======================================= " << std::endl;
      std::cout << "====> FAMILY " << funcType->c_str() << std::endl;
      std::cout << "======================================= " << std::endl;

      double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
      int order=1; int prev_order=0; int cache_order=0;

      RooAbsPdf *prev_pdf=NULL;
      RooAbsPdf *cache_pdf=NULL;
      std::vector<int> pdforders;

      std::cout << "===> F-TEST for Truth determination" << std::endl;

      int counter =0;
      while (prob<0.05 && order < 7){ 
        cout << "==> " << *funcType << " " << order << endl;
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()));
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {

          //bkgPdf->Print();
          int fitStatus = 0;
          runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7, fit_range);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
       
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          if (prev_pdf!=NULL){
            prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)),nBinsForFit,fit_range);
            std::cout << "[INFO] F-test Prob == " << prob << std::endl;
          } else {
            prob = 0;
          }
          // otherwise we get it later ...
          //if (!saveMultiPdf) plot(mass,bkgPdf,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),category_label,fitStatus,&gofProb);
          cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
          cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
        counter++;
      } // end condition for performing f-test

      // next line is commented, as we want to save only the final result (that takes into account both GOF and F-test results
      //fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

      int truthOrder = cache_order;

      // Now run loop to determine functions inside envelope
      std::cout << "===> F-TEST and GOF for ENVELOPE determination" << std::endl;
      if (saveMultiPdf){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        prob=0.;
        order=1;
        prev_order=0;
        cache_order=0;
        std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

        while (prob<upperEnvThreshold){
          cout << "==> " << *funcType << " " << order << endl;
          RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()));
          if (!bkgPdf ){
            // assume this order is not allowed
            if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
            order++;
          }
          else {
            // Fit and chi-square calculation is repeated
            //RooFitResult *fitRes;
            int fitStatus=0;
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7, fit_range);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
            //thisNll = fitRes->minNll();
            if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
            double myNll = 2.*thisNll;
            chi2 = 2.*(prevNll-thisNll);
            if (chi2<0. && order>1) chi2=0.;
            prob = TMath::Prob(chi2,order-prev_order); 

            cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
            cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
            prevNll=thisNll;
            cache_order=prev_order;
            cache_pdf=prev_pdf;

            // Calculate goodness of fit (will use toys for lowstats)
            double gofProb =0; 
            plot(mass,bkgPdf,data,Form("%s/%s%d_m_%s_cat_%s",outDir.c_str(),funcType->c_str(),order,mN_label.c_str(),catname.c_str()),category_label,fit_window_min,fit_window_max,mass_window_min,mass_window_max,nBinsForFit,nBinsForPlot,fitStatus,&gofProb);

            if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

              if (gofProb > minGofThreshold || order == truthOrder ) {  // Good looking fit or one of our regular truth functions
              //if (gofProb == -1 || gofProb > minGofThreshold) { // minimal requirement on the goodness of fit (in the case where the statistics is enough)

                std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
                  << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
                allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
                storedPdfs.add(*bkgPdf);
                pdforders.push_back(order);

                // Keep track but we shall redo this later
                if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
                  simplebestFitPdfIndex = storedPdfs.getSize()-1;
                  MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
                }
              }
            }
            prev_order=order;
            prev_pdf=bkgPdf;
            order++;
          }
        } // end while

        fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
        choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    } // end loop over families

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);

    plot(mass,pdfs,data,Form("%s/truths_m_%s_cat_%s",outDir.c_str(),mN_label.c_str(),catname.c_str()),category_label,cat,fit_window_min,fit_window_max,mass_window_min,mass_window_max,nBinsForPlot);

    if (saveMultiPdf){
      // Put selectedModels into a MultiPdf
      string catindexname;
      string catname;
      catindexname = "pdfindex_bhnl_m_" + mN_label + "_cat_" + category_label[cat].c_str();
      catname = Form("%s",category_label[cat].c_str());
      RooCategory catIndex(catindexname.c_str(),"c");
      RooMultiPdf *pdf = new RooMultiPdf("qcd_multipdf","all pdfs",catIndex,storedPdfs);
      RooRealVar nBackground("qcd_multipdf_norm","nbkg",data->sumEntries(),0,3*data->sumEntries()); //TODO do we want to keep this strategy?, make sure normalisation is correct
      //nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
      //double check the best pdf!
      int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,fit_range,!verbose);
      catIndex.setIndex(bestFitPdfIndex);
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
      std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
      storedPdfs.Print();
      std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
      std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;

      mass->setBins(nBinsForFit);
      //RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

      // Save it (also a binned version of the dataset
      outputws->import(*pdf);
      outputws->import(nBackground);
      outputws->import(catIndex);
      //outputws->import(dataBinned);
      outputws->import(*data);
      plot(mass,pdf,&catIndex,data,Form("%s/multipdf_m_%s_cat_%s",outDir.c_str(),mN_label.c_str(),catname.c_str()),category_label,cat,fit_window_min,fit_window_max,mass_window_min,mass_window_max,nBinsForPlot,fit_range,bestFitPdfIndex);

    } // end if saveMultiPdf

  } // end loop over categories 

  if (saveMultiPdf){
    outputfile->cd();
    outputws->Write();
    outputfile->Close();  
  }

  /* Save the results in a .dat file
  // Write recommended options to screen and to file
  //system(Form("mkdir -p %s",outDir.c_str()));
  //FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "[RESULT] Recommended options based on truth" << endl;

  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"\n");
  }


  cout << "[RESULT] Recommended options for envelope" << endl;
  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
      std::vector<int> ords = it->second;
      for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
        cout << "\t" << it->first << " - " << *ordit << endl;
        fprintf(dfile,"envel=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
      }
    }
    fprintf(dfile,"\n");
  }
  */

  inFile->Close();

  return 0;
}
