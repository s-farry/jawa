#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TMath.h>
#include <TObjArray.h>
#include <math.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TF1.h>
#include <Utils.h>
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TRandom3.h>
#include <TEntryList.h>
#include <TH1.h>
#include <TMinuit.h>
#include <TParameter.h>
#include <Tune2D.h>
#include <TCanvas.h>

using namespace std;

Tune2D::Tune2D() : Tune(){}
Tune2D::Tune2D(string name): Tune(name){}

Tune2D::Tune2D(string name, Tree* data, Tree* mc, Expr* tuneVar, Var2D* fVar, string cut) : Tune(name, data, mc, tuneVar, 0, cut){

  m_data_array = new TObjArray();
  m_mc_array = new TObjArray();
  m_mc_corr_array = new TObjArray();

  m_2DfVar = fVar;
  vector<double> binedges_x = Utils::GetBinEdgesX(fVar->GetHist());
  vector<double> binedges_y = Utils::GetBinEdgesY(fVar->GetHist());
  m_2Dres_sigma = new TH2F((name + "_sigma_2D").c_str(), (name+"_sigma").c_str(), binedges_x.size() - 1,
			   &binedges_x[0], binedges_y.size() - 1, &binedges_y[0]);
  m_2Dres_mean  = new TH2F((name + "_mean_2D").c_str() , (name+"_mean").c_str() , binedges_x.size() - 1,
			   &binedges_x[0], binedges_y.size() - 1, &binedges_y[0]);
  m_2Dstddev    = new TH2F((name + "_stddev_2D").c_str() , (name+"_stddev").c_str() , binedges_x.size() - 1,
			   &binedges_x[0], binedges_y.size() - 1, &binedges_y[0]);

  m_data_2Dstddev = vector< vector<double> >(binedges_x.size() - 1, vector<double>(binedges_y.size() - 1, 0.0));

}



void Tune2D::fill2DVals(){
  if (m_2DfVar){
    m_data_2Dvec = getVals(m_data, m_tuneVar, m_2DfVar, m_cut);
    m_mc_2Dvec   = getVals(m_mc  , m_tuneVar, m_2DfVar, m_cut);
  }
  for (unsigned int i = 0 ; i < m_data_2Dvec.size(); ++i){
    for (unsigned int j = 0 ; j < m_data_2Dvec[i].size(); ++j){
      double stddev = standard_deviation(m_data_2Dvec[i][j], 0.5);
      //cout<<"value: "<<m_data_2Dvec[i][j][5]<<endl;
      m_data_2Dstddev[i][j] = stddev;
      //cout<<m_data_2Dstddev[i][j]<<endl;
    }
  }
}

vector< vector< vector< pair<double, double> > > > Tune2D::getVals(Tree* tree, Expr* var, Var2D* binvar, TCut cut, vector<ReweightVar*> rwvars){
  TTree* t = tree->GetTTree();
  t->Draw(">>myList", cut , "entrylist");
  TEntryList* list = (TEntryList*)gDirectory->Get("myList");
  t->SetBranchStatus("*",0);
  tree->SetBranches(var->GetVarNames());
  tree->SetBranches(binvar->GetVar1()->GetExpr()->GetVarNames());
  tree->SetBranches(binvar->GetVar2()->GetExpr()->GetVarNames());

  //Set names for reweighted variables for data
  for (std::vector<ReweightVar*>::iterator iv = rwvars.begin(); iv!= rwvars.end();++iv){
    vector<Expr*> rwnames = (*iv)->GetExprs();
    for (unsigned int i = 0 ; i < rwnames.size() ; ++i) tree->SetBranches(rwnames.at(i)->GetVarNames());
  }
  
  Long64_t nentries = list->GetN();
  
  vector<double> bin_edges_x = binvar->GetVar1()->GetEdges();
  vector<double> bin_edges_y = binvar->GetVar2()->GetEdges();
  int xbins = bin_edges_x.size() - 1;
  int ybins = bin_edges_y.size() - 1;
  //int totbins = (bin_edges_x.size() - 1) * (bin_edges_y.size() - 1);

  vector< vector< vector< pair<double, double> > > > vals (xbins, vector< vector< pair<double, double> > >(ybins, vector< pair<double, double> >(0))); 
  
  for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
    if (jentry%10000==0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
    int entry = list->GetEntry(jentry);
    tree->GetEntry(entry);
    double xval = tree->GetVal(binvar->GetVar1()->GetExpr());
    double yval = tree->GetVal(binvar->GetVar2()->GetExpr());

    if ( xval >= bin_edges_x[0] && xval < bin_edges_x[bin_edges_x.size() - 1]
	 && yval >= bin_edges_y[0] && yval < bin_edges_y[bin_edges_y.size() - 1]){
      int xbin = binvar->GetVar1()->GetHist()->FindBin(xval);
      int ybin = binvar->GetVar2()->GetHist()->FindBin(yval);
      double val = tree->GetVal(var);
      double w = 1.0;
      if (rwvars.size() > 0) w = w * GetWeight(tree, rwvars);
      //cout<<"Filling : ("<<xbin<<","<<ybin<<") with "<<val<<" for "<<xval<<", "<<yval<<endl;
      //int bin = (ybin - 1)*(bin_edges_y.size() - 1) + (xbin - 1);
      vals[ xbin - 1 ][ ybin - 1 ].push_back(pair<double, double>(val, w));
    }

  }
  t->SetBranchStatus("*",1);
  return vals;
}

void Tune2D::tune(){

  fill2DVals();
  
  ROOT::Minuit2::Minuit2Minimizer fitter;
  ROOT::Math::Functor function(this, &Tune::metric, 2);
  fitter.SetLimitedVariable(0, "mean",  m_mean_init  , m_mean_step  , m_mean_lolimit   , m_mean_uplimit  );
  fitter.SetLimitedVariable(1, "sigma", m_sigma_init , m_sigma_step , m_sigma_lolimit  , m_sigma_uplimit );
  fitter.SetFunction(function);
  fitter.SetPrecision(m_precision);
  fitter.SetTolerance(m_tolerance);

  for (int i = 0; i < (int)m_data_2Dvec.size(); ++i){
    TObjArray* x_data_array    = new TObjArray();
    TObjArray* x_mc_array      = new TObjArray();
    TObjArray* x_mc_corr_array = new TObjArray();

    for (int j = 0; j < (int)m_data_2Dvec[i].size(); ++j){
      m_current_data   = &(m_data_2Dvec[i][j]);
      m_current_mc     = &(m_mc_2Dvec[i][j]);
      m_current_stddev = m_data_2Dstddev[i][j];

      double mean = 0.0, mean_err = 0.0, sigma = 0.0, sigma_err = 0.0;
      //double status = 0;

      if (m_current_data->size() > 0 && m_current_mc->size() > 0){
	fitter.Minimize();
	
	//status = fitter.Status();
	//if (status == 1 || status ==3){
	mean = fitter.X()[0];
	mean_err = fitter.Errors()[0];
	sigma = fitter.X()[1];
	sigma_err = fitter.Errors()[1];	
	//}
      }
      ostringstream data_name;
      data_name<<m_name<<"_data_"<<i<<"_"<<j;
      ostringstream mc_name;
      mc_name<<m_name<<"_mc_"<<i<<"_"<<j;
      ostringstream mc_corr_name;
      mc_corr_name<<m_name<<"_mc_corr_"<<i<<"_"<<j;
      
      vector< pair<double, double> > smear_mc = smear_vals(*m_current_mc, mean, sigma);
      TH1F* data    = new TH1F(data_name.str().c_str()   ,"data"     ,100, -m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      TH1F* mc      = new TH1F(mc_name.str().c_str()     ,"mc"      , 100, -m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      TH1F* mc_corr = new TH1F(mc_corr_name.str().c_str(),"mc_corr" , 100, -m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      
      m_2Dres_mean->SetBinContent( i+1, j+1, mean      );
      m_2Dres_mean->SetBinError(   i+1, j+1, mean_err  );
      m_2Dres_sigma->SetBinContent(i+1, j+1, sigma     );
      m_2Dres_sigma->SetBinError(  i+1, j+1, sigma_err );
      m_2Dstddev->SetBinContent( i + 1 , j+1, m_current_stddev);
      
      
      fillhist(data,*m_current_data);
      fillhist(mc, *m_current_mc);
      fillhist(mc_corr, smear_mc);
      
      x_data_array->Add(data);
      x_mc_array->Add(mc);
      x_mc_corr_array->Add(mc_corr);
      
      
    }
    m_data_array->Add(x_data_array);
    m_mc_array->Add(x_mc_array);
    m_mc_corr_array->Add(x_mc_corr_array);
  }
  
} 

void Tune2D::SaveToFile(){
  string fName = m_name + ".root";
  TFile f(fName.c_str(),"RECREATE");
  m_data_array->Write("data", 1);
  m_mc_array->Write("mc", 1);
  m_mc_corr_array->Write("mc_corr", 1);
  m_2Dres_mean->Write();
  m_2Dres_sigma->Write();
  m_2Dstddev->Write();
  f.Close();
  
}
