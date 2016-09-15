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
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TRandom3.h>
#include <TEntryList.h>
#include <TH1.h>
#include <TMinuit.h>
#include <TParameter.h>
#include <Tune.h>
#include <TCanvas.h>

using namespace std;

Tune::Tune() : JawaObj ("Tune"), m_r3(0) {
}
Tune::Tune(string name): JawaObj("Tune", name), m_r3(0){
  m_mean_init       =  0.00  ;
  m_mean_step       =  0.00001;
  m_mean_lolimit    = -0.1   ;
  m_mean_uplimit    =  0.1   ;

  m_sigma_init      =  0.00  ;
  m_sigma_step      =  0.00001;
  m_sigma_lolimit   =  0.0   ;
  m_sigma_uplimit   =  0.5   ;

  m_precision = 0.0001;
  m_tolerance = 0.0001;
}

Tune::Tune(string name, Tree* data, Tree* mc, Expr* tuneVar, Var* fVar, string cut) : JawaObj("Tune",name), m_r3(0) {
  m_data = data;
  m_mc = mc;
  m_tuneVar = tuneVar;
  m_cut = TCut(cut.c_str());
  m_sdfac = 3;

  m_mean_init       =  0.00  ;
  m_mean_step       =  0.0005;
  m_mean_lolimit    = -0.05  ;
  m_mean_uplimit    =  0.05  ;

  m_sigma_init      =  0.00  ;
  m_sigma_step      =  0.0005;
  m_sigma_lolimit   =  0.00  ;
  m_sigma_uplimit   =  0.05  ;

  m_precision = 0.0001;
  m_tolerance = 0.0001;

  m_data_array = new TObjArray();
  m_mc_array = new TObjArray();
  m_mc_corr_array = new TObjArray();

  if (fVar) {
    m_fVar = fVar;
    vector<double> binedges = fVar->GetEdges();
    m_res_sigma = new TH1F((name + "_sigma").c_str(), (name+"_sigma").c_str(), binedges.size() - 1, &binedges[0]);
    m_res_mean  = new TH1F((name + "_mean").c_str() , (name+"_mean").c_str() , binedges.size() - 1, &binedges[0]);
    m_stddev    = new TH1F((name + "_stddev").c_str() , (name+"_stddev").c_str() , binedges.size() - 1, &binedges[0]);
  }
  else{
    m_res_sigma = new TH1F((name + "_sigma").c_str(), (name+"_sigma").c_str(), 1, 0, 1);
    m_res_mean  = new TH1F((name + "_mean").c_str() , (name+"_mean").c_str() , 1, 0, 1);
    m_stddev    = new TH1F((name + "_stddev").c_str() , (name+"_stddev").c_str() , 1, 0, 1);
  }

}


vector<double> Tune::getCorrelatedRandoms(vector< vector<double> > corrs ){
  vector< vector<double> > chol = cholesky(corrs);
  vector<double> randoms(corrs.size(), 0);
  vector<double> corrrandoms(corrs.size(), 0);
  for (unsigned int i = 0; i < corrs.size(); ++i){
    //Get unit gaus random number
    randoms[i] = m_r3.Gaus(0, 1);
  }
  
  //Correlate them
  for (unsigned int i = 0; i < corrs.size(); ++i){
    for (unsigned int j = 0; j < corrs.size(); ++j){
      corrrandoms[i] += chol[i][j] * randoms[j];
    }
  }
  return corrrandoms;
}
/*
vector<double> Tune::getRandoms(vector< vector<double> > corrs ){
  vector<double> randoms(corrs.size(), 0);
  for (unsigned int i = 0; i < corrs.size(); ++i){
    //Get unit gaus random number
    randoms[i] = m_r3.Gaus(0, 1);
  }
  return randoms;
}
*/

vector< pair<double, double> > Tune::smear_vals(vector< pair<double, double> > & vals, double mean, double sigma){
    vector<pair<double, double > > smeared_vals;
    for (vector< pair<double, double > >::iterator iv = vals.begin(); iv != vals.end(); ++iv){
      double sm = m_r3.Gaus(mean,sigma);
      smeared_vals.push_back(pair<double, double>((*iv).first+sm, (*iv).second));
    }
    return smeared_vals;
}

void Tune::fillhist(TH1F* h, vector< pair<double, double> >& vals){
  for (unsigned int i = 0; i < vals.size(); ++i){
    h->Fill(vals.at(i).first, vals.at(i).second);
  }
  return;
}
double Tune::get_mean(vector< pair<double, double> >& vals){
  double mean = 0.0;
  double size = 0.0;
  for (unsigned int i = 0; i < vals.size(); ++i){
    mean +=vals[i].first*vals[i].second;
    size += vals[i].second;
  }
  mean = mean/vals.size();
  return mean;
}
double Tune::standard_deviation(vector< pair<double, double> >& vals, double max){
  double mean = get_mean(vals);
  double sum_deviation = 0.0;
  double size = 0.0;
  for (unsigned int i = 0; i < vals.size(); ++i){
    if (abs(vals[i].first) < max || max == -1) {
      sum_deviation+= ((vals[i].first - mean)*(vals[i].first - mean))*vals[i].second;
      size+= vals[i].second;
    }
    
  }
  return sqrt(sum_deviation/size);
}

double Tune::metric(const double *params){
  double m = params[0];
  double s = params[1];
  return getchi2(m, s);
}

double Tune::getchi2(double mean, double sigma){
  //double std_dev = standard_deviation(*m_data);
  
  TH1F* data = new TH1F("data","data",100,-(m_current_stddev) * m_sdfac, (m_current_stddev) * m_sdfac);
  TH1F* mc   = new TH1F("mc","mc", 100, -(m_current_stddev) * m_sdfac, (m_current_stddev) * m_sdfac);
  
  vector< pair<double, double> > smear_mc = Tune::smear_vals((*m_current_mc), mean, sigma);

  Tune::fillhist(data,*m_current_data);
  Tune::fillhist(mc, smear_mc);

  double chi2 = data->Chi2Test(mc, "CHI2/NDF");

  //TCanvas c1("c1");
  //data->DrawNormalized();
  //mc->SetLineColor(kRed);
  //mc->DrawNormalized("same");
  //c1.Print("test.pdf");

  data->Delete();
  mc->Delete();
  return chi2;
}
/*
void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  //Function for TMinuit to minimise
  double chisq = 0;
  chisq = getchi2(par[0], par[1]);
  f = chisq;
  return;
  }*/

void Tune::SetVals(vector< vector< pair<double, double> > > data, vector< vector< pair<double, double> > > mc){
  m_data_vec = data;
  m_mc_vec   = mc;
  for (unsigned int i = 0 ; i < m_data_vec.size(); ++i){
    double stddev = standard_deviation(m_data_vec[i], 0.5);
    m_data_stddev.push_back(stddev);
  }
}

vector< vector< pair<double, double> > > Tune::GetMCVals(){
  return m_mc_vec;
}
vector< vector< pair<double, double> > > Tune::GetDataVals(){
  return m_data_vec;
}

void Tune::fillVals(){
  //TTree* datat = m_data->GetTTree();
  //TTree* mct   = m_mc->GetTTree();
  if (m_fVar){
    m_data_vec = getVals(m_data, m_tuneVar, m_fVar, m_cut, m_data_rwvars);
    m_mc_vec   = getVals(m_mc  , m_tuneVar, m_fVar, m_cut, m_mc_rwvars);
  }
  else{
    vector< pair<double, double> > data_vec = getVals(m_data, m_tuneVar);
    vector< pair<double, double> > mc_vec   = getVals(m_mc  , m_tuneVar);
    m_data_vec = vector< vector< pair<double, double> > >(1, data_vec);
    m_mc_vec = vector< vector< pair<double, double> > >(1, mc_vec);
  }
  for (unsigned int i = 0 ; i < m_data_vec.size(); ++i){
    double stddev = standard_deviation(m_data_vec[i], 0.5);
    m_data_stddev.push_back(stddev);
  }
}
/*
double* Tune::cholesky(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
          float s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }

    return L;
    }*/

vector< vector<double> > Tune::cholesky(vector< vector<double> > A) {
  //double *L = (double*)calloc(n * n, sizeof(double));
  int n = A.size();
  vector< vector<double> > L( n, vector<double>(n, 0.0) );

    for (int i = 0; i < n; i++){
      for (int j = 0; j < (i+1); j++) {
	double s = 0;
	for (int k = 0; k < j; k++){
	  s += L[i][k] * L[j][k];
	}
	  L[i][j] = (i == j) ?
	    sqrt(A[i][i] - s) :
	    (1.0 / L[j][j] * (A[i][j] - s));
	  //}
      }
    }
    //printMatrix(L);

    return L;
}

void Tune::printMatrix(vector< vector<double> > A){
  int N = A.size();
  cout<<"---------------------"<<endl;
  for (int i = 0 ; i < N; ++i){
    for (int j = 0 ; j < N; ++j){
      cout<<A[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<"---------------------"<<endl;
}


vector< pair<double, double> > Tune::getVals(Tree* t, Expr* var, TCut cut){
  t->GetTTree()->Draw(">>myList", cut , "entrylist");
  TEntryList* list = (TEntryList*)gDirectory->Get("myList");
  t->GetTTree()->SetBranchStatus("*",0);
  t->SetBranches(var->GetVarNames());


  //t->SetBranchStatus(var.c_str(), 1);
  //double d_var;
  //t->SetBranchAddress(var.c_str(), &d_var);
 
  vector< pair<double, double> > vals; 

  Long64_t nentries = list->GetN();

  for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
    if (jentry%10000==0) info()<<"Entry "<<jentry<<" of "<<nentries<<endl;
      int entry = list->GetEntry(jentry);
      t->GetEntry(entry);
      double w = 1.0;
      vals.push_back(pair<double, double>(t->GetVal(var), w));
  }
  t->GetTTree()->SetBranchStatus("*",1);
  return vals;
}

vector< vector< pair<double, double> > > Tune::getVals(Tree* tree, Expr* var, Var* binvar, TCut cut, vector<ReweightVar*> rwvars){
  TTree* t = tree->GetTTree();
  t->Draw(">>myList", cut , "entrylist");
  TEntryList* list = (TEntryList*)gDirectory->Get("myList");
  t->SetBranchStatus("*",0);
  tree->SetBranches(var->GetVarNames());
  tree->SetBranches(binvar->GetExpr()->GetVarNames());
  
  //Set names for reweighted variables
  for (std::vector<ReweightVar*>::iterator iv = rwvars.begin(); iv!= rwvars.end();++iv){
    vector<Expr*> rwnames = (*iv)->GetExprs();
    for (unsigned int i = 0 ; i < rwnames.size(); ++i) tree->SetBranches(rwnames.at(i)->GetVarNames());
  }
  

  Long64_t nentries = list->GetN();
  
  vector<double> bin_edges = binvar->GetEdges();
  vector< vector< pair<double, double> > > vals (bin_edges.size() - 1, vector< pair<double, double> >(0)); 
  
  for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
    if (jentry%10000==0) info()<<"Entry "<<jentry<<" of "<<nentries<<endl;
    int entry = list->GetEntry(jentry);
    tree->GetEntry(entry);
    double binval = tree->GetVal(binvar->GetExpr());
    if (binval >= bin_edges[0] && binval < bin_edges[bin_edges.size() - 1]){
      int bin = binvar->GetHist()->FindBin(binval);
      double w = 1.0;
      if (rwvars.size() > 0) w = w * GetWeight(tree, rwvars);
      double val = tree->GetVal(var);
      vals[bin - 1].push_back(pair<double, double>(val, w));
    }
  }
  t->SetBranchStatus("*",1);
  return vals;
}

void Tune::tune(){
  fillVals();
  
  ROOT::Minuit2::Minuit2Minimizer fitter;
  ROOT::Math::Functor function(this, &Tune::metric, 2);

  fitter.SetLimitedVariable(0, "mean",  m_mean_init  , m_mean_step  , m_mean_lolimit   , m_mean_uplimit  );
  fitter.SetLimitedVariable(1, "sigma", m_sigma_init , m_sigma_step , m_sigma_lolimit  , m_sigma_uplimit );
  fitter.SetFunction(function);
  fitter.SetPrecision(m_precision);
  fitter.SetTolerance(m_tolerance);

  //fitter.SetLimitedVariable(0, "mean", 0.00, 0.00001, -0.5, 0.5);
  //fitter.SetLimitedVariable(1, "sigma", 0.00, 0.00001, 0, 0.2);
  //fitter.SetLimitedVariable(0, "mean", 0.0, 0.0005, -0.05, 0.05);
  //fitter.SetLimitedVariable(1, "sigma", 0.0, 0.0005, 0, 0.05);
  //fitter.SetVariableValue(0, 0.0);
  //fitter.SetVariableValue(1, 0.0);
  //fitter.SetFunction(function);
  //fitter.SetPrecision(0.0001);
  //fitter.SetTolerance(0.0001);

  for (int i = 0; i < (int)m_data_vec.size(); ++i){
      m_current_data   = &(m_data_vec[i]);
      m_current_mc     = &(m_mc_vec[i]);
      m_current_stddev = m_data_stddev[i];
      
      fitter.Minimize();
      
      double mean = 0.0, mean_err = 0.0, sigma = 0.0, sigma_err = 0.0;
      //double status = 0;

      mean = fitter.X()[0];
      mean_err = fitter.Errors()[0];
      sigma = fitter.X()[1];
      sigma_err = fitter.Errors()[1];
      //status = fitter.Status();	
      
      ostringstream data_name;
      data_name<<m_name<<"_data_"<<i;
      ostringstream mc_name;
      mc_name<<m_name<<"_mc_"<<i;
      ostringstream mc_corr_name;
      mc_corr_name<<m_name<<"_mc_corr_"<<i;

      vector< pair<double, double> > smear_mc = smear_vals(*m_current_mc, mean, sigma);
      TH1F* data    = new TH1F(data_name.str().c_str(),"data",100,-m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      TH1F* mc      = new TH1F(mc_name.str().c_str()  ,"mc", 100, -m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      TH1F* mc_corr = new TH1F(mc_corr_name.str().c_str(),"mc_corr", 100, -m_current_stddev*m_sdfac, m_current_stddev*m_sdfac);
      
      m_res_mean->SetBinContent( i+1, mean      );
      m_res_mean->SetBinError(   i+1, mean_err  );
      m_res_sigma->SetBinContent(i+1, sigma     );
      m_res_sigma->SetBinError(  i+1, sigma_err );
      m_stddev->SetBinContent( i + 1 , m_current_stddev);
      
        
      fillhist(data,*m_current_data);
      fillhist(mc, *m_current_mc);
      fillhist(mc_corr, smear_mc);
      
      m_data_array->Add(data);
      m_mc_array->Add(mc);
      m_mc_corr_array->Add(mc_corr);
      

    }
	
} 
void Tune::SaveToFile(){
  
  string fName = m_name + ".root";
  TFile f(fName.c_str(),"RECREATE");
  m_data_array->Write("data", 1);
  m_mc_array->Write("mc", 1);
  m_mc_corr_array->Write("mc_corr", 1);
  m_res_mean->Write();
  m_res_sigma->Write();
  m_stddev->Write();
  f.Close();
  
}

void Tune::SetSigmaPars(double init, double step, double lolimit, double uplimit){
  m_sigma_init = init;
  m_sigma_step = step;
  m_sigma_lolimit = lolimit;
  m_sigma_uplimit = uplimit;
}
void Tune::SetMeanPars(double init, double step, double lolimit, double uplimit){
  m_mean_init = init;
  m_mean_step = step;
  m_mean_lolimit = lolimit;
  m_mean_uplimit = uplimit;
}
void Tune::SetPrecision(double precision){
  m_precision = precision;
}
void Tune::SetTolerance(double tolerance){
  m_tolerance = tolerance;
}

void Tune::PrintSigmaPars(){
  cout<<m_sigma_init<<" "<<m_sigma_step<<" "<<m_sigma_lolimit<<" "<<m_sigma_uplimit<<endl;
}
void Tune::PrintMeanPars(){
  cout<<m_mean_init<<" "<<m_mean_step<<" "<<m_mean_lolimit<<" "<<m_mean_uplimit<<endl;
}

//For exposing to python

boost::python::list Tune::GetStdDevs_py(){
  boost::python::list ns;
  for (unsigned int i = 0 ; i < m_data_stddev.size() ; ++i){
    ns.append(m_data_stddev[i]);
  }
  return ns;
}
boost::python::list Tune::getCorrelatedRandoms_py(boost::python::list& ns){
  vector< vector<double> > corrmat(len(ns), vector<double>(len(ns), 0.0));
  for (unsigned int i = 0; i < len(ns); ++i){
    for (unsigned int j = 0; j < len(ns[i]); ++j){
    double d = boost::python::extract<double>(ns[i][j]);
    corrmat[i][j] = d;
    }
  }

  vector<double> rs = getCorrelatedRandoms(corrmat);

  boost::python::list ns2;
  for (unsigned int i = 0 ; i < rs.size() ; ++i){
    ns2.append(rs[i]);
  }
  return ns2;
}
boost::python::list Tune::cholesky_py(boost::python::list& ns){
  vector< vector<double> > corrmat(len(ns), vector<double>(len(ns), 0.0));
  for (unsigned int i = 0; i < len(ns); ++i){
    for (unsigned int j = 0; j < len(ns[i]); ++j){
    double d = boost::python::extract<double>(ns[i][j]);
    corrmat[i][j] = d;
    }
  }

  vector< vector<double> > chol = cholesky(corrmat);

  boost::python::list ns2;
  for (unsigned int i = 0 ; i < chol.size() ; ++i){
    boost::python::list ns3;
    for (unsigned int j = 0 ; j < chol.size() ; ++j){
      ns3.append(chol[i][j]);
    }
    ns2.append(ns3);
  }
  return ns2;
}

boost::python::list Tune::GetDataVec(int j){
  boost::python::list ns;
  for (unsigned int i = 0 ; i < m_data_vec[j].size() ; ++i){
    ns.append(m_data_vec[j][i].first);
  }
  return ns;
}

double Tune::standard_deviation_py(boost::python::list& ns){
  vector< pair<double, double> > vals;
  for (unsigned int i = 0; i < len(ns); ++i){
    double d = boost::python::extract<double>(ns[i]);
    vals.push_back(pair<double, double> (d, 1.0));
  }
  return standard_deviation(vals);
}

double Tune::GetWeight(Tree* tree, vector<ReweightVar*> rwvars){
  double w = 1.0;
  for (std::vector<ReweightVar*>::iterator iv = rwvars.begin(); iv!= rwvars.end();++iv){
    if ((*iv)->GetExprs().size() == 1){
      Expr* var = (*iv)->GetExpr();
      double val = tree->GetVal(var);
      w = w * ((*iv)->GetWeight(val));
    }
    else if ((*iv)->GetExprs().size() == 2){
      double val1 = tree->GetVal((*iv)->GetExprs().at(0));
      double val2 = tree->GetVal((*iv)->GetExprs().at(1));
      w = w * ((*iv)->GetWeight(val1, val2));
    }
  }
  return w;
}

void Tune::ReweightMC(ReweightVar* rwvar){
  m_mc_rwvars.push_back(rwvar);
}
void Tune::ReweightData(ReweightVar* rwvar){
  m_data_rwvars.push_back(rwvar);
}
void Tune::SetSDFactor(double s){
  m_sdfac = s;
}
