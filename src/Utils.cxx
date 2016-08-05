#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TMath.h>
#include <TObjArray.h>
#include <math.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TF1.h>
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TRandom3.h>
#include <TEntryList.h>
#include <TH1.h>
#include <TMinuit.h>
#include <TParameter.h>
#include <Utils.h>
#include <TCanvas.h>
#ifdef WITHPYTHON
#include <boost/python.hpp>
#include <boost/python/handle.hpp>
#endif

using namespace std;

namespace Utils{
  TH1F* tgraph2hist(string name, TGraphAsymmErrors* graph){
    double x = 0.0, y= 0.0;
    int N = graph->GetN();
    graph->GetPoint(0, x, y);
    double low = x - graph->GetErrorXlow(0);
    graph->GetPoint(N - 1, x, y);
    double hi  = x + graph->GetErrorXhigh(N-1);
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), N, low, hi);
    for (int i = 0; i < N; ++i){
      graph->GetPoint(i, x, y);
      hist->SetBinContent(i+1,y);
    }
    return hist;
  }

  TH1D* geteff(string name, TH3F* data, TGraphAsymmErrors* eff, int varybin, string form){
    TH1F* eff_hist = tgraph2hist("eff_hist", eff);
    Expr* f = new Expr(form);
    ostringstream title;
    title<<data->GetTitle()<<"_"<<name<<"_eff_corr_"<<varybin;
    ostringstream efftitle;
    efftitle<<data->GetTitle()<<"_"<<name<<"_eff_"<<varybin;
    cout<<efftitle.str()<<endl;
    cout<<"A: "<<eff_hist->GetBinContent(1)<<" "<<eff_hist->GetBinContent(2)<<" "<<eff_hist->GetBinContent(3)<<" "<<
      eff_hist->GetBinContent(4)<<" "<<eff_hist->GetBinContent(5)<<" "<<eff_hist->GetBinContent(6)<<" "<<
      eff_hist->GetBinContent(7)<<" "<<eff_hist->GetBinContent(8)<<" "<<eff_hist->GetBinContent(9)<<" "<<
      eff_hist->GetBinContent(10)<<endl;
    TH3F* eff_corr = (TH3F*)data->Clone(title.str().c_str());

    if (varybin == 999 ){
      for (int i = 0 ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) + eff->GetErrorYhigh(i));
      }
    }
    else if (varybin == -999){
      for (int i = 0 ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) - eff->GetErrorYlow(i));
      }
    }
    else if (varybin != 0 && varybin <= eff_hist->GetNbinsX()){
      double err = varybin > 0 ? eff->GetErrorYhigh(abs(varybin)-1) : -1*eff->GetErrorYlow(abs(varybin)-1);
      //double errlo = eff->GetErrorYlow(varybin);
      //double err = errhi > errlo ? errhi : -errlo;
      //cout<<"VARYBIN: "<<varybin<<" ERR: "<<err< " "<<eff->GetErrorYhigh(abs(varybin-1))<<" "<<-1*eff->GetErrorYlow(abs(varybin-1))<<endl;
      eff_hist->SetBinContent(abs(varybin), eff_hist->GetBinContent(abs(varybin)) + err);
    }
    cout<<"B: "<<eff_hist->GetBinContent(1)<<" "<<eff_hist->GetBinContent(2)<<" "<<eff_hist->GetBinContent(3)<<" "<<
      eff_hist->GetBinContent(4)<<" "<<eff_hist->GetBinContent(5)<<" "<<eff_hist->GetBinContent(6)<<" "<<
      eff_hist->GetBinContent(7)<<" "<<eff_hist->GetBinContent(8)<<" "<<eff_hist->GetBinContent(9)<<" "<<
      eff_hist->GetBinContent(10)<<endl;
    
    //note z axis should be final variable and x, y should be in same bins as graph

    for (int i = 0 ; i < data->GetXaxis()->GetNbins() ; ++i){
      for (int j = 0 ; j < data->GetYaxis()->GetNbins() ; ++j){
	for (int k = 0 ; k < data->GetZaxis()->GetNbins() ; ++k){
	  vector<double> effs;
	  effs.push_back(eff_hist->GetBinContent(i+1));
	  effs.push_back(eff_hist->GetBinContent(j+1));
	  if (f->GetVal(effs) < 0) cout<<"WHAT THE HELL "<<f->GetVal(effs)<<" "<<eff_hist->GetBinContent(i+1)<<" "<<
				     eff_hist->GetBinContent(j+1)<<" "<<i+1<<" "<<j+1<<" "<<k+1<<endl;
	  eff_corr->SetBinContent(i+1,j+1,k+1, eff_corr->GetBinContent(i+1,j+1,k+1) / f->GetVal(effs));	  
	}
      }
    }
    
    TH1D* final_eff = data->ProjectionZ(efftitle.str().c_str());
    TH1D* denom     = eff_corr->ProjectionZ();
    final_eff->Sumw2();
    denom->Sumw2();
    final_eff->Divide(denom);

    denom->Delete();
    eff_hist->Delete();
    eff_corr->Delete();

    return final_eff;
      
  }

  TH1D* geteff(string name, TH2F* data, TGraphAsymmErrors* eff, int varybin, string form){
    TH1F* eff_hist = tgraph2hist("eff_hist", eff);
    Expr* f = new Expr(form);
    ostringstream title;
    title<<data->GetTitle()<<"_"<<name<<"_eff_corr_"<<varybin;
    ostringstream efftitle;
    efftitle<<data->GetTitle()<<"_"<<name<<"_eff_"<<varybin;
    TH2F* eff_corr = (TH2F*)data->Clone(title.str().c_str());


    if (varybin == 999 ){
      for (int i = 0 ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) + eff->GetErrorYhigh(i));
      }
    }
    else if (varybin == -999){
      for (int i = 0 ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) - eff->GetErrorYlow(i));
      }
    }
    else if (varybin != 0 && varybin <= eff_hist->GetNbinsX()){
      double err = varybin > 0 ? eff->GetErrorYhigh(abs(varybin-1)) : -1*eff->GetErrorYlow(abs(varybin-1));
      //double errlo = eff->GetErrorYlow(varybin);
      //double err = errhi > errlo ? errhi : -errlo;
      eff_hist->SetBinContent(varybin, eff_hist->GetBinContent(varybin) + err);
    }


    //note y axis should be final variable and x should be in same bins as graph

    for (int i = 0 ; i < data->GetXaxis()->GetNbins() ; ++i){
      for (int j = 0 ; j < data->GetYaxis()->GetNbins() ; ++j){
	  vector<double> effs;
	  effs.push_back(eff_hist->GetBinContent(i+1));
	  eff_corr->SetBinContent(i+1,j+1, eff_corr->GetBinContent(i+1,j+1) / f->GetVal(effs));	  
	}
      }
    TH1D* final_eff = data->ProjectionY(efftitle.str().c_str());
    TH1D* denom     = eff_corr->ProjectionY();
    final_eff->Divide(denom);

    denom->Delete();
    eff_hist->Delete();
    eff_corr->Delete();

    return final_eff;
      
  }

  TH1D* geteff(string name, TH1F* data, TGraphAsymmErrors* eff, int varybin, string form){
    TH1F* eff_hist = tgraph2hist("eff_hist", eff);
    Expr* f = new Expr(form);
    ostringstream title;
    title<<data->GetTitle()<<"_"<<name<<"_eff_corr_"<<varybin;
    ostringstream efftitle;
    efftitle<<data->GetTitle()<<"_"<<name<<"_eff_"<<varybin;
    TH2F* eff_corr = (TH2F*)data->Clone(title.str().c_str());
    

    //correlated uncertainty - increase all efficiencies up by their uncertainties
    if (varybin == 999 ){
      for (int i = 0  ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) + eff->GetErrorYhigh(i));
      }
    }
    //correlated uncertainty - decrease all efficiencies by their uncertainties
    else if (varybin == -999){
      for (int i = 0  ; i < eff_hist->GetNbinsX() ; ++i){
	eff_hist->SetBinContent(i+1, eff_hist->GetBinContent(i+1) - eff->GetErrorYlow(i));
      }
    }
    else if (varybin != 0 && varybin <= eff_hist->GetNbinsX()){
      double err = varybin > 0 ? eff->GetErrorYhigh(abs(varybin-1)) : -1*eff->GetErrorYlow(abs(varybin-1));
      //double errlo = eff->GetErrorYlow(varybin);
      //double err = errhi > errlo ? errhi : -errlo;
      eff_hist->SetBinContent(varybin, eff_hist->GetBinContent(varybin) + err);
    }

    //note y axis should be final variable and x should be in same bins as graph

    for (int i = 0 ; i < data->GetXaxis()->GetNbins() ; ++i){
	  vector<double> effs;
	  effs.push_back(eff_hist->GetBinContent(i+1));
	  eff_corr->SetBinContent(i+1, eff_corr->GetBinContent(i+1) / f->GetVal(effs));	  
    }
    TH1D* final_eff = (TH1D*)data->Clone(efftitle.str().c_str());
    TH1D* denom     = (TH1D*)eff_corr->Clone("denom");
    final_eff->Divide(denom);
    
    denom->Delete();
    eff_hist->Delete();
    eff_corr->Delete();
    
    return final_eff;
    
  }

  pair< TH1D*, vector<TH1D*> > getEffVariations(string name, TH3F* data, TGraphAsymmErrors* eff, string form, bool correlated){
    TH1D* central = geteff(name, data, eff, 0, form);
    vector<TH1D*> errors;

    if (correlated) {
      ostringstream title;
      title<<name<<"_"<<999;
      errors.push_back(geteff(title.str(), data, eff, 999, form));
      errors.push_back(geteff(title.str(), data, eff, -999, form));
    }
    else{
      for (int i = 0 ; i < eff->GetN(); ++i){
	ostringstream title;
	title<<name<<"_"<<i;
	errors.push_back(geteff(title.str(), data, eff, i+1, form));
	errors.push_back(geteff(title.str(), data, eff, -1*(i+1), form));
      }
    }
    return std::pair< TH1D*, vector<TH1D*> >(central, errors);
  }

  pair< TH1D*, vector<TH1D*> > getEffVariations(string name, TH2F* data, TGraphAsymmErrors* eff, string form, bool correlated){
    TH1D* central = geteff(name, data, eff, 0, form);
    vector<TH1D*> errors;
    if (correlated) {
      ostringstream title;
      title<<name<<"_"<<999;
      errors.push_back(geteff(title.str(), data, eff, 999, form));
      errors.push_back(geteff(title.str(), data, eff, -999, form));
    }
    else{
      for (int i = 0 ; i < eff->GetN(); ++i){
	ostringstream title;
	title<<name<<"_"<<i;
	errors.push_back(geteff(title.str(), data, eff, i+1, form));
	errors.push_back(geteff(title.str(), data, eff, -1*(i+1), form));
      }
    }
    return pair< TH1D*, vector<TH1D*> >(central, errors);
  }

  pair< TH1D*, vector<TH1D*> > getEffVariations(string name, TH1F* data, TGraphAsymmErrors* eff, string form, bool correlated){
    TH1D* central = geteff(name, data, eff, 0, form);
    vector<TH1D*> errors;
    if (correlated) {
      ostringstream title;
      title<<name<<"_"<<999;
      errors.push_back(geteff(title.str(), data, eff, 999, form));
      errors.push_back(geteff(title.str(), data, eff, -999, form));
    }
    else{
      for (int i = 0 ; i < eff->GetN(); ++i){
	ostringstream title;
	title<<name<<"_"<<i;
	errors.push_back(geteff(title.str(), data, eff, i+1, form));
	errors.push_back(geteff(title.str(), data, eff, -1*(i+1), form));
      }
    }
    return pair< TH1D*, vector<TH1D*> >(central, errors);
  }

  matrix getEffErrMatrix( pair < TH1D*, vector<TH1D* > > p ){
    //Get matrix of % deviations by varying each bin
    TH1D* central = p.first;
    vector<TH1D*>::iterator ih;
    int bins = central->GetNbinsX();
    matrix errs;//(p.second.size(), vector<double>(bins, 0.0));
    //vector<double> total(bins, 0.0);
    for (int j = 0; j < (int)p.second.size()/2; ++j){
      vector<double> a;
      a.reserve(bins);
      for (int i = 0 ; i < bins ; ++i){
	double err = max(abs(p.second.at(j*2)->GetBinContent(i+1) - central->GetBinContent(i+1)),
			 abs(p.second.at(j*2+1)->GetBinContent(i+1) - central->GetBinContent(i+1)))
	  /central->GetBinContent(i + 1);
	a.push_back(err);
	if (err > 0.1) cout<<err<<" "<<central->GetBinContent(i+1)<<" "<<p.second.at(j*2)->GetBinContent(i+1)<<
			 " "<<p.second.at(j*2+1)->GetBinContent(i+1)<<" "<<i<<" "<<j<<endl;
	//total.at(i) =total.at(i) + err*err;
      }
      errs.push_back(a);
    }
    //for (int k = 0 ; k < bins ; ++k){
    //  total.at(k) = sqrt(total.at(k));
    //}

    //return make_pair<vector<double>, vector<vector<double> > >(total, errs);
    return errs;
    
  }

  void saveMatrix(string name, matrix A){
    ofstream ofile;
    ofile.open(name.c_str());
    for ( unsigned int i = 0 ; i < A.size() ; ++i ){
          for ( unsigned int j = 0 ; j < A[i].size() ; ++j ){
	    if ( j != 0 ) ofile<<",";
	    ofile<<A[i][j];
	  }
	  ofile<<"\n";
    }
    cout<<"wrote matrix to "<<name<<endl;
  }

  void saveTH1F(string name, TH1F* h){
    ofstream ofile;
    ofile.open(name.c_str());
    for ( int i = 0 ; i < h->GetNbinsX() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<h->GetBinContent(i+1);
    }
    ofile<<"\n";
    for ( int i = 0 ; i < h->GetNbinsX() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<h->GetBinError(i+1);
    }
    ofile<<"\n";
    cout<<"wrote th1f to "<<name<<endl;
  }
  void saveTGraph(string name, TGraph* g){
    ofstream ofile;
    ofile.open(name.c_str());
    for ( int i = 0 ; i < g->GetN() ; ++i ){
      if ( i != 0 ) ofile<<",";
      double x, y;
      g->GetPoint(i, x, y);
      ofile<<y;
    }
    ofile<<"\n";
    for ( int i = 0 ; i < g->GetN() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<g->GetErrorYhigh(i);
    }
    ofile<<"\n";
    for ( int i = 0 ; i < g->GetN() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<g->GetErrorYlow(i);
    }
    cout<<"wrote tgraph to "<<name<<endl;
  }
  void saveTGraphErrs(string name, TGraph* g){
    ofstream ofile;
    ofile.open(name.c_str());
    for ( int i = 0 ; i < g->GetN() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<g->GetErrorYhigh(i);
    }
    ofile<<"\n";
    for ( int i = 0 ; i < g->GetN() ; ++i ){
      if ( i != 0 ) ofile<<",";
      ofile<<g->GetErrorYlow(i);
    }
    ofile<<"\n";
    cout<<"wrote tgraph to "<<name<<endl;
  }

  


  vector<double> getCorrelatedRandoms(TRandom3* r3, vector< vector<double> > corrs ){
    vector< vector<double> > chol = cholesky(corrs);
    vector<double> randoms(corrs.size(), 0);
    vector<double> corrrandoms(corrs.size(), 0);
    for (unsigned int i = 0; i < corrs.size(); ++i){
      //Get unit gaus random number
      randoms[i] = r3->Gaus(0, 1);
    }
    
    //Correlate them
    for (unsigned int i = 0; i < corrs.size(); ++i){
      for (unsigned int j = 0; j < corrs.size(); ++j){
	corrrandoms[i] += chol[i][j] * randoms[j];
      }
    }
    return corrrandoms;
  }

  vector<double> getRandoms(TRandom3* r3, int n ){
    vector<double> randoms(n, 0);
    for (int i = 0; i < n; ++i){
      //Get unit gaus random number
      randoms[i] = r3->Gaus(0, 1);
    }
    return randoms;
  }
  

  void fillhist(TH1F* h, vector<double>& vals){
    for (unsigned int i = 0; i < vals.size(); ++i){
      h->Fill(vals.at(i));
    }
    return;
  }
  double get_mean(vector<double>& vals){
    double mean = 0.0;
    for (unsigned int i = 0; i < vals.size(); ++i){
      mean+=vals[i];
    }
    mean = mean/vals.size();
    return mean;
  }
  double standard_deviation(vector<double>& vals, double max){
    double mean = get_mean(vals);
    double sum_deviation = 0.0;
    for (unsigned int i = 0; i < vals.size(); ++i){
      if (abs(vals[i]) < max || max == -1) sum_deviation+= ((vals[i] - mean)*(vals[i] - mean));
    }
    return sqrt(sum_deviation/vals.size());
  }
  
  vector< vector<double> > cholesky(vector< vector<double> > A) {
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
  
  void printMatrix(vector< vector<double> > A){
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
  
  
  /*
  boost::python::list GetStdDevs_py(){
    boost::python::list ns;
    for (unsigned int i = 0 ; i < m_data_stddev.size() ; ++i){
      ns.append(m_data_stddev[i]);
    }
    return ns;
  }
  boost::python::list getCorrelatedRandoms_py(boost::python::list& ns){
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
    }*/

  vector<double> getColumn(int n, vector<vector<double> > vals){
    vector<double> col;
    col.reserve(vals.size());
    for (unsigned int i = 0 ; i < vals.size() ; ++i){
      col.push_back(vals[i][n]);
    }
    return col;
  }

  vector< vector<double> > getCorrelationMatrix(vector<vector<double> > vals){
    //nentries x nvals

    unsigned int nentries = vals.size();
    unsigned int nvals    = vals[0].size();
    vector<double> means;
    vector<double> stddevs;
    means.reserve(vals.size());
    stddevs.reserve(vals.size());

    for ( unsigned int i = 0 ; i < vals.size() ; ++ i ){
      vector<double> col = getColumn(i, vals);
      double mean = GetMean(col);
      double stddev = GetStdDev(col, mean);
      means.push_back(mean);
      stddevs.push_back(stddev);
    }
    
    vector< vector<double> > covs( nvals, vector<double>(nentries, 0.0) );

    for ( unsigned int n = 0; n < nentries ; ++n){
      for (unsigned int i = 0 ; i < nvals ; ++ i){
	for ( unsigned int j = 0 ; j < nvals ; ++ j ){
	  double xval  = vals[n][i];
	  double yval  = vals[n][j];
	  double xmean = means[i];
	  double ymean = means[j];
	  covs[i][j] += (((xval - xmean) * (yval - ymean)))/(nentries * stddevs[i] * stddevs[j]);
	}
      }
    }
    return covs;
  }

  double GetMean(vector<double> vals){
    double mean = 0;
    for (unsigned int i = 0; i < vals.size() ; ++i){
      mean += vals[i];
    }
    mean = mean/vals.size();
    return mean;
  }
  double GetStdDev(vector<double> vals, double mean){
    double sum_deviation = 0;
    for (unsigned int i = 0; i < vals.size() ; ++i){
      sum_deviation+= ((vals[i] - mean)*(vals[i] - mean));
    }
    return sqrt(sum_deviation/vals.size());
  }
  std::vector<double> GetBinEdgesX(TH2F* hist){
    int nbins = hist->GetXaxis()->GetNbins();
    std::vector<double> edges;
    edges.reserve(nbins+1);
    for (int i = 0; i <nbins+1 ; ++i){
      edges.push_back(hist->GetXaxis()->GetBinLowEdge(i+1));
    }
    return edges;
  }
  std::vector<double> GetBinEdgesY(TH2F* hist){
    int nbins = hist->GetYaxis()->GetNbins();
    std::vector<double> edges;
    edges.reserve(nbins+1);
    //double* edges = new double[nbins+1];
    //double edges[nbins+1];
    for (int i = 0; i <nbins+1 ; ++i){
      edges.push_back(hist->GetYaxis()->GetBinLowEdge(i+1));
    }
    return edges;
  }

  
  vector<double> getVals(Tree* t, Expr* var, TCut cut){
    t->GetTTree()->Draw(">>myList", cut , "entrylist");
    TEntryList* list = (TEntryList*)gDirectory->Get("myList");
    t->GetTTree()->SetBranchStatus("*",0);
    t->SetBranches(var->GetVarNames());
    //t->SetBranchStatus(var.c_str(), 1);
    //double d_var;
    //t->SetBranchAddress(var.c_str(), &d_var);
    
    vector<double> vals; 
    
    Long64_t nentries = list->GetN();
    
    vals.reserve(nentries);
    for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
      if (jentry%10000==0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
      int entry = list->GetEntry(jentry);
      t->GetEntry(entry);
      vals.push_back(t->GetVal(var));
    }
    t->GetTTree()->SetBranchStatus("*",1);
    return vals;
  }
  
  vector< vector<double> > getVals(Tree* tree, Expr* var, Var* binvar, TCut cut){
    TTree* t = tree->GetTTree();
    t->Draw(">>myList", cut , "entrylist");
    TEntryList* list = (TEntryList*)gDirectory->Get("myList");
    t->SetBranchStatus("*",0);
    tree->SetBranches(var->GetVarNames());
    tree->SetBranches(binvar->GetExpr()->GetVarNames());
    
    
    Long64_t nentries = list->GetN();
    
    vector<double> bin_edges = binvar->GetEdges();
    vector< vector<double> > vals (bin_edges.size() - 1, vector<double>(0)); 
    
    for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
      if (jentry%10000==0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
      int entry = list->GetEntry(jentry);
      tree->GetEntry(entry);
      double binval = tree->GetVal(binvar->GetExpr());
      if (binval >= bin_edges[0] && binval < bin_edges[bin_edges.size() - 1]){
	int bin = binvar->GetHist()->FindBin(binval);
	double val = tree->GetVal(var);
	vals[bin - 1].push_back(val);
      }
    }
    t->SetBranchStatus("*",1);
    return vals;
  }
  vector< vector< vector<double> > > getVals(Tree* tree, Expr* var, Var2D* binvar, TCut cut){
    TTree* t = tree->GetTTree();
    t->Draw(">>myList", cut , "entrylist");
    TEntryList* list = (TEntryList*)gDirectory->Get("myList");
    t->SetBranchStatus("*",0);
    tree->SetBranches(var->GetVarNames());
    tree->SetBranches(binvar->GetVar1()->GetExpr()->GetVarNames());
    tree->SetBranches(binvar->GetVar2()->GetExpr()->GetVarNames());
    
    Long64_t nentries = list->GetN();
    
    vector<double> bin_edges_x = binvar->GetVar1()->GetEdges();
    vector<double> bin_edges_y = binvar->GetVar2()->GetEdges();
    int xbins = bin_edges_x.size() - 1;
    int ybins = bin_edges_y.size() - 1;
    //int totbins = (bin_edges_x.size() - 1) * (bin_edges_y.size() - 1);
    
    vector< vector< vector<double> > > vals (xbins, vector< vector<double> >(ybins, vector<double>(0))); 
    
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
	//cout<<"Filling : ("<<xbin<<","<<ybin<<") with "<<val<<" for "<<xval<<", "<<yval<<endl;
	//int bin = (ybin - 1)*(bin_edges_y.size() - 1) + (xbin - 1);
	vals[ xbin - 1 ][ ybin - 1 ].push_back(val);
      }
    }
    t->SetBranchStatus("*",1);
    return vals;
  }

  TH1F* GetWeightHist(string name, TH1F* histA, TH1F* histB){
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), histA->GetNbinsX(), histA->GetBinLowEdge(1), histA->GetBinLowEdge(histA->GetXaxis()->GetLast() + 1));
    if (histA->GetNbinsX() == histB->GetNbinsX() && histA->GetBinLowEdge(0) == histB->GetBinLowEdge(0) &&
        histA->GetBinLowEdge(histA->GetXaxis()->GetLast() +1) == histB->GetBinLowEdge(histB->GetXaxis()->GetLast() + 1)){
      for (int i = 0; i < histA->GetNbinsX()+1; ++i){
	if (histB->GetBinContent(i) == 0) hist->SetBinContent(i, 1);
	else hist->SetBinContent(i, histA->GetBinContent(i)/histB->GetBinContent(i)*histB->Integral()/histA->Integral());
      }
    }
    return hist;
  }

  double GetLumi(TFile* f){
    TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
    double lumi = 0;
    double lumi_job;
    int nentries = lumit->GetEntries();
    lumit->SetBranchAddress("IntegratedLuminosity", &lumi_job);
    for (int i = 0 ; i < nentries; ++i){
      lumit->GetEntry(i);
      lumi += lumi_job;
    }
    
    //lumit->Draw("IntegratedLuminosity>>lumihist()","","goff");
    //TH1F* lumihist=(TH1F*)gDirectory->Get("lumihist");
    //double Lumi = (double)lumihist->GetMean()*lumihist->GetEntries();
    return lumi;
    
  }
  double GetLumiError(TFile* f){
    TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
    double lumierr = 0;
    double lumierr_job;
    int nentries = lumit->GetEntries();
    lumit->SetBranchAddress("IntegratedLuminosityErr", &lumierr_job);
    for (int i = 0 ; i < nentries; ++i){
      lumit->GetEntry(i);
      lumierr += lumierr_job;
    }
    /*  TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
	lumit->Draw("IntegratedLuminosityErr>>lumierrhist()","","goff");
	TH1F* lumierrhist=(TH1F*)gDirectory->Get("lumierrhist");
	double LumiError = (double)lumierrhist->GetMean()*lumierrhist->GetEntries();
    */
    return lumierr;
  }

  double GetWeightSum(TTree* t, string w, string cut){
    t->Draw(">>e", cut.c_str() , "entrylist");
    double val;
    TEntryList* l = (TEntryList*)gDirectory->Get("e");
    t->SetBranchStatus("*",0);
    t->SetBranchStatus(w.c_str(),1);
    t->SetBranchAddress(w.c_str(), &val);
    int nentries = l->GetN();
    double sum = 0;
    for (int i = 0; i < nentries; ++i){
      int entry = l->GetEntry(i);
      t->GetEntry(entry);
      sum += val;
    }
    t->SetBranchStatus("*",1);
    return sum;

  }

  vector<double> GetWeightSum(TTree* t, vector<string> weights, string cut){
    t->Draw(">>e", cut.c_str() , "entrylist");
    vector<double> vals(weights.size(), 0.0);
    vector<double> sums(weights.size(), 0.0);
    TEntryList* l = (TEntryList*)gDirectory->Get("e");
    t->SetBranchStatus("*",0);
    for (unsigned int i = 0 ; i < weights.size() ; ++i){
      t->SetBranchStatus(weights.at(i).c_str(),1);
      t->SetBranchAddress(weights.at(i).c_str(), &vals.at(i));
    }
    int nentries = l->GetN();
    for (int i = 0; i < nentries; ++i){
      int entry = l->GetEntry(i);
      t->GetEntry(entry);
      for (int j = 0 ; j < (int)weights.size() ; ++j) sums.at(j) += vals.at(j);
    }
    t->SetBranchStatus("*",1);
    return sums;

  }
  double GetSum(TTree* t, string leaf){
    double sum = 0;
    double d;
    float f;
    TClass* cl = 0; EDataType type;
    if (!t){
      cout<<"Return - null tree passed"<<endl;
      return 0;
    }
    TBranch* branch = (TBranch*)t->GetBranch(leaf.c_str());
    if (!branch){
      cout<<"Branch - "<<leaf<<" - not found"<<endl;
      return 0;
    }
    int vartype = 0;
    branch->GetExpectedType(cl, type);
    t->SetBranchStatus("*",0);
    t->SetBranchStatus(leaf.c_str(),1);
    
    if (type == 8){
      t->SetBranchAddress(leaf.c_str(), &d);
      vartype = 1;
    }
    else if (type == 5){
      t->SetBranchAddress(leaf.c_str(), &f);
      vartype = 2;
    }
    else{
      cout<<"Cannot get sum for branch of type "<<type<<endl;
      return 0;
    }
    int nentries = t->GetEntries();
    for (int i = 0 ; i < nentries; ++i){
      t->GetEntry(i);
      if (vartype == 1){
	sum += d;
      }
      else if (vartype == 2){
	sum += f;
      }
    }
    t->SetBranchStatus("*",1);
    return sum;
  }


  void RemoveErrors(TGraphAsymmErrors* graph){
    int i = 0;
    double x = 0, y = 0;
    Int_t hasPoint = graph->GetPoint(i,x,y);
    while (hasPoint == i){
      graph->SetPointEYhigh(i, 0.0);
      graph->SetPointEYlow(i, 0.0);
      i++;
      hasPoint = graph->GetPoint(i,x,y);
    }
  }

  void saveAsTree(string fileName, vector<string> varNames, string output){
    std::ifstream f(fileName.c_str());
    std::string line;
    std::vector<string> outputsS(varNames.size(), "");
    std::vector<double> outputsD(varNames.size(), 0.0);
    TFile* outputFile = new TFile(output.c_str(), "RECREATE");
    TTree* t = new TTree("tree", "tree");
    for (unsigned int i = 0 ; i < varNames.size() ; ++i){
      //t->Branch(varNames.at(i).c_str());
      t->Branch(varNames.at(i).c_str(), &outputsD.at(i));
    }


    while(std::getline(f, line)){
      std::stringstream ss(line);
      for (unsigned int i = 0 ; i < outputsS.size() ; ++i){
	ss >> outputsS.at(i);
	outputsD.at(i) = atof(outputsS.at(i).c_str());
	//cout<<outputsD.at(i)<<endl;
      }
      t->Fill();
    }
    t->Write();
    outputFile->Close();

  }


  #ifdef WITHPYTHON
  PyObject* geteff_py(string name, PyObject* data, PyObject* eff){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    PyObject* output = 0;

    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      TH1D* eff = geteff(name, data3, graph);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      TH1D* eff = geteff(name, data2, graph);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      TH1D* eff = geteff(name, data1, graph);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }
    return output;

  }
  PyObject* geteff2_py(string name, PyObject* data, PyObject* eff, int varybin){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    PyObject* output = 0;

    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      TH1D* eff = geteff(name, data3, graph, varybin);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      TH1D* eff = geteff(name, data2, graph, varybin);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      TH1D* eff = geteff(name, data1, graph, varybin);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }
    return output;

  }
  PyObject* geteff3_py(string name, PyObject* data, PyObject* eff, int varybin, string form){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    PyObject* output = 0;

    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      TH1D* eff = geteff(name, data3, graph, varybin, form);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      TH1D* eff = geteff(name, data2, graph, varybin, form);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      TH1D* eff = geteff(name, data1, graph, varybin, form);
      output = TPython::ObjectProxy_FromVoidPtr(eff, eff->ClassName());
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }
    return output;
    
  }
  boost::python::list getEffVariations_py(string name, PyObject* data, PyObject* eff){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    pair< TH1D*, vector<TH1D*> > output;


    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      output = getEffVariations(name, data3, graph);
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      output = getEffVariations(name, data2, graph);
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      output = getEffVariations(name, data1, graph);
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }
   boost::python::list ns;
   //having library trouble with boost::python::handle
   /*
   if (output.first) {
     PyObject* py = TPython::ObjectProxy_FromVoidPtr(output.first, output.first->ClassName());
     boost::python::object o(boost::python::handle<>(py));
     ns.append(o);
   }
   for (unsigned int i = 0 ; i < output.second.size() ; ++i ){
     ns.append(boost::python::object(boost::python::handle<>(TPython::ObjectProxy_FromVoidPtr(output.second.at(i), output.second.at(i)->ClassName()))));
     }*/
   
   return ns;

  }
  boost::python::list getEffErrMatrix_py(string name, PyObject* data, PyObject* eff){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    pair< TH1D*, vector<TH1D*> > output;


    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      output = getEffVariations(name, data3, graph);
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      output = getEffVariations(name, data2, graph);
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      output = getEffVariations(name, data1, graph);
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }

    matrix totMat = getEffErrMatrix(output);
    
    boost::python::list ns = mat2PyList<double>(totMat);
    return ns;

  }

  boost::python::list getEffErrMatrix2_py(string name, PyObject* data, PyObject* eff, string f){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    pair< TH1D*, vector<TH1D*> > output;


    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      output = getEffVariations(name, data3, graph, f);
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      output = getEffVariations(name, data2, graph, f);
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      output = getEffVariations(name, data1, graph, f);
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }

    matrix totMat = getEffErrMatrix(output);
    
    boost::python::list ns = mat2PyList<double>(totMat);
    return ns;

  }
  boost::python::list getEffErrMatrix3_py(string name, PyObject* data, PyObject* eff, string f, bool correlated){
    TH3F* data3 = (TH3F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH2F* data2 = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TH1F* data1 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(data));
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    pair< TH1D*, vector<TH1D*> > output;


    if (strcmp(data3->ClassName(), "TH3F") == 0) {
      output = getEffVariations(name, data3, graph, f, correlated);
    }
    else if (strcmp(data2->ClassName(), "TH2F") == 0) {
      output = getEffVariations(name, data2, graph, f, correlated);
    }
    else if (strcmp(data1->ClassName(), "TH1F") == 0) {
      output = getEffVariations(name, data1, graph, f, correlated);
    }
    else{
      cout<<"Can't process input of type :"<<data3->ClassName()<<endl;
    }

    matrix totMat = getEffErrMatrix(output);
    
    boost::python::list ns = mat2PyList<double>(totMat);
    return ns;

  }

  

  PyObject* tgraph2hist_py(string name, PyObject* eff){
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(eff));
    TH1F* o = tgraph2hist(name, graph);
    return TPython::ObjectProxy_FromVoidPtr(o, o->ClassName());
  }

  double GetLumi_py(PyObject* pyf){
    TFile* f = (TFile*)(TPython::ObjectProxy_AsVoidPtr(pyf));
    return GetLumi(f);
  }

  double GetLumiError_py(PyObject* pyf){
    TFile* f = (TFile*)(TPython::ObjectProxy_AsVoidPtr(pyf));
    return GetLumiError(f);
  }
  
  double GetWeightSum_py(PyObject* pyObj, string w, string cut){
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    return GetWeightSum(t, w, cut);
  }
  
  boost::python::list GetWeightSum2_py(PyObject* pyObj, boost::python::list& weights, string cut){
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    vector<string> ws;
    for (int i = 0 ; i < len(weights); ++i) ws.push_back(boost::python::extract<string>(weights[i]));
    vector<double> sums = GetWeightSum(t, ws, cut);
    boost::python::list lOut;
    for (int j = 0 ; j < (int)sums.size() ; ++j) lOut.append(sums.at(j));
    return lOut;
  }
  double GetSum_py(PyObject* pyf, string leaf){
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyf));
    return GetSum(t, leaf);
  }
  void RemoveErrors_py(PyObject* pyObj){
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    int i = 0;
    double x = 0, y = 0;
    Int_t hasPoint = graph->GetPoint(i,x,y);
    while (hasPoint == i){
      graph->SetPointEXhigh(i, 0.0);
      graph->SetPointEYlow(i, 0.0);
    }
  }
  boost::python::list cholesky_py(boost::python::list& ns){
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

  double standard_deviation_py(boost::python::list& ns){
    vector<double> vals;
    for (unsigned int i = 0; i < len(ns); ++i){
      double d = boost::python::extract<double>(ns[i]);
      vals.push_back(d);
    }
    return standard_deviation(vals);
  }

  void saveMatrix_py(string name, boost::python::list& ns){
    vector< vector<double> > m = pyList2Mat<double>(ns);
    saveMatrix(name, m);
  }
  void saveTH1F_py(string name, PyObject* pyObj){
    TH1F* h = Py2RootObj<TH1F>(pyObj);
    saveTH1F(name, h);
  }
  void saveTGraph_py(string name, PyObject* pyObj){
    TGraph* g = Py2RootObj<TGraph>(pyObj);
    saveTGraph(name, g);
  }
  void saveTGraphErrs_py(string name, PyObject* pyObj){
    TGraph* g = Py2RootObj<TGraph>(pyObj);
    saveTGraphErrs(name, g);
  }

  void saveAsTree_py(string fileName, boost::python::list& varNames, string output){
    vector<string> varNamesS = pyList2Vec<string>(varNames);
    saveAsTree(fileName, varNamesS, output);

  }

  
  template <typename T> boost::python::list vec2PyList(std::vector<T> vect){
    boost::python::list ns3;
    for (unsigned int j = 0 ; j < vect.size() ; ++j){
      ns3.append(vect[j]);
    }
    return ns3;
  }
  template <typename T> boost::python::list mat2PyList(std::vector< std::vector<T> > mat){
    boost::python::list ns2;
    for (unsigned int i = 0 ; i < mat.size() ; ++i){
      boost::python::list ns3;
      for (unsigned int j = 0 ; j < mat[i].size() ; ++j){
	ns3.append(mat[i][j]);
      }
      ns2.append(ns3);
    }
    return ns2;
  }
  template <typename T> vector<T> pyList2Vec(boost::python::list& ns){
    vector<T> vals;
    for (unsigned int i = 0; i < len(ns); ++i){
      T d = boost::python::extract<T>(ns[i]);
      vals.push_back(d);
    }
    return vals;
  }
  template <typename T> vector< vector<T> > pyList2Mat(boost::python::list& ns){
    vector<vector<T> > vals;
    for (unsigned int i = 0; i < len(ns); ++i){
      boost::python::list ns2 = boost::python::extract<boost::python::list>(ns[i]);
      //boost::python::list ns2 = ns[i];
      vector<T> v;
      for ( unsigned int j = 0 ; j < len(ns2); ++j){
	T d = boost::python::extract<T>(ns2[j]);
	v.push_back(d);
      }
      vals.push_back(v);
    }
    return vals;
  }
  template <typename T> T* Py2RootObj(PyObject* pyObj){
    return (T*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  }
  template <typename T> PyObject* Root2PyObj(T* cxxObj){
    return TPython::ObjectProxy_FromVoidPtr(cxxObj, cxxObj->ClassName());
  }
  
  #endif
  
}
