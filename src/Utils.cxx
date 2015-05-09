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
#include <Utils.h>
#include <TCanvas.h>

using namespace std;

namespace Utils{
  
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

  vector<double> getColumn(int n, vector<vector<double> > vals){
    vector<double> col;
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
    for (int i = 0; i <nbins+1 ; ++i){
      edges.push_back(hist->GetXaxis()->GetBinLowEdge(i+1));
    }
    return edges;
  }
  std::vector<double> GetBinEdgesY(TH2F* hist){
    int nbins = hist->GetYaxis()->GetNbins();
    std::vector<double> edges;
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
    return sum;

  }
  double GetWeightSum_py(PyObject* pyObj, string w, string cut){
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    return GetWeightSum(t, w, cut);
  }
  double GetSum(TTree* t, string leaf){
    double sum = 0;
    double e;
    t->SetBranchStatus("*",0);
    t->SetBranchStatus(leaf.c_str(),1);
    t->SetBranchAddress(leaf.c_str(), &e);
    int nentries = t->GetEntries();
    for (int i = 0 ; i < nentries; ++i){
      t->GetEntry(i);
      sum += e;
    }
    t->SetBranchStatus("*",1);
    return sum;
  }

  double GetSum_py(PyObject* pyf, string leaf){
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyf));
    return GetSum(t, leaf);
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

}
