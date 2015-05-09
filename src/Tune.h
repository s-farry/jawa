#ifndef Include_Tune_H
#define Include_Tune_H
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <TF1.h>
#include <TPython.h>
#include <boost/python.hpp>
#include <Expr.h>
#include <Var.h>
#include <Tree.h>
#include <ReweightVar.h>
#include <TRandom3.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

using namespace std;

class Tune{
 public:
  Tune();
  Tune(string name);
  Tune(string name, Tree* data, Tree* mc, Expr* tuneVar, Var* fVar = 0, string cut = "");
  void tune();
  vector< pair<double, double> > smear_vals(vector< pair<double, double> >& vals, double mean, double sigma);
  void fillhist(TH1F* h , vector< pair<double, double> >& vals);
  void fillVals();
  void SaveToFile();
  double standard_deviation_py(boost::python::list& ns);
  boost::python::list GetStdDevs_py();
  boost::python::list GetDataVec(int j);
  double cholesky(double *A, int n);
  vector< vector<double> > cholesky( vector< vector<double> > A);
  void printMatrix(vector< vector<double> > A);
  boost::python::list getCorrelatedRandoms_py(boost::python::list& ns);
  boost::python::list cholesky_py(boost::python::list& ns);
  vector<double> getCorrelatedRandoms(vector< vector<double> > corrs );
  //void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  double metric(const double* params);

  void SetSigmaPars(double init, double step, double lolimit, double uplimit);
  void SetMeanPars(double init, double step, double lolimit, double uplimit);
  void SetPrecision(double precision);
  void SetTolerance(double tolerance);
  void PrintSigmaPars();
  void PrintMeanPars();
  void SetVals(vector< vector< pair<double, double> > > data, vector< vector< pair<double, double> > > mc);
  vector< vector< pair<double, double> > > GetMCVals();
  vector< vector< pair<double, double> > > GetDataVals();

  void ReweightMC(ReweightVar* rwvar);
  void ReweightData(ReweightVar* rwvar);
  double GetWeight(Tree* tree, vector<ReweightVar>);
  void SetSDFactor(double s);

 protected:
  Expr* m_tuneVar;           // name of templates and their location
  Var* m_fVar;               // Templates to be fit

  Tree* m_data;
  Tree* m_mc;
  TCut m_cut;
  TH1F* m_res_sigma;
  TH1F* m_res_mean;
  TH1F* m_stddev;

  TH2F* m_2Dres_sigma;
  TH2F* m_2Dres_mean;
  TH2F* m_2Dstddev;

  double getchi2(double mean, double sigma);
  vector< vector< pair<double, double> > > m_data_vec;
  vector< vector< pair<double, double> > > m_mc_vec;
  vector<ReweightVar*> m_data_rwvars;
  vector<ReweightVar*> m_mc_rwvars;
  vector< double > m_data_stddev;
  vector< pair<double, double> >* m_current_data;
  vector< pair<double, double> >* m_current_mc;
  double m_current_stddev;
  double m_sdfac;

  //Pointer to current vals for the FCN class
  //vector<double>* m_current_data;
  //vector<double>* m_current_mc;
  //double m_current_stddev;
  TObjArray* m_data_array;
  TObjArray* m_mc_array;
  TObjArray* m_mc_corr_array;
  string m_name;
  TRandom3 m_r3;

  double standard_deviation(vector< pair<double, double> >& vals, double max = -1);
  double get_mean(vector< pair<double, double > >& vals);
  vector< vector< pair<double, double> > > getVals(Tree* tree, Expr* expr, Var* var, TCut cut="", vector<ReweightVar*> rwvars = vector<ReweightVar*>(0));
  vector< pair<double, double> > getVals(Tree* tree, Expr* expr, TCut cut="");
  double GetWeight(Tree* tree, vector<ReweightVar*> rwvars);

  double m_mean_init;
  double m_mean_step;
  double m_mean_lolimit;
  double m_mean_uplimit;
  double m_sigma_init;
  double m_sigma_step;
  double m_sigma_lolimit;
  double m_sigma_uplimit;
  double m_precision;
  double m_tolerance;

};

//void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//vector<double>* m_current_data;
//vector<double>* m_current_mc;
//double m_current_stddev;
#endif
