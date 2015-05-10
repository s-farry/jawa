#ifndef Include_Utils_H
#define Include_Utils_H
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
#include <Var2D.h>
#include <Tree.h>
#include <TRandom3.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

using namespace std;

namespace Utils{
  void fillhist(TH1F* h , vector<double>& vals);
  double get_mean(vector<double>& vals);
  double standard_deviation(vector<double>& vals, double max = -1);
  double cholesky(double *A, int n);
  vector< vector<double> > cholesky( vector< vector<double> > A);
  void printMatrix(vector< vector<double> > A);
  //boost::python::list getCorrelatedRandoms_py(boost::python::list& ns);
  vector<double> getCorrelatedRandoms(TRandom3* r3, vector< vector<double> > corrs );
  vector<double> getRandoms(TRandom3* r3, int n );

  vector< vector<double> > getCorrelationMatrix(vector<vector<double> > vals);
  vector<double> getColumn(int n, vector< vector<double> > vals);
  double GetMean(vector<double> vals);
  double GetStdDev(vector<double> vals, double mean);
  std::vector<double> GetBinEdgesX(TH2F* hist);
  std::vector<double> GetBinEdgesY(TH2F* hist);

  vector< vector< vector<double> > > getVals(Tree* tree, Expr* expr, Var2D* var, TCut cut="");
  vector< vector<double> > getVals(Tree* tree, Expr* var, Var* binvar, TCut cut);
  vector<double> getVals(Tree* t, Expr* var, TCut cut);
  TH1F* GetWeightHist(string name, TH1F* histA, TH1F* histB);
  double GetSum(TTree* t, string leaf);

  void RemoveErrors(TGraphAsymmErrors* graph);
  //void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  double GetWeightSum(TTree* t, string w, string cut);
  #ifdef WITHPYTHON
  double GetWeightSum_py(PyObject* py, string w, string cut);
  void RemoveErrors_py(PyObject* pyObj);
  double GetSum_py(PyObject* t, string leaf);
  double GetLumi_py(PyObject* f);
  double standard_deviation_py(boost::python::list& ns);
  boost::python::list GetStdDevs_py();
  boost::python::list cholesky_py(boost::python::list& ns);

  #endif

  
};
#endif
