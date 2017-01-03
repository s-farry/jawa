#ifndef Include_Utils_H
#define Include_Utils_H
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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
//#include "Math/Functor.h"
#include <iostream>

using namespace std;
typedef vector< vector< double > > matrix;

namespace Utils{
  enum VaryType{AllUp = 999, AllDown = -999, None = 0};

  struct weight {
    double val;
    double err;
    int bin;
  };

  void saveMatrix(string name, matrix A);
  void saveTH1F(string name, TH1F* h);
  void saveTGraph(string name, TGraph* g);
  void saveTGraphErrs(string name, TGraph* g);
  TH1D* geteff(string name, TH3F* data, TGraphAsymmErrors* eff, int varybin = 0, string form = "a*b");
  TH1D* geteff(string name, TH2F* data, TGraphAsymmErrors* eff, int varybin = 0, string form = "a");
  TH1D* geteff(string name, TH1F* data, TGraphAsymmErrors* eff, int varybin = 0, string form = "a");
  pair< TH1D* , vector<TH1D*> >  getEffVariations(string name, TH3F* data, TGraphAsymmErrors* eff, string form = "a*b", bool correlated = false);
  pair< TH1D* , vector<TH1D*> >  getEffVariations(string name, TH2F* data, TGraphAsymmErrors* eff, string form = "a", bool correlated = false);
  pair< TH1D* , vector<TH1D*> >  getEffVariations(string name, TH1F* data, TGraphAsymmErrors* eff, string form = "a", bool correlated = false); 
  matrix getEffErrMatrix( pair < TH1D*, vector<TH1D* >  > p );
  
  TH1F* tgraph2hist(string name, TGraphAsymmErrors* graph);
  void fillhist(TH1F* h , vector<double>& vals);
  double get_mean(vector<double>& vals);
  double standard_deviation(vector<double>& vals, double max = -1);
  double cholesky(double *A, int n);
  vector< vector<double> > cholesky( vector< vector<double> > A);
  void printMatrix(vector< vector<double> > A);
  //boost::python::list getCorrelatedRandoms_py(boost::python::list& ns);
  vector<double> getCorrelatedRandoms(TRandom3* r3, vector< vector<double> > corrs );
  vector<double> getRandoms(TRandom3* r3, int n );

  matrix getCorrelationMatrix(vector<vector<double> > vals);
  vector<double> getColumn(int n, vector< vector<double> > vals);
  double GetMean(vector<double> vals);
  double GetStdDev(vector<double> vals, double mean);
  std::vector<double> GetBinEdgesX(TH2F* hist);
  std::vector<double> GetBinEdgesY(TH2F* hist);

  vector< matrix > getVals(Tree* tree, Expr* expr, Var2D* var, TCut cut="");
  matrix getVals(Tree* tree, Expr* var, Var* binvar, TCut cut);
  vector<double> getVals(Tree* t, Expr* var, TCut cut);
  TH1F* GetWeightHist(string name, TH1F* histA, TH1F* histB);
  double GetSum(TTree* t, string leaf);

  void RemoveErrors(TGraphAsymmErrors* graph);
  //void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  double GetWeightSum(TTree* t, string w, string cut);
  vector<double> GetWeightSum(TTree* t, vector<string> w, string cut);
  void saveAsTree(string fileName, vector<string> varNames, string output);


  #ifdef WITHPYTHON
  PyObject* tgraph2hist_py(string name, PyObject* graph);
  PyObject* geteff_py(string name, PyObject* data, PyObject* eff);
  PyObject* geteff2_py(string name, PyObject* data, PyObject* eff, int varybin);
  PyObject* geteff3_py(string name, PyObject* data, PyObject* eff, int varybin, string form);
  boost::python::list getEffVariations_py(string name, PyObject* data, PyObject* eff);
  //boost::python::list getEffVariations_py(PyObject* data, PyObject* eff, string form);
  boost::python::list getEffErrMatrix_py(string name, PyObject* data, PyObject* eff);
  boost::python::list getEffErrMatrix2_py(string name, PyObject* data, PyObject* eff, string s);
  boost::python::list getEffErrMatrix3_py(string name, PyObject* data, PyObject* eff, string s, bool correlated);

  double GetWeightSum_py(PyObject* py, string w, string cut);
  boost::python::list GetWeightSum2_py(PyObject* py, boost::python::list& weights, string cut);
  void RemoveErrors_py(PyObject* pyObj);
  double GetSum_py(PyObject* t, string leaf);
  double GetLumi_py(PyObject* f);
  double GetLumiError_py(PyObject* f);
  double standard_deviation_py(boost::python::list& ns);
  boost::python::list GetStdDevs_py();
  boost::python::list cholesky_py(boost::python::list& ns);
  void saveMatrix_py(string name, boost::python::list& ns);
  void saveTH1F_py(string name, PyObject* pyObj);
  void saveTGraph_py(string name, PyObject* pyObj);
  void saveTGraphErrs_py(string name, PyObject* pyObj);

  template<typename T> boost::python::list vec2PyList(std::vector<T> vect);
  template<typename T> boost::python::list mat2PyList(std::vector< std::vector<T> > mat);
  template <typename T> vector<T> pyList2Vec(boost::python::list& ns);
  template <typename T> vector< vector<T> > pyList2Mat(boost::python::list& ns);
  template <typename T> T* Py2RootObj(PyObject* pyObj);
  template <typename T> PyObject* Root2PyObj(T* cxxObj);

  void saveAsTree_py(string fileName, boost::python::list& varNames, string output);
  #endif

  
};
#endif
