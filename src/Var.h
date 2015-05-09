#ifndef Include_Var_H
#define Include_Var_H
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
#include <TPython.h>
#include <Expr.h>

#ifdef WITHPYTHON
#include <boost/python.hpp>
#endif

using namespace std;

class Var{
 public:
  Var(string name , string var , int bins , double lo , double hi , string prefix="");
  Var(string name, string var, std::vector<double>& edges, string prefix="");
  Var(string name);
  Var(string name, Var* varA, Var* varB, string prefix = "");
  Var(string name, Var* v, string prefix = "");
  //~Var();
  void FillHist(double val);
  void FillHist(double val, double w);
  void FillAsymmetry(double val, double eta1, double eta2, double w);
  string GetName();
  string GetVar();
  void SetVarName(string var);
  void SmearHist();
  void ResetHist();
  void Draw(const char* options = "");

  Expr* GetExpr();
  TH1F* GetHist();
  TH1F* GetAsymmHist();
  TH1F* GetAsymmFwdHist();
  TH1F* GetAsymmBwdHist();
  void NormaliseToEvts(double evts);
  void NormaliseToMC(double xsec, double acc, double Lumi, double nEvts);
  void Scale(double scale);
  double GetLo();
  double GetHi();
  int GetBins();
  void MakeAsymmetry();
  void SetScaleFactor(double scale);
  double GetScaleFactor();
  bool IsScaled();
  bool hasMultiVar();
  void Add2ndVar(string name);
  string Get2ndVar();
  vector<double> GetEdges();

  static std::vector<double> GetBinEdges(TH1F* hist);
  static std::vector<double> GetBinEdges(TGraphAsymmErrors* graph);

  double GetVal(std::vector<double>& input);
  void SetVarExp(string varexp);
  vector<string> GetVarNames();
  int FindBin(double val);

  //For python
  #ifdef WITHPYTHON
  boost::python::list GetVarNames_py();
  double GetVal_py(boost::python::list& input);
  boost::python::list GetBinEdges_py(PyObject* pyObj);
  boost::python::list GetEdges_py();
  PyObject* GetHist_py();
  void Draw1_py();
  void Draw2_py(const char* options);
  #endif
  
 protected:
  string m_name;
  //string m_var;
  //string m_var2;
  Expr* m_expr;
  Expr* m_expr2;
  string m_prefix;
  int    m_nbins;
  double m_lo;
  double m_hi;
  TH1F* m_hist;
  std::vector<std::pair<double,double> > m_vals;
  //Forward-Backward Asymmetry hists
  TH1F* m_hist_fwd;
  TH1F* m_hist_bwd;
  TH1F* m_hist_asymm;
  double m_sf;
  bool m_scale;
  bool m_multivar;
  vector<double> m_edges;

};

std::pair< std::vector<string>, std::vector<double> > CombineBinEdges(std::vector<double> edges1, std::vector<double> edges2);
TH1F* MakeAsymmetry(TH1F* fwd, TH1F* bwd, string name = "");
bool is_number(const std::string& s);
#endif
