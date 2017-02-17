// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParameter.h>
#include <Expr.h>
#include <Var.h>
#include <Template.h>

using namespace std;

class FitAnalysis : public JawaObj {
 protected:
  string m_func;
  string m_tofit;
  Expr*  m_expr;
  std::map<string, vector<Fit*> > m_fits;
  //std::map<string, std::map<string, TH1F*> > m_hists;
  std::map<string, double> m_initvals;
  std::map<string, double> m_fixvals;
  std::map<string, std::pair<double, double> > m_initrange;
  std::map<string, int> m_paridx;
  std::map<string, vector<TH1F*> > m_hists;
  std::map<string, vector<TH1F*> > m_parhists;
  std::vector<int> m_tomean, m_toentries, m_tomax, m_torms;

  std::map<string, double> m_pars;
  std::map<string, double> m_parerrs;

  std::map<string, std::tuple<int, int, double> > m_finitvals;

  TH1F* m_hist;
  Fit* m_fit;

 public:
  FitAnalysis();
  FitAnalysis(string name, string tofit, string func);
  void SetVal(string par, double val);
  void SetVal(string var, int bin, string par, double val);
  void FixVal(string par, double val);
  void SetRange(string par, double lo, double hi);

  void SaveToFile();
  void SaveToFile(string name);
  std::vector<TH1F*> GetHists(string name);
  std::vector<TH1F*> GetParHists(string name);
  TH1F* GetHist();

  void Init(Template* a);
  void FitIt(string opt = "");
  void SetToRMS(string name);
  void SetToMean(string name);
  void SetToEntries(string name);
  void SetToMax(string name);




  //for python
  #ifdef WITHPYTHON
  void FitIt1_py();
  void FitIt2_py(string opt);
  PyObject* GetHist_py(string name, int i);
  PyObject* GetParHist_py(string name, int i);
  PyObject* GetHist2_py();
  void SaveToFile1_py();
  void SaveToFile2_py(string name);
  void SetVal1_py(string par, double val);
  void SetVal2_py(string var, int bin, string par, double val);
  #endif


};
