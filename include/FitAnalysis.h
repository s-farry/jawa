// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
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
 public:
  FitAnalysis();
  FitAnalysis(string name, string tofit, string func);
  void SetVal(string par, double val);
  void FixVal(string par, double val);
  void SetRange(string par, double lo, double hi);

  void SaveToFile();

  void Init(Template* a);
  void FitIt(string opt = "");

  //for python
  #ifdef WITHPYTHON
  void FitIt1_py();
  void FitIt2_py(string opt);
  #endif


};
