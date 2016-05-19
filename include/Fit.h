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

class Fit{
 protected:
  string m_name;
  Expr* m_expr;
  TF1* m_tf1;
  TH1F* m_hist;
 public:
  Fit();
  Fit(string name, string func, TH1F* hist = 0);
  TF1* GetTF1();
  Expr* GetExpr();
  void FitHist(TH1F* hist);
  void FitHist();
  void SetParameter(int param, double value);
  void FixParameter(int param, double value);
  double GetParameter(int param);
  double GetParError(int param);
  void SetParLimits(int param, double lo, double hi);
  void PrintParameters();
  void SetRange(double lo, double hi);

  //For python
  PyObject* GetTF1_py();
  void FitHist_py(PyObject* pyObj);
};
