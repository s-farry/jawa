#ifndef Include_ReweightVar_H
#define Include_ReweightVar_H
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
#include <Expr.h>

using namespace std;

class ReweightVar{
 public:
  ReweightVar(string name, TH1F* hist);
  ReweightVar(string name1, string name2, TH2F* hist);
  ReweightVar(string name, TF1* func);
  ReweightVar(string name1, TH1F* hist1, string name2, TH1F* hist2);
  ReweightVar(string name, std::map<int,double> map);
  ReweightVar(string weight);

  enum WeightType{Hist1D, Hist2D, Func, Leaf, Vect};
  
  //~ReweightVar();
  double GetWeight(double val);
  double GetWeight(double val1, double val2);
  double GetWeight(double val1, double val2, double val3, double val4);
  double GetWeight(int val);
  double GetWeight(float val);
  Expr* GetExpr();
  vector<Expr*> GetExprs();

  string GetWeightName();



 private:
  vector<string> m_names;
  vector<Expr*> m_exprs;
  vector<TH1F*>  m_hists;
  vector<TH2F*>  m_2dhists;
  TF1* m_func;
  std::map<int,double> m_vect;
  string m_weightname;
  WeightType m_weighttype;  
};
#endif
