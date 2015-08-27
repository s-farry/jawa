#ifndef Include_Tune2D_H
#define Include_Tune2D_H
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
#include <Tune.h>
#include <Tree.h>
#include <TRandom3.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

using namespace std;

class Tune2D : public Tune{
 public:
  Tune2D();
  Tune2D(string name);
  Tune2D(string name, Tree* data, Tree* mc, Expr* tuneVar, Var2D* fVar = 0, string cut = "");
  void tune();
  void fill2DVals();
  void SaveToFile();

 private:
  Var2D* m_2DfVar;               // Templates to be fit

  vector< vector< vector< pair<double, double> > > > m_data_2Dvec;
  vector< vector< vector< pair<double, double> > > > m_mc_2Dvec;
  vector< vector< double > > m_data_2Dstddev;
  vector< vector< vector<pair<double, double> > > > getVals(Tree* tree, Expr* expr, Var2D* var, TCut cut="", vector<ReweightVar*> rwvars = vector<ReweightVar*>(0));


  TH2F* m_2Dres_sigma;
  TH2F* m_2Dres_mean;
  TH2F* m_2Dstddev;

};

//void FCN_func(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//vector<double>* m_current_data;
//vector<double>* m_current_mc;
//double m_current_stddev;
#endif
