#include <iostream>
#ifndef Include_Var2D_H
#define Include_Var2D_H
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <THStack.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <EntryList.h>
#include <Var.h>

using namespace std;

class Var2D : public JawaObj{
 public:
  Var2D(string name);
  Var2D(string name , string var1 , int bins1 , double lo1 , double hi1 , string var2 , int bins2 , double lo2 , double hi2 , string prefix="");
  Var2D(string name, Var* var1, Var* var2, string prefix="");
  Var2D(string name, Var2D* var1, Var2D* var2, string prefix="");
  Var2D(string name, Var2D* v, string prefix="");
  //~Var2D();
  void FillHist(double val1, double val2);
  void FillHist(double val1, double val2, double w);
  string GetVarName1();
  string GetVarName2();
  string GetName1();
  string GetName2();
  //double GetLo1();
  //double GetLo2();
  //double GetHi1();
  //double GetHi2();
  TH2F* GetHist();
  TProfile* GetProfile();
  TProfile* GetProfile2();
  TH2F* GetFwdHist();
  TH2F* GetBwdHist();
  TH1F* GetCombHist();
  std::vector<double> GetBinEdges1();
  std::vector<double> GetBinEdges2();
  void NormaliseToEvts(double evts);
  void NormaliseToMC(double xsec, double acc, double Lumi, double nEvts);
  void Scale(double scale);
  void NormaliseHist(bool doX = false);
  void FillAsymmetry(double xval, double yval, double eta1, double eta2, double w);
  Var* GetVar1();
  Var* GetVar2();


  int GetBins1();
  double GetVar1Lo();
  double GetVar1Hi();
  int GetBins2();
  double GetVar2Lo();
  double GetVar2Hi();
#ifdef WITHPYTHON
  PyObject* GetHist_py();
  PyObject* GetProfile_py();
  PyObject* GetProfile2_py();
#endif
  
 protected:
  string m_varname1;
  string m_varname2;
  string m_name1;
  string m_name2;
  string m_prefix;
  int    m_nbins1;
  int    m_nbins2;
  double m_lo1;
  double m_hi1;
  double m_lo2;
  double m_hi2;
  Var* m_Var1;
  Var* m_Var2;

  TH2F* m_hist;
  TH1F* m_hist_comb;
  TH2F* m_hist_fwd;
  TH2F* m_hist_bwd;
  TProfile* m_prof;
  TProfile* m_prof2;

  std::vector<double> m_histbinedges1;
  std::vector<double> m_histbinedges2;

};
#endif
