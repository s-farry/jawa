#include <iostream>
#ifndef Include_Var3D_H
#define Include_Var3D_H
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <THStack.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <EntryList.h>
#include <Var.h>

using namespace std;

class Var3D : public JawaObj {
 public:
  Var3D(string name);
  Var3D(string name , string var1 , int bins1 , double lo1 , double hi1 , string var2 , int bins2 , double lo2 , double hi2 , string var3, int bins3, double lo3, double hi3, string prefix="");
  Var3D(string name, Var* var1, Var* var2, Var* var3, string prefix="");
  //~Var3D();
  Var3D(string name , Var3D* varA , Var3D* varB, string prefix = "");
  void FillHist(double val1, double val2, double val3);
  void FillHist(double val1, double val2, double val3, double w);
  string GetName();
  string GetVarName1();
  string GetVarName2();
  string GetVarName3();
  string GetName1();
  string GetName2();
  string GetName3();
  //double GetLo1();
  //double GetLo2();
  //double GetHi1();
  //double GetHi2();
  TH3F* GetHist();
  TProfile* GetProfile();
  TProfile* GetProfile2();
  TH3F* GetFwdHist();
  TH3F* GetBwdHist();
  TH1F* GetCombHist();
  void NormaliseToEvts(double evts);
  void NormaliseToMC(double xsec, double acc, double Lumi, double nEvts);
  void Scale(double scale);
  void NormaliseHist(bool doX = false);
  void FillAsymmetry(double xval, double yval, double eta1, double eta2, double w);
  Var* GetVar1();
  Var* GetVar2();
  Var* GetVar3();

  
  int GetBins1();
  double GetVar1Lo();
  double GetVar1Hi();
  int GetBins2();
  double GetVar2Lo();
  double GetVar2Hi();
  int GetBins3();
  double GetVar3Lo();
  double GetVar3Hi();
#ifdef WITHPYTHON
  PyObject* GetHist_py();
  PyObject* GetProfile_py();
  PyObject* GetProfile2_py();
#endif


 protected:
  string m_varname1;
  string m_varname2;
  string m_varname3;
  string m_name1;
  string m_name2;
  string m_name3;
  string m_prefix;
  int    m_nbins1;
  int    m_nbins2;
  int    m_nbins3;
  double m_lo1;
  double m_hi1;
  double m_lo2;
  double m_hi2;
  double m_lo3;
  double m_hi3;
  Var* m_Var1;
  Var* m_Var2;
  Var* m_Var3;

  TH3F* m_hist;
  TH1F* m_hist_comb;
  TH3F* m_hist_fwd;
  TH3F* m_hist_bwd;
  TProfile* m_prof;
  TProfile* m_prof2;
  TProfile* m_prof3;
  TProfile* m_prof4;
  TProfile* m_prof5;
  TProfile* m_prof6;

};
#endif
