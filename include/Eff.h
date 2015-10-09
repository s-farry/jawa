// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;

class Eff{
 protected:
  const char* m_name;
  double m_eff;
  double m_effErrHi;
  double m_effErrLo;
  int m_Ntot;
  int m_Npass;
 public:
  Eff();
  Eff(const char* name, double eff, double efferrhi, double efferrlo);
  Eff(const char* name, int Ntot, int Npass);
  Eff(const char* name, double Ntot, double Npass);
  double GetEff();
  double GetEffErrHi();
  double GetEffErrLo();
  void SetEff(double eff);
  void SetEffErrHi(double efferrhi);
  void SetEffErrLo(double efferrlo);
  void Print();
  void AddSystematic(double pc);
  void AddInvSystematic(double pc);

};
