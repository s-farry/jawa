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
//#include <AlgVar.h>
#include <Var2D.h>
#include <Var3D.h>
#include <ReweightVar.h>
#include <Tree.h>
#include <Template.h>
#ifdef WITHPYTHON
#include <boost/python.hpp>
#endif

using namespace std;

class MWTemplate : public Template {
 public:
  MWTemplate(string name);
  MWTemplate(string name, TTree* t, TCut* cut);
  MWTemplate(string name, TTree* t, TCut* cut, enum EColor color);
  MWTemplate(string name, MWTemplate* A, MWTemplate* B);

  virtual void SaveToCurrentFile();
  void SaveToFile();

  void FillVars();
  void AddWeight(string wName, string weightname);
  TH1F* GetWeightHist(string var, string wname);
  TH2F* GetWeight2DHist(string var, string wname);
  std::map<string, std::map<string, TH1F*> > GetWeightHists();
  std::map<string, std::map<string, TH2F*> > Get2DWeightHists();

  void ScaleAllWeights(double s);
  void PrintWeights();
  void ScaleWeight(string w, double s);


  #ifdef WITHPYTHON
  MWTemplate(string name, PyObject* t, PyObject* cut);
  PyObject* GetWeightHist_py(string var, string wname);
  PyObject* GetWeight2DHist_py(string var, string wname);
  #endif

 private:
  std::map<string, std::map<string, TH1F*> > m_varhists;
  std::map<string, std::map<string, TH2F*> > m_2dvarhists;
  std::map<string, ReweightVar*> m_mwreweightvars;
};
