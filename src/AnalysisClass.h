#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <Template.h>
#include <Fitter.h>
#include <boost/python.hpp>


using namespace std;

class AnalysisClass{
 public:
  AnalysisClass(string name);
  //~AnalysisClass();
  void SetSelCut(TCut cut);
  void SetSelCut_py(PyObject* cut);
  void SetData(Template* temp);
  void SetData(string name, TTree* t);
  void SetData_py(string name, PyObject* t);
  void SetData(string name, TTree* t, TCut cut);
  void AddTemplate(Template* temp);
  void AddTemplate(string name, TTree* t);
  void AddTemplate(string name, TTree* t, TCut cut);
  void AddTemplate(string name, TTree* t, enum EColor color);
  void AddTemplate(string name, TTree* t, TCut cut, enum EColor color);
  void AddTemplate1_py(string name, PyObject* t);
  void AddTemplate2_py(string name, PyObject* t, PyObject* cut);
  void AddTemplate3_py(string name, PyObject* t, int color);
  void AddTemplate4_py(string name, PyObject* t, PyObject* cut, int color);
  void AddTemplate5_py(Template* t);
  void NormaliseTemplates(double n = 1);
  void ApplyCuts();
  void FillVars();
  void AddVar(string name, string var, int bins, double lo, double hi);
  void AddVar1_py(string name, string var, int bins, double lo, double hi);
  void AddVar2_py(boost::python::list& ns);
  void AddVar3_py(string name, string var, boost::python::list& ns);
  void AddVars_py(boost::python::list& ns);
  void AddVar(string name, string var, std::vector<double> edges);
  //void AddAlgVar(string name, string varexp, int bins, double lo, double hi);
  //void AddAlgVar(string name, string varexp, std::vector<double> edges);
  void Add2DVar(string var1, string var2);
  void Add2DVars_py(boost::python::list& ns);
  void SaveToFile(string output = "");
  void SaveToFile1_py();
  void SaveToFile2_py(string output);
  Template* GetTemplate(string name);
  Template* GetData();
  Fitter* GetFitter();
  string GetName();
  THStack* MakeStack(string name);
  THStack* GetStack(string name);
  THStack* MakeCombStack(string name);
  void MakeStacks();
  void StyleTemplates();
  Var GetVar(string name);
  static double GetLumi(TFile* f);
  double GetLumi_py(PyObject* f);
  static double GetLumiError(TFile* f);
  void ConstrainRatio(string tempA, string tempB, double r);
  TFractionFitter* TFracFit(string var, double lo = 0, double hi = 0, bool combine = false);
  TFractionFitter* RedoFit(string var, double lo = 0, double hi = 0, bool combine = false);
  TFractionFitter* TFracFit(string var1, string var2);
  void AddFitter(string var, double lo = 0, double hi = 0, bool combine = false);
  void Add2DFitter(string var);
  void ClearFitter();
  void ApplyFitResults(bool combine = false);
  void Apply2DFitResults();
  void AddFitter_py(string var);
  void ApplyFitResults_py();
  TCanvas* DrawFitted();
  double GetChi2nDoF();

  PyObject* TFracFit_py(string var);
  PyObject* RedoFit_py(string var);
  PyObject* TFracFit2_py(string var, double lo, double hi);
  PyObject* TFracFit3_py(string var1, string var2);
  PyObject* GetStack_py(string name);

  TFractionFitter* CombTFracFit(string var1, string var2);
  PyObject* CombTFracFit_py(string var1, string var2);

  void Replace(string toRemove, string toAdd);
  void ReplaceInStack(string toRemove, string toAdd);
  void ReplaceInFit(string toRemove, string toAdd);

  void RemoveFromStack(string name);
  void RemoveFromFit(string name);
  void Remove(string name);
  void AddToStack(string name);
  void AddToFit(string name);
  void UnscaleTemplates();

  vector<string> GetToFit();
  void SetToFit(vector<string> toFit);
  vector<string> GetToStack();
  void SetToStack(vector<string> toStack);


  boost::python::list GetToFit_py();
  void SetToFit_py(boost::python::list& ns);

  boost::python::list GetToStack_py();
  void SetToStack_py(boost::python::list& ns);

  PyObject* DrawFitted_py();


 private:
  string m_name;
  TCut   m_selcut;
  std::vector<string> m_stackorder;
  std::vector<string> m_fitorder;
  std::map<string,Template*> m_templates;
  Template* m_data;
  std::map<string, THStack*> m_stacks;
  std::vector<string> m_vars;
  std::vector<string> m_2Dvars;
  TH1F* m_fit;
  Fitter* m_fitter;
  double m_fitchi2;
  double m_ndof;
  pair<pair<string, string>,double> m_conratio;
  bool m_fitratio;

};

