// classes example
#include <iostream>
#include <Tree.h>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <EffVar2D.h>
#include <boost/python.hpp>

using namespace std;

Double_t fitGaussLine(Double_t* x, Double_t* par);
Double_t fitCB(Double_t* x, Double_t* par);
Double_t fitCBGauss(Double_t* x, Double_t* par);

class EfficiencyBaseClass {
 public:
  
  int m_npltbins;
 protected:
  string m_name;
  string m_pltvar;
  bool m_verbose;
  std::vector<const char*> passvars;

  std::map<std::string, EffVar*> m_variables;
  std::map<std::string, EffVar2D*> m_2Dvariables;

  TCut m_selcut;

  double m_pltrangelow;
  double m_pltrangehi;
  //double countrangelow;
  //double countrangehi;

  double m_efflo, m_effhi;

  
  Eff m_toteff;

  TH1F* m_tot;
  TH1F* m_pass;
  TH1F* m_fail;

  std::vector<TEntryList*> entryListTot;
  std::vector<TEntryList*> entryListPass;

  //TTree* t;
  TCut m_passcut;

  TObjArray* effhists;

  public:
  string GetName();
  std::vector<Tree*> m_trees;
  void SetPltRange(string var, int bins, double lo, double hi);
  void AddPassVar(const char* var);
  void SetTree(TTree* tree);
  void SetTree_py(PyObject* tree);
  void SetTrees(TTree* priTree, TTree* secTree);
  void SetVerbose(bool verbose);
  bool GetVerbose();
  void AddTree(TTree* tree);
  void AddTree_py(PyObject* tree);
  void SetSelectionCut(TCut cut);
  void SetSelectionCut_py(PyObject* cut);
  std::pair<TF1*,TF1*> FitHistogram(TH1F* massplot, double lo, double hi, string opt = "Z0_CB");
  void StripTree(TCut cut);
  void StripTrees(TCut cut);
  void SetPassCut(TCut cut);
  void SetPassCut_py(PyObject* cut);
  void MakeEntryLists();
  std::map<string, EffVar*> GetVariables();
  void SetVariables(std::map<string, EffVar*>);
  void SetEffRange(double lo, double hi);
  Eff GetTotEff();
  //double as they can be weighted
  double m_Ntot;
  double m_Npass;
  
};

class EfficiencyClass: public EfficiencyBaseClass {

 public:
  EfficiencyClass(string name);
  EfficiencyClass(string name, EfficiencyClass* effA, EfficiencyClass* effB);
  void MakeHists();
  void AddVar(string name, string var, int bins, double lo, double hi);
  void AddVar(string name, string var, vector<double> edges);
  void AddVar1_py(string name, string var, int bins, double lo, double hi);
  void AddVar2_py(string name, string var, vector<double> edges);
  void AddVar3_py(string name, string var, int bins, float lo, float hi);
  void AddVar4_py(string name, string var, boost::python::list& ns);
  void AddVar5_py(boost::python::list& ns);
  void AddVars_py(boost::python::list& ns);
  void FitHists(double lo, double hi);
  void MakeEfficiencyGraph(bool fromFit=false);
  void MakeEfficiencyGraph_py();
  void LoadFromFile(const char* file="");
  void LoadFromFile_py();
  void SaveToFile_py();
  void SaveToFile(const char* file="");
  void PrintEfficiencies(string name);
  void PrintTwikiEfficiencies(string name);
  static void PrintNTwikiEfficiencies(string name, std::vector<std::pair<string, EfficiencyClass> > classes);
  //Eff GetEfficiency(const char* varname, const char* treevar, TTree* tree, TCut selcut);
  Eff GetEfficiency(const char* varname, TH1F* hist);
  EffVar* GetVar(string name);
  EffVar2D* Get2DVar(string name);
  void LoopEntries();
  void SetBranches(Tree* t);
  void PrintVars();
  void FreeBranches(Tree* t);
  double FillVars(bool pass, Tree* t);
  string GetRootName(const char* file);
  void AddFits(TObjArray* hists, TObjArray* fits);
  void AddFits(TObjArray* hists, TObjArray* fitssig, TObjArray* fitsbkg);
  void Add2DVar(string varA, string varB, string name="");
  void Add2DVar_py(string varA, string varB);
  void Add2DVars_py(boost::python::list& ns);

  void ReweightVar(string var, std::map<int, double> );
  void ReweightVar(string var, TF1* func  );
  void ReweightVar(string var, TH1F* hist );


  void Reweight1_py(string var, PyObject* pyObj);

  string m_fitopt;
  void SetFitOpts(string var);

  void Normalise(double N = 1);

  void CorrectGraphs(EfficiencyClass* classA, EfficiencyClass* classB, string opt = "M");
  TGraphAsymmErrors* DivideTGraphs(TGraphAsymmErrors* numer, TGraphAsymmErrors* denom);
  TGraphAsymmErrors* CombineTGraphs(std::vector<TGraphAsymmErrors*> graphs);
  TH2F* DivideTHists(TH2F* numer, TH2F* denom);
  TH2F* CombineTHists(std::vector<TH2F*> graphs);


  void AddSystematic(double pc);
  void AddInvSystematic(double pc);
  void AddSystematic(string name, double pc);
  void AddSystematic(string name, std::vector<double> pc);

  void AddSystematic2_py(double pc);
  void AddSystematic1_py(string name, double pc);
  void AddSystematic3_py(string name, boost::python::list& ns);

  void FillBkgHists();
  void FillMeanHists();


  int GetBin(double value, TH1F* hist);

  void RemoveErrors();

  bool VarExists(string var);

  bool   m_reweight;
  bool   m_reweight_func;
  bool   m_reweight_map;
  bool   m_reweight_bin;
  std::map<int, double> m_reweightmap;
  TF1*   m_reweight_tf1;
  TH1F*  m_reweight_hist;
  string m_reweightvar;
  bool   m_fillbkg;
  bool   m_fillmean;

  double GetCorrectedEfficiency1_py(string var, PyObject* h);
  double GetCorrectedEfficiency2_py(string var, PyObject* t, string leaf);
  boost::python::list GetCorrectedEfficiency3_py(string var, boost::python::list& hists, bool smear);

  double GetCorrectedEfficiency(string var, TH1F* h);
  double GetCorrectedEfficiency(string var, TTree* t, string name);
  std::vector<double> GetCorrectedEfficiency(string var, std::vector<TH1F*> hists, bool smear = false);


};
