#ifndef Include_Template_H
#define Include_Template_H
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
#include <Var2D.h>
#include <Var3D.h>
#include <ReweightVar.h>
#include <Tree.h>

#ifdef WITHPYTHON
#include <boost/python.hpp>
#endif

using namespace std;


class Template : public JawaObj {
 public:
  Template(string name);
  Template(string name, TTree* t, TCut* cut);
  Template(string name, TTree* t, TCut* cut, enum EColor color);
  Template(string name, vector<TTree*> trees, TCut* cut);
  Template(string name, vector<TTree*> trees, TCut* cut, enum EColor color);
  Template(string name, Template* t); // essentially a clone
  Template(string name, Template* A, Template* B);
  //Template operator+(const Template& rhs);

  ~Template();
  void Init();
  void AddTree(TTree* t);
  void AddTree(TTree* t, double w);
  void AddTree(string name, TTree* t);
  void AddTree(string name, TTree* t, double w);
  void AddTree(TTree* t, TCut* cut);
  void AddTree(TTree* t, double w, TCut* cut);
  void AddTree(string name, TTree* t, TCut* cut);
  void AddTree(string name, TTree* t, double w, TCut* cut);

  void AddTrees(vector<TTree*>& trees);
  void AddTrees(vector<TTree*>& trees, vector<TCut*>& cuts);

  void ApplyCut();
  void FillVars();
  void AddVar(string name, string var, int bins, double lo, double hi, string prefix = "");
  void AddVar(string name, string var, std::vector<double>& edges, string prefix = "");
  void SetSelCut(TCut* cut);
  void SaveToFile();
  virtual void SaveToCurrentFile();
  bool IsFixed();
  Var* GetVar(string name);
  map<string, Var*> GetVariables();
  map<string, Var2D*> Get2DVariables();
  map<string, Var3D*> Get3DVariables();
  //AlgVar* GetAlgVar(string name);
  Var2D* Get2DVar(string name);
  Var2D* Get2DVar(string name1, string name2);
  Var3D* Get3DVar(string name);
  Var3D* Get3DVar(string name1, string name2, string name3);

  void NormaliseToEvts(double evts, bool fixed = true);
  int GetEvents();
  bool IsStyled();
  void SetStyle(enum EColor fillcolor);
  void SetStyle(enum EColor fillcolor , enum EColor linecolor);
  void SetStyle(enum EColor fillcolor, enum EColor linecolor , Style_t fillstyle);
  void StyleHists();
  void Style(TH1F* hist);
  void NormaliseToMC(double xsec, double acc, double Lumi, double nEvts, bool fixed = false);
  void Scale(double scale, bool fixed = true);
  void Unscale();
  void Reweight(string var, TF1* f);
  void Reweight(string var, std::map<int,double> map);
  void Reweight(string var, TH1F* scales);
  void Reweight(string weightname);
  void Reweight(string var1, string var2, TH2F* hist);
  void SetFitFrac(double f);
  void SetFitFrac(double f, double err);
  void SetNormEvts(double evts);
  double GetFitFrac();
  double GetNormEvts();
  int GetEvts();
  void Add2DVar(string var1, string var2);
  void Add3DVar(string var1, string var2, string var3);
  void Add2DVar(string name, string var1, string var2, string prefix = "");
  void Add3DVar(string name, string var1, string var2, string var3, string prefix = "");
  void AddAsymmetry(string eta1, string eta2);
  void LoadFromFile(const char* file="");
  string GetRootName(const char* file="");
  void IsVerbose();
  TH1F* GetHist(string name);
  std::map<string, int> SetBranches(vector<double>& output_idx);
  void OutputEvts(bool output);
  Tree* GetTree();
  Tree* GetTree(string name);
  void PrintVars();
  void Run();

  // For python
  #ifdef WITHPYTHON
  Template(string name, PyObject* t, PyObject* cut);
  Template(string name, boost::python::list& ns, PyObject* cut);
  void AddVar1_py(string name, string var, int bins, double lo, double hi);
  void AddVar2_py(boost::python::list& ns);
  void AddVar3_py(string name, string var, boost::python::list& ns);
  void SetVars_py(boost::python::list& ns); // Clear variables first
  void AddVars_py(boost::python::list& ns);
  PyObject* GetHist_py(string name);
  boost::python::list GetVariables_py();
  void SetFitFrac_py(double f, double err);
  Tree* GetTree1_py();
  Tree* GetTree2_py(string name);
  void Add2DVar_py(boost::python::list& var);
  void Add2DVars_py(boost::python::list& varlist);
  void Add3DVar_py(boost::python::list& var);
  void Add3DVars_py(boost::python::list& varlist);
  void Reweight1_py(string var, PyObject* tf1);
  void Reweight2_py(string var1, string var2, PyObject* th2f);
  void Reweight3_py(string leaf);
  void SetSelCut_py(PyObject* pyObj);
  void AddTree_py(PyObject* py);
  void AddTree2_py(string name, PyObject* py);
  void AddTree3_py(PyObject* py, double w);
  void AddTree4_py(string name, PyObject* py, double w);
  void AddTree5_py(PyObject* py, PyObject* cut);
  void AddTree6_py(string name, PyObject* py , PyObject* cut);
  void AddTree7_py(PyObject* py, double w    , PyObject* cut);
  void AddTree8_py(string name, PyObject* py , PyObject* cut);
  void AddTrees_py(boost::python::list& ns);
  void AddTrees2_py(boost::python::list& ns, boost::python::list& ns2);
  Var2D* Get2DVar_py(string name1, string name2);
  Var3D* Get3DVar_py(string name1, string name2, string name3);
  void Scale1_py(double scale, bool fixed);
  void Scale2_py(double scale);
  void NormaliseToEvts1_py(double evts, bool fixed);
  void NormaliseToEvts2_py(double evts);
  void SetStyle1_py(int fillcolor);
  void SetStyle2_py(int fillcolor, int linecolor);
  void SetStyle3_py(int fillcolor, int linecolor, int fillstyle);
  void NormaliseToMC1_py(double xsec, double acc, double Lumi, double nEvts, bool fixed);
  void NormaliseToMC2_py(double xsec, double acc, double Lumi, double nEvts);
  #endif
  
 protected:
  //TTree* m_tree;
  std::vector<Tree*> m_trees;
  std::vector<TEntryList*> m_entryLists;
  std::map<string, Var*> m_variables;
  //std::map<string, AlgVar*> m_algvariables;
  std::map<string, Var2D*> m_2Dvariables;
  std::map<string, Var3D*> m_3Dvariables;
  std::vector<ReweightVar*> m_reweightvariables;
  double m_norm;
  TCut* m_selcut;
  vector<TCut*> m_selcuts;
  bool m_fixed;
  int m_evts;

  bool m_fillTree;
  std::map<int,double> m_reweightmap;
  TTree* m_tree;

  bool m_asymm;
  string m_asymmvar1;
  string m_asymmvar2;

  bool m_styled;
  string m_legendname;
  enum EColor m_fillcolor;
  Style_t m_fillstyle;
  enum EColor m_linecolor;
  double m_normN;
  double m_fitFrac;
  double m_fitFracErr;
  TH1F* m_fittedhist;
  bool m_verbose;
  bool m_outputevts;
  TFile* m_outputfile;

};

#endif
