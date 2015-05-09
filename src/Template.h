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
//#include <AlgVar.h>
#include <Var2D.h>
#include <Var3D.h>
#include <ReweightVar.h>
#include <Tree.h>
#include <boost/python.hpp>


using namespace std;


class Template{
 public:
  Template(string name);
  Template(string name, TTree* t, TCut cut);
  Template(string name, TTree* t, TCut cut, enum EColor color);
  Template(string name, PyObject* t, PyObject* cut);
  Template(string name, Template* t); // essentially a clone
  Template(string name, Template* A, Template* B);
  Template operator+(const Template& rhs);

  //~Template();
  void Init();
  void SetTree(TTree* t);
  void AddTree(string name, TTree* t);
  void SetTree_py(PyObject* py);
  void SetTree2_py(PyObject* py, double w);
  void AddTree_py(string name, PyObject* py);
  void SetTree(TTree* t, double w);
  void SetTrees(vector<TTree*>& trees);
  void SetTrees_py(boost::python::list& ns);
  void ApplyCut();
  void FillVars();
  void AddVar(string name, string var, int bins, double lo, double hi, string prefix = "");
  void AddVar(string name, string var, std::vector<double>& edges, string prefix = "");
  //void AddAlgVar(string name, string varexp, int bins, double lo, double hi, string prefix = "");
  //void AddAlgVar(string name, string varexp, std::vector<double> edges, string prefix = "");
  void SetSelCut(TCut cut);
  void SetSelCut_py(PyObject* pyObj);
  void SaveToFile();
  virtual void SaveToCurrentFile();
  bool IsFixed();
  void SetName(string name);
  string GetName() const;
  Var* GetVar(string name);
  map<string, Var*> GetVariables();
  map<string, Var2D*> Get2DVariables();
  map<string, Var3D*> Get3DVariables();
  //AlgVar* GetAlgVar(string name);
  Var2D* Get2DVar(string name);
  Var2D* Get2DVar(string name1, string name2);
  Var3D* Get3DVar(string name);
  Var3D* Get3DVar(string name1, string name2, string name3);
  Var2D* Get2DVar_py(string name1, string name2);
  Var3D* Get3DVar_py(string name1, string name2, string name3);

  void NormaliseToEvts(double evts, bool fixed = true);
  void NormaliseToEvts1_py(double evts, bool fixed);
  void NormaliseToEvts2_py(double evts);
  int GetEvents();
  bool IsStyled();
  void SetStyle(enum EColor fillcolor);
  void SetStyle(enum EColor fillcolor , enum EColor linecolor);
  void SetStyle(enum EColor fillcolor, enum EColor linecolor , Style_t fillstyle);
  void SetStyle1_py(int fillcolor);
  void SetStyle2_py(int fillcolor, int linecolor);
  void SetStyle3_py(int fillcolor, int linecolor, int fillstyle);
  void StyleHists();
  void Style(TH1F* hist);
  void NormaliseToMC(double xsec, double acc, double Lumi, double nEvts, bool fixed = false);
  void NormaliseToMC1_py(double xsec, double acc, double Lumi, double nEvts, bool fixed);
  void NormaliseToMC2_py(double xsec, double acc, double Lumi, double nEvts);
  void Scale(double scale, bool fixed = true);
  void Unscale();
  void Scale1_py(double scale, bool fixed);
  void Scale2_py(double scale);
  void Reweight(string var, TF1 f);
  void Reweight(string var, std::map<int,double> map);
  void Reweight(string var, TH1F* scales);
  void Reweight(string weightname);
  void Reweight(string var, TH1F* hist, string var2, TH1F* hist2);
  void Reweight(string var1, string var2, TH2F* hist);
  void Reweight(string var1, string var2, TH2F* hist, string var3, string var4, TH2F* hist2);
  void Reweight1_py(string var, PyObject* tf1);
  void Reweight2_py(string var, PyObject* th1f);
  void Reweight3_py(string var, PyObject* th1f, string var2, PyObject* th1f2);
  void Reweight4_py(string var1, string var2, PyObject* th2f);
  void Reweight5_py(string var1, string var2, PyObject* th2f, string var3, string var4, PyObject* th2f2);
  void Reweight6_py(string leaf);
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
  void Add2DVar_py(boost::python::list& var);
  void Add2DVars_py(boost::python::list& varlist);
  void Add3DVar_py(boost::python::list& var);
  void Add3DVars_py(boost::python::list& varlist);
  void AddAsymmetry(string eta1, string eta2);
  void LoadFromFile(const char* file="");
  string GetRootName(const char* file="");
  void IsVerbose();
  TH1F* GetHist(string name);
  std::map<string, int> SetBranches(vector<double>& output_idx);
  void OutputEvts(bool output);
  Tree* GetTree();
  Tree* GetTree(string name);
  Tree* GetTree1_py();
  Tree* GetTree2_py(string name);

  // For python
  void AddVar1_py(string name, string var, int bins, double lo, double hi);
  void AddVar2_py(boost::python::list& ns);
  void AddVar3_py(string name, string var, boost::python::list& ns);
  void SetVars_py(boost::python::list& ns); // Clear variables first
  void AddVars_py(boost::python::list& ns);
  PyObject* GetHist_py(string name);
  boost::python::list GetVariables_py();
  void PrintVars();
  void SetFitFrac_py(double f, double err);

 protected:
  //TTree* m_tree;
  std::vector<Tree*> m_trees;
  std::vector<TEntryList*> m_entryLists;
  string m_name;
  std::map<string, Var*> m_variables;
  //std::map<string, AlgVar*> m_algvariables;
  std::map<string, Var2D*> m_2Dvariables;
  std::map<string, Var3D*> m_3Dvariables;
  std::vector<ReweightVar> m_reweightvariables;
  double m_norm;
  TCut m_selcut;
  bool m_fixed;
  int m_evts;

  bool m_reweight;
  bool m_reweightfunc;
  bool m_reweight_map;
  bool m_reweight_hist;
  bool m_fillTree;
  TF1 m_reweightTF1;
  TH1F* m_reweightTH1F;
  std::map<int,double> m_reweightmap;
  string m_reweightvar;
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