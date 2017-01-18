// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Var.h>

#ifdef WITHPYTHON
#include <TPython.h>
#include <boost/python.hpp>
#endif
#include <Utils.h>


using namespace std;

class EffVar : public Var {
 public:
  EffVar();
  EffVar(string name, string var, int nbins, double lo, double hi, string prefix = "");
  EffVar(string name, string var, std::vector<double> edges, string prefix = "");
  EffVar(EffVar* varA, EffVar* varB, string prefix = "");
  EffVar(string name, TFile* f, string prefix = "");
  EffVar(string name, PyObject* f, string prefix);
  EffVar(string name, TGraphAsymmErrors* grah);
  EffVar(string name, PyObject* graph);
  
  TH1F* GetTotHist();
  TH1F* GetPassHist();
  TH1F* GetFailHist();
  TH1F* GetBkgTotHist();
  TH1F* GetBkgPassHist();
  TH1F* GetMeanTotHist();

  TObjArray* GetEffRWVaryHiGraphs();
  TObjArray* GetEffRWVaryLoGraphs();
  TObjArray* GetEffRWVaryHiPassHists();
  TObjArray* GetEffRWVaryLoPassHists();


  
  void MakeTGraph();
  void MakeEffHist(bool ClopperPearsonError = true);

  void FillVar(bool pass, double v_var, double weight = 1, double effw = 1.0);
  void FillVar(bool pass, int i_var   , double weight = 1, double effw = 1.0);
  void FillVar(bool pass, float v_var , double weight = 1, double effw = 1.0);

  void FillVar(bool pass, double v_var, Utils::weight weight, Utils::weight effw);
  //void FillVar(string type, double v_pltvar, double v_var, double efflo, double effhi, double weight = 1);
  //void FillVar(string type, double v_pltvar, float v_var, double efflo, double effhi, double weight = 1);
  //void FillVar(string type, double v_pltvar, int i_var, double efflo, double effhi, double weight = 1);

  void AddSystematic(double pc);
  void AddInvSystematic(double pc);
  void AddSystematic(std::vector<double> pc);
  void SetPrefix(string name);

  void AddEffScaleVaryHists(TH2F* scales);

  void FillBkgHists(double lo, double hi);
  void Normalise(double N);

  void AddSystematic1(double pc);
  void AddSystematic2(std::vector<double> pc);

  void RemoveErrors();

  static vector<double> GetBinEdgesX(TH2F* hist);
  static vector<double> GetBinEdgesY(TH2F* hist);
  //boost::python::list GetEdges_py();
  //boost::python::list GetBinEdges_py(PyObject* pyObj);
  TGraphAsymmErrors* GetEffGraph();
  TGraphAsymmErrors* GetSmearedEffGraph();

  TObjArray* GetTotCBFits();
  TObjArray* GetPassCBFits();
  TObjArray* GetTotBkgFits();
  TObjArray* GetPassBkgFits();
  TObjArray* GetPassHists();
  TObjArray* GetTotHists();
  TObjArray* GetFailHists();

#ifdef WITHPYTHON
  PyObject* GetTotHist_py();
  PyObject* GetPassHist_py();
  boost::python::list GetBinEdgesX_py(PyObject* pyObj);
  PyObject* GetSmearedEffGraph_py();
  PyObject* GetEffGraph_py();
#endif
  
 protected:
  string m_var;
  TObjArray* m_tothists;
  TObjArray* m_passhists;
  TObjArray* m_failhists;
  TObjArray* m_totCBFits;
  TObjArray* m_passCBFits;
  TObjArray* m_totBkgFits;
  TObjArray* m_passBkgFits;
  TH1F* m_tothist;
  TH1F* m_passhist;
  TH1F* m_passhist_effrw;
  TH1F* m_failhist;
  TH1F* m_bkgtot;
  TH1F* m_bkgpass;
  TH1F* m_meantot;
  TH1F* m_meanpass;
  TGraphAsymmErrors* m_effgraph;
  TH1F* m_effhist;
  const char* m_type;
  bool m_systematic;

  //keep errs for efficiency scaling
  TObjArray* m_reweighteffuperrs;
  TObjArray* m_reweighteffloerrs;


  TObjArray* m_rweff_varyhi_passhists;
  TObjArray* m_rweff_varylo_passhists;

  TObjArray* m_rweff_varyhi_effgraphs;
  TObjArray* m_rweff_varylo_effgraphs;
};
