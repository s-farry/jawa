// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPython.h>
#include <boost/python.hpp>
#include <Var.h>

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
  TH1F* GetMeanPassHist();

  PyObject* GetTotHist_py();
  PyObject* GetPassHist_py();
  
  void MakeTGraph();
  void MakeEffHist(bool ClopperPearsonError = true);

  void FillVar(bool pass, double v_pltvar, double v_var, double efflo, double effhi, double weight = 1);
  void FillVar(bool pass, double v_pltvar, int i_var, double efflo, double effhi, double weight = 1);
  void FillVar(bool pass, double v_pltvar, float v_var, double efflo, double effhi, double weight = 1);
  void FillVar(string type, double v_pltvar, double v_var, double efflo, double effhi, double weight = 1);
  void FillVar(string type, double v_pltvar, float v_var, double efflo, double effhi, double weight = 1);
  void FillVar(string type, double v_pltvar, int i_var, double efflo, double effhi, double weight = 1);

  void MakeHists(string name, int npltbins, double pltrangelow, double pltrangehi, bool reweight = false);

  void AddSystematic(double pc);
  void AddInvSystematic(double pc);
  void AddSystematic(std::vector<double> pc);
  void SetPrefix(string name);


  void FillBkgHists(double lo, double hi);
  void FillMeanHists();
  void Normalise(double N);

  void AddSystematic1(double pc);
  void AddSystematic2(std::vector<double> pc);

  void RemoveErrors();

  static vector<double> GetBinEdgesX(TH2F* hist);
  static vector<double> GetBinEdgesY(TH2F* hist);
  //boost::python::list GetEdges_py();
  boost::python::list GetBinEdgesX_py(PyObject* pyObj);
  //boost::python::list GetBinEdges_py(PyObject* pyObj);
  TGraphAsymmErrors* GetEffGraph();
  PyObject* GetEffGraph_py();
  TGraphAsymmErrors* GetSmearedEffGraph();
  PyObject* GetSmearedEffGraph_py();

  TObjArray* GetTotCBFits();
  TObjArray* GetPassCBFits();
  TObjArray* GetTotBkgFits();
  TObjArray* GetPassBkgFits();
  TObjArray* GetPassHists();
  TObjArray* GetTotHists();
  TObjArray* GetFailHists();

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
  TH1F* m_failhist;
  TH1F* m_bkgtot;
  TH1F* m_bkgpass;
  TH1F* m_meantot;
  TH1F* m_meanpass;
  TGraphAsymmErrors* m_effgraph;
  TH1F* m_effhist;
  const char* m_type;
  bool m_systematic;

};
