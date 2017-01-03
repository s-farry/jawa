// classes example
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <EffVar.h>
#include <Eff.h>
#include <boost/python.hpp>
#include <Var2D.h>

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

class EffVar2D: public Var2D {
 public:
  EffVar2D();
  EffVar2D(string name, string var1, int nbins1, double lo1, double hi1, string var2, int nbins2, double lo2, double hi2, const char* type1 = "D", const char* type2 = "D", string prefix = "");
  EffVar2D(EffVar2D* varA, EffVar2D* varB, string prefix = "");
  EffVar2D(string name, EffVar* varA, EffVar* varB, string prefix="");
  EffVar2D(string name, TFile* f, string prefix = "");
  EffVar2D(string name, PyObject* f, string prefix = "");
  EffVar2D(string name, string f, string prefix = "");
  EffVar2D(string name, TH2F* effgraph);

  void Init(string name, TFile* f, string prefix = "");
  void MakeTGraphs();

  //string GetName();
  //string GetName1();
  //string GetName2();
  void MakeTGraph(TH1F* pass, TH1F* tot);
  void MakeEffHist(bool ClopperPearsonError = true);
  void FillMeanHists();
  TH2F* GetTotHist();
  TH2F* GetPassHist();
  TH2F* GetFailHist();
  TH2F* Get2DEffGraph();
  PyObject* Get2DEffGraph_py();
  PyObject* GetTotHist_py();

  void FillVar(bool pass, double& v_var1, double& v_var2, double w=1.0);
  void FillVar(bool pass, int& v_var1, double& v_var2, double w);
  void FillVar(bool pass, double& v_var1, int& v_var2, double w);
  void FillVar(bool pass, int& v_var1, int& v_var2, double w);
  //void FillVar(string type, double& v_var1, double& v_var2, double w = 1.0);
  //void FillVar(bool pass, double& v_pltvar, int& i_var, double efflo, double effhi);
  void MakeHists(string name, int npltbins, double pltrangelow, double pltrangehi, bool reweight = false);

  void Normalise(double N);

  TObjArray* GetEffGraphs();
  TObjArray* GetEffRWVarHiGraphs();
  TObjArray* GetEffRWVarLoGraphs();
  TObjArray* GetEffRWVarHiPassHists();
  TObjArray* GetEffRWVarLoPassHists();


 protected:
  string m_prefix;
  TObjArray* m_totCBFits;
  TObjArray* m_passCBFits;

  TH2F* m_tothist;
  TH2F* m_passhist;
  TH2F* m_failhist;

  TObjArray* m_effgraphs;
  TObjArray* m_effhists;
  TH2F* m_2Deffgraph;

  TObjArray* m_tothists;
  TObjArray* m_passhists;
  TObjArray* m_failhists;
  TH2F* m_meantot;
  TH2F* m_meanpass;


  TObjArray* m_rweff_varyhi_passhists;
  TObjArray* m_rweff_varylo_passhists;

  TObjArray* m_rweff_varyhi_effgraphs;
  TObjArray* m_rweff_varylo_effgraphs;

  

};
