#ifndef Include_Tree_H
#define Include_Tree_H
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
#include <TPython.h>
#include <Expr.h>
#include <EntryList.h>

using namespace std;

class Data{
 public:
  Data(string type);
  ~Data();
  double* d;
  int* i;
  float* f;
  bool* b;
  unsigned int* u;
  long* l;
  Long64_t* l64;
  unsigned long* ul;
  ULong64_t* ul64;
  string type;
  string GetType();
};

class Tree : public JawaObj{
 public:
  Tree(string name , TTree* t, double w = 1.0);
  Tree(string name, PyObject* t, double w = 1.0);
  ~Tree();
  void SetBranches(std::vector<string> variables);
  void SetBranch(string name);
  double GetVal(string var);
  double GetVal(Expr* e);
  float GetFloatVal(string var);
  unsigned int GetUIntVal(string var);
  int GetIntVal(string var);
  double GetDoubleVal(string var);
  //template<typename T>  T GetVal(string var);
  Data* GetData(string var);
  double GetWeight();
  void GetEntry(int i);
  TTree* GetTTree();
  string GetBranchType(string name);
  bool isSet(string name);
  void AddBranches(Expr* e);
  void SetWeight(double w);
  EntryList* GetEntryList(Expr* e);
  int GetEntries(Expr* e);
  double GetMean(Expr* e, Expr* cut);
  double GetMean(Expr* e, EntryList* entries);
  double GetStdDev(Expr* e, double mean, Expr* cut);
  double GetStdDev(Expr* e, double mean, EntryList* entries);
  vector< vector<double> > getCorrelationMatrix(vector<Expr*> exprs, Expr* cut);


  #ifdef WITHPYTHON
  boost::python::list getCorrelationMatrix_py(boost::python::list& exprs, Expr* cut);
  //For exposing to python
  PyObject* GetTTree_py();
  double GetVal_py(string var);
  double GetVal2_py(Expr* e);
  //static bool is_number(const std::string& s);
  double GetStdDev_py(Expr* e, double mean, Expr* cut);
  double GetMean_py(Expr* e, Expr* cut);
  #endif

 private:
  TTree* m_tree;
  double m_weight;
  std::map<string, Data*> m_output;
  bool m_verbose;
};

#endif
