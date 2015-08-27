#ifndef Include_EntryList_H
#define Include_EntryList_H
#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TColor.h>
#include <stdlib.h>
#include <Expr.h>
#include <list>

using namespace std;


/*class Cut{
 public:
  Cut(string var, string op, double val);
  Cut(string cut);
  string GetVar();
  char GetOperation();
  double GetVal();
  bool Pass(double val);

 private:
  string m_var;
  string m_op;
  double m_val;

  };*/

class EntryList{
 public:
  EntryList();
  //void AddCut(string cut);
  //void AddCuts(string cuts);
  //void AddCuts(TCut cut);
  //void AddCut(string var, string op, double val);
  //void ApplyCuts();
  std::list<int>& GetList();
  void SetExpr(Expr* e);
  Expr* GetExpr();
  void SetTotEntries(int N);
  void SetPassEntries(int N);
  int GetTotEntries();
  int GetPassEntries();
  void AddEntry(int e);
  string GetCut();


 private:
  //TTree* m_tree;
  //std::vector<Cut*> m_cuts;
  std::list<int> m_entries;
  int m_passentries;
  int m_totentries;
  Expr* m_expr;

};

#endif
