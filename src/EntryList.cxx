#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TMath.h>
#include <TObjArray.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <EntryList.h>
#include <boost/algorithm/string.hpp>


using namespace std;

/*int main(){

  TFile* wmuf = new TFile("tuples/WMu.WMuNu.MC11.root");
  TTree* wmut = (TTree*)wmuf->Get("WTuple/DecayTree");
  //Cut("pt>50000");


  cout<<wmut->GetEntries()<<endl;
  EntryList el(wmut);
  el.AddCuts("muminus_PT > 60000 && muminus_cpt_0.50 < 2000 && muminus_ETA > 2");
  el.FillEntries();
  cout<<el.GetList().size()<<endl;

  }*/

EntryList::EntryList(){};

/*EntryList::EntryList(TTree* t){
  m_tree = t;
}


void EntryList::AddCut(string cut){
  m_cuts.push_back(new Cut(cut));
}

void EntryList::AddCut(string var, string op, double val){
  m_cuts.push_back(new Cut(var, op, val));
}

void EntryList::AddCuts(TCut cuts){
  string strcut(cuts.GetTitle());
  cout<<strcut<<endl;
  AddCuts(strcut);
}

void EntryList::AddCuts(string cuts){
  string delimeter = "&&";
  std::vector<string> cuts_split;
  size_t pos = cuts.find(delimeter);
  
  while (pos != std::string::npos){
    cuts_split.push_back(cuts.substr(0,pos));
    cuts.erase(0,pos+delimeter.length());
    pos = cuts.find(delimeter);
  }
  cuts_split.push_back(cuts);
  for ( unsigned int i = 0 ; i < cuts_split.size() ; ++i ){
    boost::algorithm::trim(cuts_split.at(i));
    m_cuts.push_back(new Cut(cuts_split.at(i)));
  }
}
*/
std::list<int>& EntryList::GetList(){
  return m_entries;
}

void EntryList::SetExpr(Expr* e){
  m_expr = e;
}
Expr* EntryList::GetExpr(){
  return m_expr;
}

void EntryList::SetTotEntries(int N){
  m_totentries = N;
}

void EntryList::SetPassEntries(int N){
  m_passentries = N;
}

int EntryList::GetTotEntries(){return m_totentries;}
int EntryList::GetPassEntries(){return m_passentries;}


void EntryList::AddEntry(int e){
  m_entries.push_back(e);
}

string EntryList::GetCut(){
  return m_expr->GetExpr();
}

/*
string Cut::GetVar(){
  return m_var;
}

bool Cut::Pass(double val){
  bool pass = false;
  if (m_op == "=="){
    if ( val == m_val) pass = true;
  }
  else if (m_op == ">"){
    if ( val > m_val ) pass = true;
  }
  else if (m_op == "<"){
    if ( val < m_val ) pass = true;
  }
  else if (m_op == ">="){
    if ( val >= m_val ) pass = true;
  }
  else if ( m_op == "<=" ){
    if ( val <= m_val ) pass = true;
  }
  return pass;
}

Cut::Cut(string var, string op, double val){
  m_var = var;
  m_op = op;
  m_val = val;
}

Cut::Cut(string cut){
  std::vector<string> delimeters;
  delimeters.push_back(">");
  delimeters.push_back("<");
  delimeters.push_back(">=");
  delimeters.push_back("<=");
  delimeters.push_back("==");
  int parses = 0;


  for (std::vector<string>::iterator id = delimeters.begin(); id !=delimeters.end(); ++id){
    string delimeter = (*id);
    size_t pos = 0;
    std::string token;
    //while ((pos = cut.find(delimeter)) != std::string::npos) {
    //  token = cut.substr(0,pos);
    //  std::cout <<token<<std::endl;
    pos = cut.find(delimeter);
    
    if (pos != std::string::npos){
      m_var = cut.substr(0,pos);
      boost::algorithm::trim(m_var);
      m_op = cut.substr(pos,delimeter.size());
      boost::algorithm::trim(m_op);
      const char* val = cut.substr(pos+delimeter.size(), cut.size()-1).c_str();
      m_val = atof(val);
      parses++;
    }
  }
  cout<<m_var<<" "<<m_op<<" "<<m_val<<endl;

  //check for further operations
  


  if (parses != 1) cerr<<"Error parsing: "<<cut<<endl;

}

void EntryList::ApplyCuts(){
  std::vector<double> vals(m_cuts.size(), -1.0);
  m_tree->SetBranchStatus("*",0);
  //for (std::vector<Cut*>::iterator ic = m_cuts.begin() ; ic != m_cuts.end() ; ++ic){
  for (unsigned int i = 0 ; i < m_cuts.size(); ++i){
    string var = m_cuts.at(i)->GetVar();
    m_tree->SetBranchStatus(var.c_str(), 1);
    m_tree->SetBranchAddress(var.c_str(), &vals.at(i));
  }

  int nentries = m_tree->GetEntries();

  for (int entry = 0 ; entry < nentries ; entry++ ){
    if (entry%10000==0) cout<<"Entry "<<entry<<" of "<<nentries<<endl;
    m_tree->GetEntry(entry);

    bool pass = true;
    
    unsigned int j = 0;
    while (j < m_cuts.size() && pass == true){
      Cut* selcut = m_cuts.at(j);
      pass = selcut->Pass(vals.at(j));
      j++;
    }
    if (pass) m_entries.push_back(entry);
  }
}


VarExpr::VarExpr(string expr){
  std::vector<string> delimeters;
  delimeters.push_back("+");
  delimeters.push_back("-");
  delimeters.push_back("*");
  delimeters.push_back("/");

  



}
*/
