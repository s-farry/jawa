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
#include <TEntryList.h>
#include <Tree.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>
#include "TRandom3.h"

using namespace std;

/*Tree::Tree(string name , TTree* t){
  m_name    = name;
  m_tree    = t;
  m_weight  = 1.0;
  m_verbose = false;
  }*/
Tree::Tree(string name , TTree* t, double w){
  m_name    = name;
  //t->SetBranchStatus("*",0);
  m_tree    = t;
  m_weight  = w;
  m_verbose = false;
}

Tree::Tree(string name, PyObject* t, double w){
  TTree* tree = (TTree*)(TPython::ObjectProxy_AsVoidPtr(t));
  m_name = name;
  //tree->SetBranchStatus("*",0);
  m_tree = tree;
  m_weight = w;
  m_verbose = false;
}

Tree::~Tree(){
  for (map<string, Data*>::iterator iv = m_output.begin(); iv != m_output.end(); ++iv) {
  delete (*iv).second;
  }
}

TTree* Tree::GetTTree(){
  return m_tree;
}

PyObject* Tree::GetTTree_py(){
  return TPython::ObjectProxy_FromVoidPtr(m_tree, m_tree->ClassName());
}

void Tree::SetBranches(vector<string> variables){
  for (vector<string>::iterator iv = variables.begin() ; iv != variables.end() ; ++iv ){
    string name = (*iv);
    SetBranch(name);
  }
}

void Tree::SetBranch(string name){
  if (m_tree->FindBranch(name.c_str()) == 0) {
    cout<<"No Branch - "<<name<<" - found in Tree"<<endl;
    return;
  }
  if (m_output.find(name) != m_output.end()){
    if (m_verbose) cout<<"Branch "<<name<<" already set"<<endl;
  }
  m_tree->SetBranchStatus(name.c_str(), 1);
  string type = GetBranchType(name);
  Data* da = new Data(type);
  if (type == "double"){
    m_output.insert(std::pair<string, Data*>(name, da));
    m_tree->SetBranchAddress(name.c_str() , m_output[name]->d);
  }
  else if ( type == "float" ){
    m_output.insert(std::pair<string, Data*>(name, da));
    m_tree->SetBranchAddress(name.c_str() , m_output[name]->f);
  }
  else if ( type == "int"){
    m_output.insert(std::pair<string, Data*>(name, da));
    m_tree->SetBranchAddress(name.c_str() , m_output[name]->i);
  }
  else if ( type == "bool" ){
    m_output.insert(std::pair<string, Data*>(name, da));
    m_tree->SetBranchAddress(name.c_str() , m_output[name]->b);
  }
}

Data* Tree::GetData(string var){
  if (m_output.find(var) != m_output.end()){
    return m_output[var];
  }
  else{
    cout<<"Warning: Cout not find branch "<<var<<endl;
    return 0;
  }
}

/*template<typename T>  T Tree::GetVal(string var){
  T v = -1;
  if (m_output.find(var) != m_output.end()){
    if (strcmp(m_output[var]->type.c_str(), "int") == 0){
      v = T(*(m_output[var]->i));
    }
    if (strcmp(m_output[var]->type.c_str(), "float") == 0){
      v = T(*(m_output[var]->f));
    }
    if (strcmp(m_output[var]->type.c_str(), "double") == 0){
      v = T(*(m_output[var]->f));
    }
  }
  return v;
  }*/

float Tree::GetFloatVal(string var){
  float f = -1.0;
  if (m_output.find(var) != m_output.end()){
    if ( m_output[var]->type == "float" ){
      f = *(m_output[var]->f);
    }
  }
  return f;
}
double Tree::GetDoubleVal(string var){
  double d = -1.0;
  if (m_output.find(var) != m_output.end()){
    if ( m_output[var]->type == "double" ){
      d = *(m_output[var]->d);
    }
  }
  return d;
}
int Tree::GetIntVal(string var){
  int i = -1;
  if (m_output.find(var) != m_output.end()){
    if ( m_output[var]->type == "int"){
      i = *(m_output[var]->i);
    }
  }
  return i;
}

double Tree::GetVal(string var){
  double d = -1.0;
  if (m_output.find(var) != m_output.end()){
    if ( m_output[var]->type == "double" ){
      d = *(m_output[var]->d);
    }
    else if ( m_output[var]->type == "float" ){
      d = *(m_output[var]->f);
    }
    else if ( m_output[var]->type == "int" ){
      d = *(m_output[var]->i);
    }
    else if ( m_output[var]->type == "bool" ){
      d = *(m_output[var]->b);
    }

  }
  return d;
}

void Tree::GetEntry(int i){
  m_tree->GetEntry(i);
}

double Tree::GetWeight(){
  return m_weight;
}

bool Tree::isSet(string var){
  return (m_output.find(var) != m_output.end());
}

string Tree::GetBranchType(string name){
  TClass* cl = 0;
  enum EDataType ctype;
  ((TBranch*)m_tree->GetBranch(name.c_str()))->GetExpectedType(cl, ctype);
  string type;

  if (ctype == 8) type = "double";
  else if (ctype == 5) type="float";
  else if (ctype == 18) type="bool";
  else type="int"; // ctype == 3
  return type;
}

Data::Data(string t){
  type = t;
  if ( type == "double" ){
    d = new double();
  } else if ( type == "float" ){
    f = new float();
  } else if ( type == "int" ){
    i = new int();
  } else{
    b = new bool();
  }
}

Data::~Data(){
  if ( type == "double" ){
    delete d;
  } else if ( type == "float" ){
    delete f;
  } else if ( type == "int" ){
    delete i;
  } else{
    delete b;
  }
}


string Data::GetType(){
  return type;
}

double Tree::GetVal_py(string var){
  return GetVal(var);
}
double Tree::GetVal2_py(Expr* e){
  //Expr* e = boost::python::extract<Expr*>(pyObj);
  return GetVal(e);
}


void Tree::AddBranches(Expr* e){
  vector<string> varnames = e->GetVarNames();
  for (vector<string>::iterator iv = varnames.begin(); iv != varnames.end(); ++iv){
    //cout<<"Attempting to add branch "<<(*iv)<<endl;
    //if (!(is_number((*iv)))){
      SetBranch(*iv);
      //}
  }
}

double Tree::GetVal(Expr* e){
  vector<double> input;
  vector<string> varnames = e->GetVarNames();
  for (vector<string>::iterator iv = varnames.begin(); iv != varnames.end(); ++iv){
    input.push_back(GetVal(*iv));
  }
  return e->GetVal(input);
}


EntryList* Tree::GetEntryList(Expr* e){
  m_tree->SetBranchStatus("*",0);
  EntryList* elist = new EntryList();
  int totN = m_tree->GetEntries();
  elist->SetExpr(e);
  elist->SetTotEntries(totN);

  AddBranches(e);

  for (int i = 0; i < totN; ++i){
    //if (i%10000 == 0) cout<<"Entry "<<i<<" of "<<totN<<endl;
    GetEntry(i);
    if (GetVal(e) == 1.0) elist->AddEntry(i); 
  }
  elist->SetPassEntries(elist->GetList().size());
  return elist;
}

int Tree::GetEntries(Expr* e){
  m_tree->SetBranchStatus("*",0);
  int totN = m_tree->GetEntries();
  int passN = 0;
  AddBranches(e);

  for (int i = 0; i < totN; ++i){
    //if (i%10000 == 0) cout<<"Entry "<<i<<" of "<<totN<<endl;
    GetEntry(i);
    if (GetVal(e) == 1.0) passN++; 
  }
  return passN;

}

string Tree::GetName(){
  return m_name;
}
void Tree::SetWeight(double w){
  m_weight = w;
}

bool Tree::is_number(const std::string& s){
  //cout<<"checking is number for "<<s<<endl;
  return (strspn( s.c_str(), "-.01234567890") == s.size() );
}

double Tree::GetMean(Expr* e, EntryList* l){
  list<int> entries = l->GetList();
  SetBranches(e->GetVarNames());
  double mean = 0;
  for (list<int>::iterator in = entries.begin() ; in != entries.end(); ++in ){
    GetEntry(*in);
    mean += GetVal(e);
  }
  mean = mean/entries.size();
  return mean;
}

double Tree::GetMean(Expr* e, Expr* cut){
  EntryList* l = GetEntryList(cut);
  return GetMean(e, l);
}
double Tree::GetStdDev(Expr* e, double mean, EntryList* l){
  list<int> entries = l->GetList();
  SetBranches(e->GetVarNames());
  double sum_deviation = 0;
  for (list<int>::iterator in = entries.begin() ; in != entries.end(); ++in ){
    GetEntry(*in);
    double val = GetVal(e);
    sum_deviation+= ((val - mean)*(val - mean));
  }
  return sqrt(sum_deviation/entries.size());
}
double Tree::GetStdDev(Expr* e, double mean, Expr* cut){
  EntryList* l = GetEntryList(cut);
  return GetStdDev(e, mean, l);
}

vector< vector<double> > Tree::getCorrelationMatrix(vector<Expr*> exprs, Expr* cut){


  //TMatrixD<double> m(exprs.size(), exprs.size());

  EntryList* l = GetEntryList(cut);
  std::list<int> entries = l->GetList();
  vector<double> means;
  vector<double> stddevs;
  //cout<<"Got entry list"<<endl;

  for ( unsigned int i = 0 ; i < exprs.size() ; ++ i ){
    SetBranches(exprs[i]->GetVarNames());
    //cout<<"Set Branches for "<<exprs[i]->GetExpr()<<endl;
    double mean = GetMean(exprs[i], l);
    double stddev = GetStdDev(exprs[i], mean, l);
    means.push_back(mean);
    stddevs.push_back(stddev);
  }

  //Keep track of the sums;
  vector< vector<double> > covs( exprs.size(), vector<double>(exprs.size(), 0.0) );
  int n = entries.size();

  for (list<int>::iterator in = entries.begin() ; in != entries.end(); ++in ){
    //Loop through all entries
    GetEntry(*in);

    for ( unsigned int i = 0 ; i < exprs.size() ; ++ i ){
      double xval = GetVal(exprs[i]);
      double xmean = means[i];
      for ( unsigned int j = 0 ; j < exprs.size() ; ++ j ){
	double yval = GetVal(exprs[j]);
	double ymean = means[j];
	covs[i][j] += (((xval - xmean) * (yval - ymean)))/(n * stddevs[i] * stddevs[j]);
      }
    }
  }
  return covs;
}



//For python

double Tree::GetMean_py(Expr* e, Expr* cut){return GetMean(e, cut);}
double Tree::GetStdDev_py(Expr* e, double mean, Expr* cut){return GetStdDev(e, mean, cut);}

boost::python::list Tree::getCorrelationMatrix_py(boost::python::list& ns, Expr* cut){
  vector<Expr*> exprs;
  for (int i = 0; i < len(ns); ++i){
    //boost::python::object obj = ns[i];
    Expr* e = boost::python::extract<Expr*>(ns[i]);
    exprs.push_back(e);
  }
  

  vector< vector<double> > v = getCorrelationMatrix(exprs, cut);
  cout<<"Got covariance matrix from c"<<endl;

  boost::python::list l;
  for (unsigned int i = 0; i < v.size() ; ++i){
    boost::python::list l1;
    for (unsigned int j = 0; j < v[i].size(); ++j){
      cout<<"Adding "<<v[i][j]<<endl;
      l1.append(v[i][j]);
    }
    l.append(l1);
  }
  return l;
}
