#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <Fit.h>
#include <boost/algorithm/string.hpp>

using namespace std;


Fit::Fit(  ){
  m_name = "";
  m_expr = 0;
  m_tf1 = 0;
  m_hist = 0;
}


Fit::Fit(string name, string func, TH1F* hist ){
  m_name = name;
  m_expr = new Expr(func);
  m_tf1 = new TF1((name+"_fit").c_str(), *m_expr, 20000, 70000, m_expr->GetVarNames().size() - 1);
  m_hist = hist;
  std::vector<string> vars = m_expr->GetVarNames();
  int j = 0;
  for (unsigned int i = 0 ; i < vars.size() ; ++i ){
    if ( vars.at(i) != "x" ) {
      m_tf1->SetParName(j, vars.at(i).c_str());
      j++;
    }
  }
}

TF1* Fit::GetTF1(){ return m_tf1;};


PyObject* Fit::GetTF1_py(){
  //cout<<"TF1: "<<m_tf1<<endl;
  TF1* newCxxObj = new TF1(*m_tf1);
  //cout<<"TF1: "<<newCxxObj<<endl;
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
  //return TPython::ObjectProxy_FromVoidPtr(m_tf1, m_tf1->ClassName());
}

Expr* Fit::GetExpr(){ return m_expr; };

void Fit::FitHist(){
  if (m_hist) m_hist->Fit(m_tf1);
}

void Fit::FitHist(TH1F* hist, string opt){
  m_hist = hist;
  m_hist->Fit(m_tf1, opt.c_str());
}

void Fit::FitHist_py(PyObject* pyObj){
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  FitHist(hist);
}

void Fit::SetParameter(int param, double value){
  if (m_tf1) m_tf1->SetParameter(param, value);
}
void Fit::FixParameter(int param, double value){
  if (m_tf1) m_tf1->FixParameter(param, value);
}
double Fit::GetParameter(int param){
  double par = -1;
  if (m_tf1) par = m_tf1->GetParameter(param);
  return par;
}
double Fit::GetParError(int param){
  double parerr = -1;
  if (m_tf1) parerr = m_tf1->GetParError(param);
  return parerr;
}


void Fit::SetParLimits(int param, double lo, double hi){
  if (m_tf1) m_tf1->SetParLimits(param, lo, hi);
}
void Fit::PrintParameters(){
  for (unsigned int i = 0 ; i < m_expr->GetVarNames().size(); ++i){
    cout<<m_expr->GetVarNames().at(i)<<endl;
  }
}
void Fit::SetRange(double lo, double hi){
  if (m_tf1) m_tf1->SetRange(lo, hi);
}
