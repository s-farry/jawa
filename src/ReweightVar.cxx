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
#include <ReweightVar.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

ReweightVar::ReweightVar(string name, TF1* func){
  m_exprs.push_back(new Expr(name));
  m_func = func;
  m_weighttype = ReweightVar::Func;
}

ReweightVar::ReweightVar(string name, TH1F* hist){
  m_exprs.push_back(new Expr(name));
  //m_func = 0;
  m_hists.push_back(hist);
  
  m_weighttype = ReweightVar::Hist1D;

}

ReweightVar::ReweightVar(string name1, string name2, TH2F* hist){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  //m_func = 0;
  m_2dhists.push_back(hist);
  
  m_weighttype = ReweightVar::Hist2D;
}

ReweightVar::ReweightVar(string name, std::map<int,double> map){
  m_exprs.push_back(new Expr(name));
  //m_func = 0;
  m_vect = map;
  m_weighttype = ReweightVar::Vect;

}
ReweightVar::ReweightVar(string weight){
  m_exprs.push_back(new Expr(weight));
  m_weighttype = ReweightVar::Leaf;
}

double ReweightVar::GetWeight(double val){
  //cout<<"in getweight double val"<<endl;
  //cout<<"m_weightfunc is "<<m_weightfunc<<endl;
  double w    = 1.0;
  double werr = 0.0;
  if (m_weighttype == ReweightVar::Func && m_func){
    w = w * m_func->Eval(val);
  }
  else if (m_weighttype == ReweightVar::Hist1D && m_hists.size() == 1){
    int bin = m_hists.at(0)->FindBin(val);
    if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
      w = w * m_hists.at(0)->GetBinContent(bin);
      werr = sqrt(pow(werr,2) + pow (m_hists.at(0)->GetBinError(bin), 2));
    }
  }
  else if (m_weighttype == ReweightVar::Vect){
    if ( m_vect.count(val) == 1 ) { 
      w = w * m_vect.at(val);
    }
  }
  else if (m_weighttype == ReweightVar::Leaf){
    w = w * val;
  }
  return w;
}

double ReweightVar::GetWeight(double val1, double val2){
  double w = 1.0;
  if (m_weighttype == ReweightVar::Hist2D &&
      m_2dhists.size() == 1 && m_exprs.size() == 2){
    int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
    int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
    w = w * m_2dhists.at(0)->GetBinContent(bin1,bin2);
  }
  return w;
}

double ReweightVar::GetWeight(float val){
  double d_val = double(val);
  return GetWeight(d_val);
}

double ReweightVar::GetWeight(int val){
  double d_val = double(val);
  return GetWeight(d_val);
}

vector<Expr*> ReweightVar::GetExprs(){
  return m_exprs;
}
Expr* ReweightVar::GetExpr(){
  return m_exprs.at(0);
}
