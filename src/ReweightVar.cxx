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
  if (!m_func) cout<<"Warning - Function passed is Null!"<<endl;
}

ReweightVar::ReweightVar(string name, TH1F* hist){
  m_exprs.push_back(new Expr(name));
  m_hists.push_back(hist);
  m_weighttype = ReweightVar::Hist1D;
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;
}

ReweightVar::ReweightVar(string name1, string name2, TH2F* hist){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  //m_func = 0;
  m_2dhists.push_back(hist);
  m_weighttype = ReweightVar::Hist2D;
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;

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
  double w    = 1.0;
  double werr = 0.0;
  if (m_weighttype == ReweightVar::Func){
    if (m_func) w = w * m_func->Eval(val);
    else cout<<"Function passed is null"<<endl;
  }
  else if (m_weighttype == ReweightVar::Hist1D && m_hists.size() == 1){
    if (m_hists.at(0)){
      //cout<<"Gonna get 1d weight from hist 1d"<<endl;
      int bin = m_hists.at(0)->FindBin(val);
      //cout<<"val is "<<val<<endl;
      //cout<<"bin is "<<bin<<endl;
      if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
	//cout<<"passed loops"<<endl;
	w = w * m_hists.at(0)->GetBinContent(bin);
	werr = sqrt(pow(werr,2) + pow (m_hists.at(0)->GetBinError(bin), 2));
      }
    }
    else cout<<"Hist passed is null!"<<endl;
  }
  else if (m_weighttype == ReweightVar::Vect){
    if ( m_vect.count(val) == 1 ) { 
      w = w * m_vect.at(val);
    }
  }
  else if (m_weighttype == ReweightVar::Leaf){
    w = w * val;
  }
  else{
    cout<<"Could not get weight"<<endl;
  }
  //cout<<"returning 1d weight of "<<w<<endl;
  return w;
}

double ReweightVar::GetWeight(double val1, double val2){
  double w = 1.0;
  if (m_weighttype == ReweightVar::Hist2D &&
      m_2dhists.size() == 1 && m_exprs.size() == 2){
    if (m_2dhists.at(0)){
    int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
    int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
    //cout<<"vals are: "<<val1<<" "<<val2<<endl;
    //cout<<"bins are: "<<bin1<<" "<<bin2<<endl;
    if (bin1 > 0 && bin1 <= m_2dhists.at(0)->GetNbinsX() &&
	bin2 > 0 && bin2 <= m_2dhists.at(0)->GetNbinsY()){
      w = w * m_2dhists.at(0)->GetBinContent(bin1,bin2);
    }
    }
    else cout<<"hist passed is null!"<<endl;
  }
  else{
    cout<<"Could not get 2d weight"<<endl;
  }
  //cout<<"returning 2d weight of "<<w<<endl;
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
