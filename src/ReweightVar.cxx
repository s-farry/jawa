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

ReweightVar::ReweightVar(string name, TF1* func) : JawaObj("ReweightVar", name){
  m_exprs.push_back(new Expr(name));
  m_func = func;
  m_form= 0;
  m_weighttype = ReweightVar::Func;
  if (!m_func) cout<<"Warning - Function passed is Null!"<<endl;
}

ReweightVar::ReweightVar(string name, TH1F* hist) : JawaObj("ReweightVar", name){
  m_exprs.push_back(new Expr(name));
  m_hists.push_back(hist);
  m_weighttype = ReweightVar::Hist1D;
  m_form= 0;
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;
}

ReweightVar::ReweightVar(string name1, string name2, TH2F* hist) : JawaObj("ReweightVar", name1+"_"+name2){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  //m_func = 0;
  m_2dhists.push_back(hist);
  m_weighttype = ReweightVar::Hist2D;
  m_form = 0;
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;

}

ReweightVar::ReweightVar(string name1, string name2, TH1F* hist, string form) : JawaObj("ReweightVar", name1+"_"+name2){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  //m_func = 0;
  m_hists.push_back(hist);
  m_weighttype = ReweightVar::Hist1D;
  m_form = 0;
  if (form != "") m_form = new Expr(form);
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;

}

ReweightVar::ReweightVar(string name1, string name2, string name3, string name4, TH2F* hist, string form) : JawaObj("ReweightVar", name1+"_"+name2+"_"+name3+"_"+name4){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  m_exprs.push_back(new Expr(name3));
  m_exprs.push_back(new Expr(name4));
  //m_func = 0;
  m_2dhists.push_back(hist);
  m_weighttype = ReweightVar::Hist2D;
  m_form = 0;
  if (form != "") m_form = new Expr(form);
  if (!hist) cout<<"Warning - hist passed is Null!"<<endl;

}

ReweightVar::ReweightVar(string name, std::map<int,double> map) : JawaObj("ReweightVar",name){
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
    else info()<<"Function passed is null"<<endl;
  }
  else if (m_weighttype == ReweightVar::Hist1D && m_hists.size() == 1){
    if (m_hists.at(0)){
      //info()<<"Gonna get 1d weight from hist 1d"<<endl;
      int bin = m_hists.at(0)->FindBin(val);
      //info()<<"val is "<<val<<endl;
      //info()<<"bin is "<<bin<<endl;
      if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
	//info()<<"passed loops"<<endl;
	w = w * m_hists.at(0)->GetBinContent(bin);
	werr = sqrt(pow(werr,2) + pow (m_hists.at(0)->GetBinError(bin), 2));
      }
    }
    else info()<<"Hist passed is null!"<<endl;
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
    info()<<"Could not get weight"<<endl;
  }
  //info()<<"returning 1d weight of "<<w<<endl;
  return w;
}

double ReweightVar::GetWeight(double val1, double val2){
  double w = 1.0;
  if (m_weighttype == ReweightVar::Hist2D &&
      m_2dhists.size() == 1 && m_exprs.size() == 2){
    if (m_2dhists.at(0)){
      int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
      int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
      //info()<<"vals are: "<<val1<<" "<<val2<<endl;
      //info()<<"bins are: "<<bin1<<" "<<bin2<<endl;
      if (bin1 > 0 && bin1 <= m_2dhists.at(0)->GetNbinsX() &&
	  bin2 > 0 && bin2 <= m_2dhists.at(0)->GetNbinsY()){
	w = w * m_2dhists.at(0)->GetBinContent(bin1,bin2);
      }
    }
    else info()<<"hist passed is null!"<<endl;
  }
  else if (m_weighttype == ReweightVar::Hist1D &&
	   m_hists.size() == 1 && m_exprs.size() == 2)
    {
      if(m_hists.at(0)){
	int bin1 = m_hists.at(0)->GetXaxis()->FindBin(val1);
	int bin2 = m_hists.at(0)->GetXaxis()->FindBin(val2);
	if (!m_form){
	  w = w * m_hists.at(0)->GetBinContent(bin1) *
	    m_hists.at(0)->GetBinContent(bin2);
	}
	else{
	  std::vector<double> weights;
	  weights.push_back(m_hists.at(0)->GetBinContent(bin1));
	  weights.push_back(m_hists.at(0)->GetBinContent(bin2));
	  w = w * m_form->GetVal(weights);
	}
      }
    }
  else{
    info()<<"Could not get 2d weight"<<endl;
  }
  //info()<<"returning 2d weight of "<<w<<endl;
  return w;
}

double ReweightVar::GetWeightErr(double val1, double val2){
  double err = 0.0;
  if (m_weighttype == ReweightVar::Hist2D &&
      m_2dhists.size() == 1 && m_exprs.size() == 2){
    if (m_2dhists.at(0)){
      int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
      int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
      //info()<<"vals are: "<<val1<<" "<<val2<<endl;
      //info()<<"bins are: "<<bin1<<" "<<bin2<<endl;
      if (bin1 > 0 && bin1 <= m_2dhists.at(0)->GetNbinsX() &&
	  bin2 > 0 && bin2 <= m_2dhists.at(0)->GetNbinsY()){
	err = sqrt(err*err + pow(m_2dhists.at(0)->GetBinError(bin1,bin2),2));
      }
    }
    else info()<<"hist passed is null!"<<endl;
  }
  else if (m_weighttype == ReweightVar::Hist1D &&
	   m_hists.size() == 1 && m_exprs.size() == 2)
    {
      if(m_hists.at(0)){
	int bin1 = m_hists.at(0)->GetXaxis()->FindBin(val1);
	int bin2 = m_hists.at(0)->GetXaxis()->FindBin(val2);
	if (!m_form){
	  err = sqrt(err*err +
		     pow(m_hists.at(0)->GetBinError(bin1),2) + 
		     pow(m_hists.at(0)->GetBinError(bin2),2));
	}
	else{
	  std::vector<double> weights;
	  weights.push_back(m_hists.at(0)->GetBinContent(bin1));
	  weights.push_back(m_hists.at(0)->GetBinContent(bin2));
	  err = sqrt(err*err + pow(m_form->GetVal(weights),2));
	}
      }
    }
  else{
    info()<<"Could not get 2d weight"<<endl;
  }
  //info()<<"returning 2d weight of "<<w<<endl;
  return err;
}

double ReweightVar::GetWeight(double val1, double val2, double val3, double val4){
  double w = 1.0;
  if (m_weighttype == ReweightVar::Hist2D &&
      m_2dhists.size() == 1 && m_exprs.size() == 4){
    std::vector<double> weights;
    if (m_2dhists.at(0)){
      int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
      int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
      int bin3 = m_2dhists.at(0)->GetXaxis()->FindBin(val3);
      int bin4 = m_2dhists.at(0)->GetYaxis()->FindBin(val4);
      if (bin1 > 0 && bin1 <= m_2dhists.at(0)->GetNbinsX() &&
	  bin2 > 0 && bin2 <= m_2dhists.at(0)->GetNbinsY() &&
	  bin3 > 0 && bin3 <= m_2dhists.at(0)->GetNbinsX() &&
	  bin4 > 0 && bin4 <= m_2dhists.at(0)->GetNbinsY()){
	weights.push_back(m_2dhists.at(0)->GetBinContent(bin1,bin2));
	weights.push_back(m_2dhists.at(0)->GetBinContent(bin3,bin4));
	if (!m_form) {
	  w = weights.at(0) * weights.at(1);
	}
	else{
	  w = w*m_form->GetVal(weights);
	}
      }
    }
    else info()<<"hist passed is null!"<<endl;
  }
  else{
    info()<<"Could not get 2d weight"<<endl;
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

int ReweightVar::GetBin(double val1, double val2){

  int bin1 = -1, bin2 = 0, nbins2 = 0;

  if (m_2dhists.at(0)){
    double xlo = m_2dhists.at(0)->GetXaxis()->GetBinLowEdge(1);
    double ylo = m_2dhists.at(0)->GetYaxis()->GetBinLowEdge(1);
    double xhi = m_2dhists.at(0)->GetXaxis()->GetBinUpEdge(m_2dhists.at(0)->GetXaxis()->GetNbins());
    double yhi = m_2dhists.at(0)->GetXaxis()->GetBinUpEdge(m_2dhists.at(0)->GetYaxis()->GetNbins());

    if (val1 >= xlo && val1 <= xhi && val2 >=ylo && val2 <= yhi){
      bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1) -1;
      bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2) -1;
      nbins2 = m_2dhists.at(0)->GetYaxis()->GetNbins();
    }
  }

  return (bin2*nbins2 + bin1);

}
