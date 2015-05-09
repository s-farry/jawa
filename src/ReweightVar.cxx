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

ReweightVar::ReweightVar(string name, TF1 func){
  m_names.push_back(name);
  m_exprs.push_back(new Expr(name));
  m_func = func;
  m_weightfunc = true;
  m_weighthist = false;
  m_weight2dhist = false;
  m_weightvect = false;
  m_weightleaf = false;
  m_weightname = "";
}

ReweightVar::ReweightVar(string name, TH1F* hist){
  m_exprs.push_back(new Expr(name));
  m_names.push_back(name);
  //m_func = 0;
  m_hists.push_back(hist);
  m_weightfunc = false;
  m_weighthist = true;
  m_weight2dhist = false;
  m_weightvect = false;
  m_weightleaf = false;
  m_weightname = "";

}

ReweightVar::ReweightVar(string name1, string name2, TH2F* hist){
  //m_expr = 0;
  m_exprs.push_back(new Expr(name1));
  m_exprs.push_back(new Expr(name2));
  m_names.push_back(name1);
  m_names.push_back(name2);
  //m_func = 0;
  m_2dhists.push_back(hist);
  m_weightfunc = false;
  m_weighthist = false;
  m_weight2dhist = true;
  m_weightvect = false;
  m_weightleaf = false;
  m_weightname = "";

}

ReweightVar::ReweightVar(string name1, string name2, TH2F* hist, string name3, string name4, TH2F* hist2){
  m_names.push_back(name1);
  m_names.push_back(name2);
  m_names.push_back(name3);
  m_names.push_back(name4);

  //m_func = 0;
  m_2dhists.push_back(hist);
  m_2dhists.push_back(hist2);
  m_weightfunc = false;
  m_weighthist = false;
  m_weight2dhist = true;
  m_weightvect = false;
  m_weightleaf = false;
  m_weightname = "";

}



ReweightVar::ReweightVar(string name, std::map<int,double> map){
  m_exprs.push_back(new Expr(name));
  m_names.push_back(name);
  //m_func = 0;
  m_vect = map;
  m_weightfunc = false;
  m_weighthist = false;
  m_weightvect = true;
  m_weightleaf = false;
  m_weightname = "";

}
ReweightVar::ReweightVar(string weight){
  m_exprs.push_back(new Expr(weight));
  //Take weight from another leaf in the tree
  m_names.push_back(weight);
  //m_func = 0;
  m_weightfunc = false;
  m_weighthist = false;
  m_weightvect = false;
  m_weightleaf = true;
  m_weightname = weight;
}

ReweightVar::ReweightVar(string name1, TH1F* hist1, string name2, TH1F* hist2){
  m_names.push_back(name1);
  m_names.push_back(name2);
  m_hists.push_back(hist1);
  m_hists.push_back(hist2);
  //m_func = 0; 
  m_weightfunc = false;
  m_weighthist = true;
  m_weightvect = false;
  m_weightleaf = false;
  m_weightname = "";

}



//ReweightVar::~ReweightVar(){
//  m_hist->Delete();
//  m_func->Delete();
//
//}

double ReweightVar::GetWeight(double val){
  //cout<<"in getweight double val"<<endl;
  //cout<<"m_weightfunc is "<<m_weightfunc<<endl;
  double w    = 1.0;
  double werr = 0.0;
  if (m_weightfunc){
    w = w * m_func.Eval(val);
  }
  if (m_weighthist && m_hists.size() == 1){
    int bin = m_hists.at(0)->FindBin(val);
    if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
      w = w * m_hists.at(0)->GetBinContent(bin);
      werr = sqrt(pow(werr,2) + pow (m_hists.at(0)->GetBinError(bin), 2));
    }
  }
  if (m_weightvect){
    if ( m_vect.count(val) == 1 ) { 
      w = w * m_vect.at(val);
    }
  }
  if (m_weightleaf){
    w = w * val;
  }
  return w;
}

double ReweightVar::GetWeight(double val1, double val2){
  double w = 1.0;
  if (m_weight2dhist && m_2dhists.size() == 1 && m_names.size() == 2){
    int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
    int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
    w = w * m_2dhists.at(0)->GetBinContent(bin1,bin2);
  }

  if (m_weighthist && m_hists.size() == 2 && m_names.size() == 2){
    int bin1 = m_hists.at(0)->FindBin(val1);
    int bin2 = m_hists.at(1)->FindBin(val2);

    if (bin1 > 0 && bin1 <= m_hists.at(0)->GetNbinsX()
	&& bin2 > 0 && bin2 <= m_hists.at(1)->GetNbinsX()){
      double w1 = m_hists.at(0)->GetBinContent(bin1);
      double w2 = m_hists.at(1)->GetBinContent(bin2);
      w = w * (w1*w2/(w1 + w2 - 1));
    }
  }
  return w;
}

double ReweightVar::GetWeight(double val1, double val2, double val3, double val4){
  double w = 1.0;
  if (m_weight2dhist && m_2dhists.size() == 2 && m_names.size() == 4){
    int bin1 = m_2dhists.at(0)->GetXaxis()->FindBin(val1);
    int bin2 = m_2dhists.at(0)->GetYaxis()->FindBin(val2);
    //Careful, I assume these should be 3 and 4
    int bin3 = m_2dhists.at(1)->GetXaxis()->FindBin(val3);
    int bin4 = m_2dhists.at(1)->GetYaxis()->FindBin(val4);
    double w1 = m_2dhists.at(0)->GetBinContent(bin1, bin2);
    double w2 = m_2dhists.at(1)->GetBinContent(bin3, bin4);
    w = w * (w1*w2/(w1 + w2 - 1));
  }
  return w;
}



double ReweightVar::GetWeight(float val){
  double w = 1.0;
  if (m_weightfunc){
    w = w * m_func.Eval(val);
  }
  if (m_weighthist && m_hists.size() == 1){
    int bin = m_hists.at(0)->FindBin(val);
    if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
      w = w * m_hists.at(0)->GetBinContent(bin);
    }
    
  }
  if (m_weightvect){
    if ( m_vect.count(val) == 1 ) { 
      w = w * m_vect.at(val);
    }
  }
  if (m_weightleaf){
    //cout<<val<<endl;
    w = w * val;
    //cout<<"W: "<<w<<endl;
  }

  return w;
}

double ReweightVar::GetWeight(int val){
  double w = 1.0;
  if (m_weightfunc){
    w = w * m_func.Eval(val);
  }
  if (m_weighthist && m_hists.size() == 1){
    int bin = m_hists.at(0)->FindBin(val);
    if (bin > 0 && bin <= m_hists.at(0)->GetNbinsX()){
      double val = m_hists.at(0)->GetBinContent(bin);
      w = w * val;
    }
    
  }
  if (m_weightvect){
    if ( m_vect.count(val) == 1 ) { 
      w = w * m_vect.at(val);
    }
  }
  return w;
}

string ReweightVar::GetName(){
  return m_names.at(0);
}
vector<string> ReweightVar::GetNames(){
  return m_names;
}
vector<Expr*> ReweightVar::GetExprs(){
  return m_exprs;
}
string ReweightVar::GetWeightName(){
  return m_weightname;
}


Expr* ReweightVar::GetExpr(){
  return m_exprs.at(0);
}

//bool ReweightVar::UseLeaf(){
//  return m_weightleaf;
//}
