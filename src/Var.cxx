#include <iostream>
#include <sstream>
#include <iomanip>
#include <TCut.h>
#include <TMath.h>
#include <TObjArray.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>
#include <Var.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>
#include "TRandom3.h"

using namespace std;

Var::Var(string name) : JawaObj("Var", name){}

Var::Var(string name, Var* varA, Var* varB, string prefix) : JawaObj("Var", name) {
  if ( varA && varB && varA->GetName() == varB->GetName() && 
       varA->GetHi() == varB->GetHi() && varA->GetLo() == varB->GetLo() &&
       varA->GetBins() == varB->GetBins() )
    {
      m_name               = varA->GetName();
      m_lo                 = varA->GetLo();
      m_hi                 = varA->GetHi();
      m_nbins              = varA->GetBins();
      m_edges              = varA->GetEdges();
      m_expr   = new Expr(varA->GetExpr()->GetExpr());
      m_scale  = false;
      //No extra variables
      m_multivar = false;
      m_prefix = prefix;
      
      if (prefix != "") prefix = prefix+"/";
      m_hist       = new TH1F((prefix+name).c_str() ,          name.c_str() , m_nbins, &m_edges[0]);
      m_hist_fwd   = new TH1F((prefix+name+"_fwd").c_str() ,   name.c_str() , m_nbins, &m_edges[0]);
      m_hist_bwd   = new TH1F((prefix+name+"_bwd").c_str() ,   name.c_str() , m_nbins, &m_edges[0]);
      m_hist_asymm = new TH1F((prefix+name+"_asymm").c_str() , name.c_str() , m_nbins, &m_edges[0]);
  
      m_hist->Sumw2();
      TH1F* histA = varA->GetHist();
      TH1F* histB = varB->GetHist();
      for (int i = 0 ; i < m_nbins + 2; ++i){
	m_hist->SetBinContent(i, histA->GetBinContent(i) + histB->GetBinContent(i));
	m_hist->SetBinError(i, sqrt( pow(histA->GetBinError(i), 2) + pow( histB->GetBinError(i), 2) ) );
      }

      //m_hist->Add(varA->GetHist());
      //m_hist->Add(varB->GetHist());
    }
}
Var::Var(string name, Var* v, string prefix) : JawaObj("Var", name){
  m_lo                 = v->GetLo();
  m_hi                 = v->GetHi();
  m_nbins              = v->GetBins();
  m_edges              = v->GetEdges();
  m_expr   = new Expr(v->GetExpr()->GetExpr());
  m_scale  = false;
  //No extra variables
  m_multivar = false;
  m_prefix = prefix;
  
  TH1F* hist = v->GetHist();
  
  if (prefix != "") prefix = prefix+"/";
  m_hist       = ((TH1F*)hist->Clone((prefix+name).c_str()));
  m_hist->Sumw2();
}

void Var::Add2ndVar(string name){
  m_multivar = true;
  m_expr2 = new Expr(name);
}

string Var::Get2ndVar(){
  string output = m_expr2 ? m_expr2->GetExpr() : "";
  return output;
}

bool Var::hasMultiVar(){
  return m_multivar;
}

void Var::SetVarName(string var){
  m_expr = new Expr(var);
}

void Var::SetScaleFactor(double scale){
  m_scale = true;
  m_sf = scale;
}


void Var::FillHist(double val){
  m_hist->Fill(val);
}
void Var::NormaliseToEvts(double evts){
  m_hist->Scale(evts/m_hist->Integral());
}
void Var::NormaliseToMC(double xsec, double acc, double lumi, double nEvts){
  double scale = xsec * lumi * acc / nEvts;
  m_hist->Scale(scale);

}
void Var::Scale(double scale){
  m_hist->Scale(scale);
}
void Var::FillHist(double val, double w){
  m_hist->Fill(val, w);
}

void Var::FillAsymmetry(double val, double eta1, double eta2, double w){
  if (eta2 > eta1) m_hist_fwd->Fill(val,w);
  else m_hist_bwd->Fill(val,w);
}

string Var::GetVar(){
  string output = m_expr ? m_expr->GetExpr() : "";
  return output;
}
TH1F* Var::GetHist(){
  return m_hist;
}
TH1F* Var::GetAsymmHist(){
  return m_hist_asymm;
}
TH1F* Var::GetAsymmFwdHist(){
  return m_hist_fwd;
}
TH1F* Var::GetAsymmBwdHist(){
  return m_hist_bwd;
}
double Var::GetLo(){
  return m_lo;
}
double Var::GetHi(){
  return m_hi;
}
int Var::GetBins(){
  return m_nbins;
}

void Var::MakeAsymmetry(){
  for (int i = 0; i<m_hist->GetNbinsX(); ++i){
    double fwd = m_hist_fwd->GetBinContent(i+1);
    double bwd = m_hist_bwd->GetBinContent(i+1);
    m_hist_asymm->SetBinContent(i+1,(fwd - bwd)/(fwd + bwd));
  }
}

bool Var::IsScaled(){
  return m_scale;
}

double Var::GetScaleFactor(){
  return m_sf;
}

void Var::ResetHist(){
  if (m_vals.size() > 0){
    for (int i = 0; i < m_hist->GetXaxis()->GetNbins(); ++i){
      double val = m_vals.at(i).first;
      double err = m_vals.at(i).second;
      m_hist->SetBinContent(i+1, val);
      m_hist->SetBinError(i+1, err);
    }
  }
  else{
    info("Values not backed up: Shouldn't need to be reset");
  }


}

void Var::SmearHist(){
  //Smear histogram
  //If values are not already set, set them
  if (m_vals.size() == 0){
    for (int i = 0; i < m_hist->GetXaxis()->GetNbins(); ++i){
      double val = m_hist->GetBinContent(i + 1);
      double err = m_hist->GetBinError(i + 1);
      m_vals.push_back(pair<double,double>(val,err));
      //m_hist_fix->SetBinContent(i, val);
      //m_hist_fix->SetBinError(i, err);
      //m_hist_fix->SetEntries(m_hist->GetEntries());
    }
  }


  TRandom3 r(0);
  for (int i = 0; i < m_hist->GetXaxis()->GetNbins(); ++i){
    double val = m_vals.at(i).first;
    double err = m_vals.at(i).second;
    double delta = r.Gaus(0,err);
    m_hist->SetBinContent(i+1, val + delta);
    m_hist->SetBinError(i+1, err * (val+delta)/val);
  }
}

/*
DataType::DataType(double dbl){
  d = dbl;
  i = 0;
  f = 0;
  type = "Double";
  }*/

std::pair< std::vector<string>, std::vector<double> > CombineBinEdges(std::vector<double> edges1, std::vector<double> edges2){
  std::vector<double> edges;
  std::vector<double> combEdges;
  std::vector<string> combLabels;

  double end = edges1[edges1.size()-1];
  
  for (unsigned int i = 0; i <edges1.size() ; ++i){
    combEdges.push_back(edges1[i]);
  }
  for (unsigned int j = 1; j <edges2.size() ; ++j){
    //Add end value of first variable + start value of second (in case of negatives)
    if (edges2[0] < 0){
      combEdges.push_back(edges2[j] + end - edges2[0]);
    }
    else combEdges.push_back(edges2[j] + end );
  }
  return std::pair< std::vector<string> , std::vector<double> > ( combLabels, combEdges);
}

Var::Var(string name , string varexp , int bins , double lo , double hi , string prefix) : JawaObj("Var", name){
  m_expr   = new Expr(varexp);
  m_nbins  = bins;
  m_lo     = lo;
  m_hi     = hi;
  m_prefix = prefix;
  m_scale  = false;
  //No extra variables
  m_multivar = false;
  
  if (prefix != "") prefix = prefix+"/";
  //one fixed hist
  m_hist     = new TH1F((prefix+name).c_str() , name.c_str() , bins , lo , hi);
  m_edges = GetBinEdges(m_hist);
  
  m_hist_fwd   = new TH1F((prefix+name+"_fwd").c_str() , name.c_str() , bins , lo , hi);
  m_hist_bwd   = new TH1F((prefix+name+"_bwd").c_str() , name.c_str() , bins , lo , hi);
  m_hist_asymm = new TH1F((prefix+name+"_asymm").c_str() , name.c_str() , bins , lo , hi);
  
  m_hist->Sumw2();

}

Var::Var(string name , string varexp , std::vector<double>& edges , string prefix) : JawaObj("Var", name){
  m_expr    = new Expr(varexp);
  m_nbins  = edges.size() - 1;
  m_lo     = edges[0];
  m_hi     = edges[edges.size()-1];
  m_prefix = prefix;
  m_scale = false;
  m_edges = edges;

  if (prefix != "") prefix = prefix+"/";
  m_hist = new TH1F((prefix+name).c_str() , name.c_str() , edges.size() - 1 , &edges[0]);
  m_hist_fwd   = new TH1F((prefix+name+"_fwd").c_str()   , name.c_str() , edges.size() - 1 , &edges[0]);
  m_hist_bwd   = new TH1F((prefix+name+"_bwd").c_str()   , name.c_str() , edges.size() - 1 , &edges[0]);
  m_hist_asymm = new TH1F((prefix+name+"_asymm").c_str() , name.c_str() , edges.size() - 1 , &edges[0]);
  m_hist->Sumw2();

}

Expr* Var::GetExpr(){ return m_expr;}

double Var::GetVal(std::vector<double>& input){
  return m_expr->GetVal(input);
}

int Var::FindBin(double val){
  //Will return bin - 0 is underflow i+1 is overflow
  int bin = -1;
  if (val < m_edges.at(0)) bin = 0;
  for (unsigned int i = 0; i < m_edges.size() -1; ++i){
    if (val >= m_edges.at(i) && val < m_edges.at(i+1)) bin = i+1;
  }
  if (val >= m_edges.at(m_edges.size() - 1)) bin = m_edges.size();
  return bin;
}


std::vector<string> Var::GetVarNames(){ return m_expr->GetVarNames();}

void Var::SetVarExp(string varexp){
  m_expr = new Expr(varexp);
}

std::vector<double> Var::GetEdges(){ return m_edges;}


void Var::Draw(const char* options){
  m_hist->Draw(options);
}

std::vector<double> Var::GetBinEdges(TGraphAsymmErrors* graph){
  double x = 0.0, y = 0.0;
  std::vector<double> edges;
  int i = 0;
  double xlo = 0.0, xhi = 0.0;
  Int_t hasPoint = graph->GetPoint(i,x,y);
  while (hasPoint == i){
    xlo = graph->GetErrorXlow(i);
    xhi = graph->GetErrorXhigh(i);
    edges.push_back(x - xlo);
    i++;
    hasPoint = graph->GetPoint(i,x,y);
  }
  edges.push_back(x+xhi);
  return edges;
}


bool is_number(const std::string& s){
  return (strspn( s.c_str(), "-.01234567890") == s.size() );
}
//python files
#ifdef WITHPYTHON

PyObject* Var::GetHist_py(){
  TH1F* newCxxObj = new TH1F(*m_hist);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
double Var::GetVal_py(boost::python::list& input){
  vector<double> dbl_vec;
  for (unsigned int i = 0; i < len(input); ++i){
    double d = boost::python::extract<double>(input[i]);
    dbl_vec.push_back(d);
  }
  return GetVal(dbl_vec);
}
boost::python::list Var::GetVarNames_py(){ 
  boost::python::list l;
  vector<string> varnames = GetVarNames();
  //for (vector<string>::iterator is = varnames.begin(); is != varnames.end(); ++is){
  //  l.append((*is));
  //}
  for (auto is : varnames) l.append(is);
  return l;
}
boost::python::list Var::GetEdges_py(){
  vector<double> edges = m_edges;
  boost::python::object get_iter = boost::python::iterator<std::vector<double> >();
  boost::python::object iter = get_iter(edges);
  boost::python::list l(iter);
  return l;
}
void Var::Draw1_py(){ Draw();}
void Var::Draw2_py(const char* options) {Draw(options);}

std::vector<double> Var::GetBinEdges(TH1F* hist){
  int nbins = hist->GetNbinsX();
  std::vector<double> edges;
  for (int i = 0; i <nbins+1 ; ++i){
    edges.push_back(hist->GetBinLowEdge(i+1));
  }
  return edges;
}
boost::python::list Var::GetBinEdges_py(PyObject* pyObj){
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  vector<double> edges;
  if (strcmp(hist->ClassName(), "TH1F") == 0) {
    edges = GetBinEdges(hist);
  }
  else{
    edges = GetBinEdges(graph);
  }
  boost::python::object get_iter = boost::python::iterator<std::vector<double> >();
  boost::python::object iter = get_iter(edges);
  boost::python::list l(iter);
  return l;
}


#endif
