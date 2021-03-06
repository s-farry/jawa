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
#include <Var.h>
#include <Var2D.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

Var2D::Var2D(string name) : JawaObj("Var2D", name){}

Var2D::Var2D(string name , string var1 , int bins1 , double lo1 , double hi1 ,  string var2 , int bins2 , double lo2 , double hi2, string prefix) : JawaObj("Var2D", name){
  m_name1   = var1;
  m_name2   = var2;
  m_varname1    = var1;
  m_varname2    = var2;
  m_nbins1  = bins1;
  m_nbins2  = bins2;
  m_lo1     = lo1;
  m_hi1     = hi1;
  m_lo2     = lo2;
  m_hi2     = hi2;
  m_Var1 = new Var(m_name1, m_varname1, m_nbins1, m_lo1, m_hi1, prefix);
  m_Var2 = new Var(m_name2, m_varname2, m_nbins2, m_lo2, m_hi2, prefix);

  m_prefix = prefix;

  if (prefix != "") prefix = prefix+"/";
  m_hist = new TH2F((prefix+name).c_str() , name.c_str() , bins1 , lo1 , hi1 , bins2, lo2, hi2);
  m_prof = new TProfile((prefix+name+"_prof_"+m_name1).c_str(), name.c_str(), bins1, lo1, hi1);
  m_prof2 = new TProfile((prefix+name+"_prof_"+m_name2).c_str(), name.c_str(), bins2, lo2, hi2);
  m_hist->Sumw2();
  m_prof->Sumw2();
  m_prof2->Sumw2();
}

/*Var2D::~Var2D(){
  m_hist->Delete();
  m_hist_fwd->Delete();
  m_hist_bwd->Delete();
  }*/


Var2D::Var2D(string name , Var* var1 , Var* var2, string prefix) : JawaObj("Var2D", name){
  m_name1   = var1->GetName();
  m_name2   = var2->GetName();
  m_varname1  = var1->GetVar();
  m_varname2  = var2->GetVar();
  m_nbins1  = var1->GetBins();
  m_nbins2  = var2->GetBins();
  m_lo1     = var1->GetLo();
  m_hi1     = var1->GetHi();
  m_lo2     = var2->GetLo();
  m_hi2     = var2->GetHi();

  m_Var1 = var1;
  m_Var2 = var2;

  m_prefix = prefix;

  if (prefix != "") prefix = prefix+"/";

  //double histbinedges1[m_nbins1+1];
  //double fwdbinedges1[m_nbins1+1];
  //double bwdbinedges1[m_nbins1+1];

  //double histbinedges2[m_nbins2+1];
  //double fwdbinedges2[m_nbins2+1];
  //double bwdbinedges2[m_nbins2+1];
  
  m_histbinedges1 = Var::GetBinEdges(var1->GetHist());
  m_histbinedges2 = Var::GetBinEdges(var2->GetHist());

  std::pair<std::vector<string> , std::vector<double> > histbinedges_comb = CombineBinEdges(m_histbinedges1, m_histbinedges2);

  m_hist = new TH2F((prefix+name).c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
  m_prof = new TProfile((prefix+name+"_prof_"+m_name1).c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0]);
  m_prof2 = new TProfile((prefix+name+"_prof_"+m_name2).c_str() , name.c_str() , m_nbins2, &m_histbinedges2[0]);
  m_hist_fwd = new TH2F((prefix+name+"_fwd").c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
  m_hist_bwd = new TH2F((prefix+name+"_bwd").c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
  m_hist_comb = new TH1F((prefix+name+"_comb").c_str(), name.c_str() , m_nbins1 + m_nbins2, &histbinedges_comb.second[0]);
  m_hist->Sumw2();
  m_prof->Sumw2();
  m_prof2->Sumw2();


}

Var2D::Var2D(string name , Var2D* varA , Var2D* varB, string prefix) : JawaObj("Var2D",name){
  if ( varA->GetName() == varB->GetName() && 
       varA->GetName1() == varB->GetName1() && 
       varA->GetName2() == varB->GetName2() && 
       varA->GetVar1Hi() == varB->GetVar1Hi() &&
       varA->GetVar2Hi() == varB->GetVar2Hi() &&
       varA->GetVar1Lo() == varB->GetVar1Lo() &&
       varA->GetVar2Lo() == varB->GetVar2Lo() &&
       varA->GetBins1() == varB->GetBins1() &&
       varA->GetBins2() == varB->GetBins2() )
    {
      m_name1   = varA->GetName1();
      m_name2   = varA->GetName2();
      m_varname1  = varA->GetVarName1();
      m_varname2  = varA->GetVarName2();
      m_nbins1  = varA->GetBins1();
      m_nbins2  = varA->GetBins2();
      m_lo1     = varA->GetVar1Lo();
      m_hi1     = varA->GetVar1Hi();
      m_lo2     = varA->GetVar2Lo();
      m_hi2     = varA->GetVar2Hi();
      
      //m_histbinedges1 = Var::GetBinEdges(varA->GetVar1()->GetHist());
      //m_histbinedges2 = Var::GetBinEdges(varA->GetVar2()->GetHist());
      
      m_histbinedges1 = varA->GetBinEdges1();
      m_histbinedges2 = varA->GetBinEdges2();
      
      //m_Var1 = new Var(m_name1, m_varname1, m_histbinedges1, prefix);
      //m_Var2 = new Var(m_name2, m_varname2, m_histbinedges2, prefix);
      
      m_prefix = prefix;

      if (prefix != "") prefix = prefix+"/";
      m_hist = new TH2F((prefix+name).c_str() , name.c_str() , m_nbins1 , &m_histbinedges1[0] , m_nbins2, &m_histbinedges2[0]);
      m_prof = new TProfile((prefix+name+"_prof_"+m_name1).c_str(), name.c_str(), m_nbins1, &m_histbinedges1[0]);
      m_prof2 = new TProfile((prefix+name+"_prof_"+m_name2).c_str(), name.c_str(), m_nbins2, &m_histbinedges2[0]);
      
      if (prefix != "") prefix = prefix+"/";
      
      //double m_histbinedges1[m_nbins1+1];
      //double fwdbinedges1[m_nbins1+1];
      //double bwdbinedges1[m_nbins1+1];
      
      //double m_histbinedges2[m_nbins2+1];
      //double fwdbinedges2[m_nbins2+1];
      //double bwdbinedges2[m_nbins2+1];
      
      std::pair<std::vector<string> , std::vector<double> > histbinedges_comb = CombineBinEdges(m_histbinedges1, m_histbinedges2);
      
      //m_hist = new TH2F((prefix+name).c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
      //m_prof = new TProfile((prefix+name+"_prof_"+m_name1).c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0]);
      //m_prof2 = new TProfile((prefix+name+"_prof_"+m_name2).c_str() , name.c_str() , m_nbins2, &m_histbinedges2[0]);
      //m_hist_fwd = new TH2F((prefix+name+"_fwd").c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
      //m_hist_bwd = new TH2F((prefix+name+"_bwd").c_str() , name.c_str() , m_nbins1, &m_histbinedges1[0], m_nbins2, &m_histbinedges2[0]);
      m_hist_comb = new TH1F((prefix+name+"_comb").c_str(), name.c_str() , m_nbins1 + m_nbins2, &histbinedges_comb.second[0]);
      m_hist->Sumw2();
      m_prof->Sumw2();
      m_prof2->Sumw2();

      m_hist->Add(varA->GetHist());
      m_hist->Add(varB->GetHist());
    }
}
Var2D::Var2D(string name , Var2D* v, string prefix) : JawaObj("Var2D",name){
      m_name1   = v->GetName1();
      m_name2   = v->GetName2();
      m_varname1  = v->GetVar1()->GetVar();
      m_varname2  = v->GetVar2()->GetVar();
      m_nbins1  = v->GetBins1();
      m_nbins2  = v->GetBins2();
      m_lo1     = v->GetVar1Lo();
      m_hi1     = v->GetVar1Hi();
      m_lo2     = v->GetVar2Lo();
      m_hi2     = v->GetVar2Hi();
      
      m_prefix = prefix;

      if (prefix != "") prefix = prefix+"/";
      m_hist  = (TH2F*)v->GetHist()->Clone((prefix+name).c_str());
      m_prof  = (TProfile*)v->GetProfile()->Clone((prefix+name+"_prof_"+m_name1).c_str());
      m_prof2 = (TProfile*)v->GetProfile2()->Clone((prefix+name+"_prof_"+m_name2).c_str());
      m_hist_comb = (TH1F*)v->GetCombHist()->Clone((prefix+name+"_comb").c_str());
      m_hist->Sumw2();
      m_prof->Sumw2();
      m_prof2->Sumw2();
}

void Var2D::FillHist(double val1, double val2){
  FillHist(val1, val2, 1.0);
}
void Var2D::NormaliseToEvts(double evts){
  m_hist->Scale(evts/m_hist->Integral());
  m_hist_comb->Scale(evts/m_hist_comb->Integral());
}

void Var2D::NormaliseToMC(double xsec, double acc, double lumi, double nEvts){
  double scale = xsec * lumi * acc / nEvts;
  m_hist->Scale(scale);
  m_hist_comb->Scale(scale);
}

void Var2D::Scale(double scale){
  m_hist->Scale(scale);
  m_hist_comb->Scale(scale);
}

void Var2D::FillHist(double val1, double val2, double w){
  m_hist->Fill(val1, val2, w);
  m_prof->Fill(val1, val2, w);
  m_prof2->Fill(val2, val1, w);
  if (val1 < m_hi1)  m_hist_comb->Fill(val1, w);
  if (val2 >= m_lo2)  m_hist_comb->Fill(val2 + m_hi1, w);
}

string Var2D::GetName1(){
  return m_name1;
}
string Var2D::GetName2(){
  return m_name2;
}
string Var2D::GetVarName1(){
  return m_varname1;
}
string Var2D::GetVarName2(){
  return m_varname2;
}


int Var2D::GetBins1(){
  return m_nbins1;
}

double Var2D::GetVar1Lo(){
  return m_lo1;
}
double Var2D::GetVar1Hi(){
  return m_hi1;
}
int Var2D::GetBins2(){
  return m_nbins2;
}
double Var2D::GetVar2Lo(){
  return m_lo2;
}
double Var2D::GetVar2Hi(){
  return m_hi2;
}

Var* Var2D::GetVar1(){
  return m_Var1;
}
Var* Var2D::GetVar2(){
  return m_Var2;
}


TH2F* Var2D::GetHist(){
  return m_hist;
}
TProfile* Var2D::GetProfile(){
  return m_prof;
}
TProfile* Var2D::GetProfile2(){
  return m_prof2;
}
TH1F* Var2D::GetCombHist(){
  return m_hist_comb;
}

TH2F* Var2D::GetFwdHist(){
  return m_hist_fwd;
}
TH2F* Var2D::GetBwdHist(){
  return m_hist_bwd;
}

std::vector<double> Var2D::GetBinEdges1(){
  return m_histbinedges1;
}
std::vector<double> Var2D::GetBinEdges2(){
  return m_histbinedges2;
}



void Var2D::NormaliseHist(bool doX){
  int xbins = m_hist->GetNbinsX();
  int ybins = m_hist->GetNbinsY();
  if (!doX){
    for (int i = 0; i < ybins ; ++i)
      {
	double integral = m_hist->Integral(1,xbins+1,i+1,i+1);
	    for (int j = 0; j < xbins ; ++j)
	      {
		double val = m_hist->GetBinContent(j+1,i+1);
		m_hist->SetBinContent(j+1,i+1,(double)val/integral);
	      }
      }
  }
}
void Var2D::FillAsymmetry(double xval, double yval, double eta1, double eta2, double w){
  if (eta2 > eta1) m_hist_fwd->Fill(xval,yval,w);
  else m_hist_bwd->Fill(xval,yval,w);
}

#ifdef WITHPYTHON
PyObject* Var2D::GetHist_py(){
  TH2F* newCxxObj = new TH2F(*m_hist);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

PyObject* Var2D::GetProfile_py(){
  TProfile* newCxxObj = new TProfile(*m_prof);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
PyObject* Var2D::GetProfile2_py(){
  TProfile* newCxxObj = new TProfile(*m_prof2);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
#endif
