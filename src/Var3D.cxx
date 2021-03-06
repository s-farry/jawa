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
#include <Var3D.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

Var3D::Var3D(string name) : JawaObj("Var3D", name){}

Var3D::Var3D(string name , string var1 , int bins1 , double lo1 , double hi1 ,  string var2 , int bins2 , double lo2 , double hi2, string var3, int bins3, double lo3, double hi3, string prefix) : JawaObj("Var3D", name){
  Var* fVar1 = new Var(var1, var1, bins1, lo1, hi1, prefix);
  Var* fVar2 = new Var(var2, var2, bins2, lo2, hi2, prefix);
  Var* fVar3 = new Var(var3, var3, bins3, lo3, hi3, prefix);
  Var3D(name, fVar1, fVar2, fVar3, prefix);
}

/*Var3D::~Var3D(){
  m_hist->Delete();
  m_hist_fwd->Delete();
  m_hist_bwd->Delete();
  }*/

Var3D::Var3D(string name , Var* var1 , Var* var2, Var* var3, string prefix) : JawaObj("Var3D",name) {
  m_name1   = var1->GetName();
  m_name2   = var2->GetName();
  m_name3   = var3->GetName();
  m_varname1  = var1->GetVar();
  m_varname2  = var2->GetVar();
  m_varname3  = var3->GetVar();
  m_nbins1  = var1->GetBins();
  m_nbins2  = var2->GetBins();
  m_nbins3  = var3->GetBins();
  m_lo1     = var1->GetLo();
  m_hi1     = var1->GetHi();
  m_lo2     = var2->GetLo();
  m_hi2     = var2->GetHi();
  m_lo3     = var3->GetLo();
  m_hi3     = var3->GetHi();

  m_Var1 = var1;
  m_Var2 = var2;
  m_Var3 = var3;

  m_prefix = prefix;

  if (prefix != "") prefix = prefix+"/";

  //double histbinedges1[m_nbins1+1];
  //double fwdbinedges1[m_nbins1+1];
  //double bwdbinedges1[m_nbins1+1];

  //double histbinedges2[m_nbins2+1];
  //double fwdbinedges2[m_nbins2+1];
  //double bwdbinedges2[m_nbins2+1];
  
  std::vector<double> histbinedges1 = Var::GetBinEdges(var1->GetHist());
  std::vector<double> histbinedges2 = Var::GetBinEdges(var2->GetHist());
  std::vector<double> histbinedges3 = Var::GetBinEdges(var3->GetHist());

  std::pair<std::vector<string> , std::vector<double> > histbinedges_comb = CombineBinEdges(histbinedges1, histbinedges2);

  m_hist = new TH3F((prefix+name).c_str() , name.c_str() , m_nbins1, &histbinedges1[0], m_nbins2, &histbinedges2[0], m_nbins3, &histbinedges3[0]);
  m_prof = new TProfile((prefix+name+"_prof_"+m_name1+"_"+m_name2).c_str() , name.c_str() , m_nbins1, &histbinedges1[0]);
  m_prof2 = new TProfile((prefix+name+"_prof_"+m_name2+"_"+m_name1).c_str() , name.c_str() , m_nbins2, &histbinedges2[0]);
  m_prof3 = new TProfile((prefix+name+"_prof_"+m_name1+"_"+m_name3).c_str() , name.c_str() , m_nbins1, &histbinedges1[0]);
  m_prof4 = new TProfile((prefix+name+"_prof_"+m_name2+"_"+m_name3).c_str() , name.c_str() , m_nbins2, &histbinedges2[0]);
  m_prof5 = new TProfile((prefix+name+"_prof_"+m_name3+"_"+m_name1).c_str() , name.c_str() , m_nbins3, &histbinedges3[0]);
  m_prof6 = new TProfile((prefix+name+"_prof_"+m_name3+"_"+m_name2).c_str() , name.c_str() , m_nbins3, &histbinedges3[0]);
  m_hist_fwd = new TH3F((prefix+name+"_fwd").c_str() , name.c_str() , m_nbins1, &histbinedges1[0], m_nbins2, &histbinedges2[0], m_nbins3, &histbinedges3[0]);
  m_hist_bwd = new TH3F((prefix+name+"_bwd").c_str() , name.c_str() , m_nbins1, &histbinedges1[0], m_nbins2, &histbinedges2[0], m_nbins3, &histbinedges3[0]);
  m_hist_comb = new TH1F((prefix+name+"_comb").c_str(), name.c_str() , m_nbins1 + m_nbins2, &histbinedges_comb.second[0]);

}
Var3D::Var3D(string name , Var3D* varA , Var3D* varB, string prefix) : JawaObj("Var3D",name){
  if ( varA->GetName() == varB->GetName() && 
       varA->GetName1() == varB->GetName1() && 
       varA->GetName2() == varB->GetName2() && 
       varA->GetVar1Hi() == varB->GetVar1Hi() &&
       varA->GetVar2Hi() == varB->GetVar2Hi() &&
       varA->GetVar3Hi() == varB->GetVar3Hi() &&
       varA->GetVar1Lo() == varB->GetVar1Lo() &&
       varA->GetVar2Lo() == varB->GetVar2Lo() &&
       varA->GetVar3Lo() == varB->GetVar3Lo() &&
       varA->GetBins1() == varB->GetBins1() &&
       varA->GetBins2() == varB->GetBins2() &&
       varA->GetBins3() == varB->GetBins3() )
    {
      m_name1   = varA->GetName1();
      m_name2   = varA->GetName2();
      m_name3   = varA->GetName3();
      m_varname1  = varA->GetVar1()->GetVar();
      m_varname2  = varA->GetVar2()->GetVar();
      m_varname3  = varA->GetVar3()->GetVar();
      m_nbins1  = varA->GetBins1();
      m_nbins2  = varA->GetBins2();
      m_nbins3  = varA->GetBins3();
      m_lo1     = varA->GetVar1Lo();
      m_hi1     = varA->GetVar1Hi();
      m_lo2     = varA->GetVar2Lo();
      m_hi2     = varA->GetVar2Hi();
      m_lo3     = varA->GetVar3Lo();
      m_hi3     = varA->GetVar3Hi();
      
      std::vector<double> histbinedges1 = Var::GetBinEdges(varA->GetVar1()->GetHist());
      std::vector<double> histbinedges2 = Var::GetBinEdges(varA->GetVar2()->GetHist());
      std::vector<double> histbinedges3 = Var::GetBinEdges(varA->GetVar3()->GetHist());
      
      //m_Var1 = new Var(m_name1, m_varname1, histbinedges1, prefix);
      //m_Var2 = new Var(m_name2, m_varname2, histbinedges2, prefix);
      
      m_prefix = prefix;

      if (prefix != "") prefix = prefix+"/";
      m_hist = new TH3F((prefix+name).c_str() , name.c_str() , m_nbins1 , &histbinedges1[0] , m_nbins2, &histbinedges2[0], m_nbins3, &histbinedges3[0]);
      
      if (prefix != "") prefix = prefix+"/";
      
      m_hist->Sumw2();

      m_hist->Add(varA->GetHist());
      m_hist->Add(varB->GetHist());
    }
}

void Var3D::FillHist(double val1, double val2, double val3){
  FillHist(val1, val2, val3, 1.0);
}
void Var3D::NormaliseToEvts(double evts){
  m_hist->Scale(evts/m_hist->Integral());
  m_hist_comb->Scale(evts/m_hist_comb->Integral());
}

void Var3D::NormaliseToMC(double xsec, double acc, double lumi, double nEvts){
  double scale = xsec * lumi * acc / nEvts;
  m_hist->Scale(scale);
  m_hist_comb->Scale(scale);
}

void Var3D::Scale(double scale){
  m_hist->Scale(scale);
  m_hist_comb->Scale(scale);
}

void Var3D::FillHist(double val1, double val2, double val3, double w){
  m_hist->Fill(val1, val2, val3, w);
  m_prof->Fill(val1, val2, w);
  m_prof2->Fill(val2, val1, w);
  m_prof3->Fill(val1, val3, w);
  m_prof4->Fill(val2, val3, w);
  m_prof3->Fill(val3, val1, w);
  m_prof4->Fill(val3, val2, w);
  if (val1 < m_hi1)  m_hist_comb->Fill(val1, w);
  if (val2 >= m_lo2)  m_hist_comb->Fill(val2 + m_hi1, w);
}
string Var3D::GetName(){
  return m_name;
}
string Var3D::GetName1(){
  return m_name1;
}
string Var3D::GetName2(){
  return m_name2;
}
string Var3D::GetName3(){
  return m_name3;
}
string Var3D::GetVarName1(){
  return m_varname1;
}
string Var3D::GetVarName2(){
  return m_varname2;
}
string Var3D::GetVarName3(){
  return m_varname3;
}


int Var3D::GetBins1(){
  return m_nbins1;
}

double Var3D::GetVar1Lo(){
  return m_lo1;
}
double Var3D::GetVar1Hi(){
  return m_hi1;
}
int Var3D::GetBins2(){
  return m_nbins2;
}
double Var3D::GetVar2Lo(){
  return m_lo2;
}
double Var3D::GetVar2Hi(){
  return m_hi2;
}
int Var3D::GetBins3(){
  return m_nbins3;
}
double Var3D::GetVar3Lo(){
  return m_lo3;
}
double Var3D::GetVar3Hi(){
  return m_hi3;
}

Var* Var3D::GetVar1(){
  return m_Var1;
}
Var* Var3D::GetVar2(){
  return m_Var2;
}
Var* Var3D::GetVar3(){
  return m_Var3;
}


TH3F* Var3D::GetHist(){
  return m_hist;
}
TProfile* Var3D::GetProfile(){
  return m_prof;
}
TProfile* Var3D::GetProfile2(){
  return m_prof2;
}
TH1F* Var3D::GetCombHist(){
  return m_hist_comb;
}

TH3F* Var3D::GetFwdHist(){
  return m_hist_fwd;
}
TH3F* Var3D::GetBwdHist(){
  return m_hist_bwd;
}


/*void Var3D::NormaliseHist(bool doX){
  int xbins = m_hist->GetNbinsX();
  int ybins = m_hist->GetNbinsY();
  int zbins = m_hist->GetNbinsZ();
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
  }*/
void Var3D::FillAsymmetry(double xval, double yval, double eta1, double eta2, double w){
  if (eta2 > eta1) m_hist_fwd->Fill(xval,yval,w);
  else m_hist_bwd->Fill(xval,yval,w);
}

#ifdef WITHPYTHON
PyObject* Var3D::GetHist_py(){
  TH3F* newCxxObj = new TH3F(*m_hist);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

PyObject* Var3D::GetProfile_py(){
  TProfile* newCxxObj = new TProfile(*m_prof);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
PyObject* Var3D::GetProfile2_py(){
  TProfile* newCxxObj = new TProfile(*m_prof2);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
#endif
