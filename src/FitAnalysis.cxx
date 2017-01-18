#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <TFile.h>
#include <Fit.h>
#include <FitAnalysis.h>
#include <boost/algorithm/string.hpp>

using namespace std;

FitAnalysis::FitAnalysis(string name, string tofit, string func) : JawaObj("FitAnalysis", name){
  m_func = func;
  m_tofit = tofit;
  m_expr = new Expr(func);
  int j = 0;
  for (unsigned int i = 0 ; i < m_expr->GetVarNames().size(); ++i) {
    if (m_expr->GetVarNames().at(i) != "x"){
      m_paridx[m_expr->GetVarNames().at(i)] = j;
      j++;
    }
  }
}

void FitAnalysis::SetVal(string par, double val){
  m_initvals[par] = val;
}
void FitAnalysis::FixVal(string par, double val){
  m_fixvals[par] = val;
}
void FitAnalysis::SetRange(string par, double lo, double hi){
  m_initrange[par] = std::pair<double, double>(lo, hi);
}

void FitAnalysis::Init(Template* t){
  info()<<"Initialising with template "<<t->GetName()<<endl;
  std::map<string, Var2D*> varmap = t->Get2DVariables();
  for (std::map<string, Var2D*>::iterator iv = varmap.begin(); iv != varmap.end(); ++iv){
    if (m_tofit != (*iv).second->GetName2()) continue;
    std::vector<Fit*> fits;
    TH2F* h = (TH2F*)(*iv).second->GetHist();
    m_fits[(*iv).second->GetName1()]  = vector<Fit*>(0);
    m_hists[(*iv).second->GetName1()] = vector<TH1F*>(0);
    m_parhists[(*iv).second->GetName1()] = vector<TH1F*>(0);
    
    for (int i = 0 ; i < h->GetNbinsX(); ++i){
      Fit* f = new Fit(m_name+(*iv).second->GetName1()+"_", m_func);
      
      for (std::map<string, double>::iterator ip = m_initvals.begin(); ip != m_initvals.end(); ++ip) 
	{
	  f->SetParameter(m_paridx[(*ip).first], (*ip).second);
	}
      for (std::map<string, double>::iterator ip = m_fixvals.begin(); ip != m_fixvals.end(); ++ip) 
	{
	  info()<<"fixing "<<(*ip).first<<"( "<<m_paridx[(*ip).first]<<" ) to "<<(*ip).second<<endl;;
	  f->FixParameter(m_paridx[(*ip).first], (*ip).second);
	}
      ostringstream s;
      s<<m_name<<"_"<<(*iv).second->GetName1()<<"_"<<i;
      m_fits[(*iv).second->GetName1()].push_back(f);
      m_hists[(*iv).second->GetName1()].push_back((TH1F*)h->ProjectionY(s.str().c_str(), i+1,i+1));
    }
    for (unsigned int i = 0 ; i < m_expr->GetVarNames().size(); ++i){
      if (m_expr->GetVarNames().at(i) == "x") continue;
      ostringstream s;
      s<<m_name<<"_"<<m_expr->GetVarNames().at(i)<<"_"<<(*iv).second->GetName1();
      vector<double> binedges = Utils::GetBinEdgesX(h);
      m_parhists[(*iv).second->GetName1()].push_back(new TH1F(s.str().c_str(), s.str().c_str(), binedges.size()-1, &binedges[0]));
    }
  }
}

void FitAnalysis::FitIt(string opt){
  std::map<string, vector<TH1F*> >::iterator ih;
  for (ih = m_hists.begin(); ih != m_hists.end(); ++ih){
    vector<TH1F*> hists = (*ih).second;
    for (unsigned int i = 0 ; i < hists.size() ; ++i){
      m_fits[(*ih).first].at(i)->FitHist(hists.at(i), opt);
      int k = 0;
      for (unsigned int j = 0 ; j < m_expr->GetVarNames().size(); ++j){
	if (m_expr->GetVarNames().at(j) == "x") continue;
	double par = m_fits[(*ih).first].at(i)->GetParameter(k);
	double err = m_fits[(*ih).first].at(i)->GetParError(k);
	m_parhists[(*ih).first].at(k)->SetBinContent(i+1, par);
	m_parhists[(*ih).first].at(k)->SetBinError(i+1, err);
	k++;
      }
    }
  }

}

void FitAnalysis::SaveToFile(){
  TFile* f = new TFile((m_name+".root").c_str(), "RECREATE");
  info()<<"Saving to "<<m_name<<".root"<<endl;
  for (std::map<string, std::vector<TH1F*> >::iterator iv = m_hists.begin() ; iv != m_hists.end() ; ++iv ){
    std::vector<TH1F*>::iterator ih;
    std::vector<TH1F*> hists = (*iv).second;
    std::vector<TH1F*> parhists = m_parhists[(*iv).first];
    f->mkdir((*iv).first.c_str());
    f->cd((*iv).first.c_str());
    
    for (ih = hists.begin(); ih != hists.end(); ++ih){
      (*ih)->Write();
    }
    for (ih = parhists.begin(); ih != parhists.end(); ++ih){
      (*ih)->Write();
    }
    
    f->cd();

  }
}


#ifdef WITHPYTHON
// for python
void FitAnalysis::FitIt1_py(){ FitIt() };
void FitAnalysis::FitIt2_py(string opt) { FitIt(opt) };
#endif

