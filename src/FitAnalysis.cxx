#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <TFile.h>
#include <Fit.h>
#include <FitAnalysis.h>
#include <boost/algorithm/string.hpp>
#include <Utils.h>


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
void FitAnalysis::SetVal(string var, int bin, string par, double val){
  m_finitvals[var] = std::tuple<int, int, double>();
  int paridx = m_paridx[par];
  std::get<0>(m_finitvals[var]) = bin;
  std::get<1>(m_finitvals[var]) = paridx;
  std::get<2>(m_finitvals[var]) = val;

}
void FitAnalysis::FixVal(string par, double val){
  m_fixvals[par] = val;
}
void FitAnalysis::SetRange(string par, double lo, double hi){
  m_initrange[par] = std::pair<double, double>(lo, hi);
}

void FitAnalysis::Init(Template* t){
  info()<<"Initialising with template "<<t->GetName()<<endl;
  std::map<string, Var*> vars = t->GetVariables();
  std::map<string, Var2D*> varmap = t->Get2DVariables();
  if (vars.find(m_tofit) == vars.end()) {
    info()<<"couldn't find fit variable "<<m_tofit<<endl;
    return;
  }

  m_hist = ((TH1F*)vars[m_tofit]->GetHist()->Clone((m_tofit).c_str()));
  m_fit  = new Fit(m_name, m_func);

  for (std::map<string, Var2D*>::iterator iv = varmap.begin(); iv != varmap.end(); ++iv){
    if (m_tofit != (*iv).second->GetName2()) continue;
    std::vector<Fit*> fits;
    TH2F* h = (TH2F*)(*iv).second->GetHist();
    m_fits[(*iv).second->GetName1()]  = vector<Fit*>(0);
    m_hists[(*iv).second->GetName1()] = vector<TH1F*>(0);
    m_parhists[(*iv).second->GetName1()] = vector<TH1F*>(0);
    
    for (int i = 0 ; i < h->GetNbinsX(); ++i){
      Fit* f = new Fit(m_name+(*iv).second->GetName1()+"_", m_func);
      
      /*
      for (std::map<string, double>::iterator ip = m_initvals.begin(); ip != m_initvals.end(); ++ip) 
	{
	  f->SetParameter(m_paridx[(*ip).first], (*ip).second);
	}
      for (std::map<string, double>::iterator ip = m_fixvals.begin(); ip != m_fixvals.end(); ++ip) 
	{
	  info()<<"fixing "<<(*ip).first<<"( "<<m_paridx[(*ip).first]<<" ) to "<<(*ip).second<<endl;;
	  f->FixParameter(m_paridx[(*ip).first], (*ip).second);
	  }*/
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
  // set all initial parameters
  for (std::map<string, double>::iterator ip = m_initvals.begin(); ip != m_initvals.end(); ++ip) 
    {
      m_fit->SetParameter(m_paridx[(*ip).first], (*ip).second);
    }
  for (std::map<string, double>::iterator ip = m_fixvals.begin(); ip != m_fixvals.end(); ++ip) 
    {
      m_fit->FixParameter(m_paridx[(*ip).first], (*ip).second);
    }
  
  std::vector<int>::iterator iv;
  std::map<string, vector<TH1F*> >::iterator ih;
  for (iv = m_tomean.begin(); iv != m_tomean.end(); ++iv){
    double mean = m_hist->GetMean();
    m_fit->SetParameter((*iv), mean);
  }
  for (iv = m_toentries.begin(); iv != m_toentries.end(); ++iv){
    double entries = m_hist->Integral();
    m_fit->SetParameter((*iv), entries);
  }
  for (iv = m_tomax.begin(); iv != m_tomax.end(); ++iv){
    double max = m_hist->GetBinContent(m_hist->GetMaximum());
    m_fit->SetParameter((*iv), max);
  }
  for (iv = m_torms.begin(); iv != m_torms.end(); ++iv){
    double rms = m_hist->GetRMS();
    m_fit->SetParameter((*iv), rms);
  }

  for (ih = m_hists.begin(); ih != m_hists.end(); ++ih){
    vector<TH1F*> hists = (*ih).second;
    for (unsigned int i = 0 ; i < hists.size() ; ++i){
      for (std::map<string, double>::iterator ip = m_initvals.begin(); ip != m_initvals.end(); ++ip) 
	{
	  m_fits[(*ih).first].at(i)->SetParameter(m_paridx[(*ip).first], (*ip).second);
	}
      for (std::map<string, double>::iterator ip = m_fixvals.begin(); ip != m_fixvals.end(); ++ip) 
	{
	  m_fits[(*ih).first].at(i)->FixParameter(m_paridx[(*ip).first], (*ip).second);
	}
      
      for (iv = m_tomean.begin(); iv != m_tomean.end(); ++iv){
	double mean = hists.at(i)->GetMean();
	m_fits[(*ih).first].at(i)->SetParameter((*iv), mean);
      }
      for (iv = m_toentries.begin(); iv != m_toentries.end(); ++iv){
	double entries = hists.at(i)->Integral();
	m_fits[(*ih).first].at(i)->SetParameter((*iv), entries);
      }
      for (iv = m_tomax.begin(); iv != m_tomax.end(); ++iv){
	double max = hists.at(i)->GetBinContent(hists.at(i)->GetMaximum());
	m_fits[(*ih).first].at(i)->SetParameter((*iv), max);
      }
      for (iv = m_torms.begin(); iv != m_torms.end(); ++iv){
	double rms = hists.at(i)->GetRMS();
	m_fits[(*ih).first].at(i)->SetParameter((*iv), rms);
      }
    }
  }

  for (std::map<string, std::tuple<int, int, double> >::iterator iv =  m_finitvals.begin();
       iv != m_finitvals.end(); ++iv){
    if (m_fits.find((*iv).first) != m_fits.end()){
      std::tuple<int, int, double> tup = (*iv).second;
      int bin = std::get<0>(tup);
      int par = std::get<1>(tup);
      double val = std::get<2>(tup);
      m_fits[(*iv).first].at(bin)->SetParameter(par, val);
    }
    else info()<<"couldn't add constraint for "<<(*iv).first<<" "<<endl;

  }
  m_fit->FitHist(m_hist, opt);
  int k = 0;
  for (unsigned int j = 0 ; j < m_expr->GetVarNames().size(); ++j){
    if (m_expr->GetVarNames().at(j) == "x") continue;
    m_pars[m_expr->GetVarNames().at(j)] = m_fit->GetParameter(k);
    m_parerrs[m_expr->GetVarNames().at(j)] = m_fit->GetParError(k);
    k++;
  }
  
  for (ih = m_hists.begin(); ih != m_hists.end(); ++ih){
    vector<TH1F*> hists = (*ih).second;
    for (unsigned int i = 0 ; i < hists.size() ; ++i){
      m_fits[(*ih).first].at(i)->FitHist(hists.at(i), opt);
      k = 0;
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
    
  std::map<string, double>::iterator ip;
  m_hist->Write();
  for (ip = m_pars.begin() ; ip != m_pars.end() ; ++ip){
    TParameter<double>((*ip).first.c_str(), (*ip).second).Write();
  }
  for (ip = m_parerrs.begin() ; ip != m_pars.end() ; ++ip){
    TParameter<double>(((*ip).first+"_err").c_str(), (*ip).second).Write();
  }
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
void FitAnalysis::SaveToFile(string name){
  TFile* f = new TFile(name.c_str(), "RECREATE");
  info()<<"Saving to "<<m_name<<".root"<<endl;
  std::map<string, double>::iterator ip;

  m_hist->Write();
  for (ip = m_pars.begin() ; ip != m_pars.end() ; ++ip){
    TParameter<double>((*ip).first.c_str(), (*ip).second).Write();
  }
  for (ip = m_parerrs.begin() ; ip != m_parerrs.end() ; ++ip){
    TParameter<double>(((*ip).first+"_err").c_str(), (*ip).second).Write();
  }
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


void FitAnalysis::SetToRMS(string name){
  int idx = m_paridx[name];
  m_torms.push_back(idx);
  /*
  for (std::map<string, std::vector<TH1F*> >::iterator iv = m_hists.begin() ; iv != m_hists.end() ; ++iv ){
    std::vector<TH1F*> hists = (*iv).second;
    for (int i = 0 ; i < hists.size() ;++i){
      double mean = hists.at(i)->GetMean();
      m_fits[(*iv).first].at(i)->SetParameter(idx, mean);
    }
    }*/
}
void FitAnalysis::SetToMean(string name){
  int idx = m_paridx[name];
  m_tomean.push_back(idx);
  /*
  for (std::map<string, std::vector<TH1F*> >::iterator iv = m_hists.begin() ; iv != m_hists.end() ; ++iv ){
    std::vector<TH1F*> hists = (*iv).second;
    for (int i = 0 ; i < hists.size() ;++i){
      double mean = hists.at(i)->GetMean();
      m_fits[(*iv).first].at(i)->SetParameter(idx, mean);
    }
    }*/
}
  
void FitAnalysis::SetToEntries(string name){
  int idx = m_paridx[name];
  m_toentries.push_back(idx);
  /*
  for (std::map<string, std::vector<TH1F*> >::iterator iv = m_hists.begin() ; iv != m_hists.end() ; ++iv ){
    std::vector<TH1F*> hists = (*iv).second;
    for (int i = 0 ; i < hists.size() ;++i){
      double entries = hists.at(i)->Integral();
      m_fits[(*iv).first].at(i)->SetParameter(idx, entries);
    }
  }
  */
}

void FitAnalysis::SetToMax(string name){
  int idx = m_paridx[name];
  m_tomax.push_back(idx);
  /*
  for (std::map<string, std::vector<TH1F*> >::iterator iv = m_hists.begin() ; iv != m_hists.end() ; ++iv ){
    std::vector<TH1F*> hists = (*iv).second;
    for (int i = 0 ; i < hists.size() ;++i){
      double max = hists.at(i)->GetBinContent(hists.at(i)->GetMaximum());
      m_fits[(*iv).first].at(i)->SetParameter(idx, max);
    }
  }
  */
}

std::vector<TH1F*> FitAnalysis::GetHists(string name){ return m_hists[name]; }
std::vector<TH1F*> FitAnalysis::GetParHists(string name) { return m_parhists[name];}
TH1F* FitAnalysis::GetHist(){ return m_hist; }

#ifdef WITHPYTHON
// for python
void FitAnalysis::FitIt1_py(){ FitIt(); };
void FitAnalysis::FitIt2_py(string opt) { FitIt(opt); };
PyObject* FitAnalysis::GetHist_py(string name, int i){ 
  TH1F* hist = m_hists[name].at(i);
  PyObject* pyObj =  TPython::ObjectProxy_FromVoidPtr(hist, hist->ClassName());;
  return pyObj;
}
PyObject* FitAnalysis::GetParHist_py(string name, int i) {
  TH1F* hist = m_parhists[name].at(i);
  PyObject* pyObj =  TPython::ObjectProxy_FromVoidPtr(hist, hist->ClassName());;
  return pyObj;
}
PyObject* FitAnalysis::GetHist2_py(){
  PyObject* pyObj =  TPython::ObjectProxy_FromVoidPtr(m_hist, m_hist->ClassName());;
  return pyObj;
}
void FitAnalysis::SaveToFile1_py() { SaveToFile();}
void FitAnalysis::SaveToFile2_py(string name) { SaveToFile(name); }
void FitAnalysis::SetVal1_py(string par, double val) { SetVal(par, val); }
void FitAnalysis::SetVal2_py(string var, int bin, string par, double val) { SetVal(var, bin, par, val); }
#endif
