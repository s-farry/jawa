#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <Eff.h>
#include <EffVar.h>
#include <TEfficiency.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>
#include <TRandom3.h>

using namespace std;

EffVar::EffVar(string name, string var, std::vector<double> edges, string prefix)
  : Var(name, var, edges, prefix) {

  if (prefix != "") prefix = prefix+"/";

  //Histograms of all plots
  m_tothists = new TObjArray();
  m_passhists = new TObjArray();
  m_failhists = new TObjArray();

  m_tothists->SetOwner(true);
  m_passhists->SetOwner(true);
  m_failhists->SetOwner(true);

  m_totCBFits = new TObjArray();
  m_passCBFits = new TObjArray();

  m_totBkgFits = new TObjArray();
  m_passBkgFits = new TObjArray();

  m_rweff_varyhi_passhists = new TObjArray();
  m_rweff_varylo_passhists = new TObjArray();

  //m_rweff_varyhi_effgraphs = new TObjArray();
  //m_rweff_varylo_effgraphs = new TObjArray();

  m_totCBFits->SetOwner(true);
  m_passCBFits->SetOwner(true);

  m_tothist = new TH1F((prefix+m_name+"_tot").c_str(), m_name.c_str(), m_nbins , &edges[0]);
  m_passhist = new TH1F((prefix+m_name+"_pass").c_str(), (m_name+"_pass").c_str(), m_nbins , &edges[0]);
  m_failhist = new TH1F((prefix+m_name+"_fail").c_str(), (m_name+"_fail").c_str(), m_nbins , &edges[0]);
  m_effgraph = new TGraphAsymmErrors(m_nbins);
  m_effhist  = new TH1F((prefix+m_name+"_eff").c_str(), m_name.c_str(), m_nbins , &edges[0]);

  m_bkgtot  = new TH1F((prefix+m_name+"_bkgtot").c_str() , (m_name+"_bkgtot").c_str()  , m_nbins , &edges[0]);
  m_bkgpass = new TH1F((prefix+m_name+"_bkgpass").c_str(), (m_name+"_bkgpass").c_str() , m_nbins , &edges[0]);

  m_meantot  = new TH1F((prefix+m_name+"_meantot").c_str() , (m_name+"_meantot").c_str()  , m_nbins , &edges[0]);
  m_meanpass = new TH1F((prefix+m_name+"_meanpass").c_str(), (m_name+"_meanpass").c_str() , m_nbins , &edges[0]);

  m_systematic = false;
}

EffVar::EffVar(string name, string var, int nbins, double lo, double hi, string prefix)
  : Var(name, var, nbins, lo, hi, prefix) {

  if (prefix != "") prefix = prefix+"/";

  //Histograms of all plots
  m_tothists = new TObjArray();
  m_passhists = new TObjArray();
  m_failhists = new TObjArray();

  m_tothists->SetOwner(true);
  m_passhists->SetOwner(true);
  m_failhists->SetOwner(true);

  m_totCBFits = new TObjArray();
  m_passCBFits = new TObjArray();

  m_totBkgFits = new TObjArray();
  m_passBkgFits = new TObjArray();


  m_rweff_varyhi_passhists = new TObjArray();
  m_rweff_varylo_passhists = new TObjArray();

  //m_rweff_varyhi_effgraphs = new TObjArray();
  //m_rweff_varylo_effgraphs = new TObjArray();


  m_totCBFits->SetOwner(true);
  m_passCBFits->SetOwner(true);

  m_tothist = new TH1F((prefix+m_name+"_tot").c_str(), m_name.c_str(), m_nbins, m_lo, m_hi);
  m_passhist = new TH1F((prefix+m_name+"_pass").c_str(), (m_name+"_pass").c_str(), m_nbins, m_lo, m_hi);
  m_failhist = new TH1F((prefix+m_name+"_fail").c_str(), (m_name+"_fail").c_str(), m_nbins, m_lo, m_hi);
  m_effgraph = new TGraphAsymmErrors(nbins);
  m_effhist  = new TH1F((prefix+m_name+"_eff").c_str(), m_name.c_str(), m_nbins, m_lo, m_hi);

  m_bkgtot  = new TH1F((prefix+m_name+"_bkgtot").c_str() , (m_name+"_bkgtot").c_str()  , m_nbins, m_lo, m_hi);
  m_bkgpass = new TH1F((prefix+m_name+"_bkgpass").c_str(), (m_name+"_bkgpass").c_str() , m_nbins, m_lo, m_hi);

  m_meantot  = new TH1F((prefix+m_name+"_meantot").c_str() , (m_name+"_meantot").c_str()  , m_nbins, m_lo, m_hi);
  m_meanpass = new TH1F((prefix+m_name+"_meanpass").c_str(), (m_name+"_meanpass").c_str() , m_nbins, m_lo, m_hi);

  m_edges = GetBinEdges(m_tothist);
  m_systematic = false;

}

EffVar::EffVar(string name, PyObject* pyObj): Var(name){
  TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  EffVar(name, effgraph);			   
}

EffVar::EffVar(string name, TGraphAsymmErrors* effgraph)
  :Var(name){
  m_var = name;

  m_tothists  = new TObjArray();
  m_passhists = new TObjArray();
  m_failhists = new TObjArray();

  m_totCBFits  = new TObjArray();
  m_passCBFits = new TObjArray();

  m_totBkgFits  = new TObjArray();
  m_passBkgFits = new TObjArray();


  m_rweff_varyhi_passhists = new TObjArray();
  m_rweff_varylo_passhists = new TObjArray();

  //m_rweff_varyhi_effgraphs = new TObjArray();
  //m_rweff_varylo_effgraphs = new TObjArray();

  m_tothist  = 0;
  m_passhist = 0;
  m_failhist = 0;
  
  m_nbins = (int)effgraph->GetN();

  m_effgraph = effgraph;
  m_systematic = false;


}

EffVar::EffVar(string name, PyObject* pyObj, string prefix): Var(name){
  TFile* f = (TFile*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  EffVar(name, f, prefix);
}
EffVar::EffVar(string name, TFile* f, string prefix): Var(name){
  TObjArray* totalhists = (TObjArray*)f->Get((name+"/TotalHists").c_str());
  TObjArray* passhists  = (TObjArray*)f->Get((name+"/PassHists").c_str());
  TObjArray* failhists  = (TObjArray*)f->Get((name+"/FailHists").c_str());
  
  m_totCBFits = new TObjArray();
  m_passCBFits = new TObjArray();

  m_totBkgFits = new TObjArray();
  m_passBkgFits = new TObjArray();

  m_rweff_varyhi_passhists = new TObjArray();
  m_rweff_varylo_passhists = new TObjArray();

  //m_rweff_varyhi_effgraphs = new TObjArray();
  //m_rweff_varylo_effgraphs = new TObjArray();

  m_systematic = false;

  TH1F* totalhist = (TH1F*)f->Get((name+"/TotalHist").c_str());
  TH1F* passhist  = (TH1F*)f->Get((name+"/PassHist").c_str());
  TH1F* failhist  = (TH1F*)f->Get((name+"/FailHist").c_str());
  
  TGraphAsymmErrors* effgraph = (TGraphAsymmErrors*)f->Get((name+"/EfficiencyGraph").c_str());

  m_name = name;
  m_var  = name;
  if (totalhist){
    m_nbins = totalhist->GetXaxis()->GetNbins();
    m_lo = totalhist->GetXaxis()->GetBinLowEdge(1);
    m_hi = totalhist->GetXaxis()->GetBinUpEdge(m_nbins);
  }
  if (totalhists)   m_tothists  = (TObjArray*)totalhists->Clone();
  if (passhists)    m_passhists = (TObjArray*)passhists->Clone();
  if (failhists)    m_failhists = (TObjArray*)failhists->Clone();
  if (totalhist)    {
    m_tothist   = (TH1F*)totalhist->Clone((prefix+name).c_str());
    m_edges     = GetBinEdges(totalhist);
  }
  if (passhist)     m_passhist  = (TH1F*)passhist->Clone((prefix+name+"_pass").c_str());
  if (failhist)     m_failhist  = (TH1F*)failhist->Clone((prefix+name+"_fail").c_str());
  if (effgraph)     {
    m_effgraph  = (TGraphAsymmErrors*)effgraph->Clone();
    if (m_edges.size() == 0) m_edges = GetBinEdges(effgraph);
  }
}


EffVar::EffVar(EffVar* varA, EffVar* varB, string prefix): Var(varA->GetName(), varA, varB, prefix){
  if ( varA->GetName() == varB->GetName() && 
       varA->GetHi() == varB->GetHi() && varA->GetLo() == varB->GetLo() &&
       varA->GetBins() == varB->GetBins() )
    {
      //m_name               = varA->GetName();
      m_var                = varA->GetVar();
      //m_lo                 = varA->GetLo();
      //m_hi                 = varA->GetHi();
      //m_nbins              = varA->GetBins(); 
      m_tothists           = new TObjArray();
      m_passhists          = new TObjArray();
      m_failhists          = new TObjArray();
      m_edges              = varA->GetEdges();

      m_totCBFits          = new TObjArray();
      m_passCBFits         = new TObjArray();

      m_totBkgFits         = new TObjArray();
      m_passBkgFits        = new TObjArray();
      
      m_rweff_varyhi_passhists = new TObjArray();
      m_rweff_varylo_passhists = new TObjArray();
      
      //m_rweff_varyhi_effgraphs = new TObjArray();
      //m_rweff_varylo_effgraphs = new TObjArray();

      m_prefix             = prefix;

      m_effgraph = new TGraphAsymmErrors(m_nbins);
      
      if ( m_prefix != "" ) m_prefix = m_prefix+"/";
      
      vector<double> edges = GetBinEdges(varA->GetTotHist());

      m_tothist = new TH1F((m_prefix+m_name+"_tot").c_str(), m_name.c_str(), m_nbins, &edges[0]);
      m_passhist = new TH1F((m_prefix+m_name+"_pass").c_str(), (m_name+"_pass").c_str(), m_nbins, &edges[0]);
      m_failhist = new TH1F((m_prefix+m_name+"_fail").c_str(), (m_name+"_fail").c_str(), m_nbins, &edges[0]);

      m_bkgtot  = new TH1F((m_prefix+m_name+"_bkgtot").c_str() , (m_name+"_bkgtot").c_str()  , m_nbins, &edges[0]);
      m_bkgpass = new TH1F((m_prefix+m_name+"_bkgpass").c_str(), (m_name+"_bkgpass").c_str() , m_nbins, &edges[0]);

      m_type = varA->m_type;
      m_systematic = false;

      /*

      for (int i = 0 ; i < varA->m_reweighteffuperrs->GetEntries() ; ++ i){
	TH1F* varA_up = (TH1F*)varA->m_reweighteffuperrs->At(i);
	TH1F* varA_lo = (TH1F*)varA->m_reweighteffloerrs->At(i);
	TH1F* varB_up = (TH1F*)varB->m_reweighteffuperrs->At(i);
	TH1F* varB_lo = (TH1F*)varB->m_reweighteffloerrs->At(i);


  	
      }

      */


      //Loop over tot and pass hists and add together
      for (signed int i = 0; i < m_nbins; ++i){

	double lo = m_lo + i * (m_hi - m_lo)/m_nbins;
	double hi = m_lo + (i + 1 ) * (m_hi - m_lo)/m_nbins;
	
	ostringstream loss, hiss;
	loss << std::fixed << std::setprecision(4) << lo;
	hiss << std::fixed << std::setprecision(4) << hi;
	
	string low  = loss.str();
	string high = hiss.str();
	
	string id_tot    = (m_prefix+"_"+m_name+"_"+m_var+"_"+low+"_"+high+"_tot");
	string title_tot = (m_prefix+"_"+low+"<"+m_var+"<"+high+"_tot");
	
	string id_pass    = (m_prefix+"_"+m_name+"_"+m_var+"_"+low+"_"+high+"_pass");
	string title_pass = (m_prefix+"_"+low+"<"+m_var+"<"+high+"_pass");
	
	string id_fail    = (m_prefix+"_"+m_name+"_"+m_var+"_"+low+"_"+high+"_fail");
	string title_fail = (m_prefix+"_"+low+"<"+m_var+"<"+high+"_fail");
	
	m_tothist->SetBinContent(i+1, varA->m_tothist->GetBinContent(i+1) + varB->m_tothist->GetBinContent(i+1));
	m_passhist->SetBinContent(i+1, varA->m_passhist->GetBinContent(i+1) + varB->m_passhist->GetBinContent(i+1));
	m_failhist->SetBinContent(i+1, varA->m_failhist->GetBinContent(i+1) + varB->m_failhist->GetBinContent(i+1));
	
      }

      if (varA->GetEffRWVaryHiPassHists() && varB->GetEffRWVaryHiPassHists() &&
	  varA->GetEffRWVaryHiPassHists()->GetEntries() == varB->GetEffRWVaryHiPassHists()->GetEntries()){
	for (signed int i = 0 ; i < varA->GetEffRWVaryHiPassHists()->GetEntries() ; ++ i){
	  TH1F* h  = ((TH1F*)varA->GetEffRWVaryHiPassHists()->At(i));
	  TH1F* h2 = ((TH1F*)varB->GetEffRWVaryHiPassHists()->At(i));
	  ostringstream s;
	  s<<h->GetName()<<"_combined";
	  TH1F* h3 = (TH1F*)h->Clone(s.str().c_str());
	  h3->Add(h2);
	  m_rweff_varyhi_passhists->Add(h3);
	}
      }
      if (varA->GetEffRWVaryLoPassHists() && varB->GetEffRWVaryLoPassHists() &&
	  varA->GetEffRWVaryLoPassHists()->GetEntries() == varB->GetEffRWVaryLoPassHists()->GetEntries()){
	for (signed int i = 0 ; i < varA->GetEffRWVaryLoPassHists()->GetEntries() ; ++ i){
	  TH1F* h  = ((TH1F*)varA->GetEffRWVaryLoPassHists()->At(i));
	  TH1F* h2 = ((TH1F*)varB->GetEffRWVaryLoPassHists()->At(i));
	  ostringstream s;
	  s<<h->GetName()<<"_combined";
	  TH1F* h3 = (TH1F*)h->Clone(s.str().c_str());
	  h3->Add(h2);
	  m_rweff_varylo_passhists->Add(h3);
	}
      }
      
    }
  else {
    cout<<"Warning: The variables: "<<varA->GetName()<<" and "<<varB->GetName()<<" do not match. Skipping."<<endl;
    
  }
}

TH1F* EffVar::GetTotHist(){return m_tothist;}
TH1F* EffVar::GetPassHist(){return m_passhist;}
TH1F* EffVar::GetFailHist(){return m_failhist;}
TH1F* EffVar::GetBkgTotHist(){return m_bkgtot;}
TH1F* EffVar::GetBkgPassHist(){return m_bkgpass;}


void EffVar::FillBkgHists(double lo, double hi){
  for (int i = 0 ; i<m_tothists->GetEntries(); ++i){
    TList* funcs_tot  = (TList*)(((TH1F*)m_tothists->At(i))->GetListOfFunctions());
    TList* funcs_pass = (TList*)(((TH1F*)m_passhists->At(i))->GetListOfFunctions());

    double ratio_tot = 1., ratio_pass = 1.;
    if (funcs_tot->GetEntries() == 2 && funcs_pass->GetEntries() == 2){
      TF1* totfunc_tot  = (TF1*)funcs_tot->At(0);
      TF1* bkgfunc_tot  = (TF1*)funcs_tot->At(1);

      TF1* totfunc_pass = (TF1*)funcs_pass->At(0);
      TF1* bkgfunc_pass = (TF1*)funcs_pass->At(1);
      ratio_tot  = bkgfunc_tot->Integral(lo, hi) / totfunc_tot->Integral(lo, hi);
      ratio_pass = bkgfunc_pass->Integral(lo, hi) / totfunc_pass->Integral(lo, hi);

    }
    m_bkgtot->SetBinContent( i+1 , ratio_tot );
    m_bkgpass->SetBinContent( i+1 , ratio_pass );
  }
}

void EffVar::AddEffScaleVaryHists(TH2F* scales){
  for (int i = 0 ; i < scales->GetNbinsX() ; ++i){
    for (int j = 0 ; j < scales->GetNbinsY() ; ++j ){
      ostringstream ssup, sslo;
      ssup<<m_prefix<<m_name<<"_varyup_tot_"<<i<<"_"<<j;
      sslo<<m_prefix<<m_name<<"_varydown_tot_"<<i<<"_"<<j;
      m_rweff_varyhi_passhists->Add(m_passhist->Clone(ssup.str().c_str()));
      m_rweff_varylo_passhists->Add(m_passhist->Clone(sslo.str().c_str()));
    }
  }
}

void EffVar::MakeTGraph(){
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(m_passhist,m_tothist);
  //m_effgraph = (TGraphAsymmErrors*)graph->Clone();

  for (int i = 0; i < m_nbins; ++i){
    double x, y;
    double yerrhi, yerrlo;
    double xerrhi, xerrlo;
    Int_t success = graph->GetPoint(i, x, y);
    
    if (success != -1){
      xerrhi = graph->GetErrorXhigh(i);
      xerrlo = graph->GetErrorXlow(i);
      yerrhi = graph->GetErrorYhigh(i);
      yerrlo = graph->GetErrorYlow(i);
      m_effgraph->SetPoint(i, x, y);
      m_effgraph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
      m_effgraph->GetPoint(i, x, y);
    }
  }
  //cout<<m_reweighteffuperrs<<" "<<m_reweighteffuperrs->GetEntries()<<endl;
  /*
  for (int i = 0 ; i < m_reweighteffuperrs->GetEntries() ; ++ i ){
    ((TH1F*)m_rweighteffuperrs->At(i))->Divide(m_tothist);
  }
  for (int i = 0 ; i < m_reweighteffloerrs->GetEntries() ; ++ i ){
    ((TH1F*)m_reweighteffloerrs->At(i))->Divide(m_tothist);
    }
  */
}

void EffVar::MakeEffHist(bool ClopperPearsonError){
  for (int i = 0; i < m_nbins; ++i){
    Eff eff("efficiency" , m_tothist->GetBinContent( i+1 ), m_passhist->GetBinContent( i+1 ) );

    m_effhist->SetBinContent(i+1, eff.GetEff());
    if (ClopperPearsonError) m_effhist->SetBinError(i+1, max(eff.GetEffErrHi(), eff.GetEffErrLo()));
    else {
      if (m_passhist->GetBinContent(i+1) > 0 && m_tothist->GetBinContent(i+1) > 0){
	double err = m_effhist->GetBinContent(i+1) * (sqrt(pow(m_tothist->GetBinError(i+1)/m_tothist->GetBinContent(i+1),2) + 
							   pow(m_passhist->GetBinError(i+1)/m_passhist->GetBinContent(i+1),2)));
	m_effhist->SetBinError(i+1, err);
      }
    }
  }
}

TObjArray* EffVar::GetEffRWVaryHiPassHists(){ return m_rweff_varyhi_passhists; }
TObjArray* EffVar::GetEffRWVaryLoPassHists(){ return m_rweff_varylo_passhists; }
//TObjArray* EffVar::GetEffRWVarHiGraphs(){ return m_rweff_varyhi_effgraphs; }
//TObjArray* EffVar::GetEffRWVarLoGraphs(){ return m_rweff_varylo_effgraphs; }



void EffVar::FillVar(bool pass, float v_var, double weight, double effw){
  double var = (double)v_var;
  FillVar(pass, var, weight, effw);
}

void EffVar::FillVar(bool pass, double v_var, double weight, double effw ){
  m_tothist->Fill(v_var, weight);
  if ( pass ) {
    m_passhist->Fill(v_var, weight * effw);
  } else {
    m_failhist->Fill(v_var, weight);
  }
}

void EffVar::FillVar(bool pass, double v_var, Utils::weight weight, Utils::weight effw ){
  m_tothist->Fill(v_var, weight.val);
  if ( pass ) {
    m_passhist->Fill(v_var, weight.val * effw.val);
    for (int i = 0 ; i < m_rweff_varyhi_passhists->GetEntries() ; ++ i ){
      if ( i == effw.bin ){
	((TH1F*)m_rweff_varyhi_passhists->At(i))->Fill( v_var, weight.val*(effw.val + effw.err) );
	((TH1F*)m_rweff_varylo_passhists->At(i))->Fill( v_var, weight.val*(effw.val - effw.err) );
      }
      else {
	((TH1F*)m_rweff_varyhi_passhists->At(i))->Fill( v_var, weight.val*effw.val);
	((TH1F*)m_rweff_varylo_passhists->At(i))->Fill( v_var, weight.val*effw.val);
      }
    }

  } else {
    m_failhist->Fill(v_var, weight.val);
  }
}

void EffVar::FillVar(bool pass, int i_var,  double weight, double effw){
  if (strcmp(m_type, "I") != 0) {
    cout<<"-----ERROR - Wrong Type Used"<<endl; 
    return;
  }
  m_tothist->Fill(i_var, weight);
  if ( pass ) {
    m_passhist->Fill(i_var, weight * effw);
  } else m_failhist->Fill(i_var, weight);
}

void EffVar::Normalise(double N){
  m_tothist->Scale(N/m_tothist->Integral());
  m_passhist->Scale(N/m_passhist->Integral());
  m_failhist->Scale(N/m_passhist->Integral());
}


void EffVar::AddSystematic(double pc){
  //pc is percent systematic to be added to all bins
  if (!m_systematic){
    int nentries = m_effgraph->GetN();
    for (int i = 0 ; i < nentries ; ++i){
      double x, y;
      double yerrhi, yerrlo;
      double xerrhi, xerrlo;
      Int_t s = m_effgraph->GetPoint(i, x, y);
      if (s != -1){
	xerrhi = m_effgraph->GetErrorXhigh(i);
	xerrlo = m_effgraph->GetErrorXlow(i);
	yerrhi = m_effgraph->GetErrorYhigh(i);
	yerrlo = m_effgraph->GetErrorYlow(i);
	
	yerrhi = sqrt(pow(yerrhi,2) + pow(y * pc, 2));
	yerrlo = sqrt(pow(yerrlo,2) + pow(y * pc, 2));
	
	m_effgraph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
      }
    }
    m_systematic = true;
  }  else info()<<"Systematic already added to "<<m_name<<": Skipping"<<endl;

}

void EffVar::AddInvSystematic(double pc){
  //pc is percent systematic to be added to all bins
  //This function calculates the systematic as a percentage
  // of 1-efficiency
  if (!m_systematic){
    int nentries = m_effgraph->GetN();
    for (int i = 0 ; i < nentries ; ++i){
      double x, y;
      double yerrhi, yerrlo;
      double xerrhi, xerrlo;
      Int_t s = m_effgraph->GetPoint(i, x, y);
      if (s != -1){
	xerrhi = m_effgraph->GetErrorXhigh(i);
	xerrlo = m_effgraph->GetErrorXlow(i);
	yerrhi = m_effgraph->GetErrorYhigh(i);
	yerrlo = m_effgraph->GetErrorYlow(i);
	
	yerrhi = sqrt(pow(yerrhi,2) + pow((1-y) * pc, 2));
	yerrlo = sqrt(pow(yerrlo,2) + pow((1-y) * pc, 2));
	
	m_effgraph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
      }
    }
    m_systematic = true;
  }  else cout<<"Systematic already added to "<<m_name<<" : Skipping"<<endl;

}



void EffVar::AddSystematic(std::vector<double> pc){
  //pc is percent systematic to be added to individual bins
  if (m_systematic == false){
    unsigned int nentries = m_effgraph->GetN();
    if (nentries != pc.size()){
      cout<<"Warning - Can't add systematic, incompatible sizes"<<endl;
      return;
    }
    
    for (unsigned int i = 0 ; i < nentries ; ++i){
      double x, y;
      double yerrhi, yerrlo;
      double xerrhi, xerrlo;
      Int_t s = m_effgraph->GetPoint(i, x, y);
      if (s != -1){
	xerrhi = m_effgraph->GetErrorXhigh(i);
	xerrlo = m_effgraph->GetErrorXlow(i);
	yerrhi = m_effgraph->GetErrorYhigh(i);
	yerrlo = m_effgraph->GetErrorYlow(i);
	
	yerrhi = sqrt(pow(yerrhi,2) + pow(y * pc.at(i), 2));
	yerrlo = sqrt(pow(yerrlo,2) + pow(y * pc.at(i), 2));
	
	m_effgraph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
      }
    }
  }   else {
    cout<<"Systematic already added to "<<m_name<<" : Skipping"<<endl;
  }
  m_systematic = true;
}

std::vector<double> EffVar::GetBinEdgesX(TH2F* hist){
  int nbins = hist->GetXaxis()->GetNbins();
  std::vector<double> edges;
  //double* edges = new double[nbins+1];
  //double edges[nbins+1];
  for (int i = 0; i <nbins+1 ; ++i){
    edges.push_back(hist->GetXaxis()->GetBinLowEdge(i+1));
  }
  return edges;
}
std::vector<double> EffVar::GetBinEdgesY(TH2F* hist){
  int nbins = hist->GetYaxis()->GetNbins();
  std::vector<double> edges;
  //double* edges = new double[nbins+1];
  //double edges[nbins+1];
  for (int i = 0; i <nbins+1 ; ++i){
    edges.push_back(hist->GetYaxis()->GetBinLowEdge(i+1));
  }
  return edges;
}


void EffVar::SetPrefix(string name){ m_prefix = name;}

void EffVar::RemoveErrors(){
  if (m_effgraph) Utils::RemoveErrors(m_effgraph);
}


TGraphAsymmErrors* EffVar::GetEffGraph(){
  return m_effgraph;
}

PyObject* EffVar::GetEffGraph_py(){
  TGraphAsymmErrors* newCxxObj = new TGraphAsymmErrors(*m_effgraph);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

TGraphAsymmErrors* EffVar::GetSmearedEffGraph(){
  TRandom3 r(0);
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)m_effgraph->Clone();
  for (int i = 0; i <graph->GetN(); ++i){
    double x, y;
    graph->GetPoint(i,x,y);
    double eyhi = graph->GetErrorYhigh(i);
    double eylo = graph->GetErrorYlow(i);
    y = y + r.Gaus(0, max(eyhi, eylo));
    graph->SetPoint(i, x, y);
  }
  return graph;
}

PyObject* EffVar::GetSmearedEffGraph_py(){
  TGraphAsymmErrors* newCxxObj = GetSmearedEffGraph();
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}


void EffVar::AddSystematic1(double pc) {return AddSystematic(pc);}
void EffVar::AddSystematic2(std::vector<double> pc) {return AddSystematic(pc);}

TObjArray* EffVar::GetTotCBFits(){return m_totCBFits;}
TObjArray* EffVar::GetPassCBFits(){return m_passCBFits;}
TObjArray* EffVar::GetTotBkgFits(){return m_totBkgFits;}
TObjArray* EffVar::GetPassBkgFits(){return m_passBkgFits;}
TObjArray* EffVar::GetPassHists(){return m_passhists;}
TObjArray* EffVar::GetTotHists(){return m_tothists;}
TObjArray* EffVar::GetFailHists(){return m_failhists;}

#ifdef WITHPYTHON

boost::python::list EffVar::GetBinEdgesX_py(PyObject* pyObj){
  TH2F* hist = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  vector<double> edges = GetBinEdgesX(hist);
  boost::python::object get_iter = boost::python::iterator<std::vector<double> >();
  boost::python::object iter = get_iter(edges);
  boost::python::list l(iter);
  return l;
}


PyObject* EffVar::GetTotHist_py(){
  TH1F* newCxxObj = new TH1F(*m_tothist);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

PyObject* EffVar::GetPassHist_py(){
  TH1F* newCxxObj = new TH1F(*m_passhist);
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

#endif

/*
template<typename T>
PyObject* rpwa::py::convertToPy(const T& cxxObj) {
	T* newCxxObj = new T(cxxObj);
	return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName(), true);
}

template<typename T>
T rpwa::py::convertFromPy(PyObject* pyObj) {
	TObject* TObj = (TObject*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
	T cxxObj = dynamic_cast<T>(TObj);
	return cxxObj;
}
*/
