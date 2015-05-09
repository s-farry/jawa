#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <EffVar2D.h>
#include <TEfficiency.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>

using namespace std;

EffVar2D::EffVar2D(string name, EffVar* varA, EffVar* varB, string prefix)
  : Var2D(name, varA, varB, prefix){

  //m_name1 = varA->GetName();
  //m_name2 = varB->GetName();

  //m_varname1   = varA->GetVar();
  //m_varname2   = varB->GetVar();
  
  //m_nbins1 = varA->GetBins();
  //m_lo1    = varA->GetLo();
  //m_hi1    = varA->GetHi();

  //m_nbins2 = varB->GetBins();
  //m_lo2    = varB->GetLo();
  //m_hi2    = varB->GetHi();

  vector<double> edgesA = Var::GetBinEdges(varA->GetTotHist());
  vector<double> edgesB = Var::GetBinEdges(varB->GetTotHist());

  //Histograms of all plots
  m_tothists = new TObjArray();
  m_passhists = new TObjArray();
  m_failhists = new TObjArray();

  m_tothists->SetOwner(true);
  m_passhists->SetOwner(true);
  m_failhists->SetOwner(true);

  m_totCBFits = new TObjArray();
  m_passCBFits = new TObjArray();
  
  m_totCBFits->SetOwner(true);
  m_passCBFits->SetOwner(true);

  m_effgraphs = new TObjArray();

  m_meantot  = new TH2F((prefix+m_name+"_meantot").c_str() , (m_name+"_meantot").c_str()  , edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);
  m_meanpass = new TH2F((prefix+m_name+"_meanpass").c_str(), (m_name+"_meanpass").c_str() , edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);


  m_tothist = new TH2F((prefix+m_name).c_str(), m_name.c_str(), edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);
  m_passhist = new TH2F((prefix+m_name+"_pass").c_str(), (m_name+"_pass").c_str(), edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);
  m_failhist = new TH2F((prefix+m_name+"_fail").c_str(), (m_name+"_fail").c_str(), edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);
  m_2Deffgraph = new TH2F((prefix+m_name+"_eff").c_str(), (m_name+"_eff").c_str(), edgesA.size()-1, &edgesA[0], edgesB.size()-1, &edgesB[0]);
  //m_effgraph = new TGraphAsymmErrors(nbins);

}


EffVar2D::EffVar2D(EffVar2D* varA, EffVar2D* varB, string prefix)
  :Var2D(varA->GetName()){
  if ( varA->GetName() == varB->GetName() && 
       varA->GetVar1Hi() == varB->GetVar1Hi() && varA->GetVar1Lo() == varB->GetVar1Lo() &&
       varA->GetVar2Hi() == varB->GetVar2Hi() && varA->GetVar2Lo() == varB->GetVar2Lo() &&
       varA->GetBins1() == varB->GetBins1() &&
       varA->GetBins2() == varB->GetBins2() )
    {
      m_name1               = varA->GetName1();
      m_varname1            = varA->GetVarName1();
      m_lo1                 = varA->GetVar1Lo();
      m_hi1                 = varA->GetVar1Hi();
      m_nbins1              = varA->GetBins1(); 
      
      m_name2               = varA->GetName2();
      m_varname2            = varA->GetVarName2();
      m_lo2                 = varA->GetVar2Lo();
      m_hi2                 = varA->GetVar2Hi();
      m_nbins2              = varA->GetBins2(); 

      //cout<<"Combining "<<varA->GetName()<<endl;

      m_tothists           = new TObjArray();
      m_passhists          = new TObjArray();
      m_failhists          = new TObjArray();

      m_totCBFits          = new TObjArray();
      m_passCBFits         = new TObjArray();

      m_effgraphs          = new TObjArray();

      vector<double> edgesX = EffVar::GetBinEdgesX(varA->GetTotHist());
      vector<double> edgesY = EffVar::GetBinEdgesY(varA->GetTotHist());
      
      m_prefix             = prefix;
      
      if ( prefix != "" ) prefix = prefix+"/";
      
      m_tothist = new TH2F((prefix+m_name).c_str(), m_name.c_str(), edgesX.size()-1, &edgesX[0], edgesY.size()-1, &edgesY[0]);
      m_passhist = new TH2F((prefix+m_name+"_pass").c_str(), (m_name+"_pass").c_str(), edgesX.size()-1, &edgesX[0], edgesY.size()-1, &edgesY[0]);
      m_failhist = new TH2F((prefix+m_name+"_fail").c_str(), (m_name+"_fail").c_str(), edgesX.size()-1, &edgesX[0], edgesY.size()-1, &edgesY[0]);
      m_2Deffgraph = new TH2F((prefix+m_name+"_eff").c_str(), (m_name+"_eff").c_str(), edgesX.size()-1, &edgesX[0], edgesY.size()-1, &edgesY[0]);

      //Loop over tot and pass hists and add together

      TH1F* hist = (TH1F*)varA->m_tothists->At(0);
      int npltbins = hist->GetNbinsX();
      double pltrangelow = hist->GetXaxis()->GetXmin();
      double pltrangehi  = hist->GetXaxis()->GetXmax();

      for ( int i = 0 ; i < m_nbins1 ; ++i ){

	//Get lo and hi for first variable
	
	//double lo1 = m_lo1 + i * (m_hi1 - m_lo1)/m_nbins1;
	//double hi1 = m_lo1 + (i + 1 ) * (m_hi1 - m_lo1)/m_nbins1;
	double lo1 = varA->GetTotHist()->GetXaxis()->GetBinLowEdge(i+1);
	double hi1 = varA->GetTotHist()->GetXaxis()->GetBinLowEdge(i+2);
	
	ostringstream loss1, hiss1;
	loss1 << std::fixed << std::setprecision(2) << lo1;
	hiss1 << std::fixed << std::setprecision(2) << hi1;
	
	string low1  = loss1.str();
	string high1 = hiss1.str();
	
	string graphid = m_prefix+"_"+m_varname1+"_"+low1+"_"+high1;

	m_effgraphs->Add(new TGraphAsymmErrors(m_nbins2));

	string graph_lab    = (m_prefix+"_"+m_varname1+"_"+low1+"_"+high1);

	((TGraphAsymmErrors*)m_effgraphs->At(i))->SetTitle(graph_lab.c_str());
	((TGraphAsymmErrors*)m_effgraphs->At(i))->SetName(graph_lab.c_str());

	for ( int j = 0 ; j < m_nbins2 ; ++j ) {
	  
	  //Get lo and hi for second variable
	  double lo2 = varA->GetTotHist()->GetYaxis()->GetBinLowEdge(j+1);
	  double hi2 = varA->GetTotHist()->GetYaxis()->GetBinLowEdge(j+2);
	  //double lo2 = m_lo2 + j * (m_hi2 - m_lo2)/m_nbins2;
	  //double hi2 = m_lo2 + (j + 1 ) * (m_hi2 - m_lo2)/m_nbins2;
	  
	  ostringstream loss2, hiss2;
	  loss2 << std::fixed << std::setprecision(2) << lo2;
	  hiss2 << std::fixed << std::setprecision(2) << hi2;
	  
	  string low2  = loss2.str();
	  string high2 = hiss2.str();
	  
	  string id_tot    = (m_prefix+"_"+m_varname1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_tot");
	  string title_tot    = (m_prefix+"_"+low1+"_"+m_varname1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_tot");
	  
	  string id_pass   = (m_prefix+"_"+m_varname1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_pass");
	  string title_pass  = (m_prefix+"_"+low1+"_"+m_varname1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_pass");
	  
	  string id_fail   = (m_prefix+"_"+m_varname1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_fail");
	  string title_fail  = (m_prefix+"_"+low1+"_"+m_varname1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_fail");
	  
	  TH1F* tothist  = new TH1F(id_tot.c_str(),   title_tot.c_str(),  npltbins , pltrangelow , pltrangehi );
	  TH1F* passhist = new TH1F(id_pass.c_str(),   title_pass.c_str(),  npltbins , pltrangelow , pltrangehi );
	  TH1F* failhist = new TH1F(id_fail.c_str(),   title_fail.c_str(),  npltbins , pltrangelow , pltrangehi );

	  int bin = i * m_nbins2 + j;
	  tothist->Add((TH1F*)varA->m_tothists->At(bin));
	  tothist->Add((TH1F*)varB->m_tothists->At(bin));

	  passhist->Add((TH1F*)varA->m_passhists->At(bin));
	  passhist->Add((TH1F*)varB->m_passhists->At(bin));

	  failhist->Add((TH1F*)varA->m_failhists->At(bin));
	  failhist->Add((TH1F*)varB->m_failhists->At(bin));

	  m_tothists->Add( tothist );
	  m_passhists->Add( passhist );
	  m_failhists->Add( failhist );

	  m_tothist->SetBinContent( i+1 , j+1  , varA->m_tothist->GetBinContent( i+1 , j+1 ) + varB->m_tothist->GetBinContent( i+1 , j+1 ));
	  m_passhist->SetBinContent( i+1 , j+1 , varA->m_passhist->GetBinContent( i+1 , j+1 ) + varB->m_passhist->GetBinContent( i+1 , j+1 ));
	  m_failhist->SetBinContent( i+1 , j+1 , varA->m_failhist->GetBinContent( i+1 , j+1 ) + varB->m_failhist->GetBinContent( i+1 , j+1 ));
	  
	  //m_tothist->Add(varA->m_tothist);
	  //m_tothist->Add(varB->m_tothist);

	  //m_passhist->Add(varA->m_passhist);
	  //m_passhist->Add(varB->m_passhist);

	}
      }
      //cout<<"Finished combining"<<endl;
    }
}

void EffVar2D::Init(string name, TFile* f, string prefix){
 TObjArray* totalhists = (TObjArray*)f->Get((name+"/TotalHists").c_str());
  TObjArray* passhists  = (TObjArray*)f->Get((name+"/PassHists").c_str());
  
  TH2F* totalhist = (TH2F*)f->Get((name+"/TotalHist").c_str());
  TH2F* passhist  = (TH2F*)f->Get((name+"/PassHist").c_str());
  TH2F* failhist  = (TH2F*)f->Get((name+"/FailHist").c_str());
  TH2F* effgraph  = (TH2F*)f->Get((name+"/EffGraph2D").c_str());
  
  TObjArray* effgraphs = (TObjArray*)f->Get((name+"/EffGraphs").c_str());

  std::vector<string> delim = split(name, '_');
  
  m_name1               = delim.at(0);
  m_varname1            = delim.at(0);

  if (totalhist){
    m_nbins1 = totalhist->GetXaxis()->GetNbins();
    m_lo1 = totalhist->GetXaxis()->GetBinLowEdge(1);
    m_hi1 = totalhist->GetXaxis()->GetBinUpEdge(m_nbins1);
  }
  
  m_name2               = delim.at(1);
  m_varname2            = delim.at(1);

  if (totalhist) {
    m_nbins2 = totalhist->GetYaxis()->GetNbins();
    m_lo2 = totalhist->GetYaxis()->GetBinLowEdge(1);
    m_hi2 = totalhist->GetYaxis()->GetBinUpEdge(m_nbins2);
  }


  if (totalhists) m_tothists  = (TObjArray*)totalhists->Clone();
  if (passhists)  m_passhists = (TObjArray*)passhists->Clone();
  if (totalhist)  m_tothist   = (TH2F*)totalhist->Clone((prefix+name).c_str());
  if (passhist)   m_passhist  = (TH2F*)passhist->Clone((prefix+name+"_pass").c_str());
  if (failhist)   m_failhist  = (TH2F*)passhist->Clone((prefix+name+"_fail").c_str());
  if (effgraph)   m_2Deffgraph  = (TH2F*)effgraph->Clone((prefix+name+"_eff").c_str());
  if (effgraphs)  m_effgraphs = (TObjArray*)effgraphs->Clone();

}

EffVar2D::EffVar2D(string name, TFile* f, string prefix)
  :Var2D(name){
  Init(name, f, prefix);
  /*
  TObjArray* totalhists = (TObjArray*)f->Get((name+"/TotalHists").c_str());
  TObjArray* passhists  = (TObjArray*)f->Get((name+"/PassHists").c_str());
  
  TH2F* totalhist = (TH2F*)f->Get((name+"/TotalHist").c_str());
  TH2F* passhist  = (TH2F*)f->Get((name+"/PassHist").c_str());
  TH2F* failhist  = (TH2F*)f->Get((name+"/FailHist").c_str());
  TH2F* effgraph  = (TH2F*)f->Get((name+"/EffGraph2D").c_str());
  
  TObjArray* effgraphs = (TObjArray*)f->Get((name+"/EffGraphs").c_str());

  std::vector<string> delim = split(name, '_');
  
  m_name1               = delim.at(0);
  m_varname1            = delim.at(0);

  if (totalhist){
    m_nbins1 = totalhist->GetXaxis()->GetNbins();
    m_lo1 = totalhist->GetXaxis()->GetBinLowEdge(1);
    m_hi1 = totalhist->GetXaxis()->GetBinUpEdge(m_nbins1);
  }
  
  m_name2               = delim.at(1);
  m_varname2            = delim.at(1);

  if (totalhist) {
    m_nbins2 = totalhist->GetYaxis()->GetNbins();
    m_lo2 = totalhist->GetYaxis()->GetBinLowEdge(1);
    m_hi2 = totalhist->GetYaxis()->GetBinUpEdge(m_nbins2);
  }


  if (totalhists) m_tothists  = (TObjArray*)totalhists->Clone();
  if (passhists)  m_passhists = (TObjArray*)passhists->Clone();
  if (totalhist)  m_tothist   = (TH2F*)totalhist->Clone((prefix+name).c_str());
  if (passhist)   m_passhist  = (TH2F*)passhist->Clone((prefix+name+"_pass").c_str());
  if (failhist)   m_failhist  = (TH2F*)passhist->Clone((prefix+name+"_fail").c_str());
  if (effgraph)   m_2Deffgraph  = (TH2F*)effgraph->Clone((prefix+name+"_eff").c_str());
  if (effgraphs)  m_effgraphs = (TObjArray*)effgraphs->Clone();
  cout<<"Total hist is at "<<m_tothist<<endl;
  cout<<"EffGraph is at "<<m_2Deffgraph<<endl;
  */
}


EffVar2D::EffVar2D(string name, PyObject* f, string prefix): Var2D(name)
{
  TFile* file = (TFile*)(TPython::ObjectProxy_AsVoidPtr(f));
  Init(name, file, prefix);
}

EffVar2D::EffVar2D(string name, string f, string prefix): Var2D(name)
{
  TFile* file = new TFile(f.c_str());
  EffVar2D(name, file, prefix);
}

EffVar2D::EffVar2D(string name, TH2F* effgraph) : Var2D(name){
  std::vector<string> delim = split(name, '_');
  m_name1               = delim.at(0);
  m_varname1            = delim.at(0);
  m_name2               = delim.at(1);
  m_varname2            = delim.at(1);

  m_nbins1 = effgraph->GetXaxis()->GetNbins();
  m_lo1    = effgraph->GetXaxis()->GetBinLowEdge(1);
  m_hi1    = effgraph->GetXaxis()->GetBinUpEdge(m_nbins1);
  m_nbins2 = effgraph->GetYaxis()->GetNbins();
  m_lo2    = effgraph->GetYaxis()->GetBinLowEdge(1);
  m_hi2    = effgraph->GetYaxis()->GetBinUpEdge(m_nbins2);

  m_tothists = new TObjArray();
  m_passhists = new TObjArray();
  m_failhists = new TObjArray();

  m_totCBFits = new TObjArray();
  m_passCBFits = new TObjArray();

  m_tothist = 0;
  m_passhist = 0;
  m_failhist = 0;
  
  m_2Deffgraph = effgraph;
  m_effgraphs = new TObjArray();


}

TH2F* EffVar2D::Get2DEffGraph(){ return m_2Deffgraph;}

TH2F* EffVar2D::GetTotHist(){
  return m_tothist;
}
PyObject* EffVar2D::Get2DEffGraph_py(){
  TH2F* newCxxObj = (TH2F*)m_2Deffgraph->Clone();
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}
PyObject* EffVar2D::GetTotHist_py(){
  TH2F* newCxxObj = (TH2F*)m_tothist->Clone();
  return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
}

void EffVar2D::MakeEffHist(bool ClopperPearsonError){
  for ( int i = 0 ; i < m_nbins1+2 ; ++i ){
    for ( int j = 0 ; j < m_nbins2+2 ; ++j ){
      //if (m_name == "ETA_PT") cout<<i<<" "<<j<<endl;
      Eff eff("efficiency", m_tothist->GetBinContent(i,j), m_passhist->GetBinContent( i , j));
      m_2Deffgraph->SetBinContent(i, j, eff.GetEff());

      if (ClopperPearsonError) m_2Deffgraph->SetBinError(i, j, max(eff.GetEffErrHi(), eff.GetEffErrLo()));
      else {
	double err = m_2Deffgraph->GetBinContent(i,j) * 
	  (sqrt(pow(m_tothist->GetBinError(i,j)/m_tothist->GetBinContent(i,j),2) + 
		pow(m_passhist->GetBinError(i,j)/m_passhist->GetBinContent(i,j),2)));
	m_2Deffgraph->SetBinError(i, j, err);
      }
    }
  }
}

void EffVar2D::FillMeanHists(){
  for ( int i = 0 ; i < m_nbins1 ; ++i ){
	for ( int j = 0 ; j < m_nbins2 ; ++j ) {
	  int bin = i * m_nbins2 + j;

	  TList* funcs_tot  = (TList*)(((TH1F*)m_tothists->At(bin))->GetListOfFunctions());
	  TList* funcs_pass = (TList*)(((TH1F*)m_passhists->At(bin))->GetListOfFunctions());

	  double mean_tot     = 0;
	  double mean_pass    = 0;
	  double mean_toterr  = 0;
	  double mean_passerr = 0;
	  
	  if (funcs_tot->GetEntries() == 2 && funcs_pass->GetEntries() == 2){
	    TF1* totfunc_tot  = (TF1*)funcs_tot->At(0);
	    //TF1* bkgfunc_tot  = (TF1*)funcs_tot->At(1);
	    
	    TF1* totfunc_pass = (TF1*)funcs_pass->At(0);
	    //TF1* bkgfunc_pass = (TF1*)funcs_pass->At(1);
	    
	    mean_tot     = totfunc_tot->GetParameter(3);
	    mean_pass    = totfunc_pass->GetParameter(3);
	    mean_toterr  = totfunc_tot->GetParError(3);
	    mean_passerr = totfunc_pass->GetParError(3);

	  }
	  m_meantot->SetBinContent( i+1 , j+1, mean_tot );
	  m_meanpass->SetBinContent( i+1 , j+1, mean_pass );
	  m_meantot->SetBinError( i+1 , j+1, mean_toterr );
	  m_meanpass->SetBinError( i+1 , j+1, mean_passerr );
	}
  }
}


void EffVar2D::Normalise(double N){
  m_tothist->Scale(N/m_tothist->Integral());
  m_passhist->Scale(N/m_passhist->Integral());
}

void EffVar2D::MakeHists(string name, int npltbins, double pltrangelow, double pltrangehi, bool reweight){
  if (reweight){
    m_tothist->Sumw2();
    m_passhist->Sumw2();
  }
  
  for ( int i = 0 ; i < m_nbins1 ; ++i ){
    
    //Get lo and hi for first variable
    //double lo1 = m_lo1 + i * (m_hi1 - m_lo1)/m_nbins1;
    //double hi1 = m_lo1 + (i + 1 ) * (m_hi1 - m_lo1)/m_nbins1;
    
    double lo1 = m_tothist->GetXaxis()->GetBinLowEdge(i+1);
    double hi1 = m_tothist->GetXaxis()->GetBinLowEdge(i+2);

    ostringstream loss1, hiss1;
    loss1 << std::fixed << std::setprecision(2) << lo1;
    hiss1 << std::fixed << std::setprecision(2) << hi1;
    
    string low1  = loss1.str();
    string high1 = hiss1.str();
    
    string graph_lab    = (name+"_"+m_name1+"_"+low1+"_"+high1);
    
    m_effgraphs->Add(new TGraphAsymmErrors(m_nbins2));
    
    ((TGraphAsymmErrors*)m_effgraphs->At(i))->SetTitle(graph_lab.c_str());
    ((TGraphAsymmErrors*)m_effgraphs->At(i))->SetName(graph_lab.c_str());
    
    for ( int j = 0 ; j < m_nbins2 ; ++j ) {
      
      //Get lo and hi for second variable
      double lo2 = m_tothist->GetYaxis()->GetBinLowEdge(j+1);
      double hi2 = m_tothist->GetYaxis()->GetBinLowEdge(j+2);
      //double lo2 = m_lo2 + j * (m_hi2 - m_lo2)/m_nbins2;
      //double hi2 = m_lo2 + (j + 1 ) * (m_hi2 - m_lo2)/m_nbins2;
      
      ostringstream loss2, hiss2;
      loss2 << std::fixed << std::setprecision(2) << lo2;
      hiss2 << std::fixed << std::setprecision(2) << hi2;
      
      string low2  = loss2.str();
      string high2 = hiss2.str();
      
      string id_tot    = (name+"_"+m_name1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_tot");
      string title_tot    = (name+"_"+low1+"_"+m_name1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_tot");
      
      string id_pass   = (name+"_"+m_name1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_pass");
      string title_pass  = (name+"_"+low1+"_"+m_name1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_pass");
      
      string id_fail   = (name+"_"+m_name1+"_"+low1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_fail");
      string title_fail  = (name+"_"+low1+"_"+m_name1+"_"+high1+"_"+m_name2+"_"+low2+"_"+high2+"_fail");
      
      m_tothists->Add(new TH1F(id_tot.c_str(), title_tot.c_str(), npltbins,pltrangelow,pltrangehi));
      m_passhists->Add(new TH1F(id_pass.c_str(), title_pass.c_str(), npltbins,pltrangelow,pltrangehi));
      m_failhists->Add(new TH1F(id_fail.c_str(), title_fail.c_str(), npltbins,pltrangelow,pltrangehi));
      
      
      if (reweight){
	((TH1F*)m_tothists->At(i*m_nbins2 + j))->Sumw2();
	((TH1F*)m_passhists->At(i*m_nbins2 + j))->Sumw2();
      }
      
    }
  }
}

void EffVar2D::MakeTGraphs(){
  for ( int i = 0 ; i < m_nbins1 ; ++i ){
    for ( int j = 0 ; j < m_nbins2 ; ++j ){
      double total = m_tothist->GetBinContent( i+1 , j+1 );
      double pass  = m_passhist->GetBinContent( i+1 , j+1 );

      //int bin = m_tothist->GetBin( i+1, j+1 );

      double x = m_tothist->GetYaxis()->GetBinCenter( j+1 );
      double xerrhi = m_tothist->GetYaxis()->GetBinWidth( j+1 )/2.;
      double xerrlo = xerrhi;

      Eff eff = Eff("eff", total, pass);
      double yerrhi = eff.GetEffErrHi();
      double yerrlo = eff.GetEffErrLo();
      double y = eff.GetEff();

      ((TGraphAsymmErrors*)m_effgraphs->At(i))->SetPoint(j , x , y);
      ((TGraphAsymmErrors*)m_effgraphs->At(i))->SetPointError( j , xerrlo, xerrhi, yerrlo, yerrhi );

    }
  }
}

void EffVar2D::FillVar(bool pass, double& v_pltvar, int& v_var1, double& v_var2, double efflo, double effhi, double w){
  double d = (double)v_var1;

  FillVar(pass, v_pltvar, d, v_var2, efflo, effhi, w);
}
void EffVar2D::FillVar(bool pass, double& v_pltvar, double& v_var1, int& v_var2, double efflo, double effhi, double w){
  double d = (double)v_var2;
  FillVar(pass, v_pltvar, v_var1, d, efflo, effhi, w);
}
void EffVar2D::FillVar(bool pass, double& v_pltvar, int& v_var1, int& v_var2, double efflo, double effhi, double w){
  double d1 = (double)v_var1;
  double d2 = (double)v_var2;
  FillVar(pass, v_pltvar, d1, d2, efflo, effhi, w);
}
void EffVar2D::FillVar(bool pass, double& v_pltvar, double& v_var1, double& v_var2, double efflo, double effhi, double w){
  bool inRange = (v_var1 >= m_lo1 && v_var1 < m_hi1) && (v_var2 >= m_lo2 && v_var2 < m_hi2);
  if ( inRange ){
	  
    //------------- Calculate bin and fill the individual mass histograms
    int bin1 = floor((v_var1 - m_lo1 ) * m_nbins1/( m_hi1 - m_lo1));
    int bin2 = floor((v_var2 - m_lo2 ) * m_nbins2/( m_hi2 - m_lo2));
    
    int bin = bin1 * m_nbins2 + bin2;
    
    ((TH1F*)m_tothists->At(bin))->Fill(v_pltvar, w);
    if (pass) {
      ((TH1F*)m_passhists->At(bin))->Fill(v_pltvar, w);
    }
    else {
      ((TH1F*)m_failhists->At(bin))->Fill(v_pltvar,w);
    }
  }
  //------------ Fill the total histograms within the efficiency range specified
  if (v_pltvar >= efflo && v_pltvar < effhi) 
    {
      m_tothist->Fill(v_var1, v_var2, w);
      if ( pass ) {
	m_passhist->Fill(v_var1, v_var2, w);
      } else m_failhist->Fill(v_var1, v_var2, w);
    }
}

void EffVar2D::FillVar(string type, double& v_pltvar, double& v_var1, double& v_var2, double efflo, double effhi, double w){
  if ( type != "tot" && type != "pass" ) cout<<"Please Enter either tot or pass for type - SKIPPING!"<<endl;

  bool inRange = (v_var1 >= m_lo1 && v_var1 < m_hi1) && (v_var2 >= m_lo2 && v_var2 < m_hi2);
  if ( inRange ){
    //------------ Fill the total histograms within the efficiency range specified
    if (v_pltvar >= efflo && v_pltvar < effhi) 
      {
	
	if ( type == "tot" ) m_tothist->Fill(v_var1, v_var2, w);
	else if ( type == "pass" ) m_passhist->Fill(v_var1, v_var2, w);
	
	//------------- Calculate bin and fill the individual mass histograms
	int bin1 = floor((v_var1 - m_lo1 ) * m_nbins1/( m_hi1 - m_lo1));
	int bin2 = floor((v_var2 - m_lo2 ) * m_nbins2/( m_hi2 - m_lo2));

	int bin = bin1 * m_nbins2 + bin2;
	
	if ( type == "tot" ) ((TH1F*)m_tothists->At(bin))->Fill(v_pltvar, w);
	else if (type == "pass" ) ((TH1F*)m_passhists->At(bin))->Fill(v_pltvar, w);
	
      }
  }
}



TH2F* EffVar2D::GetPassHist(){return m_passhist;}
TH2F* EffVar2D::GetFailHist(){return m_failhist;}
TH2F* EffVar2D::GetMeanTotHist(){return m_meantot;}
TH2F* EffVar2D::GetMeanPassHist(){return m_meanpass;}

TObjArray* EffVar2D::GetTotCBFits(){return m_totCBFits;}
TObjArray* EffVar2D::GetPassCBFits(){return m_passCBFits;}
TObjArray* EffVar2D::GetPassHists(){return m_passhists;}
TObjArray* EffVar2D::GetTotHists(){return m_tothists;}
TObjArray* EffVar2D::GetFailHists(){return m_failhists;}
TObjArray* EffVar2D::GetEffGraphs(){return m_effgraphs;}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems){
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim){
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
  
}
