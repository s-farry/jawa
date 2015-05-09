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
#include <MWTemplate.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

MWTemplate::MWTemplate(string name): Template(name){};

MWTemplate::MWTemplate(string name, PyObject* pyt, PyObject* pycut) : Template(name, pyt, pycut){};

MWTemplate::MWTemplate(string name, TTree* t, TCut cut): Template(name, t, cut){
};

MWTemplate::MWTemplate(string name, TTree* t, TCut cut, enum EColor color) : Template(name, t, cut, color){};

MWTemplate::MWTemplate(string name, MWTemplate* A, MWTemplate* B) : Template(name, A, B){
 
  std::map<string, std::map<string, TH1F*> > varhistsA = A->GetWeightHists();
  std::map<string, std::map<string, TH1F*> > varhistsB = B->GetWeightHists();
  std::map<string, std::map<string, TH2F*> > var2dhistsA = A->Get2DWeightHists();
  std::map<string, std::map<string, TH2F*> > var2dhistsB = B->Get2DWeightHists();

  //fill the maps with the combined histograms

  //1d
  for (std::map<string, std::map<string, TH1F*> >::iterator im = varhistsA.begin(); im != varhistsA.end(); ++im){
    std::map<string, TH1F*> varmap;
    for (std::map<string, TH1F*>::iterator imm = (*im).second.begin(); imm != (*im).second.end(); ++imm){
      TH1F* histA = (*imm).second;
      TH1F* histB = 0;
      if (varhistsB.find((*im).first) != varhistsB.end() &&
	  varhistsB.at((*im).first).find((*imm).first)!= varhistsB.at((*im).first).end()){
	histB = varhistsB.at((*im).first).at((*imm).first);
      }
      if (histA && histB){
	string title(histA->GetTitle());
	title += "_combined";
	TH1F* histC = (TH1F*)histA->Clone(title.c_str());
	histC->Add(histA);
	varmap.insert(pair<string, TH1F*>((*imm).first, histC));
      }
    }
    m_varhists.insert(pair<string, map<string,TH1F*> >((*im).first, varmap));
  }

  //2d
  for (std::map<string, std::map<string, TH2F*> >::iterator im = var2dhistsA.begin(); im != var2dhistsA.end(); ++im){
    std::map<string, TH2F*> varmap;
    for (std::map<string, TH2F*>::iterator imm = (*im).second.begin(); imm != (*im).second.end(); ++imm){
      TH2F* histA = (*imm).second;
      TH2F* histB = 0;
      if (var2dhistsB.find((*im).first) != var2dhistsB.end() &&
	  var2dhistsB.at((*im).first).find((*imm).first)!= var2dhistsB.at((*im).first).end()){
	histB = var2dhistsB.at((*im).first).at((*imm).first);
      }
      if (histA && histB){
	string title(histA->GetTitle());
	title += "_combined";
	TH2F* histC = (TH2F*)histA->Clone(title.c_str());
	histC->Add(histA);
	varmap.insert(pair<string, TH2F*>((*imm).first, histC));
      }
    }
    m_2dvarhists.insert(pair<string, map<string,TH2F*> >((*im).first, varmap));
  }

}

void MWTemplate::FillVars(){
  map<string, int> tree_idx;
  
  for (unsigned int ti = 0; ti< m_trees.size(); ++ti){
    Tree* tree = m_trees.at(ti);
    TTree* t = tree->GetTTree();
    
    t->SetBranchStatus("*" , 0 );
    
    for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
      std::vector<string> varnames = (*iv).second->GetVarNames();
      tree->SetBranches(varnames);
    }
    
    if (m_verbose) cout<<"Branches set for vars"<<endl;
    
    for (std::vector<ReweightVar>::iterator iv = m_reweightvariables.begin(); iv!=m_reweightvariables.end();++iv){
      vector<string> rwnames;
      if ((*iv).GetExpr()){
	rwnames = (*iv).GetExpr()->GetVarNames();
      }
      else{
	rwnames = (*iv).GetNames();
      }
      tree->SetBranches(rwnames);
    }
    
    for (std::map<string, ReweightVar*>::iterator im = m_mwreweightvars.begin(); im != m_mwreweightvars.end();++im){
      vector<string> rwnames;
      if ((*im).second->GetExpr()){
	rwnames = (*im).second->GetExpr()->GetVarNames();
      }
      else{
	rwnames = (*im).second->GetNames();
      }
      tree->SetBranches(rwnames);
    }
    
    if (m_verbose) cout<<"Branches set for reweighted variables"<<endl;

    //Run over entries

    if (m_entryLists.size() == 0) cerr<<"Cuts haven't been applied"<<endl;

    Long64_t nentries = m_entryLists.at(ti)->GetN();
    //Long64_t nb = 0;
    
    if (m_verbose) cout<<"Looping over events"<<endl;

    for (Long64_t jentry = 0 ; jentry < nentries ; jentry++) {
      map<string, double> weights;

      if (jentry%10000==0 && m_outputevts) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
      int entry = m_entryLists.at(ti)->GetEntry(jentry);
      t->GetEntry(entry);

      double tw = m_trees.at(ti)->GetWeight();

      // Get The Weight
      for (std::map<string, ReweightVar*>::iterator iv = m_mwreweightvars.begin(); iv!=m_mwreweightvars.end();++iv){
	double w = tw;
	if ((*iv).second->GetNames().size() == 1){
	  string var = (*iv).second->GetName();
	  Expr* e  = (*iv).second->GetExpr();
	  double val = tree->GetVal(e);
	  w = w * ((*iv).second->GetWeight(val));
	}
	weights[(*iv).first] = w;
      }
      m_normN += tw;
      
      // Loop through 1-D variables

      for (std::map<string,Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
	Var* var = (*iv).second;
	double output = tree->GetVal(var->GetExpr());
	var->FillHist(output);
	//Fill with different weights
	for (std::map<string, double>::iterator iw = weights.begin(); iw != weights.end(); ++iw){
	  m_varhists[(*iv).first][(*iw).first]->Fill(output, (*iw).second);
	}
      }
      
      // Loop through 2-D variables
      
      for (std::map<string,Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
	Var2D* var = (*iv).second;
	double xval = 0.0, yval = 0.0;
	Var* var1 = var->GetVar1();
	Var* var2 = var->GetVar2();
	xval = tree->GetVal(var1->GetExpr());
	yval = tree->GetVal(var2->GetExpr());
	(*iv).second->FillHist(xval, yval);
	for (std::map<string, double>::iterator iw = weights.begin(); iw != weights.end(); ++iw){
	  m_2dvarhists[(*iv).first][(*iw).first]->Fill(xval, yval, (*iw).second);
	}

      }

      // Loop through 3-D variables
      for (std::map<string,Var3D*>::iterator iv = m_3Dvariables.begin() ; iv != m_3Dvariables.end() ; ++iv ){
	Var3D* var = (*iv).second;
	double xval = 0.0, yval = 0.0, zval = 0.0;
	Var* var1 = var->GetVar1();
	Var* var2 = var->GetVar2();
	Var* var3 = var->GetVar3();
	xval = tree->GetVal(var1->GetExpr());
	yval = tree->GetVal(var2->GetExpr());
	zval = tree->GetVal(var3->GetExpr());
	(*iv).second->FillHist(xval, yval, zval);
	//j++;
      }
    }
    t->SetBranchStatus("*",1);
  }
  if ( m_verbose ) cout<<"Finished looping over events"<<endl;
}

void MWTemplate::AddWeight(string wName, string weightname) {
  //m_reweightvariables.push_back(ReweightVar(weightname));
  m_mwreweightvars[wName] = new ReweightVar(weightname);
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    string name = (*iv).first;
    Var* var = (*iv).second;
    TH1F* varhist0 = var->GetHist();
    string varhist0_name = varhist0->GetName();
    m_varhists[name][wName] = ((TH1F*)varhist0->Clone( (varhist0_name + wName).c_str() ));
  }
  for (std::map<string, Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
    string name = (*iv).first;
    Var2D* var = (*iv).second;
    TH2F* varhist0 = var->GetHist();
    string varhist0_name = varhist0->GetName();
    m_2dvarhists[name][wName] = ((TH2F*)varhist0->Clone( (varhist0_name + wName).c_str() ));
  }
}

void MWTemplate::SaveToCurrentFile(){
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    Var* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    map<string, TH1F*>& varwhists = m_varhists[var->GetName()];
    for (map<string, TH1F*>::iterator im = varwhists.begin(); im != varwhists.end(); ++im){
      //char wNumber[4];
      //sprintf(wNumber, "%i", i);
      //std::string wN(wNumber);
      string name = var->GetName() + "_" + (*im).first;
      m_varhists[(*iv).first][(*im).first]->Write(name.c_str());
    }
  }

  for (std::map<string, Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
    Var2D* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    map<string, TH2F*>& varwhists = m_2dvarhists[var->GetName()];
    for (map<string, TH2F*>::iterator im = varwhists.begin(); im != varwhists.end(); ++im){
      //char wNumber[4];
      //sprintf(wNumber, "%i", i);
      //std::string wN(wNumber);
      string name = var->GetName() + "_" + (*im).first;
      m_2dvarhists[(*iv).first][(*im).first]->Write(name.c_str());
    }
  }

  
  TParameter<double>* totEvts     = new TParameter<double>("",m_evts);
  TParameter<double>* normEvts    = new TParameter<double>("", m_normN);
  TParameter<double>* fitFrac     = new TParameter<double>("", m_fitFrac);
  TParameter<double>* fitFracErr  = new TParameter<double>("", m_fitFracErr);
  TObjString* cut = new TObjString(m_selcut.GetTitle());
  totEvts->Write("TotEvts");
  normEvts->Write("NormEvts");
  fitFrac->Write("FitFrac");
  fitFracErr->Write("FitFracErr");
  cut->Write("Cut");
  if (m_tree) m_tree->Write("tree");
}

void MWTemplate::SaveToFile(){
  cout<<"Name is "<<m_name<<" (MWTemplate)"<<endl;
  cout<<"Saving "<<m_name<<" to file"<<endl;
  TFile* f = new TFile((m_name+".root").c_str(),"RECREATE");
  
  SaveToCurrentFile();
  cout<<"Written"<<endl;
  //f->Write();
  gROOT->cd();
  f->Close();
}

TH1F* MWTemplate::GetWeightHist(string var, string wname){
  if (m_varhists.find(var) == m_varhists.end()){
    cout<<"Variable "<<var<<" not found - cannot get weighted histogram "<<wname<<endl;
    return 0;
  }
  if (m_varhists.at(var).find(wname) == m_varhists.at(var).end()){
    cout<<"Variable "<<var<<" found but cannot get weighted histogram "<<wname<<endl;
    return 0;
  }
  return m_varhists[var][wname];
}

TH2F* MWTemplate::GetWeight2DHist(string var, string wname){
  if (m_2dvarhists.find(var) == m_2dvarhists.end()){
    cout<<"Variable "<<var<<" not found - cannot get weighted histogram "<<wname<<endl;
    return 0;
  }
  if (m_2dvarhists.at(var).find(wname) == m_2dvarhists.at(var).end()){
    cout<<"Variable "<<var<<" found but cannot get weighted histogram "<<wname<<endl;
    return 0;
  }
  return m_2dvarhists[var][wname];
}

void MWTemplate::ScaleAllWeights(double w){
  for (std::map<string, std::map<string, TH1F*> >::iterator im = m_varhists.begin(); im != m_varhists.end(); ++im){
    for (std::map<string, TH1F*>::iterator ih = (*im).second.begin(); ih != (*im).second.end(); ++ih){
      (*ih).second->Scale(w);
    }
  }
  for (std::map<string, std::map<string, TH2F*> >::iterator im = m_2dvarhists.begin(); im != m_2dvarhists.end(); ++im){
    for (std::map<string, TH2F*>::iterator ih = (*im).second.begin(); ih != (*im).second.end() ; ++ih){
      (*ih).second->Scale(w);
    }
  }
}

void MWTemplate::ScaleWeight(string w, double s){
  for (std::map<string, std::map<string, TH1F*> >::iterator im = m_varhists.begin(); im != m_varhists.end(); ++im){
    if ((*im).second.find(w) != (*im).second.end()){
      (*im).second.at(w)->Scale(s);
    }
    else{
      cout<<"Weight "<<w<<" does not exist"<<endl;
    }
  }
  for (std::map<string, std::map<string, TH2F*> >::iterator im = m_2dvarhists.begin(); im != m_2dvarhists.end(); ++im){
    if ((*im).second.find(w) != (*im).second.end()){
      (*im).second.at(w)->Scale(s);
    }
    else{
      cout<<"Weight "<<w<<" does not exist"<<endl;
    }
  }
}


PyObject* MWTemplate::GetWeightHist_py(string var, string wname){
  TH1F* h = GetWeightHist(var, wname);
  if (h) {
    return TPython::ObjectProxy_FromVoidPtr(h, h->ClassName());
  }
  else{
    return 0;
  }
}
PyObject* MWTemplate::GetWeight2DHist_py(string var, string wname){
  TH2F* h2 = GetWeight2DHist(var, wname);
  if(h2) {
    return TPython::ObjectProxy_FromVoidPtr(h2, h2->ClassName());
  }
  else{
    return 0;
  }
}
std::map<string, std::map<string, TH1F*> > MWTemplate::GetWeightHists(){ return m_varhists;}
std::map<string, std::map<string, TH2F*> > MWTemplate::Get2DWeightHists(){return m_2dvarhists;}
