#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <EfficiencyClass.h>
#include <TEfficiency.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>
#include <TParameter.h>
#include <Utils.h>

using namespace std;

EfficiencyClass::EfficiencyClass(string name) : JawaObj("EfficiencyClass", name){
  m_verbose = false;
  m_fillbkg       = false;
  m_scaleerrs = false;

  m_tot  = 0;
  m_pass = 0;
  m_fail = 0;

  m_Npass_rweff_varyhi = new TObjArray();
  m_Npass_rweff_varylo = new TObjArray();


}

EfficiencyClass::EfficiencyClass( string name , EfficiencyClass* classA , EfficiencyClass* classB ) : JawaObj("EfficiencyClass", name){
  //cout<<"Combining efficiency classes"<<endl;
  m_verbose       = false;
  m_fillbkg       = false;
  
  std::map<string, EffVar*>   classAvars   = classA->m_variables;
  std::map<string, EffVar*>   classBvars   = classB->m_variables;
  std::map<string, EffVar2D*> classA2Dvars = classA->m_2Dvariables;
  std::map<string, EffVar2D*> classB2Dvars = classB->m_2Dvariables;
  
  m_trees.insert(m_trees.end(), classA->m_trees.begin(), classA->m_trees.end());
  m_trees.insert(m_trees.end(), classB->m_trees.begin(), classB->m_trees.end());
  
  if ((classAvars.size() != classBvars.size()) || ( classA2Dvars.size() != classB2Dvars.size() ) ) return;
  
  for ( std::map<string, EffVar*>::iterator ie = classAvars.begin(); ie != classAvars.end(); ++ie )
    {
      EffVar* varA = ie->second;
      EffVar* varB = classBvars.at(varA->GetName());
      if (varA->GetName() == varB->GetName())  {
	EffVar* combVar = new EffVar(varA, varB, name);
	m_variables.insert(std::pair<string, EffVar*>(varA->GetName(),combVar));
      }
    }
  
  for ( std::map<string, EffVar2D*>::iterator ie = classA2Dvars.begin(); ie != classA2Dvars.end(); ++ie )
    {
      EffVar2D* varA = ie->second;
      EffVar2D* varB = classB2Dvars.at(varA->GetName());
      if (varA->GetName() == varB->GetName())  {
	EffVar2D* combVar = new EffVar2D(varA, varB, name);
	m_2Dvariables.insert(std::pair<string, EffVar2D*>(varA->GetName(),combVar));
      }
    }
  
  double totN = classA->m_Ntot + classB->m_Ntot;
  double passN = classA->m_Npass + classB->m_Npass;


  m_Npass_rweff_varyhi = new TObjArray();
  m_Npass_rweff_varylo = new TObjArray();

  for (int i = 0 ; i < classA->m_Npass_rweff_varyhi->GetEntries() ; ++i){
    ostringstream a, b;
    a<<((TParameter<double>*)(classA->m_Npass_rweff_varyhi)->At(i))->GetName()<<"_combined";
    b<<((TParameter<double>*)(classA->m_Npass_rweff_varylo)->At(i))->GetName()<<"_combined";
    double val1 = ((TParameter<double>*)classA->m_Npass_rweff_varyhi->At(i))->GetVal();
    double val2 = ((TParameter<double>*)classB->m_Npass_rweff_varyhi->At(i))->GetVal();
    m_Npass_rweff_varyhi->Add(new TParameter<double>(a.str().c_str(), val1+val2));
    double val3 = ((TParameter<double>*)classA->m_Npass_rweff_varylo->At(i))->GetVal();
    double val4 = ((TParameter<double>*)classB->m_Npass_rweff_varylo->At(i))->GetVal();
    m_Npass_rweff_varylo->Add(new TParameter<double>(b.str().c_str(), val3+val4));
  }

  m_Ntot = totN;
  m_Npass = passN;
  
  m_toteff = Eff("TotalEff", totN, passN);
  
}

void EfficiencyClass::CorrectGraphs(EfficiencyClass* classA, EfficiencyClass* classB, string opt){
  std::map<string, EffVar*>   classAvars   = classA->m_variables;
  std::map<string, EffVar*>   classBvars   = classB->m_variables;
  std::map<string, EffVar2D*> classA2Dvars = classA->m_2Dvariables;
  std::map<string, EffVar2D*> classB2Dvars = classB->m_2Dvariables;

  if (opt == "M" )
    {
      double effA = classA->m_toteff.GetEff();
      double effB = classB->m_toteff.GetEff();
      double effAerrhi = classA->m_toteff.GetEffErrHi();
      double effAerrlo = classA->m_toteff.GetEffErrLo();
      double effBerrhi = classB->m_toteff.GetEffErrHi();
      double effBerrlo = classB->m_toteff.GetEffErrLo();
      double eff = effA*effB;
      double errhi =  sqrt(pow(effAerrhi/effA,2)+ pow(effBerrhi/effB,2))*eff;
      double errlo =  sqrt(pow(effAerrlo/effA,2)+ pow(effBerrlo/effB,2))*eff;
      m_toteff = Eff("TotEff", eff, errhi, errlo);
    }
  else {
      double effA = classA->m_toteff.GetEff();
      double effB = classB->m_toteff.GetEff();
      double effAerrhi = classA->m_toteff.GetEffErrHi();
      double effAerrlo = classA->m_toteff.GetEffErrLo();
      double effBerrhi = classB->m_toteff.GetEffErrHi();
      double effBerrlo = classB->m_toteff.GetEffErrLo();
      double eff = effA/effB;
      double errhi =  sqrt(pow(effAerrhi/effA,2)+ pow(effBerrhi/effB,2))*eff;
      double errlo =  sqrt(pow(effAerrlo/effA,2)+ pow(effBerrlo/effB,2))*eff;
      m_toteff = Eff("TotEff", eff, errhi, errlo);
  }

  for ( std::map<string, EffVar*>::iterator ie = classAvars.begin(); ie != classAvars.end(); ++ie )
    {
      EffVar* varA = ie->second;
      if (classBvars.find(varA->GetName()) != classBvars.end()){
	EffVar* varB = classBvars.at(varA->GetName());
	if (varA->GetName() == varB->GetName())  {
	  TGraphAsymmErrors* graph_comb;
	  if (opt == "M") {
	    std::vector<TGraphAsymmErrors*> vect;
	    vect.push_back(varA->GetEffGraph());
	    vect.push_back(varB->GetEffGraph());
	    graph_comb = CombineTGraphs(vect);
	  }
	  else {
	    graph_comb = DivideTGraphs(varA->GetEffGraph(), varB->GetEffGraph());
	  }
	  EffVar* combVar = new EffVar(varA->GetName(), graph_comb);
	  m_variables.insert(std::pair<string, EffVar*>(varA->GetName(),combVar));
	}
      }
    }
  for (std::map<string, EffVar2D*>::iterator ie = classA2Dvars.begin(); ie != classA2Dvars.end(); ++ie)
    {
      EffVar2D* varA = ie->second;
      if (classB2Dvars.find(varA->GetName()) != classB2Dvars.end()){
	EffVar2D* varB = classB2Dvars.at(varA->GetName());
	if (varA->GetName() == varB->GetName())  {
	  TH2F* graph_comb;
	  if (opt == "M") {
	    std::vector<TH2F*> vect;
	    vect.push_back(varA->Get2DEffGraph());
	    vect.push_back(varB->Get2DEffGraph());
	    graph_comb = CombineTHists(vect);
	  }
	  else {
	    info()<<"Name is: "<<varA->GetName()<<endl;
	    graph_comb = DivideTHists(varA->Get2DEffGraph(), varB->Get2DEffGraph());
	  }
	  EffVar2D* combVar = new EffVar2D(varA->GetName(), graph_comb);
	  m_2Dvariables.insert(std::pair<string, EffVar2D*>(varA->GetName(),combVar));
	}
      }
    }
}


void EfficiencyClass::SetVariables(std::map<string, EffVar*> variables){
  m_variables = variables;
  for (std::map<string, EffVar*>::iterator ie = m_variables.begin(); ie != m_variables.end(); ++ie){
    ie->second->SetPrefix(m_name);
  }
}

std::map<string, EffVar*> EfficiencyClass::GetVariables(){
  return m_variables;
}

void EfficiencyClass::SetTree(TTree* tree){
  tree->SetEntryList(0);
  tree->SetBranchStatus("*", 1);
  m_trees.push_back(new Tree("", tree, 1.0));
}

void EfficiencyClass::SetVerbose(bool verbose){
  m_verbose = verbose;
}
bool EfficiencyClass::GetVerbose(){
  return m_verbose;
}

void EfficiencyClass::SetEffRange(double lo, double hi){
  m_efflo = lo;
  m_effhi = hi;

}
void EfficiencyClass::SetTrees(TTree* priTree, TTree* secTree){
  priTree->SetEntryList(0);
  priTree->SetBranchStatus("*",1);
  m_trees.push_back(new Tree("", priTree, 1.0));
  
  secTree->SetEntryList(0);
  secTree->SetBranchStatus("*",1);
  m_trees.push_back(new Tree("", secTree, 1.0));

}

void EfficiencyClass::AddTree(TTree* tree){
  tree->SetEntryList(0);
  tree->SetBranchStatus("*",1);
  m_trees.push_back(new Tree("", tree, 1.0));
}
void EfficiencyClass::AddTrees(std::vector<TTree*> trees){
  for (unsigned int i = 0 ; i < trees.size() ; ++i) AddTree(trees.at(i));
}
Eff EfficiencyClass::GetTotEff(){
  return m_toteff;
}

void EfficiencyClass::AddVar(string name, string var, int bins, double lo, double hi){
  EffVar* variable  = new EffVar(name, var, bins, lo, hi, m_name);
  m_variables.insert(std::pair<string, EffVar*>(name, variable));
}
void EfficiencyClass::AddVar(string name, string var, vector<double> edges){
  EffVar* variable  = new EffVar(name, var, edges, m_name);
  m_variables.insert(std::pair<string, EffVar*>(name, variable));
}

void EfficiencyClass::Add2DVar(string var1, string var2, string name){
  EffVar* varA = m_variables.at(var1);
  EffVar* varB = m_variables.at(var2);
  if (name=="")  name = varA->GetName() + "_" + varB->GetName();
  verbose()<<"Adding 2D variable with name: "<<name<<endl;
  EffVar2D* variable  = new EffVar2D(name, varA, varB, m_name);
  m_2Dvariables.insert(std::pair<string, EffVar2D*>(name, variable));
  verbose()<<"Variable added"<<endl;
 
}


void EfficiencyClass::AddPassVar(const char* var){
  passvars.push_back(var);
}
void EfficiencyClass::SetPassCut(TCut cut){
  m_passcut = cut;
}

string EfficiencyClass::GetRootName(const char* file){
  string output = "";
  if ( strlen( file ) == 0 ) {
    if ( m_name.size() == 0 ) 
      {
	output = "output.root";
      }
    else 
      {
	char buffer[50];
	sprintf(buffer, "%s%s", m_name.c_str(), ".root");
	//ostringstream filess;
	//filess<<m_name<<".root";
	//output = filess.str();
	output=buffer;
      }
  }
  else{
    output = file;
  }
  return output;
}

void EfficiencyClass::AddFits(TObjArray* hists, TObjArray* fits){
  if (hists->GetEntries() != fits->GetEntries() ) return;
  if (fits->GetEntries() == 0)
  for (int i = 0; i < hists->GetEntries(); ++i){
    TH1F* hist = (TH1F*)hists->At(i);
    TF1* fit = (TF1*)fits->At(i);
    //fit->SetOwner(true);
    if (fit) hist->GetListOfFunctions()->Add(fit);
  }
}
void EfficiencyClass::AddFits(TObjArray* hists, TObjArray* fitsCB, TObjArray* fitsBkg){
  if (hists->GetEntries() != fitsCB->GetEntries() || hists->GetEntries() != fitsBkg->GetEntries() ) return;
  if (fitsCB->GetEntries() == 0)
  for (int i = 0; i < hists->GetEntries(); ++i){
    TH1F* hist = (TH1F*)hists->At(i);
    TF1* fitsig = (TF1*)fitsCB->At(i);
    TF1* fitbkg = (TF1*)fitsBkg->At(i);
    
    //fitsig->SetOwner(true);
    //fitbkg->SetOwner(true);
    if (fitsig && fitbkg) {
      hist->GetListOfFunctions()->Add(fitsig);
      hist->GetListOfFunctions()->Add(fitbkg);
    }
  }
}

void EfficiencyClass::SaveToFile(const char* file){
  string output = GetRootName(file);
  //const char* cwd = gDirectory->pwd();

  verbose()<<"Outputting to file - "<<output<<endl;
  TFile* f = new TFile(output.c_str(),"RECREATE");
  //gROOT->cd();
  verbose()<<"Writing - opened file"<<endl;
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){


    if (m_verbose) cout<<"Writing - "<<ei->first<<endl;
    EffVar* evar = ei->second;
    f->cd();
    f->mkdir(evar->GetName().c_str());
    f->cd(evar->GetName().c_str());

    if (evar->GetTotHist())  evar->GetTotHist()->Write("TotalHist");
    if (evar->GetPassHist()) evar->GetPassHist()->Write("PassHist");
    if (evar->GetFailHist()) evar->GetFailHist()->Write("FailHist");

    if (m_scaleerrs){
      if ( evar->GetEffRWVaryHiPassHists() ) {
      evar->GetEffRWVaryHiPassHists()->Write("EffRWVaryHiPassHists",1);
      }
      if ( evar->GetEffRWVaryLoPassHists() ) {
      evar->GetEffRWVaryLoPassHists()->Write("EffRWVaryLoPassHists",1);
      }
    }
    

    if (m_fillbkg && evar->GetBkgTotHist() && evar->GetBkgPassHist() ){
      evar->GetBkgTotHist()->Write("TotalBkg");
      evar->GetBkgPassHist()->Write("PassBkg");
    }


    verbose()<<"Wrote histograms"<<endl;
    
    if (evar->GetEffGraph() == 0) 
      { 
	info()<<"No Eff Graph"<<endl;
      }
    else evar->GetEffGraph()->Write("EfficiencyGraph");
    
    verbose()<<"Wrote efficiency graph"<<endl;

  }
  for (std::map<string, EffVar2D*>::iterator ei = m_2Dvariables.begin(); ei != m_2Dvariables.end(); ++ei){
    EffVar2D* evar = ei->second;
    
    f->cd();
    f->mkdir(evar->GetName().c_str());
    f->cd(evar->GetName().c_str());
    
    if (evar->GetTotHist())  evar->GetTotHist()->Write("TotalHist");
    if (evar->GetPassHist()) evar->GetPassHist()->Write("PassHist");
    if (evar->GetFailHist()) evar->GetFailHist()->Write("FailHist");
    
    if (evar->Get2DEffGraph()) evar->Get2DEffGraph()->Write("EffGraph2D");
    if (evar->GetEffGraphs()) evar->GetEffGraphs()->Write("EffGraphs",1);
    
  }
  
  TParameter<double>* toteff       = new TParameter<double>("TotalEff",      m_toteff.GetEff()      );
  TParameter<double>* toteff_errhi = new TParameter<double>("TotalEffErrHi", m_toteff.GetEffErrHi() );
  TParameter<double>* toteff_errlo = new TParameter<double>("TotalEffErrLo", m_toteff.GetEffErrLo() );

  TParameter<double>* Npass = new TParameter<double>("Npass", m_Npass);
  TParameter<double>* Ntot  = new TParameter<double>("Ntot", m_Ntot);

  f->cd();
  toteff->Write();
  toteff_errhi->Write();
  toteff_errlo->Write();
  Npass->Write();
  Ntot->Write();

  if (m_scaleerrs){
    TParameter<double>* totaleffrwerr = new TParameter<double>("TotalEffRWErr", m_toteffrw_err);
    totaleffrwerr->Write();
    if (m_Npass_rweff_varyhi) m_Npass_rweff_varyhi->Write("Npass_effrwhi",1);
    if (m_Npass_rweff_varylo) m_Npass_rweff_varylo->Write("Npass_effrwlo",1);
  }

  info()<<"Written to file - "<<output<<endl;
  gROOT->cd();
  f->Close();
}

void EfficiencyClass::FillBkgHists(){
  m_fillbkg = true;
  
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    evar->FillBkgHists(m_efflo, m_effhi);
  }
}

void EfficiencyClass::LoadFromFile(const char* file){
  string output = GetRootName(file);
  TFile* f = new TFile(output.c_str());
  gROOT->cd();
  TList* list = f->GetListOfKeys();

  //First let's get the total efficiency
  TParameter<double>* toteff = (TParameter<double>*)f->Get("TotalEff");
  TParameter<double>* toteff_errhi = (TParameter<double>*)f->Get("TotalEffErrHi");
  TParameter<double>* toteff_errlo = (TParameter<double>*)f->Get("TotalEffErrLo");

  if (toteff && toteff_errhi && toteff_errlo){
    m_toteff = Eff("TotalEff", toteff->GetVal(), toteff_errhi->GetVal(), toteff_errlo->GetVal());
  }


  //Now loop through and get variables
  for (int i = 0; i < list->GetEntries(); ++i){
    string name = list->At(i)->GetName();
    bool is2D = ((int)name.find("_")!=-1);
    
    if (name != "TotalEff" && name!= "TotalEffErrHi" && name!="TotalEffErrLo" 
	&& name != "TotalHist" && name != "PassHist" && name != "FailHist"
	&& f->cd(name.c_str())){
      if (!is2D) m_variables.insert(std::pair<string, EffVar*>(name, new EffVar(name, f, m_name)));
      else if (is2D) m_2Dvariables.insert(std::pair<string, EffVar2D*>(name, new EffVar2D(name, f, m_name)));
    }
  }
  gROOT->cd();
  //This unfortunately causes a seg fault, not sure why as the cd above should move the directory
  //f->Close();
}

void EfficiencyClass::StripTree(TCut cut){
  m_trees.at(0)->GetTTree()->Draw(">>myList",cut,"entrylist");
  TEntryList* list = (TEntryList*)gDirectory->Get("myList");
  entryListTot.push_back(list);
  //trees.at(0)->SetEntryList(list);
  //t_sel = t_tot->CopyTree(cut);
  //t=t_sel;
}

void EfficiencyClass::StripTrees(TCut cut){
  m_trees.at(0)->GetTTree()->Draw(">>myListA",cut,"entrylist");
  TEntryList* priList = (TEntryList*)gDirectory->Get("myListA");
  entryListTot.push_back(priList);

  m_trees.at(1)->GetTTree()->Draw(">>myListB",cut,"entrylist");
  TEntryList* secList = (TEntryList*)gDirectory->Get("myListB");
  //trees.at(1)->SetEntryList(secList);
  entryListTot.push_back(secList);

}

void EfficiencyClass::SetSelectionCut(TCut cut){
  m_selcut = cut;

}

void EfficiencyClass::MakeEntryLists(){
  for (unsigned int i = 0; i < m_trees.size(); ++i){
    ostringstream ss;
    ss<<">>myList"<<i;
    string label = ss.str();
    string trimlabel = ss.str();
    trimlabel.erase(0,2);
    m_trees.at(i)->GetTTree()->Draw(label.c_str(), m_selcut , "entrylist");
    const TEntryList* list = (TEntryList*)gDirectory->Get(trimlabel.c_str());
    entryListTot.push_back(new TEntryList(*list));

    if (m_passcut != ""){
      ostringstream ss_p;
      ss_p<<">>myPassLabel"<<i;
      string passlabel = ss_p.str();
      string passtrimlabel = ss_p.str();
      passtrimlabel.erase(0,2);
      m_trees.at(i)->GetTTree()->Draw(passlabel.c_str(), m_selcut && m_passcut , "entrylist");
      const TEntryList* passlist = (TEntryList*)gDirectory->Get(passtrimlabel.c_str());
      entryListPass.push_back(new TEntryList(*passlist));

    }
  }
}

Eff EfficiencyClass::GetEfficiency(const char* varname, TH1F* hist){
  EffVar* evar = m_variables.at(varname);
  ostringstream ss;
  ss<<varname<<"Eff";
  string name = ss.str();
  if (hist->GetNbinsX() == evar->GetBins() &&
      hist->GetXaxis()->GetXmin() == evar->GetLo() &&
      hist->GetXaxis()->GetXmax() == evar->GetHi())
    {
      TGraphAsymmErrors* effgraph = evar->GetEffGraph();
      double N = hist->Integral(1,hist->GetNbinsX());
      double denom = 0;
      double errnumhi = 0;
      double errnumlo = 0;
      for (int i = 0; i < hist->GetNbinsX();++i){
	double x, effi, effierrhi, effierrlo;
	effgraph->GetPoint(i, x, effi);
	effierrhi = effgraph->GetErrorYhigh(i);
	effierrlo = effgraph->GetErrorYlow(i);
	double Ni = hist->GetBinContent(i+1);
	denom += (Ni/effi);
	errnumhi += pow(N*Ni*effierrhi/pow(effi,2),2);
	errnumlo += pow(N*Ni*effierrlo/pow(effi,2),2);
      }
      double eff = N/denom;
      double efferrhi = sqrt(errnumhi)/pow(denom,2);
      double efferrlo = sqrt(errnumlo)/pow(denom,2);

      return Eff(name.c_str(), eff, efferrhi, efferrlo);
    }
  else {
    info()<<"Histogram and Variable do not match";
    return Eff(name.c_str(), 0., 0., 0.);
  }

}

void EfficiencyClass::SetBranches(Tree* t){
  t->GetTTree()->SetBranchStatus("*",0);
  
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    t->SetBranches(evar->GetVarNames());
  }
  for (std::vector<ReweightVar*>::iterator iv = m_reweightvariables.begin(); iv!=m_reweightvariables.end();++iv){
    vector<string> rwnames;
    vector<Expr*> exprs = (*iv)->GetExprs();
    for (vector<Expr*>::iterator ie = exprs.begin(); ie != exprs.end(); ++ie){
      vector<string> names = (*ie)->GetVarNames();
      rwnames.insert(rwnames.end(), names.begin(), names.end());
    } 
    t->SetBranches(rwnames);
  }
  for (std::vector<ReweightVar*>::iterator iv = m_reweighteffvariables.begin(); iv!=m_reweighteffvariables.end();++iv){
    vector<string> rwnames;
    vector<Expr*> exprs = (*iv)->GetExprs();
    for (vector<Expr*>::iterator ie = exprs.begin(); ie != exprs.end(); ++ie){
      vector<string> names = (*ie)->GetVarNames();
      rwnames.insert(rwnames.end(), names.begin(), names.end());
    } 
    t->SetBranches(rwnames);
  }
}
void EfficiencyClass::FreeBranches(Tree* t){
  //Just in a function for consistency
  t->GetTTree()->SetBranchStatus("*",1);
}

pair<Utils::weight,Utils::weight> EfficiencyClass::FillVars(bool pass, Tree* t){
  //Initialise weights to 1
  double w = 1.0;
  double w_err = 0.0;
  double effw = 1.0;
  double effw_err = 0.0;
  int effw_bin = -1;
    
  // Get The Weights from Reweight Variables
  for (std::vector<ReweightVar*>::iterator iv = m_reweightvariables.begin(); iv!=m_reweightvariables.end();++iv){
    if ((*iv)->GetExprs().size() == 1){
      Expr* e  = (*iv)->GetExprs()[0];
      double val = t->GetVal(e);
      w = w * ((*iv)->GetWeight(val));
    }
    else if ((*iv)->GetExprs().size() == 2){
      Expr* e1 = (*iv)->GetExprs().at(0);
      Expr* e2 = (*iv)->GetExprs().at(1);
      double val1 = t->GetVal(e1);
      double val2 = t->GetVal(e2);
      w = w * ((*iv)->GetWeight(val1, val2));
    }
    else if ((*iv)->GetExprs().size() == 4){
      Expr* e1 = (*iv)->GetExprs().at(0);
      Expr* e2 = (*iv)->GetExprs().at(1);
      Expr* e3 = (*iv)->GetExprs().at(2);
      Expr* e4 = (*iv)->GetExprs().at(3);
      double val1 = t->GetVal(e1);
      double val2 = t->GetVal(e2);
      double val3 = t->GetVal(e3);
      double val4 = t->GetVal(e4);
      w = w * ((*iv)->GetWeight(val1, val2, val3, val4));
      w_err = sqrt(w_err*w_err + pow((*iv)->GetWeightErr(val1,val2),2));
    }
  }

  //Get any efficiency reweight value
  for (std::vector<ReweightVar*>::iterator iv = m_reweighteffvariables.begin(); iv!=m_reweighteffvariables.end();++iv){
    if ((*iv)->GetExprs().size() == 2){
      Expr* e1 = (*iv)->GetExprs().at(0);
      Expr* e2 = (*iv)->GetExprs().at(1);
      double val1 = t->GetVal(e1);
      double val2 = t->GetVal(e2);
      effw = effw * ((*iv)->GetWeight(val1, val2));
      effw_err = sqrt(effw_err*effw_err + pow((*iv)->GetWeightErr(val1,val2),2));
      effw_bin = (*iv)->GetBin(val1, val2);
    }
  }

  Utils::weight eff_weight, scale_weight;
  eff_weight.val = effw;
  eff_weight.err = effw_err;
  eff_weight.bin = effw_bin;
  scale_weight.val = w;
  scale_weight.err = w_err;

  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    if (!m_scaleerrs){
      evar->FillVar(pass, t->GetVal(evar->GetExpr()), w, effw);
    }
    else{
      evar->FillVar(pass, t->GetVal(evar->GetExpr()), scale_weight, eff_weight );
    }
  }
  
  for (std::map<string, EffVar2D*>::iterator ei = m_2Dvariables.begin(); ei != m_2Dvariables.end(); ++ei){
    EffVar2D* evar = ei->second;
    double val1 = t->GetVal(evar->GetVar1()->GetExpr());
    double val2 = t->GetVal(evar->GetVar2()->GetExpr());
    evar->FillVar(pass, val1, val2, w);
  }
  //return the weights for use in rest of program
  pair<Utils::weight, Utils::weight> weights(scale_weight, eff_weight);
  return weights;
}

bool EfficiencyClass::VarExists(string var){
  bool exists = false;
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    if (evar->GetVar() == var || evar->GetName() == var ) exists = true;
  }
  return exists;
}

void EfficiencyClass::PrintVars(){
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    cout<<evar->GetName()<<endl;
  }
  for (std::map<string, EffVar2D*>::iterator ei = m_2Dvariables.begin(); ei != m_2Dvariables.end(); ++ei){
    EffVar2D* evar = ei->second;
    cout<<evar->GetName()<<endl;
  }
}

void EfficiencyClass::LoopEntries(){
  
  //Now fill bins - set up variables for a plt and pass var and multiple vars

  double totN = 0;
  double passN = 0;


  for(unsigned int i = 0; i < m_trees.size(); i++){
    
    verbose()<<m_name<<": Tree Number: "<<i<<endl;
    Tree* t = m_trees.at(i);
    SetBranches(t);
    
    TEntryList* totList = entryListTot.at(i);
    TEntryList* passList = entryListPass.at(i);

    Long64_t nentries = totList->GetN();
    bool pass = false;
    //Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      pass = false;
      Int_t entry = totList->GetEntry(jentry);
      t->GetEntry(entry);
      //Int_t ientry=t->LoadEntry(entry);
      //if (ientry < 0) break;
      //nb = t->GetEntry(entry);   nbytes += nb;
      if (jentry%10000==0) info()<<"Entry "<<jentry<<" of "<<nentries<<endl;
      //cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
      
      //See if entry passes
      pass = passList->Contains(entry);

      pair<Utils::weight,Utils::weight> weights = FillVars(pass, t);
      totN = totN + weights.first.val;
      //info()<<weights.first.val<<" "<<weights.second.val<<endl;

      if (pass) {
	passN = passN + weights.first.val*weights.second.val;
	if (m_scaleerrs){
	  int bin = weights.second.bin;
	  for (int i = 0 ; i < m_Npass_rweff_varyhi->GetEntries() ; ++i){
	    double valhi = ((TParameter<double>*)m_Npass_rweff_varyhi->At(i))->GetVal(); 
	    double vallo = ((TParameter<double>*)m_Npass_rweff_varylo->At(i))->GetVal(); 
	    if (bin == i ){
	      ((TParameter<double>*)m_Npass_rweff_varyhi->At(i))->SetVal(valhi + weights.first.val*(weights.second.val + weights.second.err));
	      ((TParameter<double>*)m_Npass_rweff_varylo->At(i))->SetVal(vallo + weights.first.val*(weights.second.val - weights.second.err));
	    }
	    else{
	      ((TParameter<double>*)m_Npass_rweff_varyhi->At(i))->SetVal(valhi + weights.first.val*(weights.second.val));
	      ((TParameter<double>*)m_Npass_rweff_varylo->At(i))->SetVal(vallo + weights.first.val*(weights.second.val));
	    }
	  }
	}
      }
    }
    FreeBranches( t );
  }
  
  verbose()<<"Final totN: "<<totN<<endl;
  
  m_Ntot = totN;
  m_Npass = passN;
  
  m_toteff = Eff("TotalEff", totN, passN);
  
}

void EfficiencyClass::Normalise(double N){
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    evar->Normalise(N);
  }
  
  for (std::map<string, EffVar2D*>::iterator ei = m_2Dvariables.begin(); ei != m_2Dvariables.end(); ++ei){
    EffVar2D* evar = ei->second;
    evar->Normalise(N);
  }
}


void EfficiencyClass::MakeEfficiencyGraph(){
  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    EffVar* evar = ei->second;
    evar->MakeTGraph();
    verbose()<<"Set TGraph"<<endl;
  }

  for (std::map<string, EffVar2D*>::iterator ei = m_2Dvariables.begin(); ei != m_2Dvariables.end(); ++ei){
    EffVar2D* evar = ei->second;
    evar->MakeEffHist();
    evar->MakeTGraphs();
  }

  if (m_scaleerrs){
    m_toteffrw_err = 0;
    for (int i = 0 ; i < m_Npass_rweff_varyhi->GetEntries() ; ++i ){
      double effhi = ((TParameter<double>*)m_Npass_rweff_varyhi->At(i))->GetVal() / m_Ntot;
      double efflo = ((TParameter<double>*)m_Npass_rweff_varylo->At(i))->GetVal() / m_Ntot;
      m_toteffrw_err += max (abs(effhi - m_Npass/m_Ntot), abs(efflo - m_Npass/m_Ntot));
    }
  }
  verbose()<<"leaving"<<endl;
}

std::pair<TF1*,TF1*> EfficiencyClass::FitHistogram(TH1F* massplot, double lo, double hi, string opt){
  if (massplot->GetListOfFunctions()->GetSize() > 0) massplot->GetListOfFunctions()->Clear();

  double max = massplot->GetMaximum();

  TF1* pol1 = new TF1("poly1", "pol1", lo, hi);

  massplot->Fit("poly1","QN");

  TF1* gausline = new TF1("gaussline",fitGaussLine, lo, hi, 5);
  gausline->SetParName(0, "norm"  );
  gausline->SetParName(1, "mean"  );
  gausline->SetParName(2, "sigma" );
  gausline->SetParName(3, "c"     );
  gausline->SetParName(4, "slope" );

  gausline->SetParameter(0, max);
  if ((int)opt.find("Jpsi")!=-1){
    gausline->SetParameter(1, 3100);
    gausline->SetParameter(2, 50);
  }
  else if ((int)opt.find("Z0")!=-1){
    gausline->SetParameter(1,90000);
    gausline->SetParameter(2,2000);
  }
  gausline->SetParameter(3, pol1->GetParameter(0));
  gausline->SetParameter(4, pol1->GetParameter(1));
  massplot->Fit(gausline,"QN");

  TF1* myCB;
  if ((int)opt.find("CB+G")!=-1){
    myCB = new TF1("myCB", fitCBGauss, lo , hi , 10);
    myCB->SetParName(0, "norm"   );
    myCB->SetParName(1, "alpha"  );
    myCB->SetParName(2, "n"      );
    myCB->SetParName(3, "mean"   );
    myCB->SetParName(4, "sigma"  );
    myCB->SetParName(5, "c"      );
    myCB->SetParName(6, "slope"  );
    myCB->SetParName(7, "gaussnorm" );
    myCB->SetParName(8, "gaussmean" );
    myCB->SetParName(9, "gausssigma" );
    myCB->SetParameters( gausline->GetParameter(0) , 1., 1., gausline->GetParameter(1), gausline->GetParameter(2) * 3./4., 
			 gausline->GetParameter(3), gausline->GetParameter(4), gausline->GetParameter(0), gausline->GetParameter(1),
			 gausline->GetParameter(2)/4.);
    myCB->FixParameter(2, 1);
  }
  else {
    myCB = new TF1("myCB", fitCB , lo , hi , 7);
    myCB->SetParName(0, "norm"   );
    myCB->SetParName(1, "alpha"  );
    myCB->SetParName(2, "n"      );
    myCB->SetParName(3, "mean"   );
    myCB->SetParName(4, "sigma"  );
    myCB->SetParName(5, "c"      );
    myCB->SetParName(6, "slope"  );
    
    myCB->SetParameters( gausline->GetParameter(0) , 1., 1., gausline->GetParameter(1), gausline->GetParameter(2), 
			 gausline->GetParameter(3), gausline->GetParameter(4));
    myCB->FixParameter(2, 1);
  }
  
  myCB->SetLineColor(2);

  massplot->Fit(myCB,"Q","", lo, hi);
  
  pol1->SetParameter( 0 , myCB->GetParameter(5) );
  pol1->SetParameter( 1 , myCB->GetParameter(6) );
  pol1->SetLineColor( 4 );
  pol1->SetLineStyle( 7 );

  //Return crystal ball
  massplot->GetListOfFunctions()->Add(pol1);
  std::pair<TF1*, TF1*> pair;
  pair.first = myCB;
  pair.second = pol1;

  return pair;
  //massplot->GetListOfFunctions()->Add(pol1);

}



void EfficiencyClass::Reweight(string var, TF1* f){
  ReweightVar* v = new ReweightVar(var,f);
  m_reweightvariables.push_back(v);
}
void EfficiencyClass::Reweight(string var, std::map<int,double> weights){
  m_reweightvariables.push_back(new ReweightVar(var,weights));
}
void EfficiencyClass::Reweight(string var, TH1F* scale){
  m_reweightvariables.push_back(new ReweightVar(var,scale));
}
void EfficiencyClass::Reweight(string var1, string var2, TH2F* scale){
  m_reweightvariables.push_back(new ReweightVar(var1,var2, scale));
}
void EfficiencyClass::Reweight(string var1, string var2, TH1F* scale, string form){
  m_reweightvariables.push_back(new ReweightVar(var1,var2, scale, form));
}
void EfficiencyClass::Reweight(string var1, string var2, string var3, string var4, TH2F* scale, string form){
  m_reweightvariables.push_back(new ReweightVar(var1,var2, var3, var4, scale, form));
}
void EfficiencyClass::Reweight(string weightname){
  m_reweightvariables.push_back(new ReweightVar(weightname));
}

void EfficiencyClass::ReweightEff(string var1, string var2, TH2F* scale){
  m_reweighteffvariables.push_back(new ReweightVar(var1,var2, scale));

  //if (m_scaleerrs){
    for (int i = 0 ; i < scale->GetXaxis()->GetNbins() ; ++i ){
      for ( int j = 0 ; j < scale->GetYaxis()->GetNbins() ; ++j ){
	ostringstream a,b;
	a<<"Npass_varyhi_"<<i<<"_"<<j;
	b<<"Npass_varylo_"<<i<<"_"<<j;
	m_Npass_rweff_varyhi->Add(new TParameter<double>(a.str().c_str(),0.0));
	m_Npass_rweff_varylo->Add(new TParameter<double>(b.str().c_str(),0.0));
      }
    }
    //}

  for (std::map<string, EffVar*>::iterator ei = m_variables.begin(); ei != m_variables.end(); ++ei){
    if ((*ei).second->GetEffRWVaryHiPassHists()->GetEntries() == 0 ) (*ei).second->AddEffScaleVaryHists(scale);
  }
}

void EfficiencyClass::SetScaleErr(bool scaleerr){ m_scaleerrs = scaleerr; }
bool EfficiencyClass::GetScaleErr(){ return m_scaleerrs; }


/*
void EfficiencyClass::ReweightVar(string var, std::map<int, double> map){
  m_reweight = true;
  m_reweight_map = true;
  m_reweightmap = map;
  m_reweightvar = var;
  
}

void EfficiencyClass::ReweightVar(string var, TF1* func){
  m_reweight = true;
  m_reweight_func = true;
  m_reweight_tf1 = func;
  m_reweightvar = var;

}

void EfficiencyClass::ReweightVar(string var, TH1F* hist){
  m_reweight = true;
  m_reweight_bin = true;
  m_reweight_hist = hist;
  m_reweightvar = var;

}

void EfficiencyClass::Reweight1_py(string var, PyObject* pyObj){
  
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  TF1*  f    = (TF1*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  if (strcmp(hist->ClassName(),"TH1F") == 0){
    ReweightVar(var,hist);
  }
  else if (strcmp(f->ClassName(),"TF1") == 0){
    ReweightVar(var,f);
  }
  else{
    info()<<hist->ClassName()<<" is not matched to any possibilities"<<endl;
  }
  //std::cout<<hist->ClassName()<<std::endl;
  //Reweight(var,hist);
}
*/


int EfficiencyClass::GetBin(double val, TH1F* hist){
  double lo = hist->GetXaxis()->GetXmin();
  double hi = hist->GetXaxis()->GetXmax();
  int nbins = hist->GetXaxis()->GetNbins();

  int bin = floor(( (val - lo) / (hi - lo) )*nbins) + 1;
  
  return bin;
  
}

EffVar* EfficiencyClass::GetVar(string name){
  return m_variables.at(name);
}
EffVar2D* EfficiencyClass::Get2DVar(string name){
  return m_2Dvariables.at(name);
}

void EfficiencyClass::AddSystematic(double pc){
  for (std::map<string, EffVar*>::iterator ie = m_variables.begin(); ie != m_variables.end(); ++ie){
    //EffVar evar = ie->second;
    ie->second->AddSystematic(pc);
    //cout<<"EfficiencyClass - Systematic Before: "<<evar->m_systematic<<endl;
    //evar->AddSystematic(pc);
    //cout<<"EfficiencyClass - Systematic After: "<<evar->m_systematic<<endl;
  }
  m_toteff.AddSystematic(pc);
}
void EfficiencyClass::AddInvSystematic(double pc){
  for (std::map<string, EffVar*>::iterator ie = m_variables.begin(); ie != m_variables.end(); ++ie){
    //EffVar evar = ie->second;
    //evar->AddInvSystematic(pc);
    ie->second->AddInvSystematic(pc);

  }
  m_toteff.AddInvSystematic(pc);
}

void EfficiencyClass::AddSystematic(string name, double pc){
  if (m_variables.find(name) != m_variables.end()){
    //EffVar evar = m_variables.at(name);
    //evar->AddSystematic(pc);
    //cout<<"EfficiencyClass - Systematic: "<<evar->m_systematic<<endl;
    //evar->m_systematic = true;
    m_variables.at(name)->AddSystematic(pc);
  }
  else info()<<"Could not find variable to add systematic"<<endl;
}
void EfficiencyClass::AddSystematic(string name, std::vector<double> pc){
  if (m_variables.find(name) != m_variables.end()){
    //EffVar evar = m_variables.at(name);
    //evar->AddSystematic(pc);
    //cout<<"EfficiencyClass - Systematic: "<<evar->m_systematic<<endl;
    //evar->m_systematic = true;
    m_variables.at(name)->AddSystematic(pc);
  }
  else info()<<"Could not find variable to add systematic"<<endl;
}

void EfficiencyClass::PrintEfficiencies(string name){
  EffVar* effvar = m_variables.at(name);
  cout<<"-------------------------------------------"<<endl;
  cout<<"Printing "<<m_name<<" Efficiency Values For "<<effvar->GetName()<<endl;
  cout<<"-------------------------------------------"<<endl;
  for (int i = 0; i < effvar->GetBins(); i++){
      double x, y;
      TGraphAsymmErrors* effgraph = effvar->GetEffGraph();
      effgraph->GetPoint(i, x, y);
      double lo = x - effgraph->GetErrorXlow(i);
      double hi = x + effgraph->GetErrorXhigh(i);
      double erroryhi = effgraph->GetErrorYhigh(i);
      double errorylo = effgraph->GetErrorYlow(i);
      
      cout<<"Value: "<<lo<<" - "<<hi<<" Efficiency: "<<y<<" + "<<erroryhi<<" - "<<errorylo<<endl;
      
    }
  cout<<"-------------------------------------------"<<endl;
}

void EfficiencyClass::PrintTwikiEfficiencies(string name){
  EffVar* effvar = m_variables.at(name);
  cout<<"-------------------------------------------"<<endl;
  cout<<"Printing "<<m_name<<" Twiki Efficiency Values For "<<effvar->GetName()<<endl;
  cout<<"-------------------------------------------"<<endl;
  cout<<" | *"<<name<<"* | *Efficiency* | *Uncertainty* |"<<endl;
  for (int i = 0; i < effvar->GetBins(); i++){
      double x, y;
      TGraphAsymmErrors* effgraph = effvar->GetEffGraph();
      Int_t sc = effgraph->GetPoint(i, x, y);
      if (sc == -1) break;
      double lo = x - effgraph->GetErrorXlow(i);
      double hi = x + effgraph->GetErrorXhigh(i);
      double erroryhi = effgraph->GetErrorYhigh(i);
      double errorylo = effgraph->GetErrorYlow(i);
      cout.precision(3);
      

      cout<<" | "<<lo<<" - "<<hi<<" | "<<y<<" | "<<max(erroryhi, errorylo)<<" | "<<endl;
      
    }
  cout<<"-------------------------------------------"<<endl;
}

void EfficiencyClass::PrintNTwikiEfficiencies(string name, std::vector<std::pair< string, EfficiencyClass> > classes){

  EffVar* effvar = classes.begin()->second.GetVariables().at(name);

  cout<<"-------------------------------------------"<<endl;
  cout<<"Printing Twiki Efficiency Values For "<<effvar->GetName()<<endl;
  cout<<"-------------------------------------------"<<endl;
  for (int i = -1; i < effvar->GetBins(); i++){
    for (std::vector<std::pair<string,EfficiencyClass> >::iterator ie = classes.begin() ; ie != classes.end() ; ++ie ){
      string label             = (*ie).first;
      EfficiencyClass effclass = (*ie).second;
      EffVar* effvar            = effclass.GetVariables().at(name);
      if (ie == classes.begin()){
	if (i==-1){
	  cout<<"| *"<<effvar->GetName()<<" * ";
	}
	else{
	  double x, y;
	  TGraphAsymmErrors* effgraph = effvar->GetEffGraph();
	  effgraph->GetPoint(i, x, y);
	  double lo = x - effgraph->GetErrorXlow(i);
	  double hi = x + effgraph->GetErrorXhigh(i);
	  //cout.precision(2);
	  cout.unsetf(std::ios::fixed);
	  cout<<" | "<<lo<<" - "<<hi<<" | ";
	}
      }

      if (i==-1)     {
	cout<<" | *"<<label<<"* | *Uncertainty* |";
      }
      else{
	double x, y;
	TGraphAsymmErrors* effgraph = effvar->GetEffGraph();
	effgraph->GetPoint(i, x, y);
	//double lo = x - effgraph->GetErrorXlow(i);
	//double hi = x + effgraph->GetErrorXhigh(i);
	double erroryhi = effgraph->GetErrorYhigh(i);
	double errorylo = effgraph->GetErrorYlow(i);
	cout.precision(3);
	cout<<fixed<<y<<" | "<<max(erroryhi, errorylo)<<" | ";
      
      }
    }
  cout<<" "<<endl;
  }
  cout<<"-------------------------------------------"<<endl;
}


Double_t fitGaussLine(Double_t* x, Double_t *par){
  Double_t PDF = 0.0;
  Double_t N   = par[0];
  Double_t m   = par[1];
  Double_t s   = par[2];
  Double_t c   = par[3];
  Double_t sl  = par[4];

  PDF = exp((-pow(((*x)-m),2))/(2*pow(s,2)));

  PDF = PDF * N;
  PDF = PDF + sl * (*x) + c;

  return PDF;

}

Double_t fitCBGauss(Double_t* x, Double_t* par){
  Double_t PDF = 0.0;
  Double_t N   = par[0];
  Double_t a   = par[1];
  Double_t n   = par[2];
  Double_t m   = par[3];
  Double_t s   = par[4];
  Double_t c   = par[5];
  Double_t sl  = par[6];

  Double_t N2 = par[7];
  Double_t m2 = par[8];
  Double_t s2 = par[9];


if (((*x) - m)/s > -a){
    PDF = exp((-pow(((*x)-m),2))/(2*pow(s,2)));
  }
  else{
    Double_t A = pow(n/fabs(a),n) * exp(-pow(a,2)/2);
    Double_t B = (n/fabs(a)) - fabs(a);
    PDF = A * pow(B-(((*x)-m)/s),-n);
  }
  PDF = PDF * N;
  PDF = PDF + sl * (*x) + c;
  PDF += N2*exp((-pow(((*x)-m2),2))/(2*pow(s2,2)));

  return PDF;
  

}

Double_t fitCB(Double_t* x, Double_t *par){
  Double_t PDF = 0.0;
  Double_t N   = par[0];
  Double_t a   = par[1];
  Double_t n   = par[2];
  Double_t m   = par[3];
  Double_t s   = par[4];
  Double_t c   = par[5];
  Double_t sl  = par[6];


  if (((*x) - m)/s > -a){
    PDF = exp((-pow(((*x)-m),2))/(2*pow(s,2)));
  }
  else{
    Double_t A = pow(n/fabs(a),n) * exp(-pow(a,2)/2);
    Double_t B = (n/fabs(a)) - fabs(a);
    PDF = A * pow(B-(((*x)-m)/s),-n);
  }
  PDF = PDF * N;
  PDF = PDF + sl * (*x) + c;
  
  return PDF;
}

TGraphAsymmErrors* EfficiencyClass::DivideTGraphs(TGraphAsymmErrors* numer, TGraphAsymmErrors* denom){
  int nentries = numer->GetN();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nentries);
  for (int i = 0 ; i < nentries ; ++i ){
    double x = 0.0;
    double y = 1.0;
    double yerrhi, yerrlo;
    double xerrhi, xerrlo;
    
    double xnum = 0, ynum = 0, xdenom = 0, ydenom = 0;

    Int_t s1 = numer->GetPoint(i, xnum, ynum);
    Int_t s2 = denom->GetPoint(i, xdenom, ydenom);
    if (s1 != -1 && s2 != -1){
      double xnum_errhi  = numer->GetErrorXhigh(i);
      double xnum_errlo  = numer->GetErrorXlow(i);
      //double xdenom_errhi  = denom->GetErrorXhigh(i);
      //double xdenom_errlo  = denom->GetErrorXlow(i);
      
      double ynum_errhi  = numer->GetErrorYhigh(i);
      double ynum_errlo  = numer->GetErrorYlow(i);
      double ydenom_errhi  = denom->GetErrorYhigh(i);
      double ydenom_errlo  = denom->GetErrorYlow(i);
      
      if (ydenom > 0){
	y = ynum/ydenom;
      }
      else y = 0;
      x = xnum;
      
      yerrhi = y * sqrt(pow(ynum_errhi/ynum,2) + pow(ydenom_errhi/ydenom,2)  );
      yerrlo = y * sqrt(pow(ynum_errlo/ynum,2) + pow(ydenom_errlo/ydenom,2)  );
      
      xerrlo = xnum_errlo;
      xerrhi = xnum_errhi;

      graph->SetPoint(i, x, y);
      graph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
    }
    
  }
  return graph;
}


void EfficiencyClass::RemoveErrors(){
  for (std::map<string, EffVar*>::iterator ie = m_variables.begin(); ie != m_variables.end(); ++ie){
    (*ie).second->RemoveErrors();
  }
}

TH2F* EfficiencyClass::DivideTHists(TH2F* numer, TH2F* denom){
  verbose()<<"Gonna divide "<<numer<<" by "<<denom<<endl;
  TH2F* res = (TH2F*)numer->Clone();
  verbose()<<"Cloned"<<endl;
  res->Divide(denom);
  verbose()<<"Divided"<<endl;

  return res;
}

TH2F* EfficiencyClass::CombineTHists(std::vector<TH2F*> hists){
  TH2F* res = 0;
  if (hists.size() == 0 ) return res;
  res = hists.at(0);
  if (hists.size() == 1 ) return res;
  for (unsigned int i = 1; i < hists.size() ; ++i){
    res->Multiply(hists.at(i));
  }
  return res;

}

double EfficiencyClass::GetCorrectedEfficiency(string var, TTree* t, string leaf){
  TGraphAsymmErrors* graph = m_variables.at(var)->GetEffGraph();
  TH1F* hist = m_variables.at(var)->GetTotHist();
  t->SetBranchStatus("*",0);
  t->SetBranchStatus(leaf.c_str(),1);
  double val;
  int evts = 0;
  double corrEvts = 0;
  t->SetBranchAddress(leaf.c_str(),&val);
  double x, y;
  for (int i = 0; i < t->GetEntries(); ++i){
    t->GetEntry(i);
    int bin = hist->FindBin(val) - 1 ;
    graph->GetPoint(bin, x, y);
    double evteff = y;
    corrEvts += (1/evteff);
    evts ++;
  }

  double eff = (double)evts/corrEvts;
  return eff;
}

std::vector<double> EfficiencyClass::GetCorrectedEfficiency(string var, std::vector<TH1F*> hists, bool smear){
  EffVar* v = m_variables.at(var);
  TGraphAsymmErrors* graph = v->GetEffGraph();
  if (smear){
    graph = m_variables.at(var)->GetSmearedEffGraph();
  }

  //TH1F* hist = m_variables.at(var).GetTotHist();
  double x, y;
  std::vector<double> effs;

  for (std::vector<TH1F*>::iterator it = hists.begin(); it !=hists.end(); ++it){
    int evts = 0;
    double corrEvts = 0;
    TH1F* h = (*it);
    for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i){
      double lo = h->GetXaxis()->GetBinLowEdge(i+1);
      double hi = h->GetXaxis()->GetBinLowEdge(i+2);
      //cout<<"lo: "<<lo<<" hi: "<<hi<<endl;
      //int bin1 = hist->FindBin(lo);
      //int bin2 = hist->FindBin(hi);
      int bin1 = v->FindBin(lo);
      int bin2 = v->FindBin(hi);
      if (bin1 != bin2 && hi != v->GetEdges().at(bin2 - 1 )){
	verbose()<<"lo: "<<lo<<" hi: "<<hi<<" bin1: "<<bin1<<" bin2: "<<bin2<<endl;
	verbose()<<"Precision will be lost in efficiency determination"<<endl;;
	return effs;
      }
      int p = bin1+1;
      graph->GetPoint(p, x, y);
      double evteff = y;
      double N = h->GetBinContent(bin1);
      corrEvts += (N/evteff);
      evts += N;
    }
    effs.push_back((double)evts/corrEvts);
  }

  return effs;
}


double EfficiencyClass::GetCorrectedEfficiency(string var, TH1F* h){
  TGraphAsymmErrors* graph = m_variables.at(var)->GetEffGraph();
  TH1F* hist = m_variables.at(var)->GetTotHist();
  int evts = 0;
  double corrEvts = 0;
  double x, y;
  for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i){
    double lo = h->GetXaxis()->GetBinLowEdge(i+1);
    double hi = h->GetXaxis()->GetBinLowEdge(i+2);
    int bin1 = hist->FindBin(lo);
    int bin2 = hist->FindBin(hi);
    if (bin1 != bin2 && hi != hist->GetXaxis()->GetBinLowEdge(bin2)){
      verbose()<<"Precision will be lost in efficiency determination";
      return -1;
    }
    int p = bin1+1;
    graph->GetPoint(p, x, y);
    double evteff = y;
    double N = h->GetBinContent(bin1);
    corrEvts += (N/evteff);
    evts += N;
  }

  double eff = (double)evts/corrEvts;

  return eff;
}


TGraphAsymmErrors* EfficiencyClass::CombineTGraphs(std::vector<TGraphAsymmErrors*> graphs){
  int nentries = graphs.at(0)->GetN();

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nentries);

  for (int i = 0 ; i < nentries ; ++i ){
    double x = 0.0;
    double y = 1.0;
    double yerrhi = 0, yerrlo = 0;
    double xerrhi = 0, xerrlo = 0;
    
    double sumsqerryhi = 0.;
    double sumsqerrylo = 0.;
    
    for (std::vector<TGraphAsymmErrors*>::iterator ig = graphs.begin() ; ig != graphs.end(); ++ig){
      
      TGraphAsymmErrors* igraph = (*ig);
      
      double y2;
      double yerrhi2, yerrlo2;
      
      Int_t s = igraph->GetPoint(i, x, y2);
      if (s != -1){
	xerrhi  = igraph->GetErrorXhigh(i);
	xerrlo  = igraph->GetErrorXlow(i);
	yerrhi2 = igraph->GetErrorYhigh(i);
	yerrlo2 = igraph->GetErrorYlow(i);
	
	y = y * y2;
	sumsqerryhi += pow(yerrhi2/y2,2);
	sumsqerrylo += pow(yerrlo2/y2,2);
      }
    }
    
    yerrhi = y*sqrt(sumsqerryhi);
    yerrlo = y*sqrt(sumsqerrylo);
    
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, xerrlo, xerrhi, yerrlo, yerrhi);
    
  }
  return graph;
}

void EfficiencyClass::Run(){
  MakeEntryLists();
  LoopEntries();
  MakeEfficiencyGraph();
}
#ifdef WITHPYTHON
void EfficiencyClass::Reweight1_py(string var, PyObject* pyObj){
  
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  TF1*  f    = (TF1*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  if (strcmp(hist->ClassName(),"TH1F") == 0){
    Reweight(var,hist);
  }
  else if (strcmp(f->ClassName(),"TF1") == 0){
    Reweight(var,f);
  }
  else{
    info()<<hist->ClassName()<<" is not matched to any possibilities"<<endl;
  }
}
void EfficiencyClass::Reweight2_py(string var1, string var2, PyObject* th2f){
  TH2F* hist_2d  = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  TH1F* hist_1d  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  if (strcmp(hist_1d->ClassName(),"TH1F") == 0){
    Reweight(var1, var2, hist_1d);
  }
  else{
    Reweight(var1, var2, hist_2d);
  }
}

void EfficiencyClass::Reweight3_py(string leaf){
  Reweight(leaf);
}
void EfficiencyClass::Reweight4_py(string var1, string var2, PyObject* th1f, string form){
  TH1F* hist  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(th1f));
  if (strcmp(hist->ClassName(),"TH1F") == 0){
    Reweight(var1, var2, hist, form);
  }
  else{
    info()<<"Cannot reweight for "<<var1<<" "<<var2<<endl;
  }
}
void EfficiencyClass::Reweight5_py(string var1, string var2, string var3, string var4, PyObject* th2f, string form){
  TH2F* hist  = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  if (strcmp(hist->ClassName(),"TH2F") == 0){
    Reweight(var1, var2, var3, var4, hist, form);
  }
  else{
    info()<<"Cannot reweight for "<<var1<<" "<<var2<<" "<<var3<<" "<<var4<<endl;
  }
}


void EfficiencyClass::ReweightEff_py(string var1, string var2, PyObject* th2f){
  TH2F* hist_2d  = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  //TH1F* hist_1d  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  //if (strcmp(hist_1d->ClassName(),"TH1F") == 0){
  //  ReweightEff(var1, var2, hist_1d);
  //}
  //else{
  ReweightEff(var1, var2, hist_2d);
  //}
}

double EfficiencyClass::GetTotEff_py(){
  double eff = 0.0;
  eff = m_toteff.GetEff();
  return eff;
}
PyObject* EfficiencyClass::GetTotHist_py(){
  //PyObject* pyObj = Utils::Root2PyObj<TH1F>(m_tot);
  //return pyObj;
  return TPython::ObjectProxy_FromVoidPtr(m_tot, m_tot->ClassName());
}

PyObject* EfficiencyClass::GetPassHist_py(){
  //PyObject* pyObj = Utils::Root2PyObj<TH1F>(m_pass);
  //return pyObj;
  return TPython::ObjectProxy_FromVoidPtr(m_pass, m_pass->ClassName());
}
//For python boost - deal with overloaded functions
void EfficiencyClass::AddVar1_py(string name, string var, int bins, double lo, double hi){
  AddVar(name, var, bins, lo, hi);
}
void EfficiencyClass::AddVar2_py(string name, string var, vector<double> edges){
  AddVar(name, var, edges);
}
void EfficiencyClass::AddVar3_py(string name, string var, int bins, float lo, float hi){
  AddVar(name, var, bins, lo, hi);
}
void EfficiencyClass::AddVar4_py(string name, string var, boost::python::list& ns){
  std::vector<double> edges;
  for (int i = 0; i < len(ns); ++i){
    double val = boost::python::extract<double>(ns[i]);
    edges.push_back(val);
  }
  AddVar(name, var, edges);
}
void EfficiencyClass::AddVar5_py(boost::python::list& varlist){
  if (len(varlist) == 6){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    int bins    = boost::python::extract<int>(varlist[2]);
    double lo   = boost::python::extract<double>(varlist[3]);
    double hi   = boost::python::extract<double>(varlist[4]);
    string type = boost::python::extract<string>(varlist[5]);
    AddVar(name, var, bins, lo, hi);
  }
  if (len(varlist) == 5){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    int bins    = boost::python::extract<int>(varlist[2]);
    double lo   = boost::python::extract<double>(varlist[3]);
    double hi   = boost::python::extract<double>(varlist[4]);
    AddVar(name, var, bins, lo, hi);
  }
  else if (len(varlist) == 4){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    boost::python::list ns = (boost::python::list)(varlist[2]);
    string type  = boost::python::extract<string>(varlist[3]);
    AddVar4_py(name, var, ns);
  }
  else if (len(varlist) == 3){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    boost::python::list ns = (boost::python::list)(varlist[2]);
    AddVar4_py(name, var, ns);
  }
}

void EfficiencyClass::AddVars_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    AddVar5_py(var);
  }
}

void EfficiencyClass::AddTrees_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::object obj = ns[i];
    PyObject* pyObj = obj.ptr();
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    AddTree(t);
  }
}

double EfficiencyClass::GetCorrectedEfficiency1_py(string var, PyObject* h){
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(h));
  return GetCorrectedEfficiency(var,hist);

}
double EfficiencyClass::GetCorrectedEfficiency2_py(string var, PyObject* t, string leaf){
  TTree* tree = (TTree*)(TPython::ObjectProxy_AsVoidPtr(t));
  return GetCorrectedEfficiency(var , tree , leaf);
  
}
void EfficiencyClass::AddTree_py(PyObject* pyObj){
  TTree* tree = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  AddTree(tree);
}
void EfficiencyClass::Add2DVar_py(string var1, string var2){
  Add2DVar(var1, var2);
}

void EfficiencyClass::Add2DVars_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    if (len(var) == 2){
      string var1 = boost::python::extract<string>(var[0]);
      string var2 = boost::python::extract<string>(var[1]);
      Add2DVar(var1, var2);
    }
    else if (len(var) == 3){
      string name = boost::python::extract<string>(var[0]);
      string var1 = boost::python::extract<string>(var[1]);
      string var2 = boost::python::extract<string>(var[2]);
      Add2DVar(var1, var2, name);
    }
    else{
      info()<<"Cannot Add 2D Variable - Wrong Dimensions"<<endl;
    }
  }
}
void EfficiencyClass::SetPassCut_py(PyObject* pyObj){
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  m_passcut = *cut;
}
void EfficiencyClass::SaveToFile_py(){
  SaveToFile();
}
void EfficiencyClass::LoadFromFile_py(){
  LoadFromFile();
}
void EfficiencyClass::SetSelectionCut_py(PyObject* pyObj){
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  m_selcut = *cut;
}
void EfficiencyClass::MakeEfficiencyGraph_py(){
  MakeEfficiencyGraph();
}

void EfficiencyClass::AddSystematic1_py(string name, double pc){
  AddSystematic(name, pc);
}
void EfficiencyClass::AddSystematic2_py(double pc){
  AddSystematic(pc);
}
void EfficiencyClass::AddSystematic3_py(string name, boost::python::list& ns){
  std::vector<double> pc;
  for (int i = 0; i < len(ns); ++i){
    double val = boost::python::extract<double>(ns[i]);
    pc.push_back(val);
  }
  AddSystematic(name, pc);
}
boost::python::list EfficiencyClass::GetCorrectedEfficiency3_py(string var, boost::python::list& h, bool smear){
  std::vector<TH1F*> hists;
  for (int i = 0; i < len(h); ++i){
    boost::python::object obj = h[i];
    PyObject* pyObj = obj.ptr();
    TH1F* h1f = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    hists.push_back(h1f);
  }
  std::vector<double> effs = GetCorrectedEfficiency(var, hists, smear);
  boost::python::object get_iter = boost::python::iterator<std::vector<double> >();
  boost::python::object iter = get_iter(effs);
  boost::python::list l(iter);
  return l;
}

#endif
