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
#include <AnalysisClass.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

AnalysisClass::AnalysisClass(string name){
  m_name = name;
  m_fit = 0;
  m_fitter = 0;
  m_fitchi2 = 0;
  m_ndof = 0;
  m_fitratio = false;
}

void AnalysisClass::AddTemplate(string name, TTree* t){
  m_stackorder.push_back(name);
  m_fitorder.push_back(name);
  Template* temp = new Template(name, t, m_selcut);
  m_templates.insert(std::pair<string,Template*>(temp->GetName(),temp));
}

void AnalysisClass::AddTemplate(string name, TTree* t, enum EColor color){
  m_stackorder.push_back(name);
  m_fitorder.push_back(name);
  Template* temp = new Template(name, t, m_selcut, color);
  m_templates.insert(std::pair<string,Template*>(temp->GetName(),temp));
}

void AnalysisClass::AddTemplate(string name, TTree* t, TCut cut){
  m_stackorder.push_back(name);
  m_fitorder.push_back(name);
  Template* temp = new Template(name, t, cut);
  m_templates.insert(std::pair<string,Template*>(temp->GetName(),temp));
  
}

void AnalysisClass::AddTemplate(string name, TTree* t, TCut cut, enum EColor color){
  m_stackorder.push_back(name);
  m_fitorder.push_back(name);
  Template* temp = new Template(name, t, cut, color);
  m_templates.insert(std::pair<string,Template*>(temp->GetName(),temp));
}


void AnalysisClass::AddData(string name, TTree* t){
  m_data = new Template(name,t,m_selcut);
}


void AnalysisClass::AddData(string name, TTree* t, TCut cut){
  m_data = new Template(name,t,cut);
}

void AnalysisClass::AddData(Template* temp){
  m_data = temp;
}

void AnalysisClass::AddTemplate(Template* temp){
  m_stackorder.push_back(temp->GetName());
  m_fitorder.push_back(temp->GetName());
  m_templates.insert(std::pair<string,Template*>(temp->GetName(),temp));
}

void AnalysisClass::NormaliseTemplates(double n){
  if (m_data) m_data->NormaliseToEvts(n);
  
  for (std::map<string, Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    (*it).second->NormaliseToEvts(n);
  }
}



/*
TFractionFitter* AnalysisClass::TFracFit(string var1, string var2){
  TObjArray* toFit = new TObjArray();
  string var = var1+"_"+var2;

  TH2F* data = m_data->Get2DVar(var)->GetHist();
  
  std::map<int,Template*> idx;

  int i = 0;
  
  for (std::map<string, Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    TH2F* hist = (*it).second->Get2DVar(var)->GetHist();
    toFit->Add(hist);
    idx.insert(std::pair<int,Template*>(i,(*it).second));
    i++;
    }  
  TFractionFitter* fit = new TFractionFitter(data, toFit);
  //TFractionFitter has a habit of hitting infinite loops
  //Exclude bins with < 3 events to see if this changes anything


  //Perform constraints
  for (std::map<int, Template*>::iterator it = idx.begin() ; it != idx.end() ; ++it ){
    Template* temp = (*it).second;
    if (temp->IsFixed()) {
      double constraint = temp->GetNormEvts()/m_data->GetEvents();
      std::cout<<"Template "<<temp->GetName()<<" is fixed to "<<constraint<<" Events"<<std::endl;
      fit->GetFitter()->SetParameter((*it).first, temp->GetName().c_str(),constraint, 0.0, 0.0, 0.0);
      fit->GetFitter()->FixParameter((*it).first);
    }
   else{
     fit->Constrain((*it).first + 1 , 0, 1);
     //fit->GetFitter()->SetParameter((*it).first, temp->GetName().c_str(), 0.5, 0.5, 0 , 1);
    }

  }
  
  Int_t status = fit->Fit();               // perform the fit                                                                                                
  
  if (status == 0) {
    for(int j=0;j<toFit->GetEntries();j++){
      double value, error;
      fit->GetResult(j,value,error);
      double sf = data->Integral()*value/idx.at(j)->Get2DVar(var)->GetHist()->Integral();
      //idx.at(j)->SetNormEvts(idx.at(j)->GetNormEvts() * (sf));
      idx.at(j)->Scale(sf, false);
      idx.at(j)->SetFitFrac(value, error);
      m_fit = (TH1F*)fit->GetPlot();
    }
  }
  //toFit->Delete();
  m_fitchi2 = fit->GetChisquare();
  m_ndof = fit->GetNDF();
  //m_fitter = fit->GetFitter();

  return fit;
}
*/

void AnalysisClass::ConstrainRatio(string tempA, string tempB, double r){
  pair<string, string> temps(tempA, tempB);
  m_conratio = pair< pair<string, string>, double>(temps, r);
  m_fitratio = true;
}

void AnalysisClass::AddFitter(string var, double lo, double hi, bool combine){
  //Set up fitter as a tfractionfitter

  if (!m_data) {
    cout<<"No Data Set! - Not adding Fitter"<<endl;
    return;
  }
  
  
  TObjArray* toFit = new TObjArray();
  TH1F* data;

  if (combine) data = m_data->Get2DVar(var)->GetCombHist();
  else data = m_data->GetVar(var)->GetHist(); 

  vector<string> names;
  vector< string > constrain_names;
  vector< double > constrain_vals;

  std::map<int,Template*> idx; // List of names and their index in fit

  for (std::vector<string>::iterator is = m_fitorder.begin(); is != m_fitorder.end(); ++is){
    Template* temp = m_templates.at(*is);
    TH1F* hist;
    if (combine) hist = temp->Get2DVar(var)->GetCombHist();
    else hist = temp->GetVar(var)->GetHist();
    toFit->Add(hist);
    names.push_back((*is));
    if (temp->IsFixed()) {
      double constraint = 0.0;
      //calculate constraint as a fraction of the total data
      if (combine) constraint = temp->Get2DVar(var)->GetCombHist()->Integral()/data->Integral();
      else constraint = temp->GetVar(var)->GetHist()->Integral()/data->Integral();
      constrain_names.push_back(temp->GetName());
      constrain_vals.push_back(constraint);
      //cout<<"Constraining :"<<temp->GetName()<<" "<<constraint<<endl;
    }
  }



  m_fitter = new Fitter(data, toFit, names, var);
  m_fitter->AddConstraints(0.0,1.0);
  m_fitter->AddConstraints(constrain_names, constrain_vals);
  m_fitter->SetFitRange(lo, hi);
}


void AnalysisClass::Add2DFitter(string var){
  //Set up fitter as a tfractionfitter
  if (!m_data){
    cout<<"No data set - not additing fitter"<<endl;
    return;
  }
  
  cout<<"getting 2d fitter"<<endl;
  TObjArray* toFit = new TObjArray();
  TH2F* data;

  data = m_data->Get2DVar(var)->GetHist(); 
  
  vector<string> names;
  vector< string > constrain_names;
  vector< double > constrain_vals;

  std::map<int,Template*> idx; // List of names and their index in fit

  for (std::vector<string>::iterator is = m_fitorder.begin(); is != m_fitorder.end(); ++is){
    Template* temp = m_templates.at(*is);
    TH2F* hist;
    hist = temp->Get2DVar(var)->GetHist();
    toFit->Add(hist);
    names.push_back((*is));
    if (temp->IsFixed()) {
      double constraint = temp->Get2DVar(var)->GetHist()->Integral()/data->Integral();
      constrain_names.push_back(temp->GetName());
      constrain_vals.push_back(constraint);
      //cout<<"Constraining :"<<temp->GetName()<<" "<<constraint<<endl;
    }
  }

  m_fitter = new Fitter(data, toFit, names, var);
  m_fitter->AddConstraints(0.0,1.0);
  m_fitter->AddConstraints(constrain_names, constrain_vals);
  cout<<"Got 2d fitter"<<endl;
}

void AnalysisClass::ClearFitter(){
  delete m_fitter;
  m_fitter = NULL;
}


void AnalysisClass::UnscaleTemplates(){
  for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    Template* temp = (*it).second;
    temp->Unscale();
  }
}


void AnalysisClass::Apply2DFitResults(){
  //Get the results from the fitter and scale templates appropriately
  if ((!m_data) || (!m_fitter)) {
    cout<<"No fit results to apply - not applying"<<endl;
    return;

  }
  string var = m_fitter ? m_fitter->GetVar() : "";
  
  map<string, pair<double, double> > results = m_fitter->GetResults();
  for(map<string, pair<double, double> >::iterator ij= results.begin(); ij != results.end(); ++ij){
    double value = (*ij).second.first;
    double error = (*ij).second.second;
    string name = (*ij).first;
    
    if (m_templates.find(name) != m_templates.end()){
      Template* temp = m_templates.at(name);
      temp->SetFitFrac(value, error);
      double sf = 1.0;
      sf = m_data->Get2DVar(var)->GetHist()->Integral() * value/temp->Get2DVar(var)->GetHist()->Integral();
      temp->Scale(sf, temp->IsFixed());
    }
  }

  //Save some of the key variables for later - fit and chi2/ndof
  //m_fit = (TH1F*)m_fitter->GetFitter()->GetPlot();
  m_fitchi2 = m_fitter->GetFitter()->GetChisquare();
  m_ndof = m_fitter->GetFitter()->GetNDF();


}
void AnalysisClass::ApplyFitResults(bool combine){
  //Get the results from the fitter and scale templates appropriately
  
  if ((!m_data) || (!m_fitter)) {
    cout<<"No fit results to apply - not applying"<<endl;
    return;

  }
  string var = m_fitter ? m_fitter->GetVar() : "";
  map<string, pair<double, double> > results = m_fitter->GetResults();
  for(map<string, pair<double, double> >::iterator ij= results.begin(); ij != results.end(); ++ij){
    double value = (*ij).second.first;
    double error = (*ij).second.second;
    string name = (*ij).first;
    
    if (m_templates.find(name) != m_templates.end()){
      Template* temp = m_templates.at(name);
      temp->SetFitFrac(value, error);
      double sf = 1.0;
      if (combine) sf = m_data->Get2DVar(var)->GetCombHist()->Integral() * value/temp->Get2DVar(var)->GetCombHist()->Integral();
      else sf = m_data->GetVar(var)->GetHist()->Integral() * value/temp->GetVar(var)->GetHist()->Integral();
      temp->Scale(sf, temp->IsFixed());
    }
  }

  //Save some of the key variables for later - fit and chi2/ndof
  m_fit = (TH1F*)m_fitter->GetFitter()->GetPlot();
  m_fitchi2 = m_fitter->GetFitter()->GetChisquare();
  m_ndof = m_fitter->GetFitter()->GetNDF();


}


TFractionFitter* AnalysisClass::TFracFit(string var1, string var2){
  //Helper function that performs the key TFracFit stages
  string var = var1+"_"+var2;
  cout<<"adding fitter"<<endl;
  if (!m_fitter) Add2DFitter(var);
  cout<<"2d fitted added "<<m_fitter<<endl;
  m_fitter->TFracFit();
  cout<<"fit performed"<<endl;
  Apply2DFitResults();
  cout<<"results applied"<<endl;
  return m_fitter->GetFitter();
}

TFractionFitter* AnalysisClass::TFracFit(string var, double lo, double hi, bool combine){
  //Helper function that performs the key TFracFit stages
  if (!m_fitter) AddFitter(var, lo, hi, combine);
  m_fitter->TFracFit();
  ApplyFitResults(combine);
  return m_fitter->GetFitter();

}
TFractionFitter* AnalysisClass::RedoFit(string var, double lo, double hi, bool combine){
  //Helper function that performs the key TFracFit stages
  delete m_fitter;
  AddFitter(var, lo, hi, combine);
  m_fitter->TFracFit();
  ApplyFitResults(combine);
  return m_fitter->GetFitter();
}


TFractionFitter* AnalysisClass::CombTFracFit(string var1, string var2){
  string var = var1+"_"+var2;
  return TFracFit(var, 0.0, 0.0, true);
}


Fitter* AnalysisClass::GetFitter(){
  return m_fitter;
}

Template* AnalysisClass::GetTemplate(string name){
  if (m_templates.find(name) != m_templates.end()){
    return m_templates.at(name);
  }
  else if (name == "Data") return m_data;
  cerr<<"No template found with name: "<<name<<endl;
  return 0;
}

Template* AnalysisClass::GetData(){
  return m_data;
}

void AnalysisClass::ApplyCuts(){
  if (m_data){
    cout<<"Applying cuts to "<<m_data->GetName()<<endl;
    m_data->ApplyCut();
  }
  else{
    cout<<"Note - no data to run over"<<endl;
  }
  for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    cout<<"Applying cuts to "<<(*it).first<<endl;
    Template* temp = (*it).second;
    temp->ApplyCut();
    //cout<<"Cut Applied"<<endl;
  }
}

void AnalysisClass::FillVars(){
  if (m_data) m_data->FillVars();
  for (std::map<string, Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    cout<<"Filling Variables for "<<(*it).first<<endl;
    Template* temp = (*it).second;
    temp->FillVars();
  }
}



void AnalysisClass::SaveToFile(string output){
  //cout<<"Saving to file for "<< (output == "") ? m_name : output<<endl;
  //const char* cwd = gROOT->pwd();
  string outputFile = output == "" ? m_name+".root" : output+".root";
  TFile* f = new TFile(outputFile.c_str(),"RECREATE");
  //gROOT->cd(cwd);
  if (m_data){
    f->mkdir(m_data->GetName().c_str());
    f->cd(m_data->GetName().c_str());
    m_data->SaveToCurrentFile();
  }
  f->cd();
  for (std::map<string, Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    Template* temp = (*it).second;
    f->mkdir(temp->GetName().c_str());
    f->cd(temp->GetName().c_str());
    temp->SaveToCurrentFile();
    f->cd();
    
  }
  for (std::map<string, THStack*>::iterator is = m_stacks.begin() ; is != m_stacks.end() ; ++is ){
    (*is).second->Write((*is).first.c_str());
  }
  if (m_fit) m_fit->Write("Fit");
  //if (m_fitter) m_fitter->GetFitter()->Write("TFracFit");
  TParameter<double>* fitchi2 = new TParameter<double>("FitChi2", m_fitchi2);
  TParameter<double>* ndof    = new TParameter<double>("nDoF"   , m_ndof);
  fitchi2->Write("FitChi2");
  ndof->Write("nDoF");
  //f->Write();
  gROOT->cd();
  f->Close();

  //f->Delete();
  cout<<"File saved"<<endl;
}

double AnalysisClass::GetChi2nDoF(){
  return double(m_fitchi2/m_ndof);
}



void AnalysisClass::AddVar(string name, string var, int bins, double lo, double hi){
  m_vars.push_back(name);
  if (m_templates.size() == 0) cout<<"No templates to add variables to"<<endl;

  //if (var.find("-") == std::string::npos && var.find("+") == std::string::npos)
  //  {
  if (m_data) m_data->AddVar(name , var , bins , lo , hi , m_name+"_" + m_data->GetName());
      for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
	Template* temp = (*it).second;
	temp->AddVar( name, var, bins, lo, hi, m_name+"_"+temp->GetName());
      }
      //  }
      //else {
      //m_data->AddAlgVar(name , var , bins , lo , hi , m_name+"_" + m_data->GetName());
      //for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
      //Template* temp = (*it).second;
      //temp->AddAlgVar( name, var, bins, lo, hi, m_name+"_"+temp->GetName());
      //}
      //}
}
void AnalysisClass::AddVar(string name, string var, std::vector<double> edges){
  m_vars.push_back(name);
  if (m_templates.size() == 0) cout<<"No templates to add variables to"<<endl;
  //if (var.find("-") == std::string::npos && var.find("+") == std::string::npos)
  //{
  if (m_data) m_data->AddVar(name , var , edges, m_name+"_"+m_data->GetName() );
       for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
	 Template* temp = (*it).second;
	temp->AddVar( name, var, edges, m_name+"_"+temp->GetName());
      }
       //}
       //else{
       //m_data->AddAlgVar(name , var , edges, m_name+"_"+m_data->GetName() );
       //for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
       //Template* temp = (*it).second;
       //temp->AddAlgVar( name, var, edges, m_name+"_"+temp->GetName());
       //}
       //}
}
/*
void AnalysisClass::AddAlgVar(string name, string varexp, int bins, double lo, double hi){
  m_vars.push_back(name);
  if (m_templates.size() == 0) cout<<"No templates to add variables to"<<endl;
  m_data->AddAlgVar(name , varexp , bins , lo , hi , m_name+"_" + m_data->GetName());
  for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    Template* temp = (*it).second;
    temp->AddAlgVar( name, varexp, bins, lo, hi, m_name+"_"+temp->GetName());
  }
}
void AnalysisClass::AddAlgVar(string name, string varexp, std::vector<double> edges){
  m_vars.push_back(name);
  if (m_templates.size() == 0) cout<<"No templates to add variables to"<<endl;
  m_data->AddAlgVar(name , varexp , edges, m_name+"_"+m_data->GetName() );
  for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    Template* temp = (*it).second;
    temp->AddAlgVar( name, varexp, edges, m_name+"_"+temp->GetName());
  }
}
*/
void AnalysisClass::Add2DVar(string var1, string var2){
  string name = var1+"_"+var2;
  m_2Dvars.push_back(name);
  if (m_templates.size() == 0) cout<<"No templates to add variables to"<<endl;
  if (m_data) m_data->Add2DVar(name , var1, var2 , m_name+"_"+m_data->GetName() );
  for (std::map<string,Template*>::iterator it = m_templates.begin() ; it != m_templates.end() ; ++it ){
    Template* temp = (*it).second;
    temp->Add2DVar( name, var1, var2, m_name+"_"+temp->GetName());
  }
}


void AnalysisClass::StyleTemplates(){
 for(std::map<std::string,Template*>::iterator it = m_templates.begin(); it != m_templates.end(); ++it) {
   if (it->second->IsStyled()) it->second->StyleHists();
 }
}

void AnalysisClass::MakeStacks(){
  m_stacks.clear();
  for (std::vector<string>::iterator s = m_vars.begin(); s!=m_vars.end() ; ++s){
    string name = (*s);
    THStack* stack = MakeStack(name);
    m_stacks.insert(std::pair<string,THStack*>(name, stack));
  }
  for (std::vector<string>::iterator s = m_2Dvars.begin(); s!=m_2Dvars.end() ; ++s){
    string name = (*s);
    THStack* stack = MakeCombStack(name);
    m_stacks.insert(std::pair<string,THStack*>(name, stack));
  }
};




THStack* AnalysisClass::MakeStack(string name){
  THStack *hs0 = new THStack((m_name + "_" + name+"_stack").c_str(), name.c_str());
  //Get data entries
  //double dataN = m_data->GetVar(name)->GetHist()->Integral();

  for (std::vector<string>::iterator s = m_stackorder.begin() ; s != m_stackorder.end() ; ++s){
    if (m_templates.find(*s) == m_templates.end()) {
      cout<<"Template "<<*s<<" not found"<<endl;
      return hs0;
    }

    Template* temp = m_templates.at(*s);
    Var* var = temp->GetVar(name);
    TH1F* hist = var->GetHist();

    if (temp->IsStyled()) {
      temp->Style(hist);
    }
    //Is this a dangerous bit of code?
    //hs0->Add((TH1F*)hist->Clone());
    TH1F* histClone = new TH1F(*hist);
    hs0->Add(histClone);
       
  }
  return hs0;
}

THStack* AnalysisClass::GetStack(string name){
  if (m_stacks.find(name) != m_stacks.end()) return m_stacks.at(name);
  else {
    cout<<"No stack of the name "<<name<<" present"<<endl;
    return 0;
  }

}


THStack* AnalysisClass::MakeCombStack(string name){
  THStack *hs0 = new THStack((m_name + "_" + name+"_stack").c_str(), name.c_str());
  //Get data entries
  //double dataN = m_data->Get2DVar(name)->GetHist()->Integral();
  
  for (std::vector<string>::iterator s = m_stackorder.begin() ; s != m_stackorder.end() ; ++s){
    Template* temp = m_templates.at(*s);
    Var2D* var = temp->Get2DVar(name);
    TH1F* hist = var->GetCombHist();
    //double scaleF = temp->GetNormEvts()/temp->GetEvents();
    //hist->Scale(scaleF); 
    if (temp->IsStyled()) {
      temp->Style(hist);
    }
    //hs0->Add(hist);
    
    //Is this a dangerous bit of code?
    //hs0->Add((TH1F*)hist->Clone());
    TH1F* histClone = new TH1F(*hist);
    hs0->Add(histClone);
    
    
  }
  //hs0->GetHists()->SetOwner();
  return hs0;
}


double AnalysisClass::GetLumi(TFile* f){
  TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
  double lumi = 0;
  double lumi_job;
  int nentries = lumit->GetEntries();
  lumit->SetBranchAddress("IntegratedLuminosity", &lumi_job);
  for (int i = 0 ; i < nentries; ++i){
    lumit->GetEntry(i);
    lumi += lumi_job;
  }
  


  //lumit->Draw("IntegratedLuminosity>>lumihist()","","goff");
  //TH1F* lumihist=(TH1F*)gDirectory->Get("lumihist");
  //double Lumi = (double)lumihist->GetMean()*lumihist->GetEntries();
  return lumi;

}
double AnalysisClass::GetLumiError(TFile* f){
  TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
  double lumierr = 0;
  double lumierr_job;
  int nentries = lumit->GetEntries();
  lumit->SetBranchAddress("IntegratedLuminosityErr", &lumierr_job);
  for (int i = 0 ; i < nentries; ++i){
    lumit->GetEntry(i);
    lumierr += lumierr_job;
  }
  /*  TTree* lumit = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
  lumit->Draw("IntegratedLuminosityErr>>lumierrhist()","","goff");
  TH1F* lumierrhist=(TH1F*)gDirectory->Get("lumierrhist");
  double LumiError = (double)lumierrhist->GetMean()*lumierrhist->GetEntries();
  */
  return lumierr;
}

void AnalysisClass::SetSelCut(TCut cut){
  m_selcut = cut;
}


string AnalysisClass::GetName(){return m_name;}

void AnalysisClass::Remove(string name){
  RemoveFromFit(name);
  RemoveFromStack(name);
}

void AnalysisClass::Replace(string toRemove, string toAdd){
  ReplaceInStack(toRemove, toAdd);
  ReplaceInFit(toRemove, toAdd);
}

void AnalysisClass::ReplaceInStack(string toRemove, string toAdd){
  bool found = false;
  for (unsigned int i = 0 ; i < m_stackorder.size() ; ++i){
    if (strcmp(m_stackorder.at(i).c_str(), toRemove.c_str()) == 0){
      found = true;
      m_stackorder[i] = toAdd;
    }
  }
  if (!found) cout<<"No template "<<toRemove<<" found in stack"<<endl;
}
void AnalysisClass::ReplaceInFit(string toRemove, string toAdd){
  bool found = false;
  for (unsigned int i = 0 ; i < m_fitorder.size() ; ++i){
    if (strcmp(m_fitorder.at(i).c_str(), toRemove.c_str()) == 0){
      found = true;
      m_fitorder[i] = toAdd;
    }
  }
  if (!found) cout<<"No template "<<toRemove<<" found in fit"<<endl;
}

void AnalysisClass::RemoveFromStack(string name){
  bool found = false;
  vector<string>::iterator to_erase;
  for (vector<string>::iterator is = m_stackorder.begin(); is != m_stackorder.end(); ++is){
    if (strcmp((*is).c_str(), name.c_str()) == 0) {
      to_erase = is;
      found = true;
    }
  }
  if (found)  m_stackorder.erase(to_erase);
  else cout<<"No Template: "<<name<<" found in stack for class "<<m_name<<endl;
}
void AnalysisClass::RemoveFromFit(string name){
  bool found = false;
  vector<string>::iterator to_erase;
  for (vector<string>::iterator is = m_fitorder.begin(); is != m_fitorder.end(); ++is){
    if (strcmp((*is).c_str(), name.c_str()) == 0) {
      to_erase = is;
      found = true;
    }
  }
  if (found) m_fitorder.erase(to_erase);
  else cout<<"No Template: "<<name<<" found in fit for class "<<m_name<<endl;
}

void AnalysisClass::AddToStack(string name){
  bool found = false;
  vector<string>::iterator to_erase;
  for (vector<string>::iterator is = m_stackorder.begin(); is != m_stackorder.end(); ++is){
    if (strcmp((*is).c_str(), name.c_str()) == 0) {
      to_erase = is;
      found = true;
    }
  }
  if (!found)  m_stackorder.push_back(name);
  //else cout<<"No Template: "<<name<<" found in stack for class "<<m_name<<endl;
}
void AnalysisClass::AddToFit(string name){
  bool found = false;
  vector<string>::iterator to_erase;
  for (vector<string>::iterator is = m_fitorder.begin(); is != m_fitorder.end(); ++is){
    if (strcmp((*is).c_str(), name.c_str()) == 0) {
      to_erase = is;
      found = true;
    }
  }
  if (!found) m_fitorder.push_back(name);
}

void AnalysisClass::SetToFit(vector<string> toFit){ m_fitorder = toFit;}
vector<string> AnalysisClass::GetToFit(){ return m_fitorder;}

void AnalysisClass::SetToStack(vector<string> toStack){ m_stackorder = toStack;}
vector<string> AnalysisClass::GetToStack(){ return m_stackorder;}



TCanvas* AnalysisClass::DrawFitted(){
  if ((!m_data) || (!m_fitter)) {
    cout<<"No fit to draw"<<endl;
    return 0;
  }
  string var = m_fitter ? m_fitter->GetVar(): "";
  TCanvas* c1 = new TCanvas();
  THStack* stack = GetStack(var);
  if (stack){
    stack->Draw("hist");
  }
  else cout<<"No stack found"<<endl;
  if (m_data->GetVar(var)){
    m_data->GetVar(var)->GetHist()->Draw("e1same");
  }
  if (m_fit){
    m_fit->Draw("same");
  }
  else cout <<"No variable "<<var<<" found in data template"<<endl;
  return c1;
}

void AnalysisClass::Run(){
  ApplyCuts();
  FillVars();
}

#ifdef WITHPYTHON
void AnalysisClass::AddTemplate1_py(string name, PyObject* pyt){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  AddTemplate(name, t);
}
void AnalysisClass::AddTemplate2_py(string name, PyObject* pyt, PyObject* pycut){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pycut));
  AddTemplate(name, t, *cut);
}
void AnalysisClass::AddTemplate3_py(string name, PyObject* pyt, int color){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  enum EColor col = (enum EColor)color;
  AddTemplate(name, t, col);
}

void AnalysisClass::AddTemplate4_py(string name, PyObject* pyt, PyObject* pycut, int color){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pycut));
  enum EColor col = (EColor)color;
  AddTemplate(name, t, *cut, col);
}
void AnalysisClass::AddTemplate5_py(Template* t){
  AddTemplate(t);
}
void AnalysisClass::AddData_py(string name, PyObject* pyt){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  AddData(name, t);
}
PyObject* AnalysisClass::RedoFit_py(string var){
  TFractionFitter* fit = RedoFit(var);
  return TPython::ObjectProxy_FromVoidPtr(fit, fit->ClassName());
}
PyObject* AnalysisClass::TFracFit_py(string var){
  TFractionFitter* fit = TFracFit(var);
  return TPython::ObjectProxy_FromVoidPtr(fit, fit->ClassName());
}
PyObject* AnalysisClass::TFracFit3_py(string var1, string var2){
  TFractionFitter* fit = TFracFit(var1, var2);
  return TPython::ObjectProxy_FromVoidPtr(fit, fit->ClassName());
}

PyObject* AnalysisClass::TFracFit2_py(string var, double lo, double hi){
  TFractionFitter* fit = TFracFit(var, lo, hi);
  return TPython::ObjectProxy_FromVoidPtr(fit, fit->ClassName());
}
PyObject* AnalysisClass::CombTFracFit_py(string var1, string var2){
  TFractionFitter* fit = CombTFracFit(var1, var2); 
  return TPython::ObjectProxy_FromVoidPtr(fit, fit->ClassName());
}
void AnalysisClass::AddFitter_py(string var){
  AddFitter(var);
}
void AnalysisClass::ApplyFitResults_py(){
  ApplyFitResults();
}
void AnalysisClass::SaveToFile1_py(){
  SaveToFile();
}
void AnalysisClass::SaveToFile2_py(string output){
  SaveToFile(output);
}
void AnalysisClass::AddVar1_py(string name, string var, int bins, double lo, double hi){
  return AddVar(name, var, bins, lo, hi);
}


void AnalysisClass::AddVar2_py(boost::python::list& varlist){
  if (len(varlist) == 5){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    int bins    = boost::python::extract<int>(varlist[2]);
    double lo   = boost::python::extract<double>(varlist[3]);
    double hi   = boost::python::extract<double>(varlist[4]);
    AddVar(name, var, bins, lo, hi);
  }
  else if (len(varlist) == 3){
    string name = boost::python::extract<string>(varlist[0]);
    string var  = boost::python::extract<string>(varlist[1]);
    boost::python::list ns = (boost::python::list)(varlist[2]);
    AddVar3_py(name, var, ns);
  }
}

void AnalysisClass::AddVar3_py(string name, string var, boost::python::list& ns){
  std::vector<double> edges;
  for (int i = 0; i < len(ns); ++i){
    double val = boost::python::extract<double>(ns[i]);
    edges.push_back(val);
  }
  AddVar(name, var, edges);
}

void AnalysisClass::AddVars_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    AddVar2_py(var);
  }
}

void AnalysisClass::Add2DVars_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    if (len(var) == 2){
      string var1 = boost::python::extract<string>(var[0]);
      string var2 = boost::python::extract<string>(var[1]);
      Add2DVar(var1, var2);
    }
    else{
      cout<<"Cannot Add 2D Variable - Wrong Dimensions"<<endl;
    }
  }
}

PyObject* AnalysisClass::GetStack_py(string name){
  THStack* stack = GetStack(name);
  if (stack){
    return TPython::ObjectProxy_FromVoidPtr(stack, stack->ClassName());
  }
  else {
    return 0;
  }

}
double AnalysisClass::GetLumi_py(PyObject* pyf){
  TFile* f = (TFile*)(TPython::ObjectProxy_AsVoidPtr(pyf));
  return GetLumi(f);
}
void AnalysisClass::SetSelCut_py(PyObject* pycut){

  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pycut));
  SetSelCut(*cut);

}
void AnalysisClass::SetToFit_py(boost::python::list& ns){ 
  vector<string> toFit;
  for (int i = 0; i < len(ns); ++i){
    string name = boost::python::extract<string>(ns[i]);
    toFit.push_back(name);
  }
  SetToFit(toFit);
}
boost::python::list AnalysisClass::GetToFit_py(){
  boost::python::list l;
  for (vector<string>::iterator im = m_fitorder.begin(); im != m_fitorder.end(); ++im){
    string name = (*im);
    l.append(name);
  }

  return l;
}

void AnalysisClass::SetToStack_py(boost::python::list& ns){ 
  vector<string> toStack;
  for (int i = 0; i < len(ns); ++i){
    string name = boost::python::extract<string>(ns[i]);
    toStack.push_back(name);
  }
  SetToStack(toStack);
}
boost::python::list AnalysisClass::GetToStack_py(){
  boost::python::list l;
  for (vector<string>::iterator im = m_stackorder.begin(); im != m_stackorder.end(); ++im){
    string name = (*im);
    l.append(name);
  }

  return l;
}
PyObject* AnalysisClass::DrawFitted_py(){
  TCanvas* c1 = DrawFitted();
  if (c1){
    return TPython::ObjectProxy_FromVoidPtr(c1, c1->ClassName());
  }
  else return 0;

    
}

#endif
