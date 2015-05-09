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
#include <Template.h>
#include <TParameter.h>
#include <boost/algorithm/string.hpp>

using namespace std;


void Template::Init(){

  m_reweight = false;
  m_fixed = false;
  m_asymm = false;
  m_reweightfunc = false;
  m_reweight_map = false;
  m_reweight_hist = false;
  m_verbose = false;
  m_fitFrac = 0;
  m_fitFracErr = 0;
  m_fillTree = false;
  m_tree = 0;
  m_outputevts = true;

}


void Template::OutputEvts(bool output){
  m_outputevts = output;
}

Template::Template(string name){
  m_name = name;
  Init();

}

//Template::~Template(){
  //  if (m_reweightTF1)  m_reweightTF1->Delete();
//  if (m_reweightTH1F) m_reweightTH1F->Delete();
//}

//Template::~Template(){
  //delete m_tree;
  //delete m_entryList;
//}


Template::Template(string name, PyObject* pyt, PyObject* pycut){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyt));
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pycut));
  //Template(name, t, *cut);
  m_name = name;
  m_trees.push_back(new Tree(name+"_tree", t, 1.0));
  m_selcut = (*cut);
  m_reweight = false;
  m_fixed = false;
  m_asymm = false;
  m_reweightfunc = false;
  m_reweight_map = false;
  m_reweight_hist = false;
  m_fillTree = false;
  m_verbose = false;
  m_fitFrac = 0;
  m_fitFracErr = 0;
  m_tree = 0;
  m_outputevts = true;
}

Template::Template(string name, TTree* t, TCut cut){
  Init();
  if (m_verbose) cout<<"Initialising as "<<name<<endl;
  m_name = name;
  m_trees.push_back(new Tree(name+"_tree", t, 1.0));

  m_selcut = cut;
  if (m_verbose) cout<<"Name is "<<m_name<<endl;
}

Template::Template(string name, TTree* t, TCut cut, enum EColor color){
  m_name = name;
  m_trees.push_back(new Tree(name+"_tree", t, 1.0));
  m_selcut = cut;
  SetStyle(color);
  Init();
}

Template Template::operator+(const Template& rhs){
  //Incomplete to say the least - don't use yet
  return Template(m_name+"_"+rhs.GetName());
}

map<string, Var*> Template::GetVariables(){return m_variables;}
map<string, Var2D*> Template::Get2DVariables(){return m_2Dvariables;}
map<string, Var3D*> Template::Get3DVariables(){return m_3Dvariables;}

Template::Template(string name, Template* t){
  Init();
  m_name = name;
  map<string, Var*> variables = t->GetVariables();
  map<string, Var2D*> variables2D = t->Get2DVariables();
  for (map<string, Var*>::iterator iv = variables.begin(); iv != variables.end(); ++iv){
    m_variables.insert(pair<string, Var*>((*iv).first, new Var((*iv).first, (*iv).second, name)));
  }
  for (map<string, Var2D*>::iterator iv = variables2D.begin(); iv != variables2D.end(); ++iv){
    m_2Dvariables.insert(pair<string, Var2D*>((*iv).first, new Var2D((*iv).first, (*iv).second, name)));
  }
  m_evts = t->GetEvents();
  m_normN = t->GetNormEvts();
}


Template::Template(string name, Template* A, Template* B){
  Init();
  m_name = name;
  map<string, Var*> variablesA = A->GetVariables();
  map<string, Var*> variablesB = B->GetVariables();
  map<string, Var2D*> variables2DA = A->Get2DVariables();
  map<string, Var2D*> variables2DB = B->Get2DVariables();
  map<string, Var3D*> variables3DA = A->Get3DVariables();
  map<string, Var3D*> variables3DB = B->Get3DVariables();
  for (map<string, Var*>::iterator iv = variablesA.begin(); iv != variablesA.end(); ++iv){
    if (variablesB.find((*iv).first) != variablesB.end()){
      m_variables.insert(pair<string, Var*>((*iv).first, new Var((*iv).first, (*iv).second, variablesB.at((*iv).first), name)));
    }
  }
  for (map<string, Var2D*>::iterator iv = variables2DA.begin(); iv != variables2DA.end(); ++iv){
    if (variables2DB.find((*iv).first) != variables2DB.end()){
      m_2Dvariables.insert(pair<string, Var2D*>((*iv).first, new Var2D((*iv).first, (*iv).second, variables2DB.at((*iv).first), name)));
    }
  }
  for (map<string, Var3D*>::iterator iv = variables3DA.begin(); iv != variables3DA.end(); ++iv){
    if (variables3DB.find((*iv).first) != variables3DB.end()){
      m_3Dvariables.insert(pair<string, Var3D*>((*iv).first, new Var3D((*iv).first, (*iv).second, variables3DB.at((*iv).first), name)));
    }
  }
  m_evts = A->GetEvents() + B->GetEvents();
  m_normN = A->GetNormEvts() + B->GetNormEvts();
}


Tree* Template::GetTree(){
  Tree* tree = 0;
  if (m_trees.size() > 0) tree = m_trees.at(0);
  return tree;

}

Tree* Template::GetTree(string name){
  Tree* tree = 0;
  for (unsigned int i = 0; i < m_trees.size(); ++i){
    if (strcmp((m_trees.at(i)->GetName()).c_str(), name.c_str()) == 0){
      tree = m_trees.at(i);
    } 
  }
  if (tree == 0) cout<<"Tree "<<name<<" not found in template "<<m_name<<endl;
  return tree;

}

Tree* Template::GetTree1_py(){
  return GetTree();
}
Tree* Template::GetTree2_py(string name){
  return GetTree(name);
}


void Template::SetTree_py(PyObject* pyObj){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  SetTree(t);
}

void Template::SetTree2_py(PyObject* pyObj, double w){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  SetTree(t, w);
}
void Template::AddTree_py(string name, PyObject* pyObj){
  TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  AddTree(name, t);
}

void Template::SetTrees(vector<TTree*>& trees){
  for (std::vector<TTree*>::iterator it = trees.begin(); it != trees.end(); ++it){
    m_trees.push_back(new Tree(m_name+"_tree", (*it), 1.0));
  }
}

void Template::SetTrees_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::object obj = ns[i];
    //PyObject* pyObj = boost::python::extract<PyObject*>(ns[i]);
    PyObject* pyObj = obj.ptr();
    TTree* t = (TTree*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
    SetTree(t);
  }
}

void Template::SetTree(TTree* t){
  m_trees.push_back(new Tree(m_name+"_tree", t, 1.0));
}
void Template::SetTree(TTree* t, double w){  
  m_trees.push_back(new Tree(m_name+"_tree", t, w));
}
void Template::AddTree(string name, TTree* t){
  m_trees.push_back(new Tree(name, t));
}
void Template::StyleHists(){
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    Var* var = (*iv).second;
    TH1F* hist = var->GetHist();
    Style(hist);
  }
}
void Template::Style(TH1F* hist){
  hist->SetLineColor(m_linecolor);
  hist->SetFillStyle(m_fillstyle);
  hist->SetFillColor(m_fillcolor);
}

Var* Template::GetVar(string name){
  if (m_variables.find(name) != m_variables.end()){
    return m_variables.at(name);
  }
  cerr<<"No Variable with name: "<<name<<endl;
  return 0;
}

Var2D* Template::Get2DVar(string name){
  if (m_2Dvariables.find(name) != m_2Dvariables.end()){
    return m_2Dvariables.at(name);
  }
  cerr<<"No 2D variable: "<<name<<endl;
  return 0;
}

Var2D* Template::Get2DVar(string name1, string name2){
  string name = name1+"_"+name2;
  return Get2DVar(name);
}

Var3D* Template::Get3DVar(string name1, string name2, string name3){
  string name = name1+"_"+name2+"_"+name3;
  return Get3DVar(name);
}

Var3D* Template::Get3DVar(string name){
  if (m_3Dvariables.find(name) != m_3Dvariables.end()){
    return m_3Dvariables.at(name);
  }
  cerr<<"No 3D variable: "<<name<<endl;
  return 0;
}

Var2D* Template::Get2DVar_py(string name1, string name2){
  return Get2DVar(name1, name2);
}

Var3D* Template::Get3DVar_py(string name1, string name2, string name3){
  return Get3DVar(name1, name2, name3);
}

void Template::SetFitFrac(double f, double err){
  m_fitFrac = f;
  m_fitFracErr = err;
}

void Template::SetFitFrac_py(double f, double err){
  SetFitFrac(f, err);
}

void Template::IsVerbose(){
  m_verbose=true;
}

void Template::SetSelCut(TCut cut){
  m_selcut = cut;
}

void Template::SetSelCut_py(PyObject* pyObj){
  TCut* cut = (TCut*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  SetSelCut((*cut));
}

void Template::ApplyCut(){
  m_evts = 0;
  for (unsigned int i = 0; i < m_trees.size(); ++i){
    ostringstream ss;
    ss<<">>myList"<<m_name<<"_"<<i;
    string label = ss.str();
    string trimlabel = ss.str();
    trimlabel.erase(0,2);
    m_trees.at(i)->GetTTree()->Draw(label.c_str(), m_selcut , "entrylist");
    const TEntryList* list = (TEntryList*)gDirectory->Get(trimlabel.c_str());
    m_entryLists.push_back(new TEntryList(*list));
    m_evts += m_entryLists.at(i)->GetN();
    m_trees.at(i)->GetTTree()->SetEntryList(0);
    //delete list;
    
  }
  m_normN = 0;
}

void Template::AddVar1_py(string name, string var, int bins, double lo, double hi){
  AddVar(name, var, bins, lo, hi);
}

void Template::AddVar(string name, string var, int bins, double lo, double hi, string prefix){
  if (prefix == "") prefix = m_name;

  //if (var.find("-") == std::string::npos && var.find("+") == std::string::npos && var.find("/") == std::string::npos
      //&& var.find("*") == std::string::npos)
    //{
      m_variables.insert(std::pair<string,Var*>(name, new Var(name, var, bins, lo, hi, prefix)));
      //}
  //else{
  //m_algvariables.insert(std::pair<string,AlgVar*>(name, new AlgVar(name, var, bins, lo, hi, prefix)));
    //}
}
void Template::AddVar(string name, string var, std::vector<double>& edges, string prefix){
  if (prefix == "") prefix = m_name;
  //if (var.find("-") == std::string::npos && var.find("+") == std::string::npos && var.find("/") == std::string::npos )
  //  {
      m_variables.insert(std::pair<string,Var*>(name, new Var(name, var, edges, prefix)));
  //  }
  //else{
  //  m_algvariables.insert(std::pair<string,AlgVar*>(name, new AlgVar(name, var, edges, prefix)));
    //}
}
/*
void Template::AddAlgVar(string name, string varexp, int bins, double lo, double hi, string prefix){
  if (prefix == "") prefix = m_name;
  m_algvariables.insert(std::pair<string,AlgVar*>(name, new AlgVar(name, varexp, bins, lo, hi, prefix)));
}
void Template::AddAlgVar(string name, string varexp, std::vector<double> edges, string prefix){
  if (prefix == "") prefix = m_name;
  m_algvariables.insert(std::pair<string,AlgVar*>(name, new AlgVar(name, varexp, edges, prefix)));
  }*/
void Template::Add2DVar(string var1, string var2){
  string name = var1+"_"+var2;
  Add2DVar(name, var1, var2, "");
}
void Template::Add3DVar(string var1, string var2, string var3){
  string name = var1+"_"+var2 + "_" + var3;
  Add3DVar(name, var1, var2, var3, "");
}
void Template::Add2DVar(string name, string var1, string var2, string prefix){
  if (m_variables.find(var1) == m_variables.end() ||
      m_variables.find(var2) == m_variables.end()){
    std::cout<<"Can't add 2D variable - "<<var1<<" - "<<var2<<std::endl;
    return;
  }
  Var* varn1 = m_variables.at(var1);
  Var* varn2 = m_variables.at(var2);

  if (prefix == "") prefix = m_name;

  m_2Dvariables.insert(std::pair<string,Var2D*>(name, new Var2D(name, varn1, varn2, prefix)));
}
void Template::Add3DVar(string name, string var1, string var2, string var3, string prefix){
  if (m_variables.find(var1) == m_variables.end() ||
      m_variables.find(var2) == m_variables.end() ||
      m_variables.find(var3) == m_variables.end() ){
    std::cout<<"Can't add 2D variable - "<<var1<<" - "<<var2<<std::endl;
    return;
  }
  Var* varn1 = m_variables.at(var1);
  Var* varn2 = m_variables.at(var2);
  Var* varn3 = m_variables.at(var3);

  if (prefix == "") prefix = m_name;

  m_3Dvariables.insert(std::pair<string,Var3D*>(name, new Var3D(name, varn1, varn2, varn3, prefix)));
}


void Template::Add2DVar_py(boost::python::list& var){
  if (len(var) == 2){
    string var1 = boost::python::extract<string>(var[0]);
    string var2 = boost::python::extract<string>(var[1]);
    Add2DVar(var1, var2);
  }
}
void Template::Add3DVar_py(boost::python::list& var){
  if (len(var) == 3){
    string var1 = boost::python::extract<string>(var[0]);
    string var2 = boost::python::extract<string>(var[1]);
    string var3 = boost::python::extract<string>(var[2]);
    Add3DVar(var1, var2, var3);
  }
  else{
    cout<<"Can't add 3D Var"<<endl;
  }
}

void Template::Add2DVars_py(boost::python::list& varlist){
  for (int i = 0; i < len(varlist); ++i){
    boost::python::list ns = (boost::python::list)(varlist[i]);
    Add2DVar_py(ns);
  }
}
void Template::Add3DVars_py(boost::python::list& varlist){
  for (int i = 0; i < len(varlist); ++i){
    boost::python::list ns = (boost::python::list)(varlist[i]);
    Add3DVar_py(ns);
  }
}

void Template::AddVar2_py(boost::python::list& varlist){
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

void Template::AddVar3_py(string name, string var, boost::python::list& ns){
  std::vector<double> edges;
  for (int i = 0; i < len(ns); ++i){
    double val = boost::python::extract<double>(ns[i]);
    edges.push_back(val);
  }
  AddVar(name, var, edges,"");
}

void Template::SetVars_py(boost::python::list& ns){
  m_variables.clear();
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    AddVar2_py(var);
  }
}

void Template::AddVars_py(boost::python::list& ns){
  for (int i = 0; i < len(ns); ++i){
    boost::python::list var = (boost::python::list)ns[i];
    AddVar2_py(var);
  }
}


void Template::SetName(string name){
  m_name = name;
}

string Template::GetName() const{
  return m_name;
}

int Template::GetEvents(){
  return m_evts;
}

bool Template::IsStyled(){
  return m_styled;

}

std::map<string, int> Template::SetBranches(vector<double>& output_idx){
  //Set a branch for every output variable in local m_tree
  int j = 0;
  map<string, int> name_idx;
  if (output_idx.size() < (m_variables.size() + 1)){
    cerr<<"Error in Template::SetBranches - Sizes not compatible"<<endl;
    return name_idx;
  }

  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    m_tree->Branch((*iv).second->GetName().c_str(), &output_idx.at(j), ((*iv).second->GetName()+"/D").c_str());
    name_idx[(*iv).second->GetName()] = j;
    j++;
  }
  //One branch reserved for the weight
  m_tree->Branch("w", output_idx.at(j));
  return name_idx;
}


void Template::FillVars(){
  vector<double> output_idx;
  map<string, int> tree_idx;
  if (m_fillTree){
    cout<<"Filling tree"<<endl;
    m_tree = new TTree((m_name+"_tree").c_str(),"tree");
    output_idx = vector<double>( m_variables.size() + 1, 0.0);
    tree_idx = SetBranches(output_idx);
  }
  
  //Loop over each tree in template
  for (unsigned int ti = 0; ti< m_trees.size(); ++ti){
    Tree* tree = m_trees.at(ti);
    TTree* t = tree->GetTTree();
    t->SetBranchStatus("*" , 0 );
    
    //Set the appropriate branches from the variables and reweight variables
    for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
      std::vector<string> varnames = (*iv).second->GetVarNames();
      tree->SetBranches(varnames);
    }
    if (m_verbose) cout<<"Branches set for vars"<<endl;
    for (std::vector<ReweightVar>::iterator iv = m_reweightvariables.begin(); iv!=m_reweightvariables.end();++iv){
      vector<string> rwnames;
      vector<Expr*> exprs = (*iv).GetExprs();
      for (vector<Expr*>::iterator ie = exprs.begin(); ie != exprs.end(); ++ie){
	vector<string> names = (*ie)->GetVarNames();
	rwnames.insert(rwnames.end(), names.begin(), names.end());
      } 
	  //
	  //rwnames = (*iv).GetExpr()->GetVarNames();
	  //}
	  //else{
	  //rwnames = (*iv).GetNames();
	  //}
      tree->SetBranches(rwnames);
    }
    
    if (m_verbose) cout<<"Branches set for reweighted variables"<<endl;
    if (m_entryLists.size() == 0) cerr<<"Cuts haven't been applied"<<endl;
    
    TEntryList* l = m_entryLists.at(ti);
    double tree_w = m_trees.at(ti)->GetWeight();
    int nentries = l->GetN();

    // Loop  over all entries in entrylist
    for (int jentry = 0 ; jentry < nentries ; ++jentry) {

      if (jentry%10000==0 && m_outputevts) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;

      int entry = l->GetEntry(jentry);
      t->GetEntry(entry);

      double w = tree_w;

      // Get The Weight
      for (std::vector<ReweightVar>::iterator iv = m_reweightvariables.begin(); iv!=m_reweightvariables.end();++iv){
	if ((*iv).GetNames().size() == 1){
	  string var = (*iv).GetName();
	  Expr* e  = (*iv).GetExprs()[0];
	  double val = tree->GetVal(e);
	  w = w * ((*iv).GetWeight(val));
	}
	else if ((*iv).GetNames().size() == 2){
	  Expr* e1 = (*iv).GetExprs().at(0);
	  Expr* e2 = (*iv).GetExprs().at(1);
	  double val1 = tree->GetVal(e1);
	  double val2 = tree->GetVal(e2);
	  w = w * ((*iv).GetWeight(val1, val2));
	}
	else if ((*iv).GetNames().size() == 4){
	  double val1 = tree->GetVal((*iv).GetNames().at(0));
	  double val2 = tree->GetVal((*iv).GetNames().at(1));
	  double val3 = tree->GetVal((*iv).GetNames().at(2));
	  double val4 = tree->GetVal((*iv).GetNames().at(3));
	  w = w * ((*iv).GetWeight(val1, val2, val3, val4));
	}
      }
      m_normN += w;
      
      // Loop through 1-D variables

      for (std::map<string,Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
	Var* var = (*iv).second;
	double output = tree->GetVal(var->GetExpr());
	var->FillHist(output, w);
	if (m_fillTree) output_idx.at(tree_idx[(*iv).second->GetName()]) = output;
	
      }
      
      // Loop through 2-D variables
      
      for (std::map<string,Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
	Var2D* var = (*iv).second;
	double xval = 0.0, yval = 0.0;
	Var* var1 = var->GetVar1();
	Var* var2 = var->GetVar2();
	xval = tree->GetVal(var1->GetExpr());
	yval = tree->GetVal(var2->GetExpr());
	(*iv).second->FillHist(xval, yval,w);
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
	(*iv).second->FillHist(xval, yval, zval, w);
	//j++;
      }


      if (m_fillTree) output_idx.at(tree_idx["w"]) = w;
      if (m_fillTree) m_tree->Fill();
    }
    t->SetBranchStatus("*",1);
    }
    if ( m_verbose ) cout<<"Finished looping over events"<<endl;
  }

void Template::SetNormEvts(double evts){
  m_normN = evts;
}


 
 void Template::LoadFromFile(const char* file){
   string output = GetRootName(file);
  TFile* f = new TFile(output.c_str());
  //TList* list = f->GetListOfKeys();

    
  TParameter<double>* totEvts     = (TParameter<double>*)f->Get("totEvts");
  TParameter<double>* normEvts    = (TParameter<double>*)f->Get("normEvts");
  TParameter<double>* fitFrac     = (TParameter<double>*)f->Get("fitFrac");
  TParameter<double>* fitFracErr  = (TParameter<double>*)f->Get("fitFracErr");

  m_evts = totEvts->GetVal();
  m_normN = normEvts->GetVal();
  m_fitFrac = fitFrac->GetVal();
  m_fitFracErr = fitFracErr->GetVal();
  //Now loop through and get variables
  /*
  for (int i = 0; i < list->GetEntries(); ++i){
    string name = list->At(i)->GetName();
    bool is2D = ((int)name.find("_")!=-1);
    
    if (name != "totEvts" && name!= "normEvts" && name!="fitFrac" 
	&& name != "TotalHist" && name != "PassHist" 
	&& f->cd(name.c_str())){
      if (!is2D) m_variables.insert(std::pair<string, EffVar>(name, EffVar(name, f, m_name)));
      else if (is2D) m_2Dvariables.insert(std::pair<string, EffVar2D>(name, EffVar2D(name, f, m_name)));
    }

    }*/
}

void Template::NormaliseToEvts1_py(double evts, bool fixed){
  NormaliseToEvts(evts,fixed);
}
void Template::NormaliseToEvts2_py(double evts){
  NormaliseToEvts(evts);
}
void Template::NormaliseToEvts(double evts, bool fixed){
  m_fixed = fixed;
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    (*iv).second->Scale((double)evts/m_normN);
  }
  /*for (std::map<string, AlgVar*>::iterator iv = m_algvariables.begin() ; iv != m_algvariables.end() ; ++iv ){
    (*iv).second->Scale((double)evts/m_normN);
    }*/
  for (std::map<string, Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
    (*iv).second->Scale((double)evts/m_normN);
  }
  m_normN = evts;
}
void Template::NormaliseToMC1_py(double xsec, double acc, double Lumi, double nEvts, bool fixed){
  return NormaliseToMC(xsec, acc, Lumi, nEvts, fixed);
}
void Template::NormaliseToMC2_py(double xsec, double acc, double Lumi, double nEvts){
  return NormaliseToMC(xsec, acc, Lumi, nEvts);
}

void Template::NormaliseToMC(double xsec, double acc, double Lumi, double nEvts, bool fixed){
  cout<<"Normalising to MC for: "<<m_name<<endl;
  double scale = xsec * Lumi * acc / nEvts;

  Scale(scale, fixed);
}

void Template::Unscale(){
  Scale(m_evts/m_normN);
}

void Template::Scale1_py(double scale, bool fixed){
  Scale(scale,fixed);
}
void Template::Scale2_py(double scale){
  Scale(scale);
}
void Template::Scale(double scale, bool fixed){
  m_fixed=fixed;
  m_normN = m_normN * scale;
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    (*iv).second->Scale(scale);
  }
  /*
  for (std::map<string, AlgVar*>::iterator iv = m_algvariables.begin() ; iv != m_algvariables.end() ; ++iv ){
    (*iv).second->Scale(scale);
    }*/
  for (std::map<string, Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
    (*iv).second->Scale(scale);
  }
}

void Template::Reweight1_py(string var, PyObject* tf1){
  
  TF1* f = (TF1*)(TPython::ObjectProxy_AsVoidPtr(tf1));
  Reweight(var,*f);
}
void Template::Reweight2_py(string var, PyObject* pyObj){
  
  TH1F* hist = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  TF1*  f    = (TF1*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
  if (strcmp(hist->ClassName(),"TH1F") == 0){
    Reweight(var,hist);
  }
  else if (strcmp(f->ClassName(),"TF1") == 0){
    Reweight(var,*f);
  }
  else{
    cout<<hist->ClassName()<<" is not matched to any possibilities"<<endl;
  }
  //std::cout<<hist->ClassName()<<std::endl;
  //Reweight(var,hist);
}
void Template::Reweight3_py(string var, PyObject* th1f, string var2, PyObject* th1f2){
  
  TH1F* hist  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(th1f));
  TH1F* hist2 = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(th1f2));
  Reweight(var,hist, var2, hist2);
}
void Template::Reweight4_py(string var1, string var2, PyObject* th2f){
  TH2F* hist  = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  Reweight(var1, var2, hist);
}
void Template::Reweight5_py(string var1, string var2, PyObject* th2f, string var3, string var4, PyObject* th2f2){
  TH2F* hist   = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f));
  TH2F* hist2  = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(th2f2));

  Reweight(var1, var2, hist, var3, var4, hist2);
}
void Template::Reweight6_py(string leaf){
  Reweight(leaf);
}



void Template::Reweight(string var, TF1 f){
  m_reweight = true;
  m_reweightfunc = true;
  //m_reweightTF1 = (TF1*)f->Clone();
  m_reweightTF1 = f;
  //m_reweightTF1->SetBit(kCanDelete);
  m_reweightvar = var;
  ReweightVar v = ReweightVar(var,f);
  m_reweightvariables.push_back(v);
}
void Template::Reweight(string var, std::map<int,double> weights){
  m_reweight = true;
  m_reweightvar = var;
  m_reweight_map = true;
  m_reweightmap = weights;
  m_reweightvariables.push_back(ReweightVar(var,weights));
}
void Template::Reweight(string var, TH1F* scale){
  m_reweight = true;
  m_reweightvar = var;
  m_reweight_hist = true;
  //m_reweightTH1F = (TH1F*)scale->Clone();
  m_reweightTH1F = (TH1F*)scale->Clone();
  m_reweightvariables.push_back(ReweightVar(var,scale));
}
void Template::Reweight(string var, TH1F* scale, string var2, TH1F* scale2){
  //m_reweight = true;
  //m_reweightvar = var;
  //m_reweight_hist = true;
  //m_reweightTH1F = (TH1F*)scale->Clone();
  //m_reweightTH1F = (TH1F*)scale->Clone();
  m_reweightvariables.push_back(ReweightVar(var,scale, var2, scale2));
}

void Template::Reweight(string var1, string var2, TH2F* scale){
  m_reweightvariables.push_back(ReweightVar(var1,var2, scale));
}
void Template::Reweight(string var1, string var2, TH2F* scale, string var3, string var4, TH2F* scale2){
  m_reweightvariables.push_back(ReweightVar(var1,var2, scale, var3, var4, scale2));
}

void Template::Reweight(string weightname){
  m_reweightvariables.push_back(ReweightVar(weightname));
}

void Template::SetStyle(enum EColor fillcolor){
  m_styled = true;
  m_linecolor = fillcolor;
  m_fillstyle = 1001;
  m_fillcolor = fillcolor;
  //m_linecolor = TColor(linecolor);
  //m_fillstyle = 0;
  //m_fillcolor = 0;


};

void Template::SetStyle(enum EColor fillcolor, enum EColor linecolor){
  m_styled = true;
  m_linecolor = linecolor;
  m_fillstyle = 1;
  m_fillcolor = fillcolor;


};
void Template::SetStyle(enum EColor fillcolor, enum EColor linecolor, Style_t fillstyle){
  m_styled = true;
  m_linecolor = linecolor;
  m_fillstyle = fillstyle;
  m_fillcolor = fillcolor;

  //m_linecolor = TColor(linecolor);
  //m_fillstyle = fillstyle;
  //m_fillcolor = TColor(fillcolor);


};
void Template::SetStyle1_py(int fillcolor){
  enum EColor col = (EColor)fillcolor;
  SetStyle(col);
}
void Template::SetStyle2_py(int fillcolor, int linecolor){
  enum EColor fillcol = (EColor)fillcolor;
  enum EColor linecol = (EColor)linecolor;
  SetStyle(fillcol, linecol);
}
void Template::SetStyle3_py(int fillcolor, int linecolor, int fillstyle){
  enum EColor fillcol = (EColor)fillcolor;
  enum EColor linecol = (EColor)linecolor;
  Style_t fillsty = (Style_t)fillstyle;
  SetStyle(fillcol, linecol, fillsty);
}


void Template::SaveToFile(){
  cout<<"Saving "<<m_name<<" to file"<<endl;
  TFile* f = new TFile((m_name+".root").c_str(),"RECREATE");
  
  SaveToCurrentFile();
  //f->Write();
  gROOT->cd();
  f->Close();
}

TH1F* Template::GetHist(string name){
  if (m_variables.find(name) != m_variables.end()){
    return m_variables.at(name)->GetHist();
  }
  /*else if (m_algvariables.find(name) != m_algvariables.end()){
    return m_algvariables.at(name)->GetHist();
    }*/
  else return 0;

}

PyObject* Template::GetHist_py(string name){
  TH1F* hist = GetHist(name);
  if (hist) return TPython::ObjectProxy_FromVoidPtr(hist, hist->ClassName());
  else return 0;
}

void Template::SaveToCurrentFile(){
  for (std::map<string, Var*>::iterator iv = m_variables.begin() ; iv != m_variables.end() ; ++iv ){
    Var* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    if (m_asymm) var->GetAsymmHist()->Write((var->GetName()+"_Asymm").c_str());
  }
  /*
  for (std::map<string, AlgVar*>::iterator iv = m_algvariables.begin() ; iv != m_algvariables.end() ; ++iv ){
    AlgVar* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    if (m_asymm) var->GetAsymmHist()->Write((var->GetName()+"_Asymm").c_str());
    }*/
  for (std::map<string, Var2D*>::iterator iv = m_2Dvariables.begin() ; iv != m_2Dvariables.end() ; ++iv ){
    Var2D* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    var->GetProfile()->Write((var->GetName()+"_prof_"+var->GetName1()).c_str());
    var->GetProfile2()->Write((var->GetName()+"_prof_"+var->GetName2()).c_str());
    var->GetCombHist()->Write((var->GetName()+"_comb").c_str());
  }
  for (std::map<string, Var3D*>::iterator iv = m_3Dvariables.begin() ; iv != m_3Dvariables.end() ; ++iv ){
    Var3D* var = (*iv).second;
    var->GetHist()->Write(var->GetName().c_str());
    //var->GetProfile()->Write((var->GetName()+"_prof_"+var->GetName1()).c_str());
    //var->GetProfile2()->Write((var->GetName()+"_prof_"+var->GetName2()).c_str());
    //var->GetCombHist()->Write((var->GetName()+"_comb").c_str());
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

double Template::GetFitFrac(){
  return m_fitFrac;

}

bool Template::IsFixed(){
  return m_fixed;
}

double Template::GetNormEvts(){
  return m_normN;
}

int Template::GetEvts(){
  return m_evts;
}

void Template::AddAsymmetry(string eta1, string eta2){
  m_asymm=true;
  m_asymmvar1 = eta1;
  m_asymmvar2 = eta2;

}

void Template::PrintVars(){
  for (std::map<string,Var*>::iterator it = m_variables.begin(); it != m_variables.end(); ++it){
    std::cout<< it->first <<std::endl;
  }
}

boost::python::list Template::GetVariables_py(){
  boost::python::list l;
  for (map<string, Var* >::iterator im = m_variables.begin(); im != m_variables.end(); ++im){
    l.append((*im).second);
  }
  return l;
}


string Template::GetRootName(const char* file){
  string output = "";
  if ( strlen( file ) == 0 ) {
    if ( m_name.size() == 0 ) 
      {
	output = "output.root";
      }
    else 
      {
	ostringstream filess;
	filess<<m_name<<".root";
	output = filess.str();
      }
  }
  else{
    output = file;
  }
  return output;
}
