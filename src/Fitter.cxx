#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>
#include <Fitter.h>
#include <TParameter.h>
#include <TMinuit.h>
#include <TFitter.h>
#include <boost/algorithm/string.hpp>

using namespace std;

Fitter::Fitter() : JawaObj("Fitter"){
  m_toFit = new TObjArray();
  m_status = 0;
  m_clo = 0;
  m_chi = 0;
  m_lo = 0;
  m_hi = 0;
  m_fit = 0;
  m_var = "";
  m_exclude = true;
  m_2d = false;
}

Fitter::Fitter(TH1F* hist, TObjArray* array, vector<string> names, string var) : JawaObj("Fitter"){
  m_data = hist;
  m_toFit = array;
  m_status = 0;
  m_clo = 0;
  m_chi = 0;
  m_lo = 0;
  m_hi = 0;
  m_fit = 0;
  m_exclude = true;
  m_var = var;
  m_2d = false;
  for (unsigned int i = 0; i < names.size(); ++i) m_names[names.at(i)] = i;

}

Fitter::Fitter(TH2F* hist, TObjArray* array, vector<string> names, string var) : JawaObj("Fitter"){
  m_2ddata = hist;
  m_toFit = array;
  m_status = 0;
  m_clo = 0;
  m_chi = 0;
  m_lo = 0;
  m_hi = 0;
  m_fit = 0;
  m_exclude = true;
  m_var = var;
  m_2d = true;
  for (unsigned int i = 0; i < names.size(); ++i) m_names[names.at(i)] = i;

}

/*Fitter::Fitter(PyObject* hist, boost::python::list& array, boost::python::list& ns){
  m_data  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(hist));
  m_toFit = new TObjArray();
  m_status = 0;
  for (unsigned int i = 0; i < len(ns); ++i){
    string name = boost::python::extract<string>(ns[i]);
    m_names[name] = i;
  }
  for (unsigned int i = 0; i < len(array); ++i){
    string name = boost::python::extract<string>(ns[i]);
    m_names[name] = i;
  }
  }*/

void Fitter::AddConstraints(std::vector<string> names, std::vector<double> values){
  if (names.size() != values.size()) return;
  for (unsigned int i = 0; i < names.size(); ++i){
    AddConstraint(names.at(i), values.at(i));
  }
}

void Fitter::AddConstraint(string name, double value){
  vector<double> constraint(1, value);
  m_constraints[name] = constraint;
}
void Fitter::AddConstraint(string name, double value, double pc){
  vector<double> constraint;
  constraint.push_back(value);
  constraint.push_back(value * (1 - pc));
  constraint.push_back(value * (1 + pc));
  m_constraints[name] = constraint;
}
void Fitter::AddConstraint(string name, double value, double lo, double hi){
  vector<double> constraint;
  constraint.push_back(value);
  constraint.push_back(lo);
  constraint.push_back(hi);
  m_constraints[name] = constraint;
}



/* can't use getfitter with new root version...

void Fitter::ConstrainParams(){
  //Apply all the constraints that have been added

  //Overall constraints go first - like fractions between 0 and 1
  if (m_clo != 0 || m_chi != 0 ){
    for (int i = 0 ; i < m_toFit->GetEntries(); ++i){
    m_fit->Constrain(i + 1 , m_clo , m_chi);
    }
  }

  //Perform constraints for individual templates
  for (std::map<string, vector<double> >::iterator it = m_constraints.begin() ; it != m_constraints.end() ; ++it ){
    int idx  = m_names.at((*it).first);
    //Fix to 1 value
    if ((*it).second.size() == 1){
      double c = (*it).second[0];
      m_fit->GetFitter()->SetParameter(idx, (*it).first.c_str() ,c, 0.0, 0.0, 0.0);
      m_fit->GetFitter()->FixParameter(idx);
    }
    //Fix to a value and constrain to be between higher and lower values
    else if((*it).second.size() == 3){
      double c = (*it).second[0];
      double clo = (*it).second[1];
      double chi = (*it).second[2];
      m_fit->GetFitter()->SetParameter(idx, (*it).first.c_str() ,c, abs(c - clo)/3, clo, chi);
      //m_fit->Constrain(idx + 1, clo, chi);
    }
  }

}

*/

void Fitter::AddConstraints(double lo, double hi){
  //Add overall constraints - such as fractions between 0 and 1
  m_clo = lo;
  m_chi = hi;
}


void Fitter::AddTFractionFitter(){
  m_fit = new TFractionFitter(m_data, m_toFit);
}

void Fitter::SetExclude(bool exclude){m_exclude = exclude;}
bool Fitter::GetExclude(){return m_exclude;}

void Fitter::TFracFit(){
  //Perform the fit using TFractionFitter
  if (!m_2d)  m_fit = new TFractionFitter(m_data, m_toFit);
  else {
    info()<<"preparing 2d fit"<<endl;
    m_fit = new TFractionFitter(m_2ddata, m_toFit);
    info()<<"object created "<<m_fit<<endl;
    
  }
  //TFitter* fitter = (TFitter*)m_fit->GetFitter();
  if (!m_2d){
    if ( m_lo == 0 && m_hi == 0){
      m_fit->SetRangeX(1,m_data->GetNbinsX());
    }
    else m_fit->SetRangeX(m_lo, m_hi);
  }

  info()<<"range set"<<endl;
  //Perform constraints - no can do with new root
  //ConstrainParams();

  info()<<"parameters constrained"<<endl;

  //Exclude low bins for fit stability
  if (m_exclude) ExcludeBins(3);

  info()<<"fitting"<<endl;
  //Perform the fit
  m_status = m_fit->Fit();// perform the fit
  info()<<"fitted"<<endl;
  if (m_status == 0) { // If successful
    for(map<string, int>::iterator im = m_names.begin(); im != m_names.end(); ++im){
      double value, error;
      info()<<"getting result for "<<(*im).first<<endl;
      m_fit->GetResult((*im).second,value,error);
      //Fill results with value and error
      info()<<"got result"<<endl;
      m_results[(*im).first] = pair<double, double>(value,error);
      info()<<"result set"<<endl;
    }
  }
}

void Fitter::ExcludeBins(double evts){
  //Exclude bins with < X events to see if this changes anything
  if (!m_2d){
    for (int i = 0; i <= m_data->GetNbinsX() + 2 ; ++i){
      bool exclude = false;
      if (m_data->GetBinContent(i) ==0) {
	m_fit->ExcludeBin(i);
	exclude = true;
	info()<<"Excluding bin: "<<i<<endl;
      }
      double mcsum = 0;
      for (int j = 0 ; j <m_toFit->GetEntries(); j++){
	TH1F* hist = (TH1F*)m_toFit->At(j);
	mcsum += hist->GetBinContent(i)*hist->GetEntries()/hist->Integral();
      }
      if (mcsum < evts && !exclude) {
	m_fit->ExcludeBin(i);
	exclude = true;
	info()<<"Excluding bin: "<<i<<endl;
      }
    }
  }
  else{
    for (int i = 0; i <= m_2ddata->GetNbinsX() + 2 ; ++i){
      for (int j = 0; j <= m_2ddata->GetNbinsY() + 2 ; ++j){
	bool exclude = false;
	if (m_2ddata->GetBinContent(i, j) ==0) {
	  m_fit->ExcludeBin(m_2ddata->GetBin(i, j));
	  exclude = true;
	  info()<<"Excluding bin: "<<i<<" "<<j<<endl;
	}
	double mcsum = 0;
	for (int k = 0 ; k <m_toFit->GetEntries(); k++){
	  TH2F* hist = (TH2F*)m_toFit->At(k);
	  mcsum += hist->GetBinContent(i,j)*hist->GetEntries()/hist->Integral();
	}
	if (mcsum < evts && !exclude) {
	  m_fit->ExcludeBin(m_2ddata->GetBin(i,j));
	  exclude = true;
	  info()<<"Excluding bin: "<<i<<"_"<<j<<endl;
	}
      }
    }
  }
}

void Fitter::SetFitRange(int lo, int hi){
  m_lo = lo;
  m_hi = hi;
}

map<string, pair<double, double> > Fitter::GetResults(){
  return m_results;
}

void Fitter::SetResults(map<string, pair<double, double> > results){
  m_results = results;
}

string Fitter::GetVar(){return m_var;}

map<string, vector<double> > Fitter::GetConstraints(){
  return m_constraints;
}


vector<string> Fitter::GetNames(){
  vector<string> names;
  for (map<string, int>::iterator im = m_names.begin(); im != m_names.end(); ++im){
    names.push_back((*im).first);
  }
  return names;
}
/*
TFractionFitter* Fitter::GetFitter(){
  return m_fit;
  }*/

TObjArray* Fitter::GetTemplates(){
  return m_toFit;
}
TH1F* Fitter::GetData(){
  return m_data;
}

#ifdef WITHPYTHON
Fitter::Fitter(PyObject* hist, PyObject* array, boost::python::list& ns, string var) : JawaObj("Fitter"){
  m_data  = (TH1F*)(TPython::ObjectProxy_AsVoidPtr(hist));
  m_2ddata = (TH2F*)(TPython::ObjectProxy_AsVoidPtr(hist));
  std::string cName(m_data->ClassName());
  if (cName.find("TH1") == 0){
    info()<<"Initialising 1D fitter from python: "<<m_data->ClassName()<<endl;
    m_2ddata = 0;
    m_2d = false;
  }
  else{
    info()<<"Initialising 2D fitter from python: "<<m_data->ClassName()<<endl;
    m_data = 0;
    m_2d = true;
  }

  m_toFit = (TObjArray*)(TPython::ObjectProxy_AsVoidPtr(array));
  m_status = 0;
  m_clo = 0;
  m_chi = 0;
  m_lo = 0;
  m_hi = 0;
  m_fit = 0;
  m_exclude = true;
  m_var = var;
  for (unsigned int i = 0; i < len(ns); ++i){
    string name = boost::python::extract<string>(ns[i]);
    m_names[name] = i;
  }
}

PyObject* Fitter::GetTemplates_py(){
  return TPython::ObjectProxy_FromVoidPtr(m_toFit, m_toFit->ClassName());
}
PyObject* Fitter::GetData_py(){
  return TPython::ObjectProxy_FromVoidPtr(m_data, m_data->ClassName());
}
/*
PyObject* Fitter::GetFitter_py(){
  return TPython::ObjectProxy_FromVoidPtr(m_fit, m_fit->ClassName());

}*/
void Fitter::AddConstraint_py(boost::python::list& ns){
  if (len(ns) == 2){
    string name = boost::python::extract<string>(ns[0]);
    double val  = boost::python::extract<double>(ns[1]);
    AddConstraint(name, val);
  }
  else if (len(ns) == 3){
    string name = boost::python::extract<string>(ns[0]);
    double val  = boost::python::extract<double>(ns[1]);
    double pc   = boost::python::extract<double>(ns[2]);
    AddConstraint(name, val, pc);
  }
  else if (len(ns) == 4){
    string name = boost::python::extract<string>(ns[0]);
    double val  = boost::python::extract<double>(ns[1]);
    double lo   = boost::python::extract<double>(ns[2]);
    double hi   = boost::python::extract<double>(ns[3]);
    AddConstraint(name, val, lo, hi);
  }
}
void Fitter::AddConstraints_py(boost::python::list& ns){
  for (unsigned int i = 0; i < len(ns); ++i){
    boost::python::list ns2 = (boost::python::list)ns[i];
    string name = boost::python::extract<string>(ns2[0]);
    double val  = boost::python::extract<double>(ns2[1]);
    AddConstraint(name, val);
  }
}
void Fitter::AddConstraints2_py(double lo, double hi){
  AddConstraints(lo, hi);
}
boost::python::list Fitter::GetResults_py(){
  boost::python::list l;
  for (map<string, pair<double, double> >::iterator im = m_results.begin(); im != m_results.end(); ++im){
    string name = (*im).first;
    pair<double, double> valerr = (*im).second;
    //boost::python::object get_iter = boost::python::iterator<std::vector<double> >();
    boost::python::list a;
    a.append(valerr.first);
    a.append(valerr.second);
    boost::python::list b;
    b.append(name);
    b.append(a);
    l.append(b);
  }
  return l;
}
void Fitter::SetResults_py(boost::python::list& ns){
  map<string, pair<double, double> > results;

  for (unsigned int i = 0; i < len(ns); ++i){
    if (len(ns[i]) == 3){
      string name = boost::python::extract<string>(ns[i][0]);
      double res  = boost::python::extract<double>(ns[i][1]);
      double err  = boost::python::extract<double>(ns[i][2]);
      results[name] = pair<double, double>(res, err);
    }
    else{
      info()<<"Input to Results not correct - ignoring"<<endl;
    }
  }
  SetResults(results);
}
boost::python::list Fitter::GetConstraints_py(){
  boost::python::list l;
  for (map<string, vector<double> >::iterator im = m_constraints.begin(); im != m_constraints.end(); ++im){
    string name = (*im).first;
    vector<double> constraint = (*im).second;

    boost::python::list a;
    for (unsigned int i = 0; i < constraint.size(); ++i) a.append(constraint.at(i));
    boost::python::list b;
    b.append(name);
    b.append(a);
    l.append(b);
  }
  return l;
}
boost::python::list Fitter::GetNames_py(){
  boost::python::list names;
  for (map<string, int>::iterator im = m_names.begin(); im != m_names.end(); ++im){
    names.append((*im).first);
  }
  return names;
}

#endif
