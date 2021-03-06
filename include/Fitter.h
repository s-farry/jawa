#include <iostream>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <TF1.h>
#include <JawaObj.h>
#ifdef WITHPYTHON
#include <TPython.h>
#include <boost/python.hpp>
#endif
using namespace std;

class Fitter : public JawaObj {
 public:
  Fitter();
  Fitter(TH1F* hist, TObjArray* array, vector<string> names, string var = "");
  Fitter(TH2F* hist, TObjArray* array, vector<string> names, string var = "");
  //void ConstrainParams();
  //void ConstrainParam(string name, double val, double pc);
  void TFracFit();
  //void RooFit();
  void AddConstraints(vector<string> names, vector<double> vals);
  void AddConstraint(string name, double val);
  void AddConstraint(string name, double val, double pc);
  void AddConstraint(string name, double val, double lo, double hi);
  void AddConstraints(double lo, double hi);
  void SetFitRange(int lo, int hi);
  void ExcludeBins(double evts);
  void AddTFractionFitter();
  void SetExclude(bool exclude = true);
  bool GetExclude();
  string GetVar();

  map<string, pair<double, double> > GetResults();
  void SetResults(map<string, pair<double, double> >);
  map<string, vector<double> > GetConstraints();
  vector<string> GetNames();


  TH1F* GetData();
  TObjArray* GetTemplates();
  //TFractionFitter* GetFitter();
#ifdef WITHPYTHON
  //For python
  Fitter(PyObject* hist, PyObject* array, boost::python::list& l, string var = "");
  Fitter(PyObject* hist, boost::python::list& array, boost::python::list& l, string var = "");
  PyObject* GetTemplates_py();
  PyObject* GetData_py();
  //PyObject* GetFitter_py();
  boost::python::list GetNames_py();
  boost::python::list GetConstraints_py();
  boost::python::list GetResults_py();
  void SetResults_py(boost::python::list& ns);
  void AddConstraint_py(boost::python::list& ns);
  void AddConstraints_py(boost::python::list& ns);
  void AddConstraints2_py(double lo, double hi);
#endif
  
 private:
  map<string, int> m_names;           // name of templates and their location
  TObjArray* m_toFit;                 // Templates to be fit
  TH1F* m_data;                       // Data to fit to
  TH2F* m_2ddata;                     // Data to fit to (2d)
  map<string, vector<double> > m_constraints;  // Constraints
  TFractionFitter* m_fit;
  map<string, pair<double, double> > m_results;
  int m_status;
  int m_lo;
  int m_hi;
  double m_clo;
  double m_chi;
  bool m_exclude;
  string m_var;
  bool m_2d;
  
};
