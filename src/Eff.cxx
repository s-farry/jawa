#include <iostream>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <TCut.h>
#include <TObjArray.h>
#include <Eff.h>
#include <TEfficiency.h>
#include <math.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TEntryList.h>

using namespace std;
Eff::Eff(){
  m_name = "";
  m_eff = 0.0;
  m_effErrHi = 0.0;
  m_effErrLo = 0.0;
}
Eff::Eff(const char* name, double eff, double efferrhi, double efferrlo){
  m_name = name;
  m_eff = eff;
  m_effErrHi = efferrhi;
  m_effErrLo = efferrlo;
  m_Ntot = 0;
  m_Npass = 0;
  
}

Eff::Eff(const char* name, int Ntot, int Npass){
  m_name = name;
  m_Ntot = Ntot;
  m_Npass = Npass;

  if (Ntot > 0 && Npass > 0){
    m_eff = (double)Npass/Ntot;
    m_effErrHi = TEfficiency::ClopperPearson(Ntot, Npass, 0.68, true) - m_eff;
    m_effErrLo = m_eff - TEfficiency::ClopperPearson(Ntot, Npass, 0.68, false);
  }
  else {
    m_eff = 0;
    m_effErrHi = 0;
    m_effErrLo = 0;
  }

}
Eff::Eff(const char* name, double Ntot, double Npass){
  m_name = name;
  m_Ntot = Ntot;
  m_Npass = Npass;
  if (Ntot > 0 && Npass > 0){
    m_eff = (double)Npass/Ntot;
    m_effErrHi = TEfficiency::ClopperPearson(Ntot, Npass, 0.68, true) - m_eff;
    m_effErrLo = m_eff - TEfficiency::ClopperPearson(Ntot, Npass, 0.68, false);
  }
  else{
    m_eff = 0;
    m_effErrHi = 0;
    m_effErrLo = 0;
  }

}

double Eff::GetEff(){
  return m_eff;
}
double Eff::GetEffErrHi(){
  return m_effErrHi;
}
double Eff::GetEffErrLo(){
  return m_effErrLo;
}
void Eff::SetEff(double eff){
  m_eff = eff;
}
void Eff::SetEffErrHi(double efferrhi){
  m_effErrHi = efferrhi;
}
void Eff::SetEffErrLo(double efferrlo){
  m_effErrLo = efferrlo;
}
void Eff::Print(){
  cout<<m_name<<": "<<m_eff<<" + "<<m_effErrHi<<" - "<<m_effErrLo<<endl;
  cout<<"Entries: "<<m_Npass<<"/"<<m_Ntot<<endl;

}
void Eff::AddSystematic(double pc){
  m_effErrHi = sqrt(pow(m_effErrHi,2) + pow(m_eff*pc,2));
  m_effErrLo = sqrt(pow(m_effErrLo,2) + pow(m_eff*pc,2));

}
void Eff::AddInvSystematic(double pc){
  m_effErrHi = sqrt(pow(m_effErrHi,2) + pow((1-m_eff)*pc,2));
  m_effErrLo = sqrt(pow(m_effErrLo,2) + pow((1-m_eff)*pc,2));

}

char const* greet()
{
   return "hello, world";
}
