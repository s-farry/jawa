#include <JawaObj.h>
#include <ctime>

using namespace std;

int  NullBuffer::overflow( int c ) { return c; };
NullStream::NullStream() : std::ostream(&m_sb) {};

JawaObj::JawaObj(){
  m_name = "default";
  m_class = "default";
  m_outputlevel = JawaObj::Info;
  m_null = new NullStream();
  //m_out = &std::cout;
}
JawaObj::JawaObj(string className, string name){
  m_class = className;
  m_name  = name;
  m_outputlevel = JawaObj::Info;
  m_null = new NullStream();
}
JawaObj::JawaObj(string className){
  m_class = className;
  m_name  = "default";
  m_outputlevel = JawaObj::Info;
  m_null = new NullStream();

}
void JawaObj::debug(const char* msg){
  output(msg, JawaObj::Debug);
}
void JawaObj::info(const char* msg){
  output(msg, JawaObj::Info);
}
void JawaObj::verbose(const char* msg){
  output(msg, JawaObj::Verbose);
}

void JawaObj::progress(double progress, double width){
  std::cout<< GetTime()<<" INFO [";
  int pos = width * progress;
  for (int i = 0 ; i < width; ++i){
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout<<" ";
  }
  std::cout<<" ] " << int(progress * 100.0) << "%";
  if (progress != 1) {
    std::cout<<"\r";
  }
  else{
    std::cout<<"\n";
  }
  std::cout.flush();
}

string JawaObj::GetTime(){
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  
  strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
  std::string str(buffer);
  return str;
}

ostream& JawaObj::msg(string logname, OutputLevel o){
  string logname_out(logname);
  logname_out.resize(7, ' ');
  string output(m_class);
  output.resize(15, ' ');
  string name(m_name);
  name.resize(18,' ');
  //first set to NULL
  //ostream& stream = *m_null;
  //Now check if it should be to cout
  //if 
  //ostream null_stream(&null_buffer);
  //cout<<m_null<<endl;
  //ostream& stream = (*m_null);
  //NullStream a;
  ostream& null_stream = *m_null;
  if (o > m_outputlevel) return null_stream;
  ostream& stream = cout;
  stream<<GetTime()<<" "<<logname_out<<" \t "<<output<<" \t "<<name<<"\t ";
  return stream;
}

ostream& JawaObj::info(){return msg("INFO", JawaObj::Info);}
ostream& JawaObj::verbose(){return msg("VERBOSE", JawaObj::Verbose);}
ostream& JawaObj::debug(){return msg("DEBUG", JawaObj::Debug);}

ostream& JawaObj::GetStream(){ return cout;}

//void JawaObj::debug(ostringstream& msg){
//  output(msg, JawaObj::Debug);
//}
//void JawaObj::info(ostringstream& msg){
//  output(msg, JawaObj::Info);
//}
//void JawaObj::verbose(ostringstream& msg){
//  output(msg, JawaObj::Verbose);
//}
void JawaObj::output(const char* msg, OutputLevel o){
  if (m_outputlevel < o) return;
  string output(m_class);
  output.resize(8, ' ');
  string name(m_name);
  name.resize(8,' ');
  printf("%.8s.%.8s : %s \n", output.c_str(), name.c_str(), msg);
}

/*
void JawaObj::output(ostringstream& s, OutputLevel o){
  if (m_outputlevel < o) return;
  string output(m_class);
  output.resize(8, ' ');
  string name(m_name);
  name.resize(8,' ');
  const char* msg = s.str().c_str();
  printf("%.8s : %.8s : %s \n", output.c_str(), name.c_str(), msg);
  }*/

void JawaObj::SetOutputLevel(int i){
  m_outputlevel = OutputLevel(i);
}
void JawaObj::SetOutputLevel(OutputLevel o){
  m_outputlevel = o;
}
void JawaObj::SetVerbose(bool verbose){ 
  if (verbose){
    m_outputlevel = JawaObj::Verbose;
  }
  else {
    m_outputlevel = JawaObj::Debug;
  }
}
//ostringstream JawaObj::msg(){
//  ostringstream s;
//  return s;
//}

/*
ostream JawaObj::info2(){
  ostream stream(cout.rdbuf());
  return stream;
  cout<<"hello"<<endl;
}
*/
string JawaObj::GetName(){return m_name;}
bool JawaObj::GetVerbose(){ return (m_outputlevel == JawaObj::Verbose);}

/*
JawaObj::output2::output2(){
  m_out = &std::cout;
}
template <typename T> JawaObj::output2& JawaObj::output2::operator << (T value){
  (*m_out) << value;
  return *this;
};


int JawaObj::outbuf::overflow (int c) {
  if (c != EOF){
    c = std::toupper(static_cast<char>(c),getloc());
    if (putchar(c) == EOF) {
      return EOF;
    }
  }
  return c;
}
*/
