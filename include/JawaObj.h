#ifndef Include_JawaObj
#define Include_JawaObj
// classes example
#include <iostream>
#include <streambuf>
#include <sstream>
#include <cstdio>
#include <locale>

using namespace std;

/**
Base class containing name and output functions
 */


class NullBuffer : public std::streambuf
{
  virtual int         overflow( int c ) ;
};

class NullStream : public std::ostream{
 public:
  NullStream();
 private:
  NullBuffer m_sb;
};

class JawaObj{
 protected:
  string m_name;
  string m_class;
  //std::ostream* m_out;
 public:
  JawaObj();
  JawaObj(string className);
  JawaObj(string className, string name);
  enum OutputLevel{Info=0, Debug=1, Verbose=2};
  void output (const char* msg, enum OutputLevel);
  void output (ostringstream& msg, enum OutputLevel);
  void info   (const char* msg);
  ostream& msg   (string logname, enum OutputLevel);
  ostream& info    ();
  ostream& verbose ();
  ostream& debug   ();
  ostream& GetStream();
  string GetTime();
  
  void debug  (const char* msg);
  void verbose(const char* msg);
  //void info   (ostringstream& msg);
  //void debug  (ostringstream& msg);
  //void verbose(ostringstream& msg);
  void SetOutputLevel(int i);
  void SetOutputLevel(OutputLevel o);
  void SetVerbose(bool verbose);
  bool GetVerbose();
  string GetName();
  void progress(double progress, double width);

  //ostream info2();

  /*
  class output2{
  public:
    output2();
    ostream* m_out;
    template <typename T> output2& operator << (T value);
  };
  */

  /*class outbuf : public std::streambuf
    {
      virtual int overflow(int c);
    };
  */
 public:
  NullStream* m_null;
  
  //outbuf ob;
  //ostream*  m_null;

 private:
  OutputLevel m_outputlevel;
};
#endif
