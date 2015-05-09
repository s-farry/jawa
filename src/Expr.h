#ifndef Include_Expr_H
#define Include_Expr_H
#include <iostream>
#include <iostream>  
#include <sstream>  
#include <list>  
#include <stack>  
#include <map>  
#include <string>  
#include <vector>  
#include <iterator>  
#include <stdlib.h>
#ifdef WITHPYTHON 
#include <boost/python.hpp>
#endif

using namespace std;


class Expr{
 public:
  Expr(string varexp);
  Expr operator+(const Expr& rhs);
  Expr operator&(const Expr& rhs);
  Expr operator-(const Expr& rhs);
  Expr operator*(const Expr& rhs);
  Expr operator|(const Expr& rhs);
  Expr operator>(const Expr& rhs);
  Expr operator<(const Expr& rhs);

  string GetExpr();
  string GetExpr() const;
  double GetVal();
  double GetVal(std::vector<double>& input);
  void SetVarExp(string varexp);

  //pair< vector<string>, vector<string> > GetBracketExp(string varexp);
  //boost::python::list GetBracketExp_py(string varexp);
  vector<string>& GetVarNames();

  //vector<Expr*> GetExpressions();
  //boost::python::list GetExpressions_py();
  static bool is_number(const std::string& s);
  bool isInExpr(string var);
  int VarIdx(string var);
  string FindFunc(string name);
  size_t FindNextBracket(size_t pos);
  //
  // Print iterators in a generic way  
  bool isParenthesis( const std::string& token)   ;
  bool isOperator( const std::string& token)      ;
  bool isFunction( const std::string& token)      ;
  bool isAssociative( const std::string& token, const int& type)   ;
  //template< typename T, typename InputIterator > ;
  int cmpPrecedence( const std::string& token1, const std::string& token2 )    ;
  bool infixToRPN( const std::vector<std::string>& inputTokens,     
		   const int& size,     
		   std::vector<std::string>& strArray )   ;

  std::vector<std::string> getExpressionTokens( const std::string& expression )    ; // make tokens from expression
  std::vector<std::string>& getExpressionTokens( )    ; // get them if they exist already (reference)

  std::vector<std::string> getRPN( )    ;


  //double RPNtoDouble( std::vector<std::string> tokens )          ;
  double RPNtoDouble( )          ;
  double RPNtoDouble( std::vector<double>& input )          ;

  #ifdef WITHPYTHON
  //Python functions
  boost::python::list getRPN_py(  )    ;
  boost::python::list getExpressionTokens_py( const std::string& expression )    ;
  boost::python::list getExpressionTokens2_py(  )    ;
  boost::python::list infixToRPN_py( boost::python::list& inputTokens, const int& size )   ;
  boost::python::list GetVarNames_py();
  string GetExpr_py();
  double GetVal_py(boost::python::list& input);
  double GetVal3_py();
  #endif



 private:

  //TH1F* m_dummy;
  string m_varexp;
  vector<string> m_varnames; // All the variable names
  vector<string> m_tokens;  // Tokens for shunting yard algorithm
  vector<string> m_values; // variable names + numbers
  vector<string> m_rpn;
  //vector<string> m_operations;
  //vector<string> m_delimiters;
  //Operator map is kept here - not ideal
  //map<string, pair<int,int> > opmap;

};

#endif