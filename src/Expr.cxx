#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <Expr.h>
#include <boost/algorithm/string.hpp>
#include <TMath.h>

using namespace std;


int LEFT_ASSOC  = 0;      
int RIGHT_ASSOC = 1;   

// Map the different operators: +, -, *, / etc
typedef std::map< std::string, std::pair< int,int > > OpMap;
//typedef std::vector<std::string>::const_iterator cv_iter;  
//typedef std::string::iterator s_iter;

const OpMap::value_type assocs[] =
  {  OpMap::value_type( "+", std::pair<int,int>(  0, LEFT_ASSOC ) ),      
     OpMap::value_type( "-", std::pair<int,int>(  0, LEFT_ASSOC ) ),
     OpMap::value_type( ">", std::pair<int,int>(  2, LEFT_ASSOC ) ),      
     OpMap::value_type( "<", std::pair<int,int>(  2, LEFT_ASSOC ) ),   
     OpMap::value_type( ">=", std::pair<int,int>( 2, LEFT_ASSOC ) ),      
     OpMap::value_type( "<=", std::pair<int,int>( 2, LEFT_ASSOC ) ),     
     OpMap::value_type( "==", std::pair<int,int>( 2, LEFT_ASSOC ) ),     
     OpMap::value_type( "!=", std::pair<int,int>( 2, LEFT_ASSOC ) ),
     OpMap::value_type( "*", std::pair<int,int>(  5, LEFT_ASSOC ) ),        
     OpMap::value_type( "&&", std::pair<int,int>( 1, LEFT_ASSOC ) ),       
     OpMap::value_type( "||", std::pair<int,int>( 1, LEFT_ASSOC ) ),
     OpMap::value_type( "/", std::pair<int,int>(  5, LEFT_ASSOC ) ),
     OpMap::value_type( "^", std::pair<int,int>(  6, RIGHT_ASSOC ) )};
     
const OpMap opmap( assocs, assocs + sizeof( assocs ) / sizeof( assocs[ 0 ] ) );    


Expr::Expr(string varexp) : JawaObj("Expr"){
  m_varexp = varexp;
  m_tokens = Tokenize(varexp);
  bool c = infixToRPN(m_tokens, m_rpn);
  if (!c) {
    info() << "Parse Error in " << varexp << endl;
  }
}

Expr Expr::operator+(const Expr& rhs){
  return Expr("("+m_varexp+") + (" + rhs.GetExpr()+")");
}
Expr Expr::operator-(const Expr& rhs){
  return Expr("("+m_varexp+") - (" + rhs.GetExpr()+")");
}
Expr Expr::operator&(const Expr& rhs){
  return Expr("("+m_varexp+") && (" + rhs.GetExpr()+")");
}
Expr Expr::operator*(const Expr& rhs){
  return Expr("("+m_varexp+") * (" + rhs.GetExpr()+")");
}
Expr Expr::operator|(const Expr& rhs){
  return Expr("("+m_varexp+") || (" + rhs.GetExpr()+")");
}
Expr Expr::operator>(const Expr& rhs){
  return Expr("("+m_varexp+") > (" + rhs.GetExpr()+")");
}
Expr Expr::operator<(const Expr& rhs){
  return Expr("("+m_varexp+") < (" + rhs.GetExpr()+")");
}

double Expr::GetVal(){
  if (m_varnames.size() == 0){
    return RPNtoDouble();
  }
  else{
    info()<<"Parse Error : Value(s) ";
    for (unsigned int i = 0 ; i < m_varnames.size(); ++i){
      info()<<m_varnames.at(i);
    }
    info()<<" not found"<<endl;
    return -1.0;
  }
}

double Expr::GetVal(std::vector<double>& input){
  if (input.size() != m_varnames.size() || input.size() == 0) return -1;
  if (m_rpn.size() == 1 && input.size() == 1) return input.at(0);
  return RPNtoDouble(input);
}

double Expr::operator()(double* x, double* p){
  //cout<<x[0]<<" "<<p[0]<<endl;
  std::vector<double> input;
  int j = 0;
  for (unsigned int i = 0 ; i < m_varnames.size(); ++i){
    if ("x" == m_varnames.at(i)) {
      //cout<<"x: "<<x[0]<<endl;
      input.push_back(x[0]);
    }
    else {
      //cout<<j<<" "<<p[j]<<endl;
      input.push_back(p[j]);
      j++;
    }
  }
  return GetVal(input);
}


string Expr::GetExpr() const{ return m_varexp;}
string Expr::GetExpr() {return m_varexp;}
std::vector<string>& Expr::GetVarNames(){ return m_varnames;}

bool Expr::is_number(const std::string& s){
  //cout<<"checking is number for "<<s<<endl;
  //return (strspn( s.c_str(), "-.01234567890") == s.size() );
  //std::string::const_iterator it = s.begin();
  //while (it != s.end() && std::isdigit(*it)) ++it;
  //return !s.empty() && it == s.end();
  //std::string::const_iterator it = s.begin();
  if (!s.empty() && s.find_first_not_of("-.0123456789") == std::string::npos) return true;
  else return false;

}

// Test if token is an pathensesis    
bool Expr::isParenthesis( const std::string& token)          
{          
  return token == "(" || token == ")";        
}        
      
// Test if token is an operator          
bool Expr::isOperator( const std::string& token)          
{          
  return token == "+" || token == "-" ||        
  token == "*" || token == "/" || token == "^"||
  token == ">" || token == "<" || token == "&&"||
  token == ">=" || token == "<=" ||
  token == "|"  || token == "==" || token == "!=";
}        

bool Expr::isFunction( const std::string& token)          
{
  return token == "sqrt" || token == "abs" || token == "fabs" ||     
    token == "sin" || token == "cos" ||
    token == "tan" || token == "log" || token == "log10" ||
    token == "min" || token == "max" || token == "cosh" ||
    token == "sinh" || token == "tanh" || token == "exp" ||
    token == "CB" || token == "Voigtian" || token == "Gaussian";
}

      
// Test associativity of operator token          
bool Expr::isAssociative( const std::string& token, const int& type)          
{              
  //cout<<token<<" is associative?"<<endl;
  //cout<<"printing.... "<<opmap["+"].first<<endl;
  //bool inmap = (opmap.find(token) != opmap.end());
  //cout<<token<<" is in map? : "<<inmap<<endl;
  const std::pair<int,int> p = opmap.find( token )->second;          
  //cout<<"looking to return"<<endl;
  return p.second == type;          
}          
      
// Compare precedence of operators.          
int Expr::cmpPrecedence( const std::string& token1, const std::string& token2 )          
{          
  const std::pair<int,int> p1 = opmap.find( token1 )->second;        
  const std::pair<int,int> p2 = opmap.find( token2 )->second;        
  
  return p1.first - p2.first;          
}          
/*
boost::python::list Expr::infixToRPN_py(boost::python::list& inputTokens,   const int& size){
  std::vector<string> output;
  std::list<std::string> input;
  for (unsigned int i = 0; i < len(inputTokens); ++i){
    string s = boost::python::extract<string>(inputTokens[i]);
    input.push_back(s);
  }

  infixToRPN(input, size, output);
  boost::python::list outputlist;
  for (vector<string>::iterator is = output.begin(); is != output.end(); ++is){
    outputlist.append(*is);
  }
  return outputlist;

  }*/

// Convert infix expression format into reverse Polish notation          
bool Expr::infixToRPN( const std::list<std::string>& inputTokens,     
		       //		       const int& size,     
		       std::vector<std::string>& strArray )          
{     
  bool success = true;      
  
  //const int LEFT_ASSOC  = 0;      
  //const int RIGHT_ASSOC = 1;   
  std::list<std::string> out;
  std::stack<std::string> stack;
  
  // While there are tokens to be read      
  //for ( int i = 0; i < size; i++ )          
    for ( std::list<std::string>::const_iterator it = inputTokens.begin(); it !=  inputTokens.end(); it++ )          
    {        
      // Read the token      
      //const std::string token = inputTokens.at(i);
      const std::string token = *it;
      if (isFunction( token ) ){
	const std::string f = token;          
	stack.push(f);
      }
      // If token is an operator          
      else if ( isOperator( token ) )    
	{                      
	  // While there is an operator token, o2, at the top of the stack AND      
	  // either o1 is left-associative AND its precedence is equal to that of o2,      
	  // OR o1 has precedence less than that of o2,      
	  const std::string o1 = token;          
	  
	  if ( !stack.empty() )      
	    {      
	      std::string o2 = stack.top();      
	      
	      while ( isOperator( o2 ) &&          
		      ( ( isAssociative( o1, LEFT_ASSOC ) &&  cmpPrecedence( o1, o2 ) == 0 ) ||      
			( cmpPrecedence( o1, o2 ) < 0 ) ) )      
		{          
		  // pop o2 off the stack, onto the output queue;                          
		  stack.pop();
		  out.push_back( o2 );
		  
		  if ( !stack.empty() )       
		    o2 = stack.top();        
		  else      
		    break;      
		}           
	    }      
	  
	  // push o1 onto the stack.       
	  stack.push( o1 );          
	}          
      // If the token is a left parenthesis, then push it onto the stack.      
      else if ( token == "(" )          
	{          
	  // Push token to top of the stack        
	  stack.push( token );          
	}          
      // If token is a right bracket ')'          
      else if ( token == ")" )          
	{        
	  // Until the token at the top of the stack is a left parenthesis,       
	  // pop operators off the stack onto the output queue.      
	  std::string topToken  = stack.top();    
	  
	  while ( topToken != "(" )          
	    {               
	      out.push_back(topToken );          
	      stack.pop();      

	      if ( stack.empty() ) break;    
	      topToken = stack.top();    
	    }

	  // Pop the left parenthesis from the stack, but not onto the output queue.                                
	  if ( !stack.empty() ) stack.pop();                               

	  // If the stack runs out without finding a left parenthesis,       
	  // then there are mismatched parentheses.                  
	  if ( topToken != "(" )      
	    {
	      return false;      
	    }
	  //Check for function
	  if (!stack.empty()){
	    topToken = stack.top();
	    if (isFunction(topToken)){
	      out.push_back(topToken);
	      stack.pop();
	    }
	  }

	}
      // If the token is a number, then add it to the output queue.        
      else          
	{
	  out.push_back( token );          
	}
    }        
  
  // While there are still operator tokens in the stack:      
  while ( !stack.empty() )          
    {        
      const std::string stackToken = stack.top();      
      
      // If the operator token on the top of the stack is a parenthesis,       
      // then there are mismatched parentheses.      
      if ( isParenthesis( stackToken )   )    
	{      
	  return false;      
	}      
      
      // Pop the operator onto the output queue./      
      out.push_back( stackToken );                        
      stack.pop();        
    }          
  
  strArray.assign( out.begin(), out.end() );    
  //strArray.insert(strArray.begin(), out.begin(), out.end());

  return success;      
}          


//double Expr::RPNtoDouble( std::vector<std::string> tokens )      
double Expr::RPNtoDouble(){
  vector<double> noinput;
  return RPNtoDouble(noinput);
}          
    
double Expr::RPNtoDouble( std::vector<double>& input )          
{

  //Speed things up a bit
  std::stack<double> st;
  
  // For each token          
  for ( int i = 0; i < (int) m_rpn.size(); ++i )          
    {
      const std::string token = m_rpn[ i ];         
      // If the token is a value push it onto the stack          
      if ( !isOperator(token) && !isFunction(token) )         
	{
	  double d = 0.0;
	  if ( !is_number(token) ){
	    int idx = VarIdx(token);
	    if (idx > -1){
	      d = input.at(idx);
	    }
	    else info()<<"Parsing Error with "<<token<<endl;
	  }
	  else d = strtod( token.c_str() , NULL ) ;
	  st.push(d);
	}
      else if ( isOperator(token) )   
	{
	  double result =  0.0;    
	  // Token is an operator: pop top two entries          
	  double d2 = st.top();        
	  st.pop();
	  
	  if ( !st.empty() )
	    {  
	      const double d1 = st.top();        
	      st.pop();
	      //Get the result
	      result = token == "+" ? d1 + d2 :
		token == "-" ? d1 - d2 :
		token == "*" ? d1 * d2 :
		token == "^" ? pow(d1, d2) :
		token == "&&" ? d1 * d2 :
		token == "||" ? max(1.0, d1 + d2) :
		token == ">" ? d1 > d2 :
		token == "<" ? d1 < d2 :
		token == ">=" ? d1 >= d2 :
		token == "<=" ? d1 <= d2 :
		token == "==" ? d1 == d2 :
		token == "!=" ? d1 != d2 :
		d1 / d2;
	    }  
	  else    
	    {    
	      if ( token == "-" )    
		result = d2 * -1;    
	      else     
		result = d2;
	    }    
	  
	  // Push result onto stack         
	  //std::ostringstream s;        
	  //s << result;  
	  st.push( result );          
	}
      else if (isFunction(token))
	{
	  double result = 0.0;
	  double d = st.top();
	  st.pop();
	  //const double d = strtod( val.c_str(), NULL);
	  if (token == "sqrt") result = sqrt(d);
	  else if (token == "abs") result = abs(d);
	  else if (token == "fabs") result =  fabs(d); 
	  else if (token == "sin") result =  sin(d);
	  else if (token == "cos") result =  cos(d);
	  else if (token == "tan") result =  tan(d);
	  else if (token == "sinh") result =  sinh(d);
	  else if (token == "cosh") result =  cosh(d);
	  else if (token == "tanh") result =  tanh(d);
	  else if (token == "log") result =  log(d);
	  else if (token == "log10") result =  log10(d);
	  else if (token == "exp") result = exp(d);
	  else if (token == "min") {
	    double d2 = st.top();
	    st.pop();
	    result = min(d, d2);
	  }
	  else if (token == "max") {
	    double d2 = st.top();
	    st.pop();
	    result = max(d, d2);
	  }
	  else if (token == "CB") {
	    double d2 = st.top();
	    st.pop();
	    double d3 = st.top();
	    st.pop();
	    double d4 = st.top();
	    st.pop();
	    double d5 = st.top();
	    st.pop();
	    double d6 = st.top();
	    st.pop();
	    double d7 = st.top();
	    st.pop();
	    double d8 = st.top();
	    st.pop();
	    result = CB(d8, d7, d6, d5, d4, d3, d2, d);
	  }
	  else if (token == "Voigtian") {
	    double d2 = st.top();
	    st.pop();
	    double d3 = st.top();
	    st.pop();
	    double d4 = st.top();
	    st.pop();
	    double d5 = st.top();
	    st.pop();
	    result = Voigtian(d5, d4, d3, d2, d);
	  }
	  else if (token == "Gaussian") {
	    double d2 = st.top();
	    st.pop();
	    double d3 = st.top();
	    st.pop();
	    double d4 = st.top();
	    st.pop();
	    result = Gaussian(d4, d3, d2, d);
	  }
	  else result = -1.0;
	  // Push result onto stack         
	  //std::ostringstream s;        
	  //s << result;        
	  st.push( result );          
	  
	}
    }                  
  return st.top();        
}          


/*
double Expr::RPNtoDouble( std::vector<double>& input )          
{          
  std::stack<std::string> st;
  
  // For each token          
  for ( int i = 0; i < (int) m_tokens.size(); ++i )          
    {
      const std::string token = m_tokens[ i ];        
      
      // If the token is a value push it onto the stack          
      if ( !isOperator(token) && !isFunction(token) )          
	{
	  st.push(token);
	}
      else if ( isOperator(token) )   
	{
	  double result =  0.0;    
	  // Token is an operator: pop top two entries          
	  const std::string val2 = st.top();        
	  st.pop();
	  double d2;
	  if ( !is_number(val2) ){
	    cout<<val2<<" is not a number"<<endl;
	    int idx = VarIdx(val2);
	    if (idx > -1){
	      d2 = input.at(idx);
	    }
	    else cout<<"Parsing Error"<<endl;
	    cout<<"Found "<<val2<<" at "<<idx<<" and returned "<<d2<<endl;

	  }
	  else d2 = strtod( val2.c_str() , NULL ) ;
	  
	  if ( !st.empty() )
	    {    
	      const std::string val1 = st.top();        
	      st.pop();
	      double d1;
	      if ( !is_number(val2) ){
		int idx = VarIdx(val1);
		if (idx > -1){
		  d1 = input.at(idx);
		}
		else cout <<"Parsing Error"<<endl;
	      }
	      else d1 = strtod( val1.c_str() , NULL ) ;
	      //const double d1 = strtod( val1.c_str(), NULL );
	      
	      //Get the result
	      result = token == "+" ? d1 + d2 :
		token == "-" ? d1 - d2 :
		token == "*" ? d1 * d2 :
		token == "^" ? pow(d1, d2) :
		token == "&" ? d1 * d2 :
		token == "|" ? max(1.0, d1 + d2) :
		token == ">" ? d1 > d2 :
		token == "<" ? d1 < d2 :
		d1 / d2;          
	    }  
	  else    
	    {    
	      if ( token == "-" )    
		result = d2 * -1;    
	      else     
		result = d2;
	    }    
	  
	  // Push result onto stack         
	  std::ostringstream s;        
	  s << result;        
	  st.push( s.str() );          
	}          

      else if (isFunction(token))
	{
	  double result = 0.0;
	  const std::string val = st.top();
	  st.pop();
	  const double d = strtod( val.c_str(), NULL);
	  result = token == "sqrt" ? sqrt(d) :
	    token == "abs" ? abs(d) : 
	    token == "sin" ? sin(d) :
	    token == "cos" ? cos(d) :
	    token == "tan" ? tan(d) : -1.0;
	  // Push result onto stack         
	  std::ostringstream s;        
	  s << result;        
	  st.push( s.str() );          
	  
	}
    }                  
  
  return strtod( st.top().c_str(), NULL );        
}          
*/
std::list<std::string>& Expr::getTokens(){
  return m_tokens;
}

int Expr::VarIdx(string var){
  int idx = -1 ;
  for ( unsigned int i = 0; i < m_varnames.size(); ++i){
    if (var == m_varnames.at(i)) idx = i;
  }
  return idx;

}

bool Expr::isInExpr(string var){
  return (VarIdx(var) > -1);
}

std::list<std::string> Expr::Tokenize( const std::string& expression )    
{    
  std::list<std::string> tokens;          
  std::string str = "";
  
  for ( int i = 0; i < (int) expression.length(); ++i )      
    {      
      const std::string token( 1, expression[ i ] );      
      
      //If operator or paranthesis
      if ( isOperator( token ) || isParenthesis( token ) )      
	{  
	  if ( !str.empty() )  
	    {
	      tokens.push_back( str ) ;
	      if (!is_number(str) && !isFunction(str) && !isInExpr(str) ) m_varnames.push_back(str);
	    }
	  str = "";
	  tokens.push_back( token );
	}
      else if (i == ((int)expression.length() -1) && !token.empty() && token != " " && token !=","){
	str.append(token);
	if (!str.empty() && str != " ")
	  {
	    tokens.push_back(str);
	    if (!is_number(str) && !isFunction(str) && !isInExpr(str)) m_varnames.push_back(str);
	    str = "";
	  }
      }
      //If number, variable or function
      else      
	{         
	  // Append the numbers      
	  if ( !token.empty() && token != " " && token!= "," )      
	    {                     
	      str.append( token );

	      //Cases where the operator has more than 1 character
	      if (isOperator(str) ){
		tokens.push_back(str);
		str = "";
	      }
	      //Cases where it has more than one character and is not separated by a space
	      if (str.length() >= 2 && isOperator( str.substr( str.length() - 2 ) )){
		string substr = str.substr(0,str.length() - 2);
		tokens.push_back(substr);
		if (!is_number(substr) && !isFunction(substr) && !isInExpr(substr)) m_varnames.push_back(substr);
		tokens.push_back(str.substr(str.length() - 2));
		str = "";
	      }
	      
	    }     
	  //For the last variable/number where token is empty
	  else       
	    {      
	      if ( str != "" )      
		{      
		  tokens.push_back( str );  
		  if (!is_number(str) &&!isFunction(str) &&!isInExpr(str) ) m_varnames.push_back( str );
		  str = "";      
		}      
	    }                             
	}      
    }          

  
  //deal with start token being a minus sign
  std::string firstToken = tokens.front();
  if (firstToken == "-")
    {
      std::list<string>::iterator it = tokens.begin();
      it++;
      if ( it == tokens.end() ){
	return tokens;
      }
      std::string nextToken = *(it);
      
      if (is_number( nextToken)){
	  tokens.pop_front();
	  tokens.front() = firstToken + nextToken;
      }
      else if (nextToken == "(" || isFunction( nextToken) || isOperator(nextToken) ){
	tokens.front() = firstToken + "1";
	tokens.insert(it, "*");
      }
      else if (nextToken == "-" && firstToken == "-"){
	std::list<std::string>::iterator nit = it;
	std::advance(nit, -1);
	tokens.erase(it);
	tokens.erase(nit);

      }
    }

  typedef std::list<std::string>::iterator t_iter;
  std::string prevToken = "";
  for (t_iter it = tokens.begin(); it != boost::prior(tokens.end()); it++)
    {
      std::string token = *it;
      t_iter nit = it;
      std::advance(nit, 1);

      if ( nit == tokens.end() ){
	break;
      }

      std::string ntoken = *nit;
      
      if (token == "-" && prevToken == "(")
	{
	  if (is_number(ntoken)){
	    tokens.erase(nit);
	    *it = "-"+ntoken;
	    token = *it;
	  }
	}
      else if ( token == "-" &&
		(isOperator (prevToken) || isFunction(prevToken)))
	{
	  if (token == "-" && prevToken == "-")
	    {
	      t_iter nit = it;
	      t_iter nnit = nit;
	      nnit++;
	      std::advance(nit, -1);
	      tokens.erase(it);
	      *nit = "+";

	      t_iter pnit = nit;
	      std::advance(pnit, -1);
	      if (isOperator(*pnit) || *pnit == "(")
		{
		  tokens.erase(nit);
		}

	      token = *nnit;
	      it = nnit;

	      if ( it == boost::prior(tokens.end()) )
		{
		  break;
		}

	    }
	  else if (is_number(ntoken) || isFunction(ntoken) || isOperator(ntoken) )
	    {
	      bool exit = false;
	      if (nit == boost::prior(tokens.end()) )
		{
		  exit = true;
		}
	      tokens.erase(nit);
	      *it = "-" + ntoken;
	      token = *it;
	      if (exit) break;
	    }
	  else if ( ntoken == "(")
	    {
	      *it = "-1";
	      token = *it;
	      tokens.insert(nit, "*");
	    }
	}
      prevToken = token;
    }
  
  prevToken = "";
  t_iter prevIt;

  for (t_iter it = tokens.begin(); it != tokens.end(); it++){
    std::string token = *it;
    
    if (token == "(" && prevToken == "-")
      {
	tokens.insert(it, "1");
	tokens.insert(it, "*");
      }
    prevToken = token;
    prevIt = it;
  }
  
  return tokens;    
}  

double Expr::CB(double x, double N, double a, double n, double m, double s, double sl, double c){
  double PDF = 0.0;
  if (((x) - m)/s > -a){
    PDF = exp((-pow(((x)-m),2))/(2*pow(s,2)));
  }
  else{
    double A = pow(n/fabs(a),n) * exp(-pow(a,2)/2);
    double B = (n/fabs(a)) - fabs(a);
    PDF = A * pow(B-(((x)-m)/s),-n);
  }
  PDF = PDF * N;
  PDF = PDF + sl * x + c;
  
  return PDF;
    
}

double Expr::Voigtian(double x, double N, double m, double s, double a){
  return N*TMath::Voigt(x - m, s, a, 4);
}
double Expr::Gaussian(double x, double N, double m, double s){
  return N*TMath::Gaus(x, m, s);
}

//python files
#ifdef WITHPYTHON


boost::python::list Expr::Tokenize_py( const std::string& expression ){
  list<string> tokens = Tokenize(expression);
  boost::python::list s;
  for (list<string>::iterator is = tokens.begin(); is != tokens.end(); ++is){
    s.append((*is));
  }
  return s;
}
boost::python::list Expr::getTokens_py( ){
  list<string> tokens = getTokens();
  boost::python::list s;
  for (list<string>::iterator is = tokens.begin(); is != tokens.end(); ++is){
    s.append((*is));
  }
  return s;
}

vector<string> Expr::getRPN(){
  return m_rpn;
}

boost::python::list Expr::getRPN_py( ){
  vector<string> tokens = getRPN();
  boost::python::list s;
  for (vector<string>::iterator is = tokens.begin(); is != tokens.end(); ++is){
    s.append((*is));
  }
  return s;
}



double Expr::GetVal_py(boost::python::list& input){
  vector<double> dbl_vec;
  dbl_vec.reserve(len(input));
  for (unsigned int i = 0; i < len(input); ++i){
    double d = boost::python::extract<double>(input[i]);
    dbl_vec.push_back(d);
  }
  return GetVal(dbl_vec);
}

double Expr::GetVal3_py(){
  return GetVal();
}
boost::python::list Expr::GetVarNames_py(){ 
  boost::python::list l;
    for (vector<string>::iterator is = m_varnames.begin(); is != m_varnames.end(); ++is){
    l.append((*is));
  }

  return l;
}

string Expr::GetExpr_py() {return m_varexp;}


boost::python::list Expr::infixToRPN_py(boost::python::list& inputTokens){
  std::vector<string> output;
  std::list<std::string> input;
  for (unsigned int i = 0; i < len(inputTokens); ++i){
    string s = boost::python::extract<string>(inputTokens[i]);
    input.push_back(s);
  }

  infixToRPN(input, output);
  boost::python::list outputlist;
  for (vector<string>::iterator is = output.begin(); is != output.end(); ++is){
    outputlist.append(*is);
  }
  return outputlist;

}

#endif
