#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <Expr.h>
#include <boost/algorithm/string.hpp>

using namespace std;


const int LEFT_ASSOC  = 0;      
const int RIGHT_ASSOC = 1;   

// Map the different operators: +, -, *, / etc
typedef std::map< std::string, std::pair< int,int > > OpMap;
//typedef std::vector<std::string>::const_iterator cv_iter;  
//typedef std::string::iterator s_iter;

const OpMap::value_type assocs[] =     
  {  OpMap::value_type( "+", std::make_pair<int,int>( 0, LEFT_ASSOC ) ),      
     OpMap::value_type( "-", std::make_pair<int,int>( 0, LEFT_ASSOC ) ),
     OpMap::value_type( ">", std::make_pair<int,int>( 2, LEFT_ASSOC ) ),      
     OpMap::value_type( "<", std::make_pair<int,int>( 2, LEFT_ASSOC ) ),      
     OpMap::value_type( "==", std::make_pair<int,int>( 2, LEFT_ASSOC ) ),     
     OpMap::value_type( "!=", std::make_pair<int,int>( 2, LEFT_ASSOC ) ),
     OpMap::value_type( "*", std::make_pair<int,int>( 5, LEFT_ASSOC ) ),        
     OpMap::value_type( "&&", std::make_pair<int,int>( 1, LEFT_ASSOC ) ),       
     OpMap::value_type( "||", std::make_pair<int,int>( 1, LEFT_ASSOC ) ),
     OpMap::value_type( "/", std::make_pair<int,int>( 5, LEFT_ASSOC ) ),
     OpMap::value_type( "^", std::make_pair<int,int>( 6, RIGHT_ASSOC ) )};
     
const OpMap opmap( assocs, assocs + sizeof( assocs ) / sizeof( assocs[ 0 ] ) );    


Expr::Expr(string varexp){
  m_varexp = varexp;
  m_tokens = getExpressionTokens(varexp);

  bool c = infixToRPN(m_tokens, m_tokens.size(), m_rpn);
  if (!c) {
    cout << "Parse Error in " << varexp << endl;
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

/*
pair<vector<string>, vector<string> > Expr::GetBracketExp(string varexp){
  pair<vector<string>, vector<string> > output;
  vector<size_t> left;
  vector<size_t> right;
  vector<string> strings;
  vector<string> brkdelims;
  vector<pair<size_t, char> > brackets;
  //vector<pair<size_t, char> > delimiter;
  size_t pos = 0;
  size_t brpos_left  = varexp.find( "(" , 0  );
  size_t brpos_right = varexp.find( ")" , 0  );

  while (brpos_left != std::string::npos || brpos_right != std::string::npos){
    if (brpos_left < brpos_right){
      left.push_back(brpos_left);
      brackets.push_back(pair<size_t,char>(brpos_left,'('));
      pos = brpos_left;
    }
    else {
      right.push_back(brpos_right);
      brackets.push_back(pair<size_t,char>(brpos_right,')'));
      pos = brpos_right;
    }
    brpos_left  = varexp.find( "(" , pos+1 , 1 );
    brpos_right = varexp.find( ")" , pos+1 , 1 );
  }

  //Nor we have got all the brackets


  if (left.size() == right.size()){
    int nleft = 0;
    int nright = 0;
    size_t l = 0, r = 0;
    for (std::vector<pair<size_t,char> >::iterator im = brackets.begin(); im != brackets.end(); ++im){
      if ((*im).second == '('){
	if (nleft == 0 || (r>0 && l < r && (*im).first > r)){
	  l = (*im).first;
	}
	nleft++;
      }
      if ((*im).second == ')'){
	nright++;
	if (nleft == nright){
	  r = (*im).first;
	  strings.push_back(varexp.substr(l+1, r-(l+1)));
	  //This closes the bracket, so we look in space before next open bracket

	  size_t nxt = varexp.find("(", r+1, 1);
	  string rst;
	  if (nxt != std::string::npos)	  rst = varexp.substr(r+1, nxt - (r+1));
	  else {
	    rst = varexp.substr(r+1, std::string::npos);
	    cout<<"string is "<<rst<<endl;
	  }
	  for (std::vector<string>::iterator id = m_delimiters.begin(); id !=m_delimiters.end(); ++id){
	    if (rst.find(*id) !=std::string::npos) {
	      brkdelims.push_back((*id));
	      rst.erase(rst.find(*id),(*id).size());
	      cout<<"Found "<<(*id)<<" so string is now "<<rst<<endl;
	    }
	  }
	  boost::algorithm::trim(rst);
	  cout<<"After trimming we have "<<rst<<endl;
	  if (rst.size() > 0 ) strings.push_back(rst);
	}
      }
    }
  }
  else cout<<"Brackets don't match"<<endl;
  output.first = strings;
  output.second = brkdelims;
  return output;
}

boost::python::list Expr::GetBracketExp_py(string varexp){
  pair<vector<string>, vector<string> > o = GetBracketExp(varexp);
  std::vector<string> strings = o.first;
  std::vector<string> chars = o.second;
  boost::python::list s;
  boost::python::list c;
  boost::python::list l;
  for (vector<string>::iterator is = strings.begin(); is != strings.end(); ++is){
    s.append((*is));
  }
  for (vector<string>::iterator ic = chars.begin(); ic != chars.end(); ++ic){
    c.append((*ic));
  }
  l.append(s);
  l.append(c);
  return l;
}


double Expr::GetVal(Tree* t){
  if (t == 0) {
    cout<<"No valid tree passed"<<endl;
    return -1.0;
  }
  double output;
  if (m_expressions.size() == 0){
    vector<double> algvals;
    for (std::vector<string>::iterator is = m_varnames.begin() ; is != m_varnames.end() ; ++is){
      double value = 0.0;
      if (is_number((*is))){
	value = atof((*is).c_str());
      }
      else{
	value = t->GetVal(*is);
      }
      
      algvals.push_back(value);
    }
    output = GetVal(algvals);
  }
  else{
    vector<double> algvals;
    for (std::vector<Expr*>::iterator ie = m_expressions.begin() ; ie != m_expressions.end() ; ++ie){
      double value = (*ie)->GetVal(t);
      algvals.push_back(value);
      cout<<"Pushing back "<<value<<endl;
    }
    output = GetVal(algvals);
  }
  
  return output;
}
*/

double Expr::GetVal(){
  if (m_varnames.size() == 0){
    return RPNtoDouble();
  }
  else{
    cout<<"Parse Error : Value(s) ";
    for (unsigned int i = 0 ; i < m_varnames.size(); ++i){
      cout<<m_varnames.at(i);
    }
    cout<<" not found"<<endl;
    return -1.0;
  }
}

double Expr::GetVal(std::vector<double>& input){
  if (input.size() != m_varnames.size() || input.size() == 0) return -1;
  if (m_rpn.size() == 1 && input.size() == 1) return input.at(0);
  return RPNtoDouble(input);
}

/*
double Expr::GetVal(std::vector<double>& input){
  if ((input.size() != m_varnames.size() && input.size() !=  m_expressions.size()) || input.size() == 0) return -1;
  double val = 1.0;
  double eval = 0.0;

  if (input.size() == m_expressions.size() || m_expressions.size() == 0){
    int j = 0;
    if (is_number(m_values.at(0))) val = atof(m_values.at(0).c_str());
    else {
      val = input[0];
      j++;
    }
    for (unsigned int i = 0; i < m_operations.size(); ++i){
      if (is_number(m_values.at(i+1))) {
	eval = atof(m_values.at(i+1).c_str());
      }
      else {
	eval = input.at(j);
	j++;
      }

      if (m_operations.at(i) == "+")  val = val + eval;
      if (m_operations.at(i) == "-")  val = val - eval;
      if (m_operations.at(i) == "/")  val = val/eval;
      if (m_operations.at(i) == "*")  val = val*eval;
      if (m_operations.at(i) == "==") val = (val == eval ? 1.0 : 0.0);
      if (m_operations.at(i) == ">")  val = (val > eval ? 1.0 : 0.0);
      if (m_operations.at(i) == "<")  val = (val < eval ? 1.0 : 0.0);
      if (m_operations.at(i) == ">=") val = (val >= eval ? 1.0 : 0.0);
      if (m_operations.at(i) == "<=") val = (val <= eval ? 1.0 : 0.0);
      if (m_operations.at(i) == "&&") val = val * eval;
      if (m_operations.at(i) == "||") val = ((val == 1 || eval == 1) ? 1.0 : 0.0);
    }
  }
  else if (input.size() == m_varnames.size() && m_expressions.size() > 0){
    int s = 0;
    vector<double> vals;
    int n = m_expressions[0]->GetVarNames().size();
    vals.insert(vals.end(), input.begin() + s, input.begin() + s + n); 
    val = m_expressions[0]->GetVal(vals);
    s = s + n;
    vals.clear();
    double output;
    for (unsigned int i = 0; i < m_operations.size(); ++i){
      n = m_expressions[i]->GetVarNames().size();
      vals.insert(vals.end(), input.begin() + s, input.begin() + s + n); 
      output = m_expressions[i]->GetVal(vals);
      if (m_operations.at(i) == "+") val = val + output;
      if (m_operations.at(i) == "-") val = val - output;
      if (m_operations.at(i) == "/") val = val/output;
      if (m_operations.at(i) == "*") val = val*output;
      if (m_operations.at(i) == "==") val = (val == output ? 1.0 : 0.0);
      if (m_operations.at(i) == ">") val = (val > output ? 1.0 : 0.0);
      if (m_operations.at(i) == "<") val = (val < output ? 1.0 : 0.0);
      if (m_operations.at(i) == ">=") val = (val >= output ? 1.0 : 0.0);
      if (m_operations.at(i) == "<=") val = (val <= output ? 1.0 : 0.0);
      if (m_operations.at(i) == "&&") val = val * output;
      if (m_operations.at(i) == "||") val = ((val == 1 || output == 1) ? 1.0 : 0.0);
      s = s+n;
      vals.clear();
    }
  }
  else cerr<<"Error Calculating Value - Check inputs"<<endl;
  return val;
}
*/
double Expr::GetVal_py(boost::python::list& input){
  vector<double> dbl_vec;
  for (unsigned int i = 0; i < len(input); ++i){
    double d = boost::python::extract<double>(input[i]);
    dbl_vec.push_back(d);
  }
  return GetVal(dbl_vec);
}

double Expr::GetVal3_py(){
  return GetVal();
}
/*
double Expr::GetVal2_py(PyObject* tree){
  Tree* t = boost::python::extract<Tree*>(tree);
  return GetVal(t);
}

void Expr::SetTree(Tree* t){
  for (vector<string>::iterator iv = m_varnames.begin(); iv != m_varnames.end(); ++iv){
    t->SetBranch(*iv);
  }

}

void Expr::SetTree_py(PyObject* t){
  Tree* tree = boost::python::extract<Tree*>(t);
  SetTree(tree);
}
*/
std::vector<string>& Expr::GetVarNames(){ return m_varnames;}
boost::python::list Expr::GetVarNames_py(){ 
  boost::python::list l;
    for (vector<string>::iterator is = m_varnames.begin(); is != m_varnames.end(); ++is){
    l.append((*is));
  }

  return l;
}
/*
std::vector<Expr*> Expr::GetExpressions(){ return m_expressions;}
boost::python::list Expr::GetExpressions_py(){ 
  boost::python::list l;
    for (vector<Expr*>::iterator is = m_expressions.begin(); is != m_expressions.end(); ++is){
      //boost::python::object o(boost::python::handle<>(boost::python::borrowed(*is)));
      //boost::python::reference_existing_object::apply<Expr*>::type converter;
      //PyObject* o = converter( *is );
      Expr* e = (*is);
      l.append(e);
    }
    return l;
}
*/
string Expr::GetExpr() const{ return m_varexp;}
string Expr::GetExpr() {return m_varexp;}
string Expr::GetExpr_py() {return m_varexp;}


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

string Expr::FindFunc(string name){
  size_t pos = m_varexp.find(name);
  //brpos_left  = varexp.find( "(" , pos+1 , 1 );
  //brpos_right = varexp.find( ")" , pos+1 , 1 );

  if (pos != std::string::npos){
    cout<<m_varexp[pos + name.size()]<<endl;
  }
  return name;

}

size_t Expr::FindNextBracket(size_t pos){
  size_t brpos_left  = m_varexp.find( "(" , pos+1 , 1 );
  size_t brpos_right = m_varexp.find( ")" , pos+1 , 1 );
  
  int nleft = 1;
  int nright = 0;

  while ((brpos_left != std::string::npos || brpos_right != std::string::npos) && nleft > nright){
    if (brpos_left < brpos_right){
      nleft ++;
      pos = brpos_left;
    }
    else {
      nright++;
      pos = brpos_right;
    }
    brpos_left  = m_varexp.find( "(" , pos+1 , 1 );
    brpos_right = m_varexp.find( ")" , pos+1 , 1 );
  }
  if (nleft != nright) {
    cout<<"Problem with brackets"<<endl;
    return std::string::npos;
  }
  else return pos;

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
  token == "|"  || token == "==" || token == "!=";
}        

bool Expr::isFunction( const std::string& token)          
{
  return token == "sqrt" || token == "abs" ||     
    token == "sin" || token == "cos" ||
    token == "tan" || token == "log" || token == "log10" ||
    token == "min" || token == "max";
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

boost::python::list Expr::infixToRPN_py(boost::python::list& inputTokens,   const int& size){
  vector<string> output;
  std::vector<std::string> input;
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

}

// Convert infix expression format into reverse Polish notation          
bool Expr::infixToRPN( const std::vector<std::string>& inputTokens,     
		       const int& size,     
		       std::vector<std::string>& strArray )          
{       
  bool success = true;      
  
  //const int LEFT_ASSOC  = 0;      
  //const int RIGHT_ASSOC = 1;   
  std::list<std::string> out;
  std::stack<std::string> stack;
  
  // While there are tokens to be read      
  for ( int i = 0; i < size; i++ )          
    {        
      // Read the token      
      const std::string token = inputTokens[ i ];        
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
	    else cout<<"Parsing Error with "<<token<<endl;
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
	  else if (token == "abs") result =  abs(d); 
	  else if (token == "sin") result =  sin(d);
	  else if (token == "cos") result =  cos(d);
	  else if (token == "tan") result =  tan(d);
	  else if (token == "log") result =  log(d);
	  else if (token == "log10") result =  log10(d);
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
std::vector<std::string>& Expr::getExpressionTokens(){
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

std::vector<std::string> Expr::getExpressionTokens( const std::string& expression )    
{    
  std::vector<std::string> tokens;          
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
  
  return tokens;    
}  

boost::python::list Expr::getExpressionTokens_py( const std::string& expression ){
  vector<string> tokens = getExpressionTokens(expression);
  boost::python::list s;
  for (vector<string>::iterator is = tokens.begin(); is != tokens.end(); ++is){
    s.append((*is));
  }
  return s;
}
boost::python::list Expr::getExpressionTokens2_py( ){
  vector<string> tokens = getExpressionTokens();
  boost::python::list s;
  for (vector<string>::iterator is = tokens.begin(); is != tokens.end(); ++is){
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





/*void Expr::Print( const std::string& message,  
		  const InputIterator& itbegin,     
		  const InputIterator& itend,     
		  const std::string& delimiter)    
{    
  std::cout << message << std::endl;  
  
  std::copy(itbegin,     
	    itend,     
	    std::ostream_iterator<T>(std::cout, delimiter.c_str()));    
  
  std::cout << std::endl;  
  } 

double Expr::GetVal2_py(string s){
  return GetVal(s);
}

double Expr::GetVal(string s)          
{
  // Tokenize input expression          
  std::vector<std::string> tokens = getExpressionTokens( s );                             
  
  // Evaluate feasible expressions    
  std::vector<std::string> rpn;          
  if ( infixToRPN( tokens, tokens.size(), rpn ) )        
    {
      double d = RPNtoDouble( rpn );     
      std::cout << "Result = " << d << std::endl;        
    }        
  else        
    {        
      std::cout << "Mis-match in parentheses" << std::endl;        
    }           
  
  return 0;          
}    


double Expr::GetVal(vector<string>& tokens)          
{
  // Evaluate feasible expressions    
  std::vector<std::string> rpn;          
  double d = -1.0;
  if ( infixToRPN( tokens, tokens.size(), rpn ) )        
    {
      d = RPNtoDouble( rpn );        
      //std::cout << "Result = " << d << std::endl;        
    }        
  else        
    {        
      std::cout << "Mis-match in parentheses" << std::endl;        
    }           
  
  return d;
}    
*/
