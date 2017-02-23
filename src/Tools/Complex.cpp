/*
  TODO : Copyright
*/

#include "Complex.hpp"
using namespace std;

ostream& operator << (ostream& out, const complex & c)  
{  
	//out << std::scientific;
	out << c.re << endl;
	//out << " ";
	//if (c.im >= 0) out << "+";
	out << c.im;
	return out;
}
