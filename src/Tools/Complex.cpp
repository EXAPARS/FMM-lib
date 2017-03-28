/*
  TODO : Copyright
*/

#include "Complex.hpp"
using namespace std;

ostream& operator << (ostream& out, const complex & c)  
{  
	//out << std::scientific;
	out << c.re << " " << c.im;
	return out;
}
