/**
 * @file mystring.cpp Example of extending GiNaC: writing new classes
 */
#include <iostream>
#include <string>   
#include <stdexcept>
using namespace std;

#include <ginac/ginac.h>
using namespace GiNaC;

class mystring : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(mystring, basic)
public:
	mystring(const string &s);
	ex eval(int level) const;
private:
	string str;

protected:
	void do_print(const print_context &c, unsigned level = 0) const;

};

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(mystring, basic,
  print_func<print_context>(&mystring::do_print))

// ctors
mystring::mystring() { }
mystring::mystring(const string &s) :  str(s) { }

// comparison
int mystring::compare_same_type(const basic &other) const
{
	const mystring &o = static_cast<const mystring &>(other);
	int cmpval = str.compare(o.str);
	if (cmpval == 0)
		return 0;
	else if (cmpval < 0)
		return -1;
	else
		return 1;
}

// printing
void mystring::do_print(const print_context &c, unsigned level) const
{
	// print_context::s is a reference to an ostream
	c.s << '\"' << str << '\"';
}

/**
 * evaluation: all strings automatically converted to lowercase with
 * non-alphabetic characters stripped, and empty strings removed
 */
ex mystring::eval(int level) const
{
	string new_str;
	for (size_t i=0; i<str.length(); i++) {
		char c = str[i];
		if (c >= 'A' && c <= 'Z') 
			new_str += tolower(c);
		else if (c >= 'a' && c <= 'z')
			new_str += c;
	}

	if (new_str.length() == 0)
		return 0;
	else
		return mystring(new_str).hold();
}

int main(int argc, char** argv)
{
	ex e = mystring("Hello, world!");
	cout << is_a<mystring>(e) << endl;
	cout << ex_to<basic>(e).class_name() << endl;
	cout << e << endl;
	ex another = pow(mystring("One string"), 2*sin(Pi-mystring("Another string")));
	cout << another << endl;
	return 0;
}
