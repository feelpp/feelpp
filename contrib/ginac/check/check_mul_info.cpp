#include "ginac.h"
#include <iostream>
using namespace GiNaC;

int main(int argc, char** argv)
{
	symbol x("x"), y("y");
	ex e = x*y;
	if (!e.info(info_flags::indefinite)) {
		std::cerr << "eek, product of two symbols is NOT indefinite";
		return 1;
	}
	return 0;
}
