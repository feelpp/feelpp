/** @file exam_hashmap.cpp
 *
 *  Regression tests for the exhashmap<> container. */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "ginac.h"
using namespace GiNaC;

#include <iostream>
using namespace std;

unsigned exam_hashmap()
{
	unsigned result = 0;
	unsigned N = 100;

	cout << "examining hash maps" << flush;

	// Create empty container
	exhashmap<unsigned> M1;
	if (M1.size() != 0) {
		clog << "Newly constructed container has size() != 0" << endl;
		++result;
	}
	if (!M1.empty()) {
		clog << "Newly constructed container is not empty" << endl;
		++result;
	}

	cout << '.' << flush;

	// Insert elements
	for (unsigned i = 0; i < N; ++i)
		M1.insert(make_pair(i, i));

	if (M1.size() != N) {
		clog << "After " << N << " insertions, size() returns " << M1.size() << " instead of " << N << endl;
		++result;
	}

	for (unsigned i = 0; i < N; ++i) {
		if (M1[i] != i) {
			clog << "Lookup of key " << i << " in M1 didn't return correct value" << endl;
			++result;
			break;
		}
	}

	if (M1.size() != N) {
		clog << "After " << N << " lookups, size() returns " << M1.size() << " instead of " << N << endl;
		++result;
	}

	cout << '.' << flush;

	// Test constructor from two iterators and operator==
	exhashmap<unsigned> M2(M1.begin(), M1.end());
	if (M2.size() != M1.size()) {
		clog << "Constructor from two iterators: size of destination container (" << M2.size() << ") not equal to size of source (" << M1.size() << ")" << endl;
		++result;
	}

	for (unsigned i = 0; i < N; ++i) {
		if (M2[i] != i) {
			clog << "Lookup of key " << i << " in M2 didn't return correct value" << endl;
			++result;
			break;
		}
	}

	if (M2 != M1) {
		clog << "Copied container not equal to source" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test assignment operator
	exhashmap<unsigned> M3;
	M3 = M1;
	if (M3.size() != N) {
		clog << "Assignment operator: size of assigned container not equal to size of original" << endl;
		++result;
	}

	for (unsigned i = 0; i < N; ++i) {
		if (M3[i] != i) {
			clog << "Lookup of key " << i << " in M3 didn't return correct value" << endl;
			++result;
			break;
		}
	}

	cout << '.' << flush;

	// Test insert(it, it)
	exhashmap<unsigned> M4;
	M4.insert(M1.begin(), M1.end());

	if (M4.size() != M1.size()) {
		clog << "insert(it, it): size of destination container not equal to size of source" << endl;
		++result;
	}

	for (unsigned i = 0; i < N; ++i) {
		if (M4[i] != i) {
			clog << "Lookup of key " << i << " in M4 didn't return correct value" << endl;
			++result;
			break;
		}
	}

	cout << '.' << flush;

	// Test insert()/find()
	symbol x("x"), y("y");
	exhashmap<unsigned> M5;
	M5.insert(make_pair(x-2, 1));
	M5.insert(make_pair(sin(x+y), 2));
	M5.insert(make_pair(Pi, 3));
	M5.insert(make_pair(0, 4));
	M5.insert(make_pair(4*pow(x, y), 5));
	if (M5.size() != 5) {
		clog << "After 5 insertions, size() returns " << M5.size() << " instead of 5" << endl;
		++result;
	}

	exhashmap<unsigned>::const_iterator cit = M5.find(sin(x+y));
	if (cit == M5.end()) {
		clog << "Lookup of sin(x+y) didn't find anything" << endl;
		++result;
	}
	if (!cit->first.is_equal(sin(x+y))) {
		clog << "Lookup of sin(x+y) returned an incorrect iterator" << endl;
		++result;
	}
	if (cit->second != 2) {
		clog << "Lookup of sin(x+y) returned wrong value" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test re-inserting insert()
	pair<exhashmap<unsigned>::iterator, bool> pit = M5.insert(make_pair(sin(x+y), 42));
	if (pit.second) {
		clog << "Reinsertion of sin(x+y) inserted a new value" << endl;
		++result;
	}
	if (!pit.first->first.is_equal(sin(x+y))) {
		clog << "Reinsertion of sin(x+y) returned an incorrect iterator" << endl;
		++result;
	}
	if (pit.first->second != 2) {
		clog << "Reinsertion of sin(x+y) changed the value" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test operator[]
	unsigned v = M5[sin(x+y)];
	if (M5.size() != 5) {
		clog << "operator[] with an existing key changed the container size" << endl;
		++result;
	}
	if (v != 2) {
		clog << "operator[] with an existing key returned the wrong value" << endl;
		++result;
	}

	v = M5[y+1];
	if (M5.size() != 6) {
		clog << "operator[] with a new key didn't insert a new value" << endl;
		++result;
	}
	if (v != 0) {
		clog << "operator[] with a new key returned the wrong value" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test erase()
	exhashmap<unsigned>::iterator it = M5.find(y+1);
	if (it == M5.end()) {
		clog << "Key y+1 wasn't found" << endl;
		++result;
	}
	if (!it->first.is_equal(y+1)) {
		clog << "find() returned an incorrect iterator" << endl;
		++result;
	}
	if (it->second != 0) {
		clog << "find() returned an incorrect value" << endl;
		++result;
	}

	M5.erase(it);
	if (M5.size() != 5) {
		clog << "erase(it) didn't reduce the size of the container" << endl;
		++result;
	}

	it = M5.find(y+1);
	if (it != M5.end()) {
		clog << "Key was still found after erase()" << endl;
		++result;
	}

	exhashmap<unsigned>::size_type n = M5.erase(Pi);
	if (n != 1) {
		clog << "erase(Pi) returned " << n << " instead of 1" << endl;
		++result;
	}
	if (M5.size() != 4) {
		clog << "erase(Pi) didn't reduce the size of the container" << endl;
		++result;
	}

	n = M5.erase(42);
	if (n != 0) {
		clog << "erase(42) returned " << n << " instead of 0" << endl;
		++result;
	}
	if (M5.size() != 4) {
		clog << "erase(42) reduced the size of the container" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test swap()
	exhashmap<unsigned> M6;
	M6.swap(M1);
	if (M6.size() != N) {
		clog << "After swap, size() returns " << M6.size() << " instead of " << N << endl;
		++result;
	}
	if (M1.size() != 0) {
		clog << "After swap with empty container, size() returns " << M1.size() << " instead of 0" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test clear()
	M2.clear();
	if (M2.size() != 0) {
		clog << "Size of cleared container is " << M5.size() << " instead of 0" << endl;
		++result;
	}

	cout << '.' << flush;

	// Test count()
	n = M5.count(Pi);
	if (n != 0) {
		clog << "count(Pi) returns " << n << " instead of 0" << endl;
		++result;
	}

	n = M5.count(4*pow(x, y));
	if (n != 1) {
		clog << "count(4*x^y) returns " << n << " instead of 1" << endl;
		++result;
	}
	cout << '.' << flush;

	return result;
}

int main(int argc, char** argv)
{
	return exam_hashmap();
}
