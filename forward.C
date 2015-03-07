#include <iostream>
#include <cmath>
using namespace std;

int f1(int x)	;

int f2(int x) {
	return f1(x) + 5;
}

int f1(int x) {

	return x + 5;
}

int main() {

	int x = 5;
	cout << f2(x) << endl;

	return 0 ;
}