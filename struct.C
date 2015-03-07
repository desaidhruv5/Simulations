#include <iostream>
#include <cmath>
using namespace std;

struct Ad {

	int adnumber;	//#of ads shown to users
	float percent;	//percent of users who clicked on ads
	double earned;	//amt earned on avg for each ad clicked
};

void profitperday(Ad d) {
	cout << "Total ads shown: " << d.adnumber << endl;
	cout << "Percentage of total users who clicked on ads: " << d.percent << endl;
	cout << "Amount of money earned on average, per ad: " << d.earned << endl;
	cout << "You made $" << d.adnumber * d.percent * d.earned * .01 << "." << endl;
}

int main() {
	cout << "Let's find out how much you made today!" << endl;
	Ad d1;
	cout << "Enter the total number of ads: ";
	cin >> d1.adnumber;
	cout << "Enter the percentage of total users who clicked on ads: ";
	cin >> d1.percent;
	cout << "Enter the amount of money earned on average, per ad: ";
	cin >> d1.earned;
	profitperday(d1);
	return 0;
}