#include <mpi.h>
#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <string>
//function approximation, least squares

/*[INPUT FILE GUIDE]
0 row - degree of approximating polynom
1 st row - lower nd upper bound
2 nd - number of dots
the rest - dots
*/

using namespace std;

const string filename = "Input.txt";
//function for printing the elements in a list
void showlist(list <pair <double, double> > _list) {
	list < pair <double, double> > ::iterator it;
	for (it = _list.begin(); it != _list.end(); ++it)
		cout << '\t' << "x : " << (*it).first << " y : " << (*it).second << endl;
	cout << endl;
}
// get lower and upper bounds
void getNums(double &lbound, double &ubound, string tmp) {
	string tmp1 = "";
	bool flag = false;
	for (int i = 0; i < tmp.length(); i++) {
		if (47 < tmp[i] < 58) {
			tmp1 += tmp[i];
		}

		if (tmp[i] == ' '  && flag == false) {
			lbound = stoi(tmp1);
			tmp1 = "";
			flag = true;
		}
		if (flag == true && i == tmp.length() - 1) {
			ubound = stoi(tmp1);
		}

	}
}

//read input file
void readInput(list< pair<double, double> > &dots, int &dotNums, double &lbound, double &ubound, int &deg) {

	ifstream fin(filename, ios_base::in);
	string tmp;
	double a, b;
	if (fin.is_open()) {
		cout << "success" << endl;
		getline(fin, tmp);
		deg = stoi(tmp);
		getline(fin, tmp);
		getNums(lbound, ubound, tmp);
		getline(fin, tmp);
		dotNums = 0;
		while (!fin.eof()) { 
			fin >> a;
			fin >> b;
			dots.push_back(make_pair(a,b));
			dotNums++;
		}
		fin.close();
	}
	else {
		cout << "error!" << endl;
	}
}



int main(int argc, char **argv)
{


	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		list< pair<double, double> > dots;
		int dotNums, deg;
		double lbound, ubound;
		readInput(dots, dotNums, lbound, ubound, deg);
		cout << "deg = " << deg << endl;
		cout << "lbound = " << lbound << " ubound = " << ubound << endl;
		cout << "dotNums = " << dotNums << endl;
		showlist(dots);
	}








	MPI_Finalize();
	return 0;
}