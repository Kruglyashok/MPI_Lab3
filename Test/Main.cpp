#include <mpi.h>
#include <iostream>
#include <fstream>
#include <list>
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

void readInput(list< pair<double, double> > &dots, int &dotNums, double &lbound, double &ubound, int &deg) {
	ifstream fin(filename);
	if (fin.is_open()) {
	
	}
	else {
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
	}











	MPI_Finalize();
	return 0;
}