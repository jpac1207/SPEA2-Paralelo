#include <iostream>
#include<stdlib.h>
#include<stdio.h>
//#include<cstdlib>
#include<time.h>
#include <sstream>

#include"mpi.h"

#include"Individual.h"
#include"Population.h"
#include"PopulationReader.h"
#include"HyperVolumeCalculator.h"
#include"SPEA2.h"
#include"TestUtil.h"
#include <time.h>


using namespace std;

void single() {
	
	for (int i = 0; i < 31; i++) {
		clock_t tStart = clock();
		SPEA2* spea2 = new SPEA2();
		spea2->runSingle();
		cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;
		delete spea2;
	}
	std::getchar();
}

int main(int argc, char* argv[])
{
	/*
		Implementado por João Pedro Augusto Costa
		Instituto Federal do Maranhão
		Campus São Luís, Monte Castelo
		Departamento Acadêmico de Informatica
	*/

	/*TestUtil t;
	t.run();
	getchar();*/


	srand((unsigned int)time(NULL));

	int rank = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < 31; i++) {
		clock_t tStart = clock();
		SPEA2* spea2 = new SPEA2(rank);
		spea2->run();
		if (rank == 0) {
			cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;
		}
		delete spea2;
	}

	MPI_Finalize();
	if (rank == 0) {
		cout << "FIM" << endl;
		std::getchar();
	}
	
	/*single();*/
	return 0;
}


