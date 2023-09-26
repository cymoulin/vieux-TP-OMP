// Hello,
//
// I haven't found the time to finish this TP,
// and if I had I think I would have spent this time for the current MPI TP.
// I probably spent too much time in the beginning 
// hesitating between eclipse and geany and switching from C to C++.
// So, I am sorry you have only code to read and not much comments.
// 
// Cyrille Moulin


#include <iostream>
using namespace std;

#include "omp.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>


// EXERCICES PRELIMINAIRES


void HelloWorld(){
	//int num_threads, rank;
	#pragma omp parallel //private (rank)
	{
		int num_threads = omp_get_num_threads();
		int rank = omp_get_thread_num();
		int max = omp_get_max_threads();

		cout << "Hello from thread" << rank << " out of " << num_threads;
		cout << "max = "<<  max;
	}
	}

void QuiFaitQuoi(){
	int tab[100];
	int count[4];

	for (int i = 0; i<4; i++) { count[i] = 0;}

	#pragma omp parallel
	{
		int rank = omp_get_thread_num();
		#pragma omp for schedule(dynamic, 10)
		for (int i = 0; i<100; i++){
			tab[i] = rank;
			}
	}

	for (int i = 0; i<100; i++){
			cout << tab[i] ;
			count[tab[i]] = count[tab[i]] + 1;
			}
	cout << endl;

	for (int i = 0; i<4; i++) {
		cout << i <<" --> " << count[i];
		}
	}



// EXERCICES A RENDRE

// Exercice 1

#define SIZE 300

void fillMatWithRand(double mat[SIZE][SIZE], int row, int col) {
	int i, j;
	for (i=0; i<SIZE; i++) {
		for (j=0; j<SIZE; j++) {
				mat[i][j] = rand();
			}
		}
	}

void multiplyMat(){
	double M1[SIZE][SIZE];
	double M2[SIZE][SIZE];
	double M3[SIZE][SIZE];

	fillMatWithRand(M1, SIZE, SIZE);
	fillMatWithRand(M2, SIZE, SIZE);

	double t0 = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic, 100) // slower with chunksize given
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			for (int k=0; k<SIZE; k++) {
				M3[i][j] += M1[i][k]*M2[k][j];
				}
			}
		}
	double t1 = omp_get_wtime();
	cout << "Time taken: " << t1-t0;

	}


// Exercice 2

int is_prime(int n) {
	int i;
	for (i = 2; i < n; i++) {
		if (n%i == 0) return 0;
		}
	return 1;
	}

int num_of_primes(int n){
	int count=0;
	#pragma omp parallel for schedule(dynamic) //reduction(+:count)
	for (int i = 2; i <= n; i++) {
		if (is_prime(i)==1) {
			//#pragma omp atomic
			count += 1;
			}
		}
	return count;
	}

void exercice2(int n) {
	//int n = 100000;
	double t0 = omp_get_wtime();
	cout << "The number of primes inferior to " << n << " is " << num_of_primes(n) << endl;
	double t1 = omp_get_wtime();
	cout << "Time taken: " << t1-t0;
	}

// Exercice 3

void histo(int N) {
	double tab[1000000];
	int histo[10];
	int count = 0;

	double t0 = omp_get_wtime();

	#pragma omp parallel for schedule(dynamic) // schedule(static)
	for (int i=0; i<N; i++) {
		tab[i] = 10* ((double) rand())/RAND_MAX; // to get floats between 0 and 10.
		}

	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
	#pragma omp atomic
		histo[(int) floor(tab[i])]++;
		}

	double t1 = omp_get_wtime();
	cout << "Time taken: " << t1-t0 << endl;

	for (int i = 0; i < 10; i++) {
			cout << i << " --> " << histo[i] << endl;
			count += histo[i];
			}
	cout << endl;
	cout << "count = " << count << endl;
	}

// Exercice 4:  Tri par transposition impaire-paire.

void exercice4() {
	int n;
	do {
	cout << "Entrez la taille paire du tableau: ";
	cin >> n;
	cout << endl;
	} while (n % 2 != 0);


	// création d'un tableau de 1000000 entiers au hasard.
	int tab[1000000];
	for (int i=0; i<n; i++) {
		tab[i] = rand();
		}
	// tri
	double t0 = omp_get_wtime();
	for (int phase = 0; phase < n; phase += 2) {
		// phase paire

		//omp_set_num_threads(2);
		#pragma omp parallel for
		for (int i = 0; i < n; i +=2) {
			if (tab[i] > tab[i+1]) {
				int temp = tab[i];
				tab[i] = tab[i+1];
				tab[i+1] = temp;
			}
		}
		// phase impaire
		#pragma omp parallel for
		for (int i = 1; i < n-1; i +=2) {
			if (tab[i] > tab[i+1]) {
				int temp = tab[i];
				tab[i] = tab[i+1];
				tab[i+1] = temp;
			}
		}
	}
	double t1 = omp_get_wtime();
	cout << "Time taken: " << t1-t0 << endl;

	bool sorted = true;
	int i = 0;
	while (sorted && (i < n-1)) {
		sorted = sorted && (tab[i] <= tab[i+1]);
		i++;
	}
	if (sorted) { cout << "Trié";} else {
		cout << "erreur de tri!";}
}



int main() {
	cout << "Number of procs = " << omp_get_num_procs() << endl ;
	multiplyMat();
	//exercice4();
	//histo(1000000);

	return 0;
	}








