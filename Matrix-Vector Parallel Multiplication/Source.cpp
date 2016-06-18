#include <iostream>
#include <windows.h>
#include <ppl.h>
#include <random>
#include <time.h>
#include <fstream>
#include <omp.h>
#include <stdio.h>

using namespace std;
using namespace concurrency;



template <class Generator>
double** initialize_matrix(double** m, size_t size, Generator& gen);
template <class Generator>
double* initialize_array(double* m, size_t size, Generator& gen);

void mat_vec_parallel_mul(
	double** m1, double* m2, double* result, size_t size);
void mat_vec_sequental(
	double** m1, double* m2, double* result, size_t size);

void equal_result(double* r1, double* r2, double* r3, size_t size);

void MATVEC(double** MATRIX, double* VECTOR, double* RESULT, int SIZE_I, int SIZE_J)
{
	INT I, J;

	#pragma OMP PARALLEL SHARED(MATRIX, RESULT, VECTOR) PRIVATE(I, J)
	{
		#pragma omp for schedule(static)
		for(I = 0; I<SIZE_I; I++){
			RESULT[I] = 0.;
			for(J = 0; J<SIZE_J; J++){
				RESULT[I] += ((MATRIX[I][J])*(VECTOR[J]));
			}
		}
	}
}


int main(){

	omp_set_num_threads(4);
	cout << omp_get_num_threads()<< endl;
	ofstream ofs("results.txt");
//	ofs << "[";
	for (size_t s = 350; s <= 6300; s += 350){
		
	//	const size_t size = 2500;
		mt19937 gen(42);

		double** m1 = new double*[s];
		for (size_t i = 0; i < s; ++i)
		{
			m1[i] = new double[s];
		}


		double* m2 = new double[s];

		m1 = initialize_matrix(m1, s, gen);
		m2 = initialize_array(m2, s, gen);

		double* r1 = new double[s];
		double* r2 = new double[s];
		double* r3 = new double[s];

		time_t start;
		time_t end;

		start = clock();
		mat_vec_parallel_mul(m1, m2, r1, s);
		end = clock();

		ofs << "SIZE = " << s << endl;

		long diff = end - start;
		ofs << "PARALLEL PPL: " << diff << endl;

		start = clock();
		mat_vec_sequental(m1, m2, r2, s);
		end = clock();

		diff = end - start;
		ofs << "SEQUENTAL : " << diff << endl;

		start = clock();
		MATVEC(m1, m2, r3, s,s);
		end = clock();

		diff = end - start;
		ofs << "PARALLEL OMP : " << diff << endl;

		equal_result(r1, r2,r3, s);
	

		ofs << endl << endl;

		delete[] m2;
		delete[] m1;
		delete[] r1;
		delete[] r2;
		delete[] r3;
	}
	 /*
	 ofs << size << endl;
	  for (int i = 0; i < size; i++){
		 for (int j = 0; j < size; j++){
			 ofs << m1[i][j] << "\t";
		 }
		 ofs << endl;
	 }
	  
	 ofs << endl;
	 ofs << endl;
	 for (int i = 0; i < size; i++){
		
			 ofs << m2[i] << "\t";
	
	 }


	 ofs << endl;
	 ofs << endl;
	 for (int i = 0; i < size; i++){

		 ofs << r1[i] << "\t";

	 }
	 */

	cout << omp_get_num_threads() << endl;
	return 0;
}



void equal_result(double* r1, double* r2, double* r3, size_t size){

	for (int i = 0; i < size; i++){
		if (r1[i] != r2[i] || r1[i]  != r3[i] || r2[i] != r3[i]){
			cout << "FAIL ON "<< size << endl;
			break;
		}
	}

	cout << "SUCCESS ON " <<size << endl;
}

template <class Generator>
double** initialize_matrix(double** m, size_t size, Generator& gen)
{
	for (size_t i = 0; i < size; ++i)
	{
		for (size_t j = 0; j < size; ++j)
		{
			m[i][j] = static_cast<double>(gen());
		}
	}
	return m;
}

template <class Generator>
double* initialize_array(double* m, size_t size, Generator& gen)
{
	
		for (size_t j = 0; j < size; ++j)
		{
			m[j] = static_cast<double>(gen());
		}

	return m;
}

void mat_vec_parallel_mul(
	double** m1, double* m2, double* result, size_t size)
{
	parallel_for(size_t(0), size, [&](size_t i)
	{
			result[i] = 0;
			for (int k = 0; k < size; k++)
			{
			//	temp += m1[i][k] * m2[k][j];
				result[i] += m1[i][k] * m2[k];

			}
	});
}


void mat_vec_sequental(
	double** m1, double* m2, double* result, size_t size)
{
	for (int i = 0; i < size; i++)
	{
			result[i] = 0;
			for (int k = 0; k < size; k++)
			{
				result[i] += m1[i][k] * m2[k];
			}
	}
}