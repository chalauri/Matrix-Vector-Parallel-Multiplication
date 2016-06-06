#include <iostream>
#include <windows.h>
#include <ppl.h>
#include <random>
#include <time.h>


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

void equal_result(double* r1, double* r2, size_t size);

int main(){

	const size_t size = 2500;
	mt19937 gen(42);
	
	double** m1 = new double*[size];
	for (size_t i = 0; i < size; ++i)
	{
		m1[i] = new double[size];
	}


	double* m2 = new double[size];
	
	 m1 = initialize_matrix(m1, size, gen);
	 m2 = initialize_array(m2, size, gen);

	 double* r1 = new double[size];
	 double* r2 = new double[size];

	 time_t start;
	 time_t end;

	 start = clock();
	 mat_vec_parallel_mul(m1, m2, r1, size);
	 end = clock();

	 long diff = end - start;
	 cout << "PARALLEL : " << diff << endl;

	 start = clock();
	 mat_vec_sequental(m1, m2, r2, size);
	 end = clock();

	 diff = end - start;
	 cout << "SEQUENTAL : " << diff << endl;

	 equal_result(r1, r2, size);
	 /*
	 cout << size << endl;
	  for (int i = 0; i < size; i++){
		 for (int j = 0; j < size; j++){
			 cout << m1[i][j] << "\t";
		 }
		 cout << endl;
	 }
	  
	 cout << endl;
	 cout << endl;
	 for (int i = 0; i < size; i++){
		
			 cout << m2[i] << "\t";
	
	 }


	 cout << endl;
	 cout << endl;
	 for (int i = 0; i < size; i++){

		 cout << r1[i] << "\t";

	 }
	 */
	return 0;
}



void equal_result(double* r1, double* r2, size_t size){

	for (int i = 0; i < size; i++){
		if (r1[i] != r2[i]){
			cout << "FAIL" << endl;
			break;
		}
	}

	cout << "SUCCESS" << endl;
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
		for (size_t j = 0; j < size; j++)
		{
			result[i] = 0;
			for (int k = 0; k < size; k++)
			{
			//	temp += m1[i][k] * m2[k][j];
				result[i] += m1[i][k] * m2[k];

			}
		
		}
	});
}


void mat_vec_sequental(
	double** m1, double* m2, double* result, size_t size)
{
	for (int i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			result[i] = 0;
			for (int k = 0; k < size; k++)
			{
				result[i] += m1[i][k] * m2[k];
			}

		}
	}
}