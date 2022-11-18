#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <Math.h>
#include <stdio.h>

using namespace std;


typedef double(*TestFunctTemp2)(double*&, int&);

template <typename T>
void fillArr(T*& array, int& N) {
    for (int i = 0; i < N; i++)
        array[i] = pow(i, 2.0 / 3.0) * sin(i / 2.0);
}
template <typename T>
void printArr(T*& array, int& N) {
	for (int i = 0; i < N; i++)
		cout << array[i] << endl;
}

template <typename T>
void compare_exchange(T& array1, T& array2) {
	T temp;
	if (array1 > array2) {
		temp = array2;
		array2 = array1;
		array1 = temp;
	}
}
template <typename T>
double bubbleSort(T*& array, int& N) {
	double time_start = omp_get_wtime();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N - i - 1; j++) {
			compare_exchange(array[j], array[j + 1]);
        }
    }

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

template <typename T>
double oddEvenSorting(T*& array, int& N) {
	double time_start = omp_get_wtime();

	T temp;
	int upper_bound;
	if (N % 2 == 0)
		upper_bound = N / 2 - 1;
	else
		upper_bound = N / 2;

	for (int i = 0; i < N; i++) {
		if (i % 2 == 0)
			for (int j = 0; j < N / 2; j++)
				compare_exchange(array[2 * j], array[2 * j + 1]);
		else
			for (int j = 0; j < upper_bound; j++)
				compare_exchange(array[2 * j + 1], array[2 * j + 2]);
	}

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

template <typename T>
double oddEvenSortingParallel(T*& array, int& N) {
	double time_start = omp_get_wtime();
	T temp;
	int upper_bound;
	if (N % 2 == 0)
		upper_bound = N / 2 - 1;
	else
		upper_bound = N / 2;

	for (int i = 0; i < N; i++) {
		if (i % 2 == 0)
#pragma omp parallel for
			for (int j = 0; j < N / 2; j++)
				compare_exchange(array[2 * j], array[2 * j + 1]);
		else
#pragma omp parallel for
			for (int j = 0; j < upper_bound; j++)
				compare_exchange(array[2 * j + 1], array[2 * j + 2]);
	}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}
template <typename T>
void InsertSort(T* arr, int i, int length, int half) {
	T temp = 0;
	int j = 0;

	for (int f = half + i; f < length; f = f + half)
	{
		j = f;
		while (j > i && arr[j - half] > arr[j])
		{
			temp = arr[j];
			arr[j] = arr[j - half];
			arr[j - half] = temp;
			j = j - half;
		}
	}
}
template <typename T>
double shellSort(T*& array, int& N) {
	double time_start = omp_get_wtime();

	int i, j, step;
	T tmp;
	for (step = N / 2; step > 0; step /= 2)
		for (i = 0; i < step; i++)
		{
			InsertSort(array, i, N, step);
		}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

template <typename T>
double shellSortParallel(T*& array, int& N) {
	double time_start = omp_get_wtime();

	int i, j, step;
	T tmp;
	for (step = N / 2; step > 0; step /= 2)
		#pragma omp parallel for shared( array, N, step, i)  default(none)
			for (i = 0; i < step; i++)
			{
				InsertSort(array, i, N, step);
			}
	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

template <typename T>
double quickSort(T* array, int N) {
	double time_start = omp_get_wtime();

	long i = 0, j = N;
	T temp, p;

	p = array[N >> 1];               

	do {
		while (array[i] < p) i++;
		while (array[j] > p) j--;

		if (i <= j) {
			temp = array[i]; array[i] = array[j]; array[j] = temp;
			i++; j--;
		}
	} while (i <= j);

	if (j > 0) quickSort(array, j);
	if (N > i) quickSort(array + i, N - i);

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}
template <typename T>
double quickSortParallel(T* array, int N) {
	double time_start = omp_get_wtime();

	long i = 0, j = N;
	T temp, p;

	p = array[N >> 1];

	do {
		while (array[i] < p) i++;
		while (array[j] > p) j--;

		if (i <= j) {
			temp = array[i]; array[i] = array[j]; array[j] = temp;
			i++; j--;
		}
	} while (i <= j);

#pragma omp task shared(array)
	if (j > 0) quickSortParallel(array, j);
#pragma omp task shared(array)
	if (i < N) quickSortParallel(array + i, N - i);
#pragma omp taskwait

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

template <typename T>
double testBubbleSort(T*& array, int& N) {
	return bubbleSort(array, N);
}
template <typename T>
double testOddEvenSorting(T*& array, int& N) {
	return oddEvenSorting(array, N);
}
template <typename T>
double testShellSort(T*& array, int& N) {
	return shellSort(array, N);
}
template <typename T>
double testQuickSort(T*& array, int& N) {
	return quickSort(array, N);
}

template <typename T>
double testOddEvenSortingParallel(T*& array, int& N) {
	return oddEvenSortingParallel(array, N);
}
template <typename T>
double testShellSortParallel(T*& array, int& N) {
	return shellSortParallel(array, N);
}
template <typename T>
double testQuickSortParallel(T*& array, int& N) {
	return quickSortParallel(array, N);
}


double AvgTrustedInterval(double& avg, double*& times, int& cnt)
{
	double sd = 0, newAVg = 0;
	int newCnt = 0;
	for (int i = 0; i < cnt; i++)
	{
		sd += (times[i] - avg) * (times[i] - avg);
	}
	sd /= (cnt - 1.0);
	sd = sqrt(sd);
	for (int i = 0; i < cnt; i++)
	{
		if (avg - sd <= times[i] && times[i] <= avg + sd)
		{
			newAVg += times[i];
			newCnt++;
		}
	}
	if (newCnt == 0) newCnt = 1;
	return newAVg / newCnt;
}

template <typename T>
double TestIter(void* Funct, T* array, int size, int iterations)
{
	typedef double(*TestFunctTempl)(T*&, int&);

	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];
	cout << endl;

	for (int i = 0; i < iterations; i++)
	{
		T* arr = new T[size];
		memcpy(arr, array, sizeof(int) * size);

		curtime = ((*(TestFunctTempl)Funct)(arr, size)) * 1000;

		Times[i] = curtime;
		avgTime += curtime;
		cout << "+";
		//printArr(arr, size);
		delete[] arr;
	}
	cout << endl;

	avgTime /= iterations;
	cout << "AvgTime:" << avgTime << endl;

	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "AvgTimeTrusted:" << avgTimeT << endl;

	return avgTimeT;
}

template <typename T>
void test_functions(void** Functions, string(&function_names)[7])
{
	int iters = 2;
	int nd = 0;
	double times[4][7][3];
	T* arr;
	for (int size = 20000; size <= 50000; size += 10000)
	{
		for (int threads = 1; threads <= 4; threads++)
		{
			omp_set_num_threads(threads);
			//перебор алгоритмов по условиям
			for (int alg = 0; alg <= 6; alg++)
			{
				arr = new T[size];
				fillArr(arr, size);
				if (threads == 1)
				{
					if (alg == 0 || alg == 1 || alg == 3 || alg == 5) {
						times[nd][alg][0] = TestIter<T>(Functions[alg], arr, size, iters);
						// iters - кол-во запусков алгоритма
						times[nd][alg][1] = times[nd][alg][0];
						times[nd][alg][2] = times[nd][alg][0];
						delete[] arr;
					}
				}
				else
				{
					if (alg != 0 && alg != 1 && alg != 3 && alg != 5)
					{
						times[nd][alg][threads - 2] = TestIter<T>(Functions[alg], arr, size, iters);
						delete[] arr;
					}
				}
			}
		}
		nd++;
	}
	ofstream fout("output.txt");
	fout.imbue(locale("Russian"));
	for (int ND = 0; ND < 4; ND++)
	{
		switch (ND)
		{
		case 0:
			cout << "\n----------20000 элементов массива----------" << endl;
			break;
		case 1:
			cout << "\n----------30000 элементов массива----------" << endl;
			break;
		case 2:
			cout << "\n----------40000 элементов массива----------" << endl;
			break;
		case 3:
			cout << "\n----------50000 элементов массива----------" << endl;
			break;
		default:
			break;
		}
		for (int alg = 0; alg <= 6; alg++)
		{
			for (int threads = 1; threads <= 4; threads++)
			{
			cout << "Поток " << threads << " --------------" << endl;
				if (threads == 1)
				{
					if (alg == 0 || alg == 1 || alg == 3 || alg == 5) {
						cout << function_names[alg] << "\t" << times[ND][alg][0] << " ms." << endl;
						fout << times[ND][alg][0] << endl;
					}
				}
				else
				{
					if (alg != 0 && alg != 1 && alg != 3 && alg != 5)
					{
						cout << function_names[alg] << "\t" << times[ND][alg][threads - 2] << " ms." << endl;
						fout << times[ND][alg][threads - 2] << endl;
					}
				}
			}
		}
	}
	fout.close();
}

int main() {
	setlocale(LC_ALL, "RUS");

	void** FunctionsINT = new void* [7]{ testBubbleSort<int>, testOddEvenSorting<int>, testOddEvenSortingParallel<int>, 
		testShellSort<int>, testShellSortParallel<int>, testQuickSort<int>, testQuickSortParallel<int> };
	void** FunctionsDOUBLE = new void* [7]{ testBubbleSort<double>, testOddEvenSorting<double>, testOddEvenSortingParallel<double>,
		testShellSort<double>, testShellSortParallel<double>, testQuickSort<double>, testQuickSortParallel<double> };
	string function_names[7]{ "Bubble sort(посл.)", "OddEven sort(посл.)", "OddEven sort(парал.)", "Shell sort(посл.)", "Shell sort(парал.)",
							"Quick sort(посл.)", "Quick sort(парал.)" };
	/*cout << "---------------------------ТИП INT---------------------------" << endl;
	test_functions<int>(FunctionsINT, function_names);*/
	cout << "---------------------------ТИП DOUBLE---------------------------" << endl;
	test_functions<double> (FunctionsDOUBLE, function_names);

	return 0;
}