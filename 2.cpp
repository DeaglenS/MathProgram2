#include "PolStr.h"
#include <iostream>
#include <windows.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>
#include <string>
using namespace std;

void MatrixProcessor(int m, int n, double* array, double* b);
double Determinant(int n, double* array);
double* SolveSLAAWithDecompositionMethod(int n, double* array, double* b);
double* InvertibleMatrix(int n, double* array);
double* SolveSLAAWithSimpleIterationMethod(int n, double* array, double* b);
#include <fstream>

int main()
{
	int m;
	long int n;

	FILE* file;
	fopen_s(&file, "input.txt", "r");
	if (file != 0)
	{
		fscanf_s(file, "%d", &m);
		fscanf_s(file, "%d", &n);

		double x;
		double* array;
		double* b = new double[n];

		if (m == 1 || m == 4)
		{

			array = new double[n * (n + 1)];
			for (int k = 0; k < n; k++)
			{
				for (int i = 0; i < n; i++)
				{
					fscanf_s(file, "%lf", &x);
					array[k * n + i] = x;
				}
				fscanf_s(file, "%lf", &x);
				b[k] = x;
			}
		}
		else
		{

			array = new double[n * n];
			for (int k = 0; k < n * n; k++)
			{
				fscanf_s(file, "%lf", &x);
				array[k] = x;
			}
		}





		cout << m << endl;
		cout << n << endl;



		for (int k = 0; k < n; k++)
		{
			for (int i = 0; i < n; i++)
			{
				cout << array[k * n + i] << " ";
			}
			if (m == 1 || m == 4)
			{
				cout << b[k] << " ";
			}
			cout << endl;
		}


		cout << endl << endl;
		MatrixProcessor(m, n, array, b);


	}


}

/// <summary>
/// 1. Решение СЛАУ методом Декомпозиции (точный метод).
///	2. Поиск определителя матрицы.
///	3. Поиск обратной матрицы.
/// 4. Решение СЛАУ методом Простой итерации (итерационный метод).
/// </summary>
/// <param name="m">Тип задачи</param>
/// <param name="n">Порядок матрицы</param>
/// <param name="array">Матрица</param>
/// <param name="b">Вектор свободных коэфицентов при m=1</param>
void MatrixProcessor(int m, int n, double* array, double* b = NULL)
{
	double* newArray1 = new double[n];
	double answer;
	double* newArray3 = new double[n * n];
	double* newArray4 = new double[n];
	switch (m)
	{
	case 1:

		newArray1 = SolveSLAAWithDecompositionMethod(n, array, b);
		cout << "Solution is ";
		for (int k = 0; k < n; k++)
		{
			cout << newArray1[k] << " ";
		}
		cout << endl;
		return;
		break;
	case 2:
		answer = Determinant(n, array);
		cout << "Determinant is " << answer << endl;
		cout << endl;
		return;
		break;
	case 3:
		newArray3 = InvertibleMatrix(n, array);
		for (int k = 0; k < n; k++)
		{
			cout << newArray3[k] << " ";
		}
		return;
		break;
	case 4:
		newArray4 = SolveSLAAWithSimpleIterationMethod(n, array, b);
		for (int k = 0; k < n; k++)
		{
			cout << newArray4[k] << " ";
		}
		return;
		break;
	default:
		cout << "Ошибка размерности переменной, отвечающей за тип задачи";
		return;
		break;
	}
}
/// <summary>
/// Решение СЛАУ(Система линейно-алгебраических принадлежностей) методом декомпозиции(метод Халецкого).
/// </summary>
/// <param name="n">Порядок матрицы</param>
/// <param name="array">Матрица</param>
/// <param name="b">Вектор свободных коэфицентов</param>
double* SolveSLAAWithDecompositionMethod(int n, double* array, double* b)
{
	/*
	1
	5
	122 33 - 73 1 316 3
	1 1 - 24 3 236 5
	5 52 41 2 136 6
	2 32 3 3 3 4
	1 - 4 6 5 4 1
	*/
	double* B = new double[n * n];
	double* C = new double[n * n];
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			B[k * n + i] = 0;
			if (k == i)
			{
				C[k * n + i] = 1;
			}
			else
			{
				C[k * n + i] = 0;
			}
		}
	}


	for (int k = 0; k < n; k++)
	{
		for (int i = k; i < n; i++) // для каждой строчки этого столбца 
		{
			double minus = 0;
			for (int z = 0; z < k; z++)
			{
				minus += B[i * n + z] * C[z * n + k];
				//cout << "!!" << B[i * n + z] << "!*!" << C[z * n + k] << "!!" << endl;
			}
			B[i * n + k] = array[i * n + k] - minus;
			//	cout << "B" << i << k << " = " << array[i * n + k] << " - " << minus << " = " << B[i * n + k] << endl;
		}


		for (int r = k; r < n; r++) // для каждого столбца строчки
		{
			double minus = 0;
			for (int z = 0; z < k; z++)
			{
				minus += B[k * n + z] * C[z * n + r];
			}
			C[k * n + r] = 1 / B[k * n + k] * (array[k * n + r] - minus);
			//	cout << "C" << k << r << " = " << 1 / B[k * n + k] << " * " << " ( " << array[k * n + r] << " - " << minus << " ) " << " = " << C[k * n + r] << endl;
		}
	}


	////////////
	cout << "B = " << endl;
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			cout << B[k * n + i] << " ";
		}
		cout << " ";

		cout << endl;
	}

	cout << "C = " << endl;
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			cout << C[k * n + i] << " ";
		}
		cout << " ";

		cout << endl;
	}





	double* Y = new double[n];
	double* X = new double[n];


	for (int k = 0; k < n; k++)
	{
		X[k] = 0;
		Y[k] = 0;
	}


	Y[0] = b[0] / B[0];

	for (int i = 1; i < n; i++)
	{
		double minus = 0;
		for (int k = 0; k < i; k++)
		{
			minus += B[i * n + k] * Y[k];
		}
		Y[i] = 1 / B[i * n + i] * (b[i] - minus);




	}

	X[n - 1] = Y[n - 1];

	for (int i = n - 1; i >= 0; i--)
	{

		double minus = 0;
		for (int k = i + 1; k < n; k++)
		{
			minus += C[i * n + k] * X[k];
		}
		X[i] = Y[i] - minus;
	}









	cout << "Y = ";
	for (int k = 0; k < n; k++)
	{
		cout << Y[k] << " ";
	}
	cout << endl;





	/*
	for (int k = 0; k < n; k++)
	{
		cout << Y[k] << " ";
	}
	cout << endl << endl << endl << endl;
	for (int k = 0; k < n; k++)
	{
		cout << X[k] << " ";
	}*/

	/*
	double* B = new double[n * n];


	for (int j = 1; j <= n; j++)
	{

		for (int i = j; j <= n; i++)
		{
			double MinusForB = 0;

			for (int k = 1; k <= j - 1; k++)
			{
				MinusForB += B[i * n + k] *
			}


			B[i * n + j] = array[i * n + j] - MinusForB;
		}
	}*/



	FILE* file;

	fopen_s(&file, "output.txt", "wb");
	if (file != 0)
	{

		for (int k = 0; k < n; k++)
		{
			fprintf_s(file, "%G ", X[k]);
		}

		fprintf_s(file, "\n");



		double* Xold = new double[n];
		double* e = new double[n];
		for (int k = 0; k < n; k++)
		{
			Xold[k] = X[k];
		}

		for (int k = 0; k < n; k++)
		{
			double sum = 0;
			for (int l = 0; l < n; l++)
			{
				sum += array[k * n + l] * Xold[l];
			}
			e[k] = sum;
		}

		for (int k = 0; k < n; k++)
		{
			e[k] = e[k] - b[k];
		}
		for (int k = 0; k < n; k++)
		{
			fprintf_s(file, "%G \n", e[k]);
		}
		//e = Ax - b;

		fprintf_s(file, "\n");
		double eAbs = 0;
		for (int k = 0; k < n; k++)
		{
			eAbs += e[k] * e[k];
		}
		eAbs = sqrt(eAbs);
		fprintf_s(file, "%G \n", eAbs);
	}
	fclose(file);

	return X;
}

double Determinant(int n, double* array)
{
	/*
	2
	5
	1 1 2 3 4
	3 5 6 4 0
	4 34 7 63 -8
	6 -8 5 3 1
	1 5 4 3 2
	*/
	double* B = new double[n * n];
	double* C = new double[n * n];
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			B[k * n + i] = 0;
			if (k == i)
			{
				C[k * n + i] = 1;
			}
			else
			{
				C[k * n + i] = 0;
			}
		}
	}


	for (int k = 0; k < n; k++)
	{
		for (int i = k; i < n; i++) // для каждой строчки этого столбца 
		{
			double minus = 0;
			for (int z = 0; z < k; z++)
			{
				minus += B[i * n + z] * C[z * n + k];
				//cout << "!!" << B[i * n + z] << "!*!" << C[z * n + k] << "!!" << endl;
			}
			B[i * n + k] = array[i * n + k] - minus;
			//	cout << "B" << i << k << " = " << array[i * n + k] << " - " << minus << " = " << B[i * n + k] << endl;
		}


		for (int r = k; r < n; r++) // для каждого столбца строчки
		{
			double minus = 0;
			for (int z = 0; z < k; z++)
			{
				minus += B[k * n + z] * C[z * n + r];
			}
			C[k * n + r] = 1 / B[k * n + k] * (array[k * n + r] - minus);
			//	cout << "C" << k << r << " = " << 1 / B[k * n + k] << " * " << " ( " << array[k * n + r] << " - " << minus << " ) " << " = " << C[k * n + r] << endl;
		}
	}



	double answer = 1;

	for (int i = 0; i < n; i++)
	{
		answer *= B[i * n + i];
	}



	FILE* file;

	fopen_s(&file, "output.txt", "wb");
	if (file != 0)
	{
		fprintf_s(file, "%G", answer);
	}
	fclose(file);
	return answer;
}


double* InvertibleMatrix(int n, double* array)
{
	double* answer = new double[n * n];
	double* E = new double[n * n];
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			if (k == i)
			{
				E[k * n + i] = 1;
			}
			else
			{
				E[k * n + i] = 0;
			}
		}
	}

	/*
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			cout << E[k * n + i] << " ";
		}
		cout << endl;
	}
	*/


	double* b = new double[n];

	double* tmp = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n; k++)
		{
			b[k] = E[k * n + i];
		}

		tmp = SolveSLAAWithDecompositionMethod(n, array, b);

		for (int k = 0; k < n; k++)
		{
			answer[k * n + i] = tmp[k];
		}
	}

	FILE* file;

	fopen_s(&file, "output.txt", "wb");
	if (file != 0)
	{

		for (int k = 0; k < n; k++)
		{
			for (int i = 0; i < n; i++)
			{
				fprintf_s(file, "%G ", answer[k * n + i]);
			}
			fprintf_s(file, "\n");
		}

		fprintf_s(file, "\n");
		double* e = new double[n * n];
		for (int k = 0; k < n; k++)
		{

			for (int i = 0; i < n; i++)
			{
				e[k * n + i] = 0;
				for (int j = 0; j < n; j++)
				{
					e[k * n + i] += array[k * n + j] * answer[j * n + i];
				}//ax
			}
		}

		for (int k = 0; k < n; k++)
		{
			for (int i = 0; i < n; i++)
			{
				e[k * n + i] -= E[k * n + i];
				fprintf_s(file, "%G ", e[k * n + i]);
			}
			fprintf_s(file, "\n");
		}

		double eAbs = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				eAbs += e[i * n + j] * e[i * n + j];
			}
		}
		eAbs = sqrt(eAbs);

		fprintf_s(file, "%G ", eAbs);
	}
	fclose(file);
	return answer;
}

double* SolveSLAAWithSimpleIterationMethod(int n, double* array, double* b)
{
	/*
	4
	3
	5 -1 3 5
	1 -4 2 20
	2 -1 5 10
	*/
	double eps = 0.0001;
	double* a = new double[n * n];
	double* B = new double[n];
	double* answer = new double[n];
	double* X = new double[n];

	for (int i = 0; i < n; i++)
	{
		B[i] = b[i] / array[i * n + i];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				a[i * n + j] = 0;
			}
			else
			{
				a[i * n + j] = -array[i * n + j] / array[i * n + i];
			}
		}
	}

	cout << "B = ";
	for (int i = 0; i < n; i++)
	{
		cout << B[i] << " ";
	}

	cout << endl;
	cout << "a = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << a[i * n + j] << " ";
		}
		cout << endl;
	}

	//double* X = new double[n];
	for (int k = 0; k < n; k++)
	{
		X[k] = B[k];
	}

	double* Xeps = new double[n];
	double* Xbefore = new double[n];
	for (int k = 0; k < n; k++)
	{
		Xbefore[k] = 0;
	}
	while (true)
	{


		/*for (int f = 0; f < n; f++)
		{
			sum += a[k * n + f] * X[f];
		}*/
		double* Xold = new double[n];
		for (int k = 0; k < n; k++)
		{
			Xold[k] = X[k];
		}
		for (int k = 0; k < n; k++)
		{
			double sum = 0;
			for (int l = 0; l < n; l++)
			{
				sum += a[k * n + l] * Xold[l];
			}
			X[k] = sum;
		}


		for (int r = 0; r < n; r++)
		{
			X[r] += B[r];
		}

		//X = B + a * X;


		/*
		cout << endl << endl;
		for (int k = 0; k < n; k++)
		{
			cout << X[k] << " ";
		}
		*/


		// проверка условия выхода



		bool stop = true;
		for (int k = 0; k < n; k++)
		{
			double eps0 = abs(X[k] - Xbefore[k]) / abs(Xbefore[k]);
			if (eps0 > eps)
			{
				stop = false;
			}
		}
		if (stop == true)
		{
			break;
		}
		/*
		double* Xleft = new double[n];
		for (int k = 0; k < n; k++)
		{
			Xleft[k] = X[k] - Xbefore[k];
		}


		if (Determinant(n, Xleft) < ((1 - Determinant(n, a)) / Determinant(n, a) * eps))
		{
			break;
		}
		*/
		for (int k = 0; k < n; k++)
		{
			Xbefore[k] = X[k];
		}


	}
	/*
	cout << endl << endl;
	for (int k = 0; k < n; k++)
	{
		cout << X[k] << " ";
	}*/
	/*

	for (int k = 0; k < n; k++)
	{
		double x;
		X = B;
		while ()
		{
			x = B + a * x;
			/// a --- уже получен
			// B ---
		}
		X[k] = x;
	}
	*/

	cout << endl;



	FILE* file;

	fopen_s(&file, "output.txt", "wb");
	if (file != 0)
	{

		for (int k = 0; k < n; k++)
		{
			fprintf_s(file, "%G ", X[k]);
		}

		fprintf_s(file, "\n");

		//double* e = new double[n];

		double* Xold = new double[n];
		double* e = new double[n];
		for (int k = 0; k < n; k++)
		{
			Xold[k] = X[k];
		}

		for (int k = 0; k < n; k++)
		{
			double sum = 0;
			for (int l = 0; l < n; l++)
			{
				sum += array[k * n + l] * Xold[l];
			}
			e[k] = sum;
		}

		for (int k = 0; k < n; k++)
		{
			e[k] = e[k] - b[k];
		}
		for (int k = 0; k < n; k++)
		{
			fprintf_s(file, "%G \n", e[k]);
		}
		//e = Ax - b;

		fprintf_s(file, "\n");
		double eAbs = 0;
		for (int k = 0; k < n; k++)
		{
			eAbs += e[k] * e[k];
		}
		eAbs = sqrt(eAbs);
		fprintf_s(file, "%G \n", eAbs);
	}
	fclose(file);

	return X;
}