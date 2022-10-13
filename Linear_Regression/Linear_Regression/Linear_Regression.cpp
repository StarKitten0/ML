#include <iostream>
#include <fstream>
using namespace std;

int** X;
int** Y;
double** koef;

int N, M;

void ReadFile(string name)
{
	ifstream in(name);
	in >> M;
	in >> N;
	
	 //не забыть про delete [] p_darr;
	X = new int* [M];
	Y = new int*[1];
	Y[0] = new int[M];
	for (int i=0;i < M;i++)
	{
		X[i] = new int[N];
		for (int j = 0;j < N;j++)
		{
			in >> X[i][j];
		}
	}
	for (int i = 0;i < M;i++)
	{
		in >> Y[0][i];
	}
	
}

int** Transpose(int M,int N,int** a)
{
	int** array = new int*[N];
	for (int i = 0;i < N;i++)
	{
		array[i] = new int[M];
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
			array[j][i] = a[i][j];
	}
	return array;
}

int** MultipMatrix(int row1,int col1,int row2,int col2,int** a, int** b)
{
	int** c = new int* [row1];
	for (int i = 0; i < row1; i++)
	{
		c[i] = new int[col2];
		for (int j = 0; j < col2; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < col1; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return c;
}

double** MultipDoubMatrix(int row1, int col1, int row2, int col2, double** a, int** b)
{
	double** c = new double* [row1];
	for (int i = 0; i < row1; i++)
	{
		c[i] = new double[col2];
		for (int j = 0; j < col2; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < col1; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return c;
}


double** MultipDoubDoubMatrix(int row1, int col1, int row2, int col2, double** a, double** b)
{
	double** c = new double* [row1];
	for (int i = 0; i < row1; i++)
	{
		c[i] = new double[col2];
		for (int j = 0; j < col2; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < col1; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return c;
}

double** ConvertToDouble(int N,int M,int** a)
{
	double** c = new double*[N];
	for (int i = 0; i < N;i++)
	{
		c[i] = new double[M];
		for (int j = 0;j < M;j++)
		{
			c[i][j] = (double)a[i][j];
		}
	}
	return c;
}

double** inversion(double** A, int N)
{
	double temp;

	double** E = new double* [N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;

	return A;
}

void PrintMatrix(int M,int N,int** a)
{
	for (int i = 0;i < M;i++)
	{
		for (int j = 0;j < N;j++)
			cout << a[i][j] << " ";

		cout << endl;
	}

	cout << endl;
}
void PrintDoubMatrix(int M, int N, double** a)
{
	for (int i = 0;i < M;i++)
	{
		for (int j = 0;j < N;j++)
			cout << a[i][j] << " ";

		cout << endl;
	}
	cout << endl;
}

void fit()
{
	
	

	int** XT = Transpose(M, N, X);

	//X^T * X

	int** XTX = MultipMatrix(N, M, M, N, XT, X);
	//находим обратную для всего этого матрицу

	double** inv = inversion(ConvertToDouble(N, N, XTX), N);

	//сноваумножаем на транспонированную X

	double** invXT = MultipDoubMatrix(N, N, N, M, inv, XT);
	koef = MultipDoubMatrix(N, M, M, 1, invXT, Transpose(1,M,Y));
// koef и есть матрица нужных нам величин

	/*PrintMatrix(M, N, X);
	PrintMatrix(N, M, XT);
	PrintMatrix(N, N, XTX);
	PrintDoubMatrix(N, N, inv);*/

	
}

void rmse(double** pred, int** Y)
{
	double sum = 0;
	for (int i = 0;i < M;i++)
	{
		sum += (Y[i][0] - pred[i][0]) * (Y[i][0] - pred[i][0]);
	}

	cout << "rmse = " << sum;
}

void predict()
{
	ReadFile("test.txt");
	double** predict = MultipDoubDoubMatrix(M, N, N, 1, ConvertToDouble(M,N,X), koef);

	//PrintDoubMatrix(M, 1, predict);

	rmse(predict, Transpose(1,M,Y));
}



int main()
{
	ReadFile("train.txt");
	fit();
	predict();
	
	

}

