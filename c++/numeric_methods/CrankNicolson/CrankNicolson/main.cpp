#include <iostream>
#include <iomanip>
#include <random>
#include <math.h>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <vector>
#include <numeric>
#include <climits>
#include <cstring>

#include <omp.h>

using namespace std;

#define M_PI 3.14159265

class heat_task {
public:
	double T;
	double L;
	int n;
	int m;
	explicit heat_task(double _T, double _L, int _n, int _m) :
		T(_T), L(_L), n(_n), m(_m) { }
	double initial_condition(double x) { return sin(x * M_PI) + sin(3 * M_PI * x); }
	double left_condition(double t) { return 0; }
	double right_condition(double t) { return 0; }
	double f(double x, double t) { return 0; }
};

void cyclic_reduction(int size, double **A, double *F)
{
	int i, j, k;
	int index1, index2, offset;
	double alpha, gamma;
		
	for (i = 0; i < log2(size + 1) - 1; i++) {
		for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1)) {
			offset = pow(2, i);
			index1 = j - offset;
			index2 = j + offset;
			alpha = A[j][index1] / A[index1][index1];
			gamma = A[j][index2] / A[index2][index2];
			for (k = 0; k < size; k++) {
				A[j][k] -= (alpha*A[index1][k] + gamma * A[index2][k]);
			}
			F[j] -= (alpha*F[index1] + gamma * F[index2]);
		}
	}
}

void cyclic_reduction_back_substitution(int size, double **A, double *F, double *x)
{
	int index1, index2, offset;
	int index = (size - 1) / 2;
	int i, j, k;

	x[index] = F[index] / A[index][index];
	for (i = log2(size + 1) - 2; i >= 0; i--) {
		for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1)) {
			offset = pow(2, i);
			index1 = j - offset;
			index2 = j + offset;
			x[index1] = F[index1];
			x[index2] = F[index2];
			for (k = 0; k<size; k++) {
				if (k != index1)
					x[index1] -= A[index1][k] * x[k];
				if (k != index2)
					x[index2] -= A[index2][k] * x[k];
			}
			x[index1] = x[index1] / A[index1][index1];
			x[index2] = x[index2] / A[index2][index2];
		}
	}
}

void heat_equation_crank_nicolson(heat_task task, double *v)
{
	double h = task.L / task.n;
	double tao = task.T / task.m;

	double **G = new double*[task.m + 1];
	for (int i = 0; i <= task.m; i++) {
		G[i] = new double[task.n + 1];
	}

	for (int i = 0; i <= task.n; ++i) {
		G[0][i] = task.initial_condition(i * h);
	}
	for (int j = 1; j <= task.m; ++j) {
		G[j][0] = task.left_condition(j * tao);
		G[j][task.n] = task.right_condition(j * tao);
	}

	double ac_val;
	double right_side_1_val, right_side_2_val;
	double b_val;

	double *a;
	double *right;

	if (task.n == 1) {
		goto done;
	}

	ac_val = (tao / (2 * h * h));
	right_side_1_val = (tao / (h * h) - 1);
	right_side_2_val = (tao / (2 * h * h));
	b_val = -(1 + tao / (h * h));

	right = new double[task.n - 1];

	for (int j = 1; j <= task.m; ++j) {
		int size = task.n - 1;
		double *x = new double[size];
		for (int i = 0; i<size; i++)
			x[i] = 0.0;
		double **A = new double*[size];
		for (int i = 0; i<size; i++) {
			A[i] = new double[size];
			for (int k = 0; k < size; k++)
				A[i][k] = 0.;
		}
		A[0][0] = b_val; A[0][1] = ac_val;
		A[size - 1][size - 2] = ac_val; A[size - 1][size - 1] = b_val;
		for (int i = 1; i < size - 1; i++) {
			A[i][i] = b_val;
			A[i][i - 1] = ac_val;
			A[i][i + 1] = ac_val;
		}

		for (int i = 1; i < task.n; ++i) {

			right[i - 1] = right_side_1_val * G[j - 1][i] -
				right_side_2_val * (G[j - 1][i - 1] + G[j - 1][i + 1]) +
				tao * task.f((i - 1) * h, ((j - 1) + 0.5) * tao);

			/*a[i - 1] = c[i - 1] = (tao / (2 * h * h));
			b[i - 1] = -(1 + tao / (h * h));*/

		}

		cyclic_reduction(task.n - 1, A, right);
		cyclic_reduction_back_substitution(task.n - 1, A, right, x);
		//double *a_cr = cyclic_reduction_factor(task.n - 1, a);
		//double *x = cyclic_reduction_solve(task.n - 1, a_cr, 1, right);

		for (int i = 1; i < task.n; i++	)
			G[j][i] = x[i - 1];
		//delete[] a_cr, x;
		for (int i = 0; i<size; i++)
			delete A[i];
		delete[] A;
	}
done:
	for (int i = 0; i <= task.n; i++)
		v[i] = G[task.m][i];
}

int main(int argc, char **argv)
{
	heat_task task(10000, 10000, 1024, 1024);
	double *v = new double[1001];
	heat_equation_crank_nicolson(task, v);

	for (int i = 0; i < task.n; i++) {
		cout << v[i] << " ";
	}

	getchar();

	return 0;
}