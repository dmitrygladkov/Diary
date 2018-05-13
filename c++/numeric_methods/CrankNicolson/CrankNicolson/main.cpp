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

typedef std::vector<std::vector<double>> Grid;

static inline
void cycle_reduction_method(
	double *x, double a_0, double a_1, double a_2,
	double *a, double *b, double *c, double *f,
	int n, int q)
{
	a[0] = a_0;
	b[0] = a_1;
	c[0] = a_2;
	f[0] = 0;
	f[n] = 0;
	x[0] = 0;
	x[n] = 0;
	int start = 2, elementsNum = n, step = 1;

	for (int j = 0; j <= q - 1; j++) {
		double alpha = - a[j] / b[j];
		double beta = - c[j] / b[j];
		a[j + 1] = alpha * a[j];
		b[j + 1] = b[j] + 2 * alpha * c[j];
		c[j + 1] = beta * c[j];
		elementsNum = (elementsNum - 1) / 2;
		for (int i = 0; i <= elementsNum; i++) {
			int k = start * (i + 1);
			f[k] = alpha * f[k - step] + f[k] + beta * f[k + step];
		}
		start = 2 * start;
		step = 2 * step;
	}

	start = n / 2;
	step = start;
	elementsNum = 1;
	for (int j = q - 1; j <= 0; j++) {
		double alpha = -a[j] / b[j];
		double beta = -c[j] / b[j];
		for (int i = 0; i <= elementsNum; i++) {
			int k = start * (2 * i + 1);
			x[k] = f[k] / b[j] +
				alpha * x[k - step] +
				beta * x[k + step];
		}
		start = start / 2;
		step = start;
		elementsNum = elementsNum * 2;
	}
}

void heat_equation_crank_nicolson(heat_task task, double *v)
{
	double h = task.L / task.n;
	double tao = task.T / task.m;

	double **G = new double*[task.m + 1], **func = new double*[task.n + 1];
	for (int i = 0; i <= task.m; i++) {
		G[i] = new double[task.n + 1];
		func[i] = new double[task.n + 1];
	}

	for (int i = 0; i <= task.n; ++i) {
		G[0][i] = task.initial_condition(i * h);
	}
	for (int j = 0; j <= task.m; ++j) {
		G[j][0] = task.left_condition(j * tao);
		G[j][task.n] = task.right_condition(j * tao);
	}
	for (int i = 0; i <= task.m; ++i) {
		for (int j = 0; j <= task.n; ++j)
			func[i][j] = tao * task.f(i * h, (j + 0.5) * tao);
	}

	for (int i = 1; i < task.m; ++i) {
		double *a = new double[task.n + 1];
		double *b = new double[task.n + 1];
		double *c = new double[task.n + 1];
		double *right = new double[task.n + 1];

		right[0] = 0;
		for (int j = 0; j < task.n; j++) {
			right[j] = (1 - tao / (h * h)) * G[i][j] + (tao / (2 * h * h)) * (G[i - 1][j] + G[i + 1][j]) + func[i][j];
		}

		cycle_reduction_method(G[i + 1], (1 + tao / (h * h)), tao / (2 * h * h),
				       tao / (2 * h * h), a, b, c, func[i], task.n, sqrt(task.n));


		delete[] a, b, c, right;
	}

	for (int i = 0; i <= task.n; i++)
		v[i] = G[task.m][i];
}

int main(int argc, char **argv)
{
	heat_task task(0.1, 1, 8, 16);
	double *v = new double[11];
	heat_equation_crank_nicolson(task, v);
}