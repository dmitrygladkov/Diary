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

void solve(double* a, double* b, double* c, double* d, int n) {
	/*
	// n is the number of unknowns

	|b0 c0 0 ||x0| |d0|
	|a1 b1 c1||x1|=|d1|
	|0  a2 b2||x2| |d2|

	1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

	x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

	2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
	from 1st it.: -| a1x0 + a1g0x1        = a1r0
	-----------------------------
	(b1 - a1g0)x1 + c1x2 = d1 - a1r0

	x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

	3rd iteration:      | a2x1 + b2x2   = d2
	from 2nd it. : -| a2x1 + a2g1x2 = a2r2
	-----------------------
	(b2 - a2g1)x2 = d2 - a2r2
	x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
	Finally we have a triangular matrix:
	|1  g0 0 ||x0| |r0|
	|0  1  g1||x1|=|r1|
	|0  0  1 ||x2| |r2|

	Condition: ||bi|| > ||ai|| + ||ci||

	in this version the c matrix reused instead of g
	and             the d matrix reused instead of r and x matrices to report results
	Written by Keivan Moradi, 2014
	*/
	n--; // since we start from x0 (not x1)
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {
		d[i] -= c[i] * d[i + 1];
	}
}

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
	for (int j = 1; j <= task.m; ++j) {
		G[j][0] = task.left_condition(j * tao);
		G[j][task.n] = task.right_condition(j * tao);
	}
	for (int j = 0; j <= task.m; ++j) {
		for (int i = 0; i <= task.n; ++i)
			func[j][i] = tao * task.f(i * h, (j + 0.5) * tao);
	}

	for (int j = 0; j < task.m; ++j) {
		double *a = new double[task.n + 1];
		double *b = new double[task.n + 1];
		double *c = new double[task.n + 1];
		double *right = new double[task.n + 1];

		right[0] = func[j][0];
		for (int i = 0; i <= task.n; ++i) {
			right[i] = (1 - tao / (h * h)) * G[j][i] + (tao / (2 * h * h)) * (G[j][i - 1] + G[j][i + 1]) + func[j][i];

			a[i] = (tao / (2 * h * h));
			b[i] = (1 - tao / (h * h));
			c[i] = (tao / (2 * h * h));
		}
		a[0] = 0;
		c[task.n] = 0;
		solve(a, b, c, right, task.n);
		
		/*cycle_reduction_method(G[i + 1], (1 + tao / (h * h)), tao / (2 * h * h),
				       tao / (2 * h * h), a, b, c, func[i], task.n, sqrt(task.n));*/
		for (int i = 0; i <= task.n; i++)
			G[j + 1][i] = right[i];

		delete[] a, b, c, right;
	}

	for (int i = 0; i <= task.n; i++)
		v[i] = G[task.m][i];
}

int main(int argc, char **argv)
{
	heat_task task(10000, 10000, 4, 2);
	double *v = new double[11];
	heat_equation_crank_nicolson(task, v);

	for (int i = 0; i < task.n; i++) {
		cout << v[i] << " ";
	}
	getchar();

	return 0;
}