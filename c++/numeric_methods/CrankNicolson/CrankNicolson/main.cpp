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

double *cyclic_reduction_factor(int n, double a[])
{
	double *a_cr;
	int iful;
	int ifulp;
	int ihaf;
	int il;
	int ilp;
	int inc;
	int incr;
	int ipnt;
	int ipntp;
	int j;

	a_cr = new double[3 * (2 * n + 1)];

	if (n == 1) {
		a_cr[0 + 0 * 3] = 0.0;
		a_cr[0 + 1 * 3] = 0.0;
		a_cr[0 + 2 * 3] = 0.0;
		a_cr[1 + 0 * 3] = 0.0;
		a_cr[1 + 1 * 3] = 1.0 / a[1 + 1 * 3];
		a_cr[1 + 2 * 3] = 0.0;
		a_cr[2 + 0 * 3] = 0.0;
		a_cr[2 + 1 * 3] = 0.0;
		a_cr[2 + 2 * 3] = 0.0;
		return a_cr;
	}
	/*
	Zero out the workspace entries.
	*/
	a_cr[0 + 0 * 3] = 0.0;
	for (j = 1; j <= n - 1; j++) {
		a_cr[0 + j * 3] = a[0 + j * 3];
	}
	for (j = n; j <= 2 * n; j++) {
		a_cr[0 + j * 3] = 0.0;
	}

	a_cr[1 + 0 * 3] = 0.0;
	for (j = 1; j <= n; j++) {
		a_cr[1 + j * 3] = a[1 + (j - 1) * 3];
	}
	for (j = n + 1; j <= 2 * n; j++) {
		a_cr[1 + j * 3] = 0.0;
	}
	a_cr[2 + 0 * 3] = 0.0;
	for (j = 1; j <= n - 1; j++) {
		a_cr[2 + j * 3] = a[2 + (j - 1) * 3];
	}
	for (j = n; j <= 2 * n; j++) {
		a_cr[2 + j * 3] = 0.0;
	}

	il = n;
	ipntp = 0;

	while (1 < il) {
		ipnt = ipntp;
		ipntp = ipntp + il;
		if ((il % 2) == 1) {
			inc = il + 1;
		} else {
			inc = il;
		}

		incr = inc / 2;
		il = il / 2;
		ihaf = ipntp + incr + 1;
		ifulp = ipnt + inc + 2;

		for (ilp = incr; 1 <= ilp; ilp--) {
			ifulp = ifulp - 2;
			iful = ifulp - 1;
			ihaf = ihaf - 1;
			a_cr[1 + iful * 3] = 1.0 / a_cr[1 + iful * 3];
			a_cr[2 + iful * 3] = a_cr[2 + iful * 3] * a_cr[1 + iful * 3];
			a_cr[0 + ifulp * 3] = a_cr[0 + ifulp * 3] * a_cr[1 + (ifulp + 1) * 3];
			a_cr[1 + ihaf * 3] = a_cr[1 + ifulp * 3] - a_cr[0 + iful * 3] * a_cr[2 + iful * 3]
				- a_cr[0 + ifulp * 3] * a_cr[2 + ifulp * 3];
			a_cr[2 + ihaf * 3] = -a_cr[2 + ifulp * 3] * a_cr[2 + (ifulp + 1) * 3];
			a_cr[0 + ihaf * 3] = -a_cr[0 + ifulp * 3] * a_cr[0 + (ifulp + 1) * 3];
		}
	}

	a_cr[1 + (ipntp + 1) * 3] = 1.0 / a_cr[1 + (ipntp + 1) * 3];

	return a_cr;
}

double *cyclic_reduction_solve(int n, double a_cr[], int nb, double b[])
{
	int i;
	int iful;
	int ifulm;
	int ihaf;
	int il;
	int ipnt;
	int ipntp;
	int j;
	int ndiv;
	double *rhs;
	double *x;

	if (n == 1) {
		x = new double[nb * n];
		for (j = 0; j < nb; j++) {
			x[0 + j * n] = a_cr[1 + 1 * 3] * b[0 + j * n];
		}
		return x;
	}
	//
	//  Set up RHS.
	//
	rhs = new double[(2 * n + 1) * nb];

	for (j = 0; j < nb; j++) {
		rhs[0 + j * (2 * n + 1)] = 0.0;
		for (i = 1; i <= n; i++) {
			rhs[i + j * (2 * n + 1)] = b[i - 1 + j * n];
		}
		for (i = n + 1; i <= 2 * n; i++) {
			rhs[i + j * (2 * n + 1)] = 0.0;
		}
	}

	il = n;
	ndiv = 1;
	ipntp = 0;

	while (1 < il) {
		ipnt = ipntp;
		ipntp = ipntp + il;
		il = il / 2;
		ndiv = ndiv * 2;

		for (j = 0; j < nb; j++) {
			ihaf = ipntp;
			for (iful = ipnt + 2; iful <= ipntp; iful = iful + 2) {
				ihaf = ihaf + 1;
				rhs[ihaf + j * (2 * n + 1)] = rhs[iful + j * (2 * n + 1)]
					- a_cr[2 + (iful - 1) * 3] * rhs[iful - 1 + j * (2 * n + 1)]
					- a_cr[0 + iful * 3] * rhs[iful + 1 + j * (2 * n + 1)];
			}
		}
	}

	for (j = 0; j < nb; j++) {
		rhs[ihaf + j * (2 * n + 1)] = rhs[ihaf + j * (2 * n + 1)] * a_cr[1 + ihaf * 3];
	}

	ipnt = ipntp;

	while (0 < ipnt) {
		ipntp = ipnt;
		ndiv = ndiv / 2;
		il = n / ndiv;
		ipnt = ipnt - il;

		for (j = 0; j < nb; j++) {
			ihaf = ipntp;
			for (ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2) {
				iful = ifulm + 1;
				ihaf = ihaf + 1;
				rhs[iful + j * (2 * n + 1)] = rhs[ihaf + j * (2 * n + 1)];
				rhs[ifulm + j * (2 * n + 1)] = a_cr[1 + ifulm * 3] * (
					rhs[ifulm + j * (2 * n + 1)]
					- a_cr[2 + (ifulm - 1) * 3] * rhs[ifulm - 1 + j * (2 * n + 1)]
					- a_cr[0 + ifulm * 3] * rhs[iful + j * (2 * n + 1)]);
			}
		}
	}

	x = new double[n * nb];

	for (j = 0; j < nb; j++) {
		for (i = 0; i < n; i++) {
			x[i + j * n] = rhs[i + 1 + j * (2 * n + 1)];
		}
	}

	free(rhs);

	return x;
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

	a = new double[3 * task.n - 1];
	right = new double[task.n - 1];

	for (int j = 1; j <= task.m; ++j) {

		a[0 + 0 * 3] = 0.0;
		a[1 + 0 * 3] = b_val;

		for (int i = 1; i < task.n; ++i) {
			right[i - 1] = right_side_1_val * G[j - 1][i] -
				right_side_2_val * (G[j - 1][i - 1] + G[j - 1][i + 1]) +
				tao * task.f((i - 1) * h, ((j - 1) + 0.5) * tao);

			/*a[i - 1] = c[i - 1] = (tao / (2 * h * h));
			b[i - 1] = -(1 + tao / (h * h));*/
			a[0 + (i - 1) * 3] = a[2 + (i - 1) * 3] = ac_val;
			a[1 + (i - 1) * 3] = b_val;
		}

		a[2 + (task.n - 2) * 3] = 0.0;

		double *a_cr = cyclic_reduction_factor(task.n - 1, a);
		double *x = cyclic_reduction_solve(task.n - 1, a_cr, 1, right);

		for (int i = 1; i < task.n; i++)
			G[j][i] = x[i - 1];
		delete[] a_cr, x;
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