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
	for (int j = 0; j < q; j++) {
		double alpha = -a[j] / b[j];
		double beta = -c[j] / b[j];
		a[j + 1] = alpha * a[j];
		b[j + 1] = b[j] + alpha * a[j] + beta * c[j];
		c[j + 1] = beta * c[j];

		elementsNum = (elementsNum - 1) / 2;
#pragma omp parallel for
		for (int i = 0; i < elementsNum; i++) {
			int k = start * (i + 1);
			f[k] = alpha * f[k - step] + f[k] + beta * f[k + step];
		}
		start = 2 * start;
		step = 2 * step;
	}

	start = n / 2;
	step = start;
	elementsNum = 1;
	for (int j = q - 1; j >= 0; j--) {
		double alpha = -a[j] / b[j];
		double beta = -c[j] / b[j];
#pragma omp parallel for
		for (int i = 0; i < elementsNum; i++) {
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

	double *G = new double[task.n + 1];
#pragma omp parallel for
	for (int i = 0; i <= task.n; ++i) {
		G[i] = task.initial_condition(i * h);
	}

	double ac_val;
	double right_side_1_val, right_side_2_val;
	double b_val;

	ac_val = (tao / (2 * h * h));
	right_side_1_val = (tao / (h * h) - 1);
	right_side_2_val = (tao / (2 * h * h));
	b_val = -(1 + tao / (h * h));

	{
		double *a = new double[log2(task.n) + 1];
		double *b = new double[log2(task.n) + 1];
		double *c = new double[log2(task.n) + 1];
		double *right = new double[task.n + 1];

		for (int j = 1; j <= task.m; ++j) {
#pragma omp parallel for
			for (int i = 2; i <= task.n; ++i) {
				right[i - 1] = right_side_1_val * G[(i - 1)] -
					right_side_2_val * (G[(i - 1) - 1] + G[(i - 1) + 1]) -
					tao * task.f((i - 1) * h,
						     ((j - 1) + 0.5) * tao);
			}
			right[1] -= task.left_condition((j) * tao) * (tao / (h * h * 2));
			right[task.n - 1] -= task.right_condition((j) * tao) * (tao / (h * h * 2));
			/* re-use G */
			cycle_reduction_method(G, ac_val,
				b_val, ac_val, a, b, c, right, task.n, log2(task.n));
			G[0] = task.left_condition(j * tao);
			G[task.n] = task.right_condition(j * tao);
		}
#pragma omp parallel for
		for (int i = 0; i <= task.n; i++)
			v[i] = G[i];

		delete[] a, b, c, right, G;
	}
}

int main(int argc, char **argv)
{
	/*heat_task task(10000, 10000, 8, 8);
	double *v = new double[1001];
	heat_equation_crank_nicolson(task, v);
	for (int i = 0; i <= task.n; i++) {
		cout << v[i] << " ";
	}*/
	int n = 8;
	double *a = new double[log2(n) + 1];
	double *b = new double[log2(n) + 1];
	double *c = new double[log2(n) + 1];
	double *right = new double[n + 1];
	double *v = new double[n + 1];

	for (int i = 1; i < n; i++)
		right[i] = i;

	cycle_reduction_method(v, -1,
		2, 1, a, b, c, right, n, log2(n));

	for (int i = 0; i <= n; i++) {
		cout << v[i] << " ";
	}

	getchar();

	return 0;
}