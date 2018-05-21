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
	double *x, double a_0, double a_1,
	double *ac, double *b, double *f,
	int n, int q)
{
	ac[0] = a_0;
	b[0] = a_1;
	/*f[0] = 0;*/
	f[n] = 0;
	/*
	x[0] = 0;
	x[n] = 0;*/
	int start = 2, elementsNum = n, step = 1;

	for (int j = 0; j < q - 1; j++) {
		double alpha = -ac[j] / b[j];
		ac[j + 1] = alpha * ac[j];
		b[j + 1] = b[j] + 2 * alpha * ac[j];

		elementsNum = (elementsNum - 1) / 2;
		for (int i = 0; i < elementsNum; i++) {
			int k = start * (i + 1);
			f[k] = alpha * f[k - step] + f[k] + alpha * f[k + step];
		}
		start = 2 * start;
		step = 2 * step;
	}

	start = n / 2;
	step = start;
	elementsNum = 1;
	for (int j = q - 1 - 1; j >= 0; j--) {
		double alpha = -ac[j] / b[j];
		for (int i = 0; i < elementsNum; i++) {
			int k = start * (2 * i + 1);
			x[k] = f[k] / b[j] +
				alpha * x[k - step] +
				alpha * x[k + step];
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
		double *ac = new double[task.n + 1];
		double *b = new double[task.n + 1];
		double *right = new double[task.n + 1];

		for (int j = 1; j <= task.m; ++j) {

			for (int i = 1; i < task.n; ++i) {
				right[i - 1] = right_side_1_val * G[i] -
					right_side_2_val * (G[i - 1] + G[i + 1]) +
					tao * task.f((i - 1) * h, ((j - 1) + 0.5) * tao);
				ac[i - 1] = 0;
				b[i - 1] = 0;
			}
			ac[task.n] = ac[task.n - 1] = 0;
			b[task.n] = b[task.n - 1] = 0;
			right[task.n - 1] = 0;

			/* re-use G */
			G[0] = task.left_condition(j * tao);
			G[task.n] = task.right_condition(j * tao);

			cycle_reduction_method(G, ac_val,
				b_val, ac, b, right, task.n, log2(task.n));
		}

		for (int i = 0; i <= task.n; i++)
			v[i] = G[i];

		delete[] ac, b, right, G;
	}
}

int main(int argc, char **argv)
{
	heat_task task(10000, 10000, 2048, 8192);
	double *v = new double[1001];
	heat_equation_crank_nicolson(task, v);

	for (int i = 0; i <= task.n; i++) {
		cout << v[i] << " ";
	}

	getchar();

	return 0;
}