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
#include <math.h>
#include <unistd.h>

#include <omp.h>

using namespace std;

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
		double *a = new double[(int)log2(task.n) + 1];
		double *b = new double[(int)log2(task.n) + 1];
		double *c = new double[(int)log2(task.n) + 1];
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
	int max_omp_threads = omp_get_max_threads();
	double *start_time = new double[max_omp_threads + 1];
	double *diff_time = new double[max_omp_threads + 1];
	double *speedup = new double[max_omp_threads + 1];
	double *efficiency = new double[max_omp_threads + 1];

	heat_task task(10000, 10000, 62500, 10000);
	double *v = new double[62501];
	
	omp_set_num_threads(1);
	start_time[0] = omp_get_wtime();
	heat_equation_crank_nicolson(task, v);
	diff_time[0] = omp_get_wtime() - start_time[0];
	speedup[0] = 1;
	efficiency[0] = speedup[0] / 1;

	cout << "Run time - " << diff_time[0] << " ";
	cout << "Speedup - " << speedup[0] << " ";
	cout << "Efficiency - " << efficiency[0] << endl;

	for (int i = 1; i < max_omp_threads + 1; i++) {
		omp_set_num_threads(i);
		//cout << "Result (parallel " << i <<  " threads):" << endl;

		start_time[i] = omp_get_wtime();
		heat_equation_crank_nicolson(task, v);
		diff_time[i] = omp_get_wtime() - start_time[i];
		speedup[i] = diff_time[0] / diff_time[i];
		efficiency[i] = speedup[i] / i;
		cout << diff_time[i] << ", ";
		cout << speedup[i] << ", ";
		cout << efficiency[i];
		cout << endl;
	}

	return 0;
}
