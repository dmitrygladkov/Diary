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

#include <omp.h>

using namespace std;

struct CRSMatrix {
	int n;			// number of rows in a matrix
	int m;			// number of columns in a matrix
	int nz;			// number of non-zero elements in a sparse matrix
	vector<double> val;	// array of values of matrix per string
	vector<int> colIndex;	// array of column indices
	vector<int> rowPtr;	// array of row start indices
};

static inline double getCRSElem(const CRSMatrix &crs, int i, int j)
{
	for (int k = crs.rowPtr[i]; k < crs.rowPtr[i + 1]; ++k) {
		if (crs.colIndex[k] == j)
			return crs.val[k];
	}
	return 0.0;
}

static inline void showCRSMatrix(const CRSMatrix &crs)
{
	for (int i = 0; i < crs.m; i++) {
		for (int j = 0; j < crs.n; j++) {
			cout << setw(10) << getCRSElem(crs, i, j) << " ";
		}
		cout << endl;
	}
}

class CRS_Matrix {
	vector<int> non_zero_row_elem_num;

	inline int getRowOffset(int i) {
		return std::accumulate(non_zero_row_elem_num.begin(), non_zero_row_elem_num.begin() + i, 0);
	}

	void generateCRSMatrix() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, UCHAR_MAX);

		int f, c;

		srand(time(NULL));

		for (int i = 0; i < crs.n; i++) {
			// Form column indices in row `i`
			for (int j = 0; j < non_zero_row_elem_num[i]; j++) {
				do {
					crs.colIndex[getRowOffset(i) + j] = rand() % crs.n;
					f = 0;
					for (int k = 0; k < j; k++) {
						if (crs.colIndex[getRowOffset(i) + j] ==
							crs.colIndex[getRowOffset(i) + k])
							f = 1;
					}
				} while (f == 1);
			}
			// Sort column indices in row `i`
			for (int j = 0; j < non_zero_row_elem_num[i] - 1; j++) {
				for (int k = 0; k < non_zero_row_elem_num[i] - 1; k++) {
					if (crs.colIndex[getRowOffset(i) + k] >
						crs.colIndex[getRowOffset(i) + k + 1]) {
						std::swap(crs.colIndex[getRowOffset(i) + k],
							crs.colIndex[getRowOffset(i) + k + 1]);
					}
				}
			}
		}
		// Fill array of values
		for (int i = 0; i < crs.nz; i++)
			crs.val[i] = dis(gen);
		// Fill array of rows indices
		c = 0;
		for (int i = 0; i < crs.n; i++) {
			crs.rowPtr[i] = c;
			c += non_zero_row_elem_num[i];
		}
		crs.rowPtr[crs.n] = c;
	}
public:
	CRSMatrix crs;
	CRS_Matrix(int dim = 10, vector<int> non_zero_row_elem_num = { 0, 1, 0, 1, 0, 0, 2, 0, 0, 3 }/*{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }*/)
	{
		crs.n = dim;
		crs.m = dim;

		this->non_zero_row_elem_num = non_zero_row_elem_num;
		while (non_zero_row_elem_num.size() < dim)
			this->non_zero_row_elem_num.push_back(0);

		crs.nz = std::accumulate(non_zero_row_elem_num.begin(), non_zero_row_elem_num.end(), 0);

		crs.val.resize(crs.nz);
		crs.colIndex.resize(crs.nz);
		crs.rowPtr.resize(dim + 1);

		generateCRSMatrix();
	}

	~CRS_Matrix() {
		crs.n = crs.m = crs.nz = 0;
		crs.val.clear();
		crs.colIndex.clear();
		crs.rowPtr.clear();
	}

	inline double getElem(int i, int j)
	{
		return getCRSElem(crs, i, j);
	}

	inline int getDim()
	{
		return crs.n;
	}

	inline int getNonZeroElemNumber()
	{
		return crs.nz;
	}
private:
	friend ostream& operator<<(ostream& os, const CRS_Matrix& obj)
	{
		for (int i = 0; i < obj.crs.m; i++) {
			for (int j = 0; j < obj.crs.n; j++) {
				os << setw(10) << getCRSElem(obj.crs, i, j) << " ";
			}
			os << endl;
		}
		return os;
	}
};

static inline void initCRSMatrix(CRSMatrix &crs_matrix, int size, int non_zero)
{
	crs_matrix.n = crs_matrix.m = size;
	crs_matrix.nz = non_zero;
	crs_matrix.val.resize(non_zero);
	crs_matrix.colIndex.resize(non_zero);
	crs_matrix.rowPtr.resize(size + 1);
}

static inline void initVector(double *vector, int size, double range)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, range);

	for (int i = 0; i < size; i++)
		vector[i] = dis(gen);
}

static inline void transposeCRSMatrix(const CRSMatrix &crs_matrix, CRSMatrix &crs_matrix_t) {
	std::fill(crs_matrix_t.rowPtr.begin(), crs_matrix_t.rowPtr.end(), 0);
	for (int i = 0; i < crs_matrix.nz; i++)
		crs_matrix_t.rowPtr[crs_matrix.colIndex[i] + 1]++;

	int S = 0, tmp = 0;
	for (int i = 1; i <= crs_matrix.n; i++) {
		tmp = crs_matrix_t.rowPtr[i];
		crs_matrix_t.rowPtr[i] = S;
		S += tmp;
	}

	for (int i = 0; i < crs_matrix_t.n; i++) {
		for (int j = crs_matrix.rowPtr[i]; j < crs_matrix.rowPtr[i + 1]; j++) {
			crs_matrix_t.val[crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]] = crs_matrix.val[j];
			crs_matrix_t.colIndex[crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]] = i;
			crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]++;
		}
	}
}

static inline void multiplicationCRSMatrixVector(const CRSMatrix &crs_matrix, const double *vector, double *result_vector)
{
	for (int i = 0; i < crs_matrix.n; i++) {
		result_vector[i] = 0;
		for (int j = crs_matrix.rowPtr[i]; j < crs_matrix.rowPtr[i + 1]; j++)
			result_vector[i] += crs_matrix.val[j] * vector[crs_matrix.colIndex[j]];
	}
}

static inline double multiplicationVectorVector(double *vector1, double *vector2, int size_vector)
{
	double sum = 0.0;
#pragma omp parallel for simd schedule (static,16) reduction(+:sum) 
	for (int i = 0; i < size_vector; i++)
		sum += vector1[i] * vector2[i];
	return sum;
}

int main(int argc, char **argv)
{
	CRS_Matrix crs_matrix;
	CRSMatrix crs_matrix_t;
	double *b = new double[crs_matrix.getDim()];
	double *x = new double[crs_matrix.getDim()];

	initVector(b, crs_matrix.getDim(), UCHAR_MAX);

	initCRSMatrix(crs_matrix_t, crs_matrix.getDim(), crs_matrix.getNonZeroElemNumber());
	transposeCRSMatrix(crs_matrix.crs, crs_matrix_t);

	cout << crs_matrix;
	cout << endl;
	showCRSMatrix(crs_matrix_t);

	getchar();
}