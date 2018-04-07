#include <iostream>
#include <iomanip>
#include <random>
#include <math.h>
#include <ctime>

using namespace std;

#define MATRIX(array, i, j, dim)	\
	(*((array) + (j) + (i) * (dim)))

template <typename T>
class Matrix {
public:
	Matrix(size_t dimension = 10, int gen_method = 1) : dimension(dimension), gen_method(gen_method) {
		matrix = new T*[dimension];
		array = new T[dimension * dimension];
		for (size_t i = 0; i < dimension; i++)
			matrix[i] = new T[dimension];
		if (gen_method == 1)
			generate_sym_pos_def_matrix_method1();
		else
			generate_sym_pos_def_matrix_method2();
		if (check_sym_pos_def_matrix())
			cerr << "Error! Matrix isn't symmetric and/or not positive definite" << endl;
	}
	size_t get_dimension() {
		return dimension;
	}
	T **get_matrix() {
		return matrix;
	}
	T *get_1d_array() {
		return array;
	}
	void set_matrix(T **matrix) {
		for (size_t i = 0; i < dimension; i++) {
			for (size_t j = 0; j < dimension; j++)
				MATRIX(this->array, i, j, dimension) =
					this->matrix[i][j] = matrix[i][j];
		}
	}
	void set_1d_array(T *array) {
		for (size_t i = 0; i < dimension; i++) {
			for (size_t j = 0; j < dimension; j++)
				MATRIX(this->array, i, j, dimension) =
					this->matrix[i][j] = MATRIX(array, i, j, dimension);
		}
	}
	~Matrix() {
		for (size_t i = 0; i < dimension; i++)
			delete[] matrix[i];
		delete[] matrix;
		delete[] array;
	}
private:
	T **matrix;
	T *array;
	size_t dimension;
	int gen_method;
	void generate_sym_pos_def_matrix_method1(void);
	void generate_sym_pos_def_matrix_method2(void);
	int check_sym_pos_def_matrix(void);
	void matrix_transposition(T *matrix, size_t dim);
	void matrix_multiplication(T *res_matrix, T *left_matrix, T *right_matrix, int dim);
	void get_cofactor(T **mat, T **temp, size_t p, size_t q, size_t n);
	T determinant_matrix(T **mat, size_t n);
	friend ostream& operator<<(ostream& os, const Matrix& obj)
	{
		for (size_t i = 0; i < obj.dimension; i++) {
			for (size_t j = 0; j < obj.dimension; j++) {
				os << setw(10) << obj.matrix[i][j] << " ";
			}
			os << endl;
		}
		return os;
	}
};

template <typename T>
void Matrix<T>::generate_sym_pos_def_matrix_method1(void)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, UCHAR_MAX);

	for (size_t i = 0; i < dimension; i++) {
		MATRIX(array, i, i, dimension) = matrix[i][i] = 0;
		for (size_t j = 0; j < i; j++)
			MATRIX(array, i, i, dimension) = matrix[i][i] += matrix[j][i];
		for (size_t j = i + 1; j < dimension; j++)
			MATRIX(array, i, i, dimension) = matrix[i][i] +=
			MATRIX(array, j, i, dimension) = matrix[j][i] = 
			MATRIX(array, i, j, dimension) = matrix[i][j] = dis(gen);
	}
}

template <typename T>
void Matrix<T>::generate_sym_pos_def_matrix_method2(void)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, UCHAR_MAX);
	T *tr_matrix = new T[dimension * dimension];
	T *tmp_matrix = new T[dimension * dimension];

	for (size_t i = 0; i < dimension; i++) {
		for (size_t j = 0; j < dimension; j++)
			MATRIX(tmp_matrix, i, j, dimension) =
			MATRIX(tr_matrix, i, j, dimension) =
			MATRIX(array, i, j, dimension) =
			matrix[i][j] = dis(gen);
	}
	matrix_transposition(tr_matrix, dimension);
	matrix_multiplication(array, tmp_matrix, tr_matrix, dimension);
	for (size_t i = 0; i < dimension; i++) {
		for (size_t j = 0; j < dimension; j++)
			matrix[i][j] = MATRIX(array, i, j, dimension);
	}
	delete[] tr_matrix;
	delete[] tmp_matrix;
}

template <typename T>
void Matrix<T>::matrix_multiplication(T *res_matrix, T *left_matrix, T *right_matrix, int dim)
{
	size_t i, j, k;

	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			MATRIX(res_matrix, i, j, dim) = 0;
			for (k = 0; k < dim; k++)
				MATRIX(res_matrix, i, j, dim) +=
					(MATRIX(left_matrix, i, k, dim) * MATRIX(right_matrix, k, j, dim));
		}
	}
}

template <typename T>
void Matrix<T>::matrix_transposition(T *matrix, size_t dim) {
	size_t i, j;
	T tmp = 0;

	for (i = 0; i < dim; i++) {
		for (j = i; j < dim; j++) {
			if (i != j) {
				tmp = MATRIX(matrix, i, j, dim);
				MATRIX(matrix, i, j, dim) = MATRIX(matrix, j, i, dim);
				MATRIX(matrix, j, i, dim) = tmp;
			}

		}
	}
}

// Function to get cofactor of mat[p][q] in temp[][]. n is current
// dimension of mat[][]
template <typename T>
void Matrix<T>::get_cofactor(T **mat, T **temp, size_t p, size_t q, size_t n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				temp[i][j++] = mat[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
n is current dimension of mat[][]. */
template <typename T>
T Matrix<T>::determinant_matrix(T **mat, size_t n)
{
	T D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1)
		return mat[0][0];

	// To store cofactors
	T **temp = new T*[n];
	for (size_t i = 0; i < n; i++)
		temp[i] = new T[n];

	int sign = 1; // To store sign multiplier

	// Iterate for each element of first row
	for (size_t f = 0; f < n; f++) {
		// Getting Cofactor of mat[0][f]
		get_cofactor(mat, temp, 0, f, n);
		D += sign * mat[0][f] * determinant_matrix(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

template <typename T>
int Matrix<T>::check_sym_pos_def_matrix(void)
{
	for (size_t i = 0; i < dimension; i++) {
		// Check for symmetric
		double sum = 0;
		for (size_t j = 0; j < i; j++)
			sum += matrix[i][j];
		for (size_t j = i + 1; j < dimension; j++) {
			if (MATRIX(array, j, i, dimension) != MATRIX(array, i, j, dimension) ||
				(matrix[j][i] != matrix[i][j]))
				return 1;
			sum += matrix[i][j];
		}
		if ((gen_method == 1) && (matrix[i][i] != sum))
			return 1;
		// Check for postive definite
		if (determinant_matrix(matrix, i + 1) <= 0)
			return 1;
	}
	return 0;
}

void Cholesky_Decomposition_line(double *A, double *L, int n)
{
	double sum = 0;
	memset(L, 0, (n * n) * sizeof(double));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++)
				sum += (MATRIX(L, i, k, n) * MATRIX(L, j, k, n));
			MATRIX(L, i, j, n) = (MATRIX(A, i, j, n) - sum) / MATRIX(L, j, j, n);
		}
		for (int k = 0; k < i; k++)
			sum += pow(MATRIX(L, i, k, n), 2);
		MATRIX(L, i, i, n) = sqrt(MATRIX(A, i, i, n) - sum);
	}
}

int main(char **argv, int argc)
{
	Matrix<double> matrix_obj(5), matrix_res(5);
	size_t dim = matrix_obj.get_dimension();
	double *matrix = matrix_obj.get_1d_array();
	double *result = new double[dim * dim];

	cout << matrix_obj << endl;

	Cholesky_Decomposition_line(matrix, result, dim);

	matrix_res.set_1d_array(result);

	cout << matrix_res << endl;
	getchar();

	delete[] result;
}