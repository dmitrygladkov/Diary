#include <iostream>
#include <iomanip>
#include <random>
#include <math.h>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <cassert>

using namespace std;

#define MATRIX(array, i, j, dim)	\
	(*((array) + (j) + (i) * (dim)))

#define BLOCK_MATRIX(array, i, j, begin_i, begin_j, total_dim)		\
	(*((array) + ((begin_j) + (j)) + (((begin_i) + (i)) * (total_dim))))

#define SYSTEM(system, i, dim)		\
	((system) + (i) * (dim))

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
	void matrix_multiplication(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim);
	void matrix_subtraction(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim);
	void matrix_addition(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim);
	inline void copy_from(T *matrix, size_t dim)
	{
		T *arr = array;
		for (size_t i = 0; i < dim; i++) {
			for (size_t j = 0; j < dim; j++)
				MATRIX(arr, i, j, dim) = this->matrix[i][j] = MATRIX(matrix, i, j, dim);
		}
	}
	void get_cofactor(T **mat, T **temp, size_t p, size_t q, size_t n);
	T determinant(T **mat, size_t n);
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
void Matrix<T>::matrix_multiplication(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim)
{
	Matrix_Multiplication(res_matrix, left_matrix, right_matrix, dim);
	copy_from(res_matrix, dim);
}

template <typename T>
static inline
void Matrix_Multiplication(T *res_matrix, T *left_matrix, T *right_matrix, int dim)
{
	Matrix_Multiplication_Rect(res_matrix, left_matrix, right_matrix, dim, dim, dim, dim);
}

template <typename T>
static inline
void Matrix_Multiplication_Rect(T *res_matrix, T *left_matrix, T *right_matrix, int left_dim1, int left_dim2, int right_dim1, int right_dim2)
{
	assert(left_dim2 == right_dim1);

	for (size_t i = 0; i < left_dim1; i++) {
		for (size_t j = 0; j < right_dim2; j++) {
			MATRIX(res_matrix, i, j, right_dim2) = 0;
			for (size_t k = 0; k < left_dim2; k++)
				MATRIX(res_matrix, i, j, right_dim2) +=
				(MATRIX(left_matrix, i, k, left_dim2) * MATRIX(right_matrix, k, j, right_dim2));
		}
	}
}

template <typename T>
static inline
void Matrix_Multiplication_Rect_left_block(T *res_matrix, T *left_matrix, T *right_matrix,
					int left_dim1, int left_dim2, int right_dim1, int right_dim2,
					int begin_i_left, int begin_j_left, int total_dim_left)
{
	assert(left_dim2 == right_dim1);

	for (size_t i = 0; i < left_dim1; i++) {
		for (size_t j = 0; j < right_dim2; j++) {
			MATRIX(res_matrix, i, j, right_dim2) = 0;
			for (size_t k = 0; k < left_dim2; k++) {
				MATRIX(res_matrix, i, j, right_dim2) +=
					(BLOCK_MATRIX(left_matrix, i, k, begin_i_left, begin_j_left, total_dim_left) *
					 MATRIX(right_matrix, k, j, right_dim2));
			}
		}
	}
}

template <typename T>
static inline
void Matrix_Multiplication_Rect_res_left_block(T *res_matrix, T *left_matrix, T *right_matrix, int left_dim1, int left_dim2, int right_dim1, int right_dim2,
						int begin_i_left, int begin_j_left, int total_dim_left,
						int begin_i_res, int begin_j_res, int total_dim_res)
{
	assert(left_dim2 == right_dim1);

	for (size_t i = 0; i < left_dim1; i++) {
		for (size_t j = 0; j < right_dim2; j++) {
			BLOCK_MATRIX(res_matrix, i, j, begin_i_res, begin_j_res, total_dim_res) = 0;

			for (size_t k = 0; k < left_dim2; k++) {
				BLOCK_MATRIX(res_matrix, i, j, begin_i_res, begin_j_res, total_dim_res) +=
					(BLOCK_MATRIX(left_matrix, i, k, begin_i_left, begin_j_left, total_dim_left) * MATRIX(right_matrix, k, j, right_dim2));
			}
		}
	}
}

template <typename T>
void Matrix<T>::matrix_transposition(T *matrix, size_t dim)
{
	Matrix_Transposition(matrix, dim);
	copy_from(matrix, dim);
}

template <typename T>
static inline
void Matrix_Transposition(T *matrix, size_t dim)
{
	T tmp;

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = i + 1; j < dim; j++) {
			tmp = MATRIX(matrix, i, j, dim);
			MATRIX(matrix, i, j, dim) = MATRIX(matrix, j, i, dim);
			MATRIX(matrix, j, i, dim) = tmp;
		}
	}
}

template <typename T>
static inline
void Matrix_Transposition_Rect(T *matrix_res, T *matrix, size_t dim1, size_t dim2)
{
	for (size_t i = 0; i < dim1; i++) {
		for (size_t j = 0; j < dim2; j++) {
			MATRIX(matrix_res, j, i, dim1) = MATRIX(matrix, i, j, dim2);
		}
	}
}

template <typename T>
static inline
void Matrix_Transposition_Rect_block(T *matrix_res, T *matrix, size_t dim1, size_t dim2, size_t begin_i, size_t begin_j, size_t total_dim)
{
	for (size_t i = 0; i < dim1; i++) {
		for (size_t j = 0; j < dim2; j++) {
			MATRIX(matrix_res, j, i, dim1) = BLOCK_MATRIX(matrix, i, j, begin_i, begin_j, total_dim);
		}
	}
}

template <typename T>
void Matrix<T>::matrix_subtraction(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim)
{
	Matrix_Subtraction(res_matrix, left_matrix, right_matrix, dim);
	copy_from(matrix.dim);
}

template <typename T>
static inline
void Matrix_Subtraction(T *matrix_res, T *matrix_left, T *matrix_right, size_t dim)
{
	Matrix_Subtraction_Rect(matrix_res, matrix_left, matrix_right, dim, dim);
}

template <typename T>
static inline
void Matrix_Subtraction_Rect_block(T *matrix_res, T *matrix_left, T *matrix_right, size_t dim1, size_t dim2,
				int begin_i, int begin_j, int total_dim)
{
	for (size_t i = 0; i < dim1; i++) {
		for (size_t j = 0; j < dim2; j++) {
			MATRIX(matrix_res, i, j, dim1) =
				BLOCK_MATRIX(matrix_left, i, j, begin_i, begin_j, total_dim) - MATRIX(matrix_right, i, j, dim1);
		}
	}
}

template <typename T>
static inline
void Matrix_Subtraction_Rect(T *matrix_res, T *matrix_left, T *matrix_right, size_t dim1, size_t dim2)
{
	for (size_t i = 0; i < dim1; i++) {
		for (size_t j = 0; j < dim2; j++) {
			MATRIX(matrix_res, i, j, dim1) =
				MATRIX(matrix_left, i, j, dim1) - MATRIX(matrix_right, i, j, dim1);
		}
	}
}

template <typename T>
void Matrix<T>::matrix_addition(T *res_matrix, T *left_matrix, T *right_matrix, size_t dim)
{
	Matrix_Addition(res_matrix, left_matrix, right_matrix, dim);
	copy_from(matrix.dim);
}

template <typename T>
static inline
void Matrix_Addition(T *matrix_res, T *matrix_left, T *matrix_right, size_t dim)
{
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			MATRIX(matrix_res, i, j, dim) =
				MATRIX(matrix_left, i, j, dim) + MATRIX(matrix_right, i, j, dim);
		}
	}
}

template <typename T>
static inline
void solve_lower_triangular_matrix(size_t dim, T *a, T *b, T *x)
{
	T s;

	for (size_t i = 0; i < dim; i++) {
		s = 0;
		for (size_t j = 0; j < i; j++)
			s += MATRIX(a, i, j, dim) * x[j];
		x[i] = (b[i] - s) / MATRIX(a, i, i, dim);
	}
}

template <typename T>
static inline
void solve_higher_triangular_matrix(size_t dim, T *a, T *b, T *x)
{
	T s;

	for (int64_t i = dim - 1; i >= 0; i--) {
		s = 0;
		for (int64_t j = dim - 1; j > i; j--)
			s += MATRIX(a, i, j, dim) * x[j];
		x[i] = (b[i] - s) / MATRIX(a, i, i, dim);
	}
}

template <typename T>
static inline
void solve_lower_triangular_matrix_system(size_t system_num, size_t dim, T *a, T *b_system, T *x_system)
{
	for (int64_t i = 0; i < system_num; i++)
		solve_lower_triangular_matrix(dim, a, SYSTEM(b_system, i, system_num), SYSTEM(x_system, i, system_num));
}

template <typename T>
static inline
void solve_higher_triangular_matrix_system(size_t system_num, size_t dim, T *a, T *b_system, T *x_system)
{
	for (int64_t i = 0; i < system_num; i++)
		solve_higher_triangular_matrix(dim, a, SYSTEM(b_system, i, system_num), SYSTEM(x_system, i, system_num));
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

template <typename T>
static inline
void get_cofactor_matrix(T *mat, T *temp, size_t p, size_t q, size_t n)
{
	size_t i = 0, j = 0;

	// Looping for each element of the matrix
	for (size_t row = 0; row < n; row++) {
		for (size_t col = 0; col < n; col++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				MATRIX(temp, i, j++, n - 1) = MATRIX(mat, row, col, n);

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
T Matrix<T>::determinant(T **mat, size_t n)
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
		D += sign * mat[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

template <typename T>
static inline
T determinant_matrix(T *mat, size_t n)
{
	T D = 0; // Initialize result

		 //  Base case : if matrix contains single element
	if (n == 1)
		return MATRIX(mat, 0, 0, n);

	// To store cofactors
	T *temp = new T[n*n];

	int sign = 1; // To store sign multiplier

		      // Iterate for each element of first row
	for (size_t f = 0; f < n; f++) {
		// Getting Cofactor of mat[0][f]
		get_cofactor_matrix(mat, temp, 0, f, n);

		D += sign * MATRIX(mat, 0, f, n) * determinant_matrix(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	delete[] temp;
	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
template <typename T>
static inline void adjoint_matrix(T *A, T *adj, size_t dim)
{
	if (dim == 1) {
		MATRIX(adj, 0, 0, dim) = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	T *temp = new T[dim*dim];

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			// Get cofactor of A[i][j]
			get_cofactor_matrix(A, temp, i, j, dim);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			MATRIX(adj, j, i, dim) = (sign)*(determinant_matrix(temp, dim - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
template <typename T>
static inline bool inverse_matrix(T *A, T *inverse, size_t dim)
{
	// Find determinant of A[][]
	T det = determinant_matrix(A, dim);
	if (det == 0) {
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	T * adj = new T[dim*dim];;
	adjoint_matrix(A, adj, dim);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (size_t i = 0; i < dim; i++)
		for (size_t j = 0; j < dim; j++)
			MATRIX(inverse, i, j, dim) = MATRIX(adj, i, j, dim) / det;

	return true;
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
		if (determinant(matrix, i + 1) <= 0)
			return 1;
	}
	return 0;
}

// L21 * L11T = A21
//	||
//	\/
// L21 = A21 * L11T_inverse
//
// This fucntion calculates L21
// NOTE:
//	A21			- (n - r) x r matrix
//	L11			- r x r matrix
//	L21			- (n - r) x r matrix
//	A21 * L11T_inverse 	- (n - r) x r matrix - tre result of the routine (i.e. L21)
static inline
void Cholesky_Solve_Second_Iteration(double *A21, double *L11, double *L21, int n, int r)
{
	double *L11T = new double[r * r], *L11T_inverse = new double[r * r];
	Matrix_Transposition_Rect(L11T, L11, r, r);
	bool res = inverse_matrix(L11T, L11T_inverse, r);
	assert(res);
	Matrix_Multiplication_Rect(L21, A21, L11T_inverse, (n - r), r, r, r);

	delete[] L11T_inverse, L11T;
}

// NOTE: int begin_i_A21, int begin_j_A21, int total_len_A21 are valid for L21
static inline
void Cholesky_Solve_Second_Iteration_block(double *A21, double *L11, double *L21, int n, int r,
	int begin_i_L11, int begin_j_L11, int total_len_L11,
	int begin_i_A21, int begin_j_A21, int total_len_A21)
{
	double *L11T = new double[r * r], *L11T_inverse = new double[r * r];
	Matrix_Transposition_Rect_block(L11T, L11, r, r, begin_i_L11, begin_j_L11, total_len_L11);
	bool res = inverse_matrix(L11T, L11T_inverse, r);
	assert(res);
	Matrix_Multiplication_Rect_res_left_block(L21, A21, L11T_inverse, (n - r), r, r, r,
		begin_i_A21, begin_j_A21, total_len_A21,
		begin_i_A21, begin_j_A21, total_len_A21);

	delete[] L11T_inverse, L11T;
}

// A22_red = A22 - L21 * L21T
//
// This fucntion calculates A22_red
// NOTE:
//	A22			- (n - r) x (n - r) matrix
//	L21			- (n - r) x r matrix
//	L21T			- r x (n - r) matrix
//	A22 - L21 * L21T 	- (n - r) x (n - r) matrix - the result of the routine (i.e. A22_red)
static inline
void Cholesky_Find_Reduced_Matrix(double *A22_red, double *A22, double *L21, int n, int r)
{
	double *L21T = new double[r * (n - r)], *L21_L21T = new double[(n - r) * (n - r)];

	Matrix_Transposition_Rect(L21T, L21, n - r, r);
	Matrix_Multiplication_Rect(L21_L21T, L21, L21T, (n - r), r, r, (n - r));
	Matrix_Subtraction_Rect(A22_red, A22, L21_L21T, (n - r), (n - r));

	delete[] L21T, L21_L21T;
}

static inline
void Cholesky_Find_Reduced_Matrix_block(double *A22_red, double *A22, double *L21, int n, int r,
					int begin_i_L21, int begin_j_L21, int total_len_L21,
					int begin_i_A22, int begin_j_A22, int total_len_A22)
{
	double *L21T = new double[r * (n - r)], *L21_L21T = new double[(n - r) * (n - r)];
	Matrix_Transposition_Rect_block(L21T, L21, n - r, r, begin_i_L21, begin_j_L21, total_len_L21);
	Matrix_Multiplication_Rect_left_block(L21_L21T, L21, L21T, (n - r), r, r, (n - r), begin_i_L21, begin_j_L21, total_len_L21);
	Matrix_Subtraction_Rect_block(A22_red, A22, L21_L21T, (n - r), (n - r), begin_i_A22, begin_j_A22, total_len_A22);

	delete[] L21T, L21_L21T;
}

static inline
void Cholesky_Decomposition_line(double *A, double *L, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			sum = 0;
			for (int k = 0; k < j; k++)
				sum += (MATRIX(L, i, k, n) * MATRIX(L, j, k, n));
			MATRIX(L, i, j, n) = (MATRIX(A, i, j, n) - sum) / MATRIX(L, j, j, n);
			sum = 0;
		}
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += pow(MATRIX(L, i, k, n), 2);
		MATRIX(L, i, i, n) = sqrt(MATRIX(A, i, i, n) - sum);
	}
}

static inline
void Cholesky_Decomposition_line_block(double *A, double *L, int n,
					int begin_i, int begin_j, int total_n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			sum = 0;
			for (int k = 0; k < j; k++)
				sum += (BLOCK_MATRIX(L, i, k, begin_i, begin_j, total_n) *
					BLOCK_MATRIX(L, j, k, begin_i, begin_j, total_n));
			BLOCK_MATRIX(L, i, j, begin_i, begin_j, total_n) =
				(BLOCK_MATRIX(A, i, j, begin_i, begin_j, total_n) - sum) /
				 BLOCK_MATRIX(L, j, j, begin_i, begin_j, total_n);
		}
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += pow(BLOCK_MATRIX(L, i, k, begin_i, begin_j, total_n), 2);
		BLOCK_MATRIX(L, i, i, begin_i, begin_j, total_n) =
			sqrt(BLOCK_MATRIX(A, i, i, begin_i, begin_j, total_n) - sum);
		sum = 0;
	}
}


static int BLOCK_SIZE = 4;

void Cholesky_Decomposition(double *A, double *L, int n)
{
	//BLOCK_SIZE = n / 4; // 25% of initial matrix
	//if (BLOCK_SIZE == 1) {
	//	Cholesky_Decomposition_line(A, L, n);
	//	return;
	//}
	if (n <= BLOCK_SIZE ) {
			//Cholesky_Decomposition_line(A, L, n);
			Cholesky_Decomposition_line_block(A, L, n, 0, 0, n);
			return;
	}

	// FIRST
	Cholesky_Decomposition_line_block(A, L, BLOCK_SIZE, 0, 0, n);


	// SECOND
	Cholesky_Solve_Second_Iteration_block(A, L, L, n, BLOCK_SIZE, 
						0 ,0, n,		// L11
						BLOCK_SIZE, 0, n	// L21, A21
						);


	double *A22_red = new double[(n - BLOCK_SIZE) * (n - BLOCK_SIZE)];
	/*double *A22 = new double[(n - BLOCK_SIZE) * (n - BLOCK_SIZE)];
	double *L21 = new double[(n - BLOCK_SIZE) * (BLOCK_SIZE)];
	for (size_t i = 0; i < n - BLOCK_SIZE; i++) {
		for (size_t j = 0; j < BLOCK_SIZE; j++) {
			MATRIX(L21, i, j, BLOCK_SIZE) = BLOCK_MATRIX(L, i, j, BLOCK_SIZE, 0, n);
		}
	}
	for (size_t i = 0; i < n - BLOCK_SIZE; i++) {
		for (size_t j = 0; j < n - BLOCK_SIZE; j++) {
			MATRIX(A22, i, j, n - BLOCK_SIZE) = BLOCK_MATRIX(A, i, j, BLOCK_SIZE, BLOCK_SIZE, n);
			cout << MATRIX(A22, i, j, n - BLOCK_SIZE) << ", ";
		}
		cout << endl;
	}
	Cholesky_Find_Reduced_Matrix(A22_red, A22, L21, n, BLOCK_SIZE);*/

	// THIRD
	Cholesky_Find_Reduced_Matrix_block(A22_red, A, L, n, BLOCK_SIZE,
						BLOCK_SIZE, 0, n,		// L21
						BLOCK_SIZE, BLOCK_SIZE, n	// A22
						);

	double *L22_red = new double[(n - BLOCK_SIZE) * (n - BLOCK_SIZE)];
	memset(L22_red, 0, sizeof(double) * (n - BLOCK_SIZE) * (n - BLOCK_SIZE));
	Cholesky_Decomposition(A22_red, L22_red, n - BLOCK_SIZE);
	size_t k = 0, l = 0;
	for (size_t i = BLOCK_SIZE; i < n; i++) {
		for (size_t j = BLOCK_SIZE; j < n; j++) {
			MATRIX(L, i, j, n) = MATRIX(L22_red, l, k, n - BLOCK_SIZE);
			k++;
		}
		k = 0;
		l++;
	}

	delete[] L22_red, A22_red;
}

int main(char **argv, int argc)
{
	Matrix<double> matrix_obj(6), matrix_res(6), matrix_check(6);
	size_t dim = matrix_obj.get_dimension();
	double *matrix = matrix_obj.get_1d_array();
	double *result = new double[dim * dim], *result_check = new double[dim * dim], *result_t = new double[dim * dim];
	memset(result, 0, sizeof(double) * dim * dim);

	cout << matrix_obj << endl;

	Cholesky_Decomposition(matrix, result, dim);

	matrix_res.set_1d_array(result);

	Matrix_Transposition_Rect(result_t, result, dim, dim);
	Matrix_Multiplication(result_check, result, result_t, dim);

	matrix_check.set_1d_array(result_check);

	cout << "Result:" << endl;
	cout << matrix_res << endl;

	cout << "Check:" << endl;
	cout << matrix_check << endl;

	getchar();

	delete[] result;

	/*std::copy(L21, L21 + 4, std::ostream_iterator<float>(std::cout, ","));
	cout << endl;
	std::copy(A21_red, A21_red + 4, std::ostream_iterator<float>(std::cout, ","));

	getchar();*/
}