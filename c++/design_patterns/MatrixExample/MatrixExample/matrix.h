#pragma once
#include <vector>
#include <memory>
#include <list>
#include "drawer.h"
#include <functional>

template<class T>
class IMatrix {
public:
	virtual size_t Size() const = 0;
	virtual T Get(size_t i, size_t j) const = 0;
	virtual void Set(size_t i, size_t j, const T& value) = 0;
};

template<class T>
class IDrawableMatrix : public IDrawable, public IMatrix<T> {
public:
	void Draw(IDrawer* drawer) const override {
		DrawMatrix(*this, drawer);
	}
	virtual void DrawMatrix(const IMatrix<T>& matrix, IDrawer* drawer) const = 0;
};

template<class T>
class DenseMatrix : public IDrawableMatrix<T> {
public:
	DenseMatrix(size_t size) : size(size), elements(size, std::vector<T>(size)) {}
	size_t Size() const override { return size; }
	T Get(size_t i, size_t j) const override { return elements[i][j]; }
	void Set(size_t i, size_t j, const T& value) override { elements[i][j] = value; }
private:
	void DrawMatrix(const IMatrix<T>& matrix, IDrawer* drawer) const override {
		cout << "DenseMatrix::DrawMatrix" << endl;
		drawer->DrawBorder(size, size);
		drawer->DrawValues(size, size, true);
	};
	const size_t size;
	std::vector<std::vector<T>> elements;
};

struct CRSMatrix {
	int n;			// number of rows in a matrix
	int m;			// number of columns in a matrix
	int nz;			// number of non-zero elements in a sparse matrix
	vector<double> val;	// array of values of matrix per string
	vector<int> colIndex;	// array of column indices
	vector<int> rowPtr;	// array of row start indices
};

template<class T>
class SparseMatrix : public IDrawableMatrix<T> {
public:
	SparseMatrix(size_t size) : size(size) {
		matrix.n = matrix.m = size;
		matrix.nz = 0;
		matrix.val.resize(size);
		matrix.colIndex.resize(size);
		matrix.rowPtr.resize(size + 1);
	}
	size_t Size() const override { return size; }
	T Get(size_t i, size_t j) const override {
		return 0;
	}
	void Set(size_t i, size_t j, const T& value) override {
		
	}
private:
	void DrawMatrix(const IMatrix<T>& matrix, IDrawer* drawer) const override {
		cout << "SparseMatrix::DrawMatrix" << endl;
		drawer->DrawBorder(size, size);
		drawer->DrawValues(size, size, false);
	}
	CRSMatrix matrix;
	const size_t size;
};

template<class T>
class Transpose : public IDrawableMatrix<T> {
public:
	Transpose(IDrawableMatrix<T> *matrix) : size(matrix->Size()), matrix(matrix) {}
	size_t Size() const override { return size; }
	T Get(size_t i, size_t j) const override { return matrix->Get(j, i); }
	void Set(size_t i, size_t j, const T& value) override { matrix->Set(j, i, value); }
private:
	void DrawMatrix(const IMatrix<T>& _matrix, IDrawer* drawer) const override {
		cout << "Transpose::DrawMatrix" << endl;
		matrix->DrawMatrix(_matrix, drawer);
	}
	size_t size;
	IDrawableMatrix<T> *matrix;
};

template<class T>
class NullAboveMainDiagonal : public IDrawableMatrix<T> {
public:
	NullAboveMainDiagonal(IDrawableMatrix<T> *matrix) : size(matrix->Size()), matrix(matrix) { }
	size_t Size() const override { return size; }
	T Get(size_t i, size_t j) const override {
		if (i < j) {
			return 0;
		}
		return matrix->Get(i, j);
	}
	void Set(size_t i, size_t j, const T& value) override {
		matrix->Set(i, j, value);
	}
private:
	void DrawMatrix(const IMatrix<T>& _matrix, IDrawer* drawer) const override {
		cout << "NullAbovemainDiagonal::DrawMatrix" << endl;
		matrix->DrawMatrix(_matrix, drawer);
	}
	size_t size;
	IDrawableMatrix<T> *matrix;
};

template<class T>
class Minor : public IDrawableMatrix<T> {
public:
	Minor(IDrawableMatrix<T> *matrix,
			const std::vector<size_t>& rows,
			const std::vector<size_t>& columns) :
		size(rows.size()),
		matrix(matrix),
		rows(rows),
		columns(columns) {}
	size_t Size() const override { return size; }
	T Get(size_t i, size_t j) const override { return matrix->Get(rows[i], columns[j]); }
	void Set(size_t i, size_t j, const T& value) override { matrix->Set(rows[i], columns[j], value); }
private:
	void DrawMatrix(const IMatrix<T>& _matrix, IDrawer* drawer) const override {
		cout << "Minor::DrawMatrix" << endl;
		matrix->DrawMatrix(_matrix, drawer);
	}
	size_t size;
	IDrawableMatrix<T> *matrix;
	std::vector<size_t> rows;
	std::vector<size_t> columns;
};
