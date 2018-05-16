#pragma once

#include "IMatrixDrawer.h"

class IMatrix {
public:
	IMatrix(IMatrixDrawer *drawer)
	{
		cout << "constructor IMatrix" << endl;
		this->drawer = drawer;
	}
	IMatrix()
	{
		cout << "constructor IMatrix" << endl;
	}
	void SetDrawer(IMatrixDrawer *drawer)
	{
		this->drawer = drawer;
	}
	IMatrixDrawer *drawer;
	virtual void Draw(int weight, int height) = 0;
	virtual ~IMatrix()
	{
		cout << "destructor IMatrix" << endl;
		delete drawer;
	}
};

class DenseMatrix : public IMatrix {
public:
	DenseMatrix(IMatrixDrawer *drawer) : IMatrix(drawer)
	{
		cout << "constructor DenseMatrix" << endl;
	}
	void Draw(int weight, int height)
	{
		cout << "DenseMatrix::Draw" << endl;
		drawer->DrawOutBorder(weight, height);
		drawer->DrawCellBorder(0, 0);
		drawer->DrawCell(0, 0, 10);
	}
	virtual ~DenseMatrix()
	{
		cout << "destructor DenseMatrix" << endl;
	}
};

class SparseMatrix : public IMatrix {
public:
	SparseMatrix(IMatrixDrawer *drawer) : IMatrix(drawer)
	{
		cout << "constructor SparseMatrix" << endl;
	}
	void Draw(int weight, int height)
	{
		cout << "SparseMatrix::Draw" << endl;
		drawer->DrawOutBorder(weight, height);
		drawer->DrawCell(0, 0, 10);
	}
	virtual ~SparseMatrix()
	{
		cout << "destructor SparseMatrix" << endl;
	}
};

class MinorMatrix : public IMatrix {
public:
	IMatrix *super;
	MinorMatrix(IMatrix *matrix_decorator) : IMatrix(matrix_decorator->drawer)
	{
		super = matrix_decorator;
		cout << "constructor MinorMatrix" << endl;
	}
	void Draw(int weight, int height) {
		cout << "MinorMatrix::Draw" << endl;
		super->Draw(weight, height);
	}
	virtual ~MinorMatrix()
	{
		cout << "destructor MinorMatrix" << endl;
	}
};

class TransposeMatrix : public IMatrix {
public:
	IMatrix * super;
	TransposeMatrix(IMatrix *matrix_decorator) : IMatrix(matrix_decorator->drawer)
	{
		super = matrix_decorator;
		cout << "constructor TransposeMatrix" << endl;
	}
	void Draw(int weight, int height) {
		cout << "TransposeMatrix::Draw" << endl;
		super->Draw(weight, height);
	}
	virtual ~TransposeMatrix()
	{
		cout << "destructor TransposeMatrix" << endl;
	}
};

class NullAboveMainDiagonalMatrix : public IMatrix {
public:
	IMatrix * super;
	NullAboveMainDiagonalMatrix(IMatrix *matrix_decorator) : IMatrix(matrix_decorator->drawer)
	{
		super = matrix_decorator;
		cout << "constructor NullAboveMainDiagonalMatrix" << endl;
	}
	void Draw(int weight, int height) {
		cout << "NullAboveMainDiagonalMatrix::Draw" << endl;
		super->Draw(weight, height);
	}
	virtual ~NullAboveMainDiagonalMatrix()
	{
		cout << "destructor NullAboveMainDiagonalMatrix" << endl;
	}
};
