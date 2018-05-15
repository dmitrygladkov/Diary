#pragma once

#include "IMatrixDrawer.h"

class IMatrix
{
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
