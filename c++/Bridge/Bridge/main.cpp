#include "IMatrixDrawer.h"
#include "IMatrix.h"

#include <conio.h>

int main(int argc, char **argv)
{
	{
		cout << "=== Testing Bridge ===" << endl;
		cout << "=== INIT ===" << endl;
		IMatrix *DenseConsole = new DenseMatrix(new ConsoleMatrixDrawer),
			*DenseGraph = new DenseMatrix(new GraphMatrixDrawer),
			*DenseHTML = new DenseMatrix(new HTMLMatrixDrawer);
		IMatrix *SparseConsole = new SparseMatrix(new ConsoleMatrixDrawer),
			*SparseGraph = new SparseMatrix(new GraphMatrixDrawer),
			*SparseHTML = new SparseMatrix(new HTMLMatrixDrawer);
		cout << "=== ~INIT~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== DRAW ===" << endl;
		DenseConsole->Draw(100, 100);
		DenseGraph->Draw(100, 100);
		DenseHTML->Draw(100, 100);
		SparseConsole->Draw(100, 100);
		SparseGraph->Draw(100, 100);
		SparseHTML->Draw(100, 100);
		cout << "=== ~DRAW~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== FINALIZE ===" << endl;
		delete DenseConsole, DenseGraph, DenseHTML;
		delete SparseConsole, SparseGraph, SparseHTML;
		cout << "=== ~FINALIZE~ ===" << endl;
		cout << "=== ~Testing Bridge~ ===" << endl;
	}
	_getch();
	system("cls");
	{
		cout << "=== Testing Decorator ===" << endl;
		cout << "=== INIT ===" << endl;
		IMatrix *DenseConsole = new DenseMatrix(new ConsoleMatrixDrawer);
		MinorMatrix *MinorDenseConsole = new MinorMatrix(DenseConsole);
		TransposeMatrix *TransposeMinorMatrixConsole = new TransposeMatrix(MinorDenseConsole);
		NullAboveMainDiagonalMatrix *NullTransposeMinorMatrixConsole = new NullAboveMainDiagonalMatrix(TransposeMinorMatrixConsole);
		cout << "=== ~INIT~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== DRAW ===" << endl;
		DenseConsole->Draw(100, 100);
		MinorDenseConsole->Draw(100, 100);
		TransposeMinorMatrixConsole->Draw(100, 100);
		NullTransposeMinorMatrixConsole->Draw(100, 100);
		cout << "=== ~DRAW~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== FINALIZE ===" << endl;
		delete DenseConsole, MinorDenseConsole, TransposeMinorMatrixConsole, NullTransposeMinorMatrixConsole;
		cout << "=== ~FINALIZE~ ===" << endl;
		cout << "=== ~Testing Bridge~ ===" << endl;
	}
	_getch();
	return 0;
}