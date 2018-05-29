#include "matrix.h"
#include "drawer.h"

#include <conio.h>

int main(int argc, char **argv)
{
	{
		cout << "=== Testing Bridge ===" << endl;
		cout << "=== INIT ===" << endl;
		IDrawer *ConsoleMatrixDrawer = new ConsoleDrawer();
		IDrawer *HTMLMatrixDrawer = new HTMLDrawer();
		IDrawer *GraphicsMatrixDrawer = new GraphicsDrawer();
		IDrawableMatrix<double> *Dense = new DenseMatrix<double>(100);
		IDrawableMatrix<double> *Sparse = new SparseMatrix<double>(100);
		cout << "=== ~INIT~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== DRAW ===" << endl;
		Dense->Draw(ConsoleMatrixDrawer);
		cout << endl;
		Dense->Draw(HTMLMatrixDrawer);
		cout << endl;
		Dense->Draw(GraphicsMatrixDrawer);
		cout << endl;
		Sparse->Draw(ConsoleMatrixDrawer);
		cout << endl;
		Sparse->Draw(HTMLMatrixDrawer);
		cout << endl;
		Sparse->Draw(GraphicsMatrixDrawer);
		cout << "=== ~DRAW~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== FINALIZE ===" << endl;
		delete ConsoleMatrixDrawer, HTMLMatrixDrawer, GraphicsMatrixDrawer;
		delete Dense, Sparse;
		cout << "=== ~FINALIZE~ ===" << endl;
		cout << "=== ~Testing Bridge~ ===" << endl;
	}
	_getch();
	system("cls");
	{
		cout << "=== Testing Decorator ===" << endl;
		cout << "=== INIT ===" << endl;
		IDrawer *ConsoleMatrixDrawer = new ConsoleDrawer();

		IDrawableMatrix<double> *Dense = new DenseMatrix<double>(100);
		IDrawableMatrix<double> *MinorDense = new Minor<double>(Dense, { 0, 1, 2, 3 }, { 0, 1, 2, 3 });
		IDrawableMatrix<double> *TransposeMinorDense = new Transpose<double>(MinorDense);
		IDrawableMatrix<double> *NullTransposeMinorDense = new NullAboveMainDiagonal<double>(TransposeMinorDense);
		cout << "=== ~INIT~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== DRAW ===" << endl;
		Dense->Draw(ConsoleMatrixDrawer);
		cout << endl;
		MinorDense->Draw(ConsoleMatrixDrawer);
		cout << endl;
		TransposeMinorDense->Draw(ConsoleMatrixDrawer);
		cout << endl;
		NullTransposeMinorDense->Draw(ConsoleMatrixDrawer);
		cout << "=== ~DRAW~ ===" << endl;

		cout << endl << endl << endl;

		cout << "=== FINALIZE ===" << endl;
		delete Dense, MinorDense, TransposeMinorDense, NullTransposeMinorDense;
		cout << "=== ~FINALIZE~ ===" << endl;
		cout << "=== ~Testing Decorator~ ===" << endl;
	}
	_getch();
	return 0;
}