#include "IMatrixDrawer.h"
#include "IMatrix.h"

int main(int argc, char **argv)
{
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

	getchar();

	return 0;
}