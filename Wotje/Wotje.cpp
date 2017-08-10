// Wotje.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include "MatrixTemplate.h"
#include "BladeCSD.h"

void Mat2Arr(myTYPE **A, Matrix2<myTYPE> &B) {
	// B is col major
	for (int i = 0; i < B.NI; ++i) {
		cout << &B.v_p[i] << endl;
		cout << i*B.NJ << endl;
		cout << &B.v_p[i] + i*B.NJ << endl;
		*A = &B.v_p[i] + i*B.NJ;
		for (int j = 0; j < B.NJ; ++j) {
			*(*A + j) = B(i, j);
		}
	}
}

int main()
{
	Blade blade;
	Matrix2<myTYPE> B(3, 4);
	for (int i = 0; i < B.NI; ++i) {
		for (int j = 0; j < B.NI; ++j) {
			B(i, j) = i + j + 1;
		}
	}
	B.output();

	//myTYPE A[3];
	//myTYPE **A = NULL;
	//*A = new myTYPE[3];

	////myTYPE A[3][4] = { 0.0 };
	//Mat2Arr(A, B);
	//for (int i = 0; i < B.NI; ++i) {
	//	for (int j = 0; j < B.NJ; ++j) {
	//		cout << A[i][j] << "/t";
	//	}
	//}
	
	blade.functest();

	system("pause");
	return 0;
}

