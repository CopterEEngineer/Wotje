// Wotje.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include "BladeCSD.h"
#include "BeamCSD.h"
//#include "ModelCase.h"

int main()
{
	//Blade blade;

	Beam_1D beam1d;
	ModelCase cant;

	cant.GetModel();
	cant.GetProb();
	beam1d.InitBeam1D(cant);
	beam1d.AssembleKM();
	beam1d.SetBdryCond();
	beam1d.StaticSolution();
	//beam1d.ModeSolution();

	//beam1d._eleKMgenrt_N2(me, ke);
	//Matrix3<double> mm, aa;
	//beam1d._eleKMgenrt_N2(mm, aa);

	//blade.functest();
	//blade.CantileverModel();
	//blade.AssembleKM();
	//blade.SetBdryCond();
	//blade.StaticSolving();
	//blade.ModeSolving();
	//int a[6] = { 0, 1, 2, 3, 4, 5 };
	//Matrix1<FreeDomType> a(2);
	//Matrix1<int> _a(2);
	//a(0) = DispPlaneOut;
	//a(1) = BendPlaneOut;
	/*for (int i = 0; i < 6; i++)
		a(i) = i;*/
	//FreeDomType v = DispPlaneOut;
	//int _v;
	//int b = 0, e = 1;
	//int n = 2;
	//int asub[6] = { 0 };
	//int vp[2] = { 2.5, 4.5 };
	//cout<<static_cast<int>(FreeDomType::DispAxis) << endl;
	//for (int i = 0; i < 2; ++i)
	//	_a(i) = static_cast<int>(a(i));
	//_v = static_cast<int>(v);
	//cout<<static_cast<int>(a.v_p[0]) << endl;
	//cout << BiSearch(_a.v_p, b, e, v) << endl;
	//cout << BiSearchRound(a, b, e, v) << endl;
	//cout << BiSearchFloor(a, b, e, v) << endl;
	//cout << BiSearchCeil(a, b, e, v) << endl;
	//BiSearchRange(a, b, e, vp, n, asub);
	/*for(int i=0;i<n;++i)
		cout << asub[i] << endl;

	int nEle = 10, nNodeEle = 3, nNodePart = nEle*nNodeEle-(nEle-1);
	Matrix2<int> nodeID(nEle, nNodeEle);
	for (int i = 0; i < nNodeEle; ++i)
		nodeID(0, i) = i + 1;
	for (int i = 1; i < nEle; ++i)
		for (int j = 0; j < nNodeEle; ++j)
		{
			if (j > 0)
				nodeID(i, j) = nodeID(i, j - 1) + 1;
			else
				nodeID(i, j) = nodeID(i - 1, nNodeEle - 1);
		}
	nodeID(3, 1) = 99, nodeID(3, 2) = 98, nodeID(4, 0) = 98;
	nodeID.output();

	Matrix1<int> nodeID_seq(nNodePart);
	int k;
	for (int i = 0; i<nNodeEle; ++i)
	{
		k = i;
		for (int j = 0; j < nEle; ++j)
		{
			nodeID_seq(k) = nodeID(j, i);
			k += (nNodeEle - 1);
		}
	}
	nodeID_seq.output();*/

    system("pause");
	return 0;
}

