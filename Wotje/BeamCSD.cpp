#include "stdafx.h"
#include "BeamCSD.h"


void ProbDefinition::GetProb(void)
{
	dof = 2, nbdry = 2;
	frdtype.allocate(dof);
	bdrytype.allocate(dof, nbdry);
	bdryvalue.allocate(dof, nbdry);
	bdryposition.allocate(nbdry, 2);

	probtype = Mode;
	frdtype(0) = DispPlaneIn;
	frdtype(1) = BendPlaneIn;
	for (int i = 0; i < dof; ++i)
		bdrytype(i, 0) = Fix;
	bdrytype(0, 1) = ConstForc;
	bdryvalue(0, 1) = 1.0;
	// each row index boundary, coloum 0: start point, column 1: end point
	// [0, 1]
	bdryposition(0, 0) = 0.0;
	bdryposition(0, 1) = 0.0;
	bdryposition(1, 0) = 1.0;
	bdryposition(1, 1) = 1.0;
}

void ModelCase::CantileverModel(void)
{
	myTYPE _rho = 2.0, _E = 2.0, _iflap = 1, _A = sqrt(4*PI);

	length = 1.0, ni = 6;
	ristation.allocate(ni);
	for (int i = 0; i < ni; ++i)
		ristation(i) = length * i / (ni - 1);

	rho_b.allocate(ni);
	rho_b.setvalue(_rho);

	iflap.allocate(ni);
	iflap.setvalue(_iflap*_rho/_A);

	eiflap.allocate(ni);
	eiflap.setvalue(_E*_iflap);
}

void ModelCase::GetModel(void)
{
	CantileverModel();
}

void ModelCase::GetProb(void)
{
	probdef.GetProb();
	MM.allocate(probdef.dof, probdef.dof, ni);
	DM.allocate(probdef.dof, probdef.dof, ni);
	for (int i = 0; i < ni; ++i) {
		MM(0, 0, i) = rho_b(i);
		MM(1, 1, i) = iflap(i);
		DM(1, 1, i) = eiflap(i);
	}
}

void Beam_1D::InitBeam1D(ModelCase &M)
{
	nEle = 23, nNodeperEle = 2, nNodeperPart = nEle*nNodeperEle - (nEle - 1);
	nGauss = 2;
	nodeID.allocate(nEle, nNodeperEle);
	nodePosition.allocate(nNodeperPart);

	model = M;
	dof = M.probdef.dof, nBdry = M.probdef.nbdry, ptype = M.probdef.probtype;
	bdrytype = M.probdef.bdrytype, bdryvalue = M.probdef.bdryvalue;
	frdtype = M.probdef.frdtype;
	length = M.length, ristation = M.ristation;
	MM = M.MM, DM = M.DM;

	q.allocate(dof*nNodeperPart), p0.allocate(nNodeperPart, dof), p1.allocate(nNodeperPart, dof);

	for (int i = 0; i < nNodeperEle; ++i)
		nodeID(0, i) = i + 1;
	for (int i = 1; i < nEle; ++i)
		for (int j = 0; j < nNodeperEle; ++j)
		{
			if (j > 0)
				nodeID(i, j) = nodeID(i, j - 1) + 1;
			else
				nodeID(i, j) = nodeID(i - 1, nNodeperEle - 1);
			
			nodePosition(nodeID(i, j)-1) = (double)(i*length / nEle + j*length / nEle / (nNodeperEle - 1));

		}

	/*for (int i = 0; i < nEle; ++i)
		for (for j = 0; j < nNodeperEle; ++j)
			nodePosition(nodeID(i, j)) = (double)(i*length / nEle + j*length / nEle/ (nNodeperEle - 1));*/

	//for (int i = 0; i < nNodeperPart - 1; ++i)
	//	nodePosition(i) = (double)(length*i) / (nNodeperPart - 1);
	int k;
	nodeID_seq.allocate(nNodeperPart);
	for (int i = 0; i<nNodeperEle; ++i)
	{
		k = i;
		for (int j = 0; j < nEle; ++j)
		{
			nodeID_seq(k) = nodeID(j, i);
			k += (nNodeperEle - 1);
		}
	}

	nodeID_bdry.allocate(2, nBdry); // start nodeID , end nodeID
	for (int i = 0; i < nBdry; ++i)
	{
		nodeID_bdry(0, i) = 1 + BiSearchRound(nodePosition.v_p, 0, nNodeperPart - 1, M.probdef.bdryposition(i, 0)*M.length);
		nodeID_bdry(1, i) = 1 + BiSearchRound(nodePosition.v_p, 0, nNodeperPart - 1, M.probdef.bdryposition(i, 1)*M.length);
	}

	Fp.allocate(dof*nNodeperPart);
	q.allocate(nNodeperPart*dof);
	p0.allocate(nNodeperPart, dof);
	p1.allocate(nNodeperPart, dof);

	printf("nodeID: \n");
	nodeID.output();
	printf("\n");
	printf("nodeID_seq: \n");
	nodeID_seq.output();
	printf("\n");
	printf("nodeID_bdry: \n");
	nodeID_bdry.output();
	printf("\n");
}

void Beam_1D::AssembleKM(void)
{
	AssembleKM(Ma, Ka);
}

bool Beam_1D::SetBdryCond(void)
{
	Matrix1<int> asub(nNodeperPart), _nodeID_bdry(2);
	int n;
	bool flg;

	Fp.setvalue(0);
	for (int i = 0; i < nBdry; ++i)
	{
		_nodeID_bdry = nodeID_bdry(step(0, 1), i); // i boundary start nodeID and end nodeID
		_nodeID_bdry.output();
		n = 2;
		flg = BiSearchRange(nodeID_seq.v_p, 0, nNodeperPart - 1, _nodeID_bdry.v_p, n, asub.v_p);
		for (int j = 1; j <= dof; ++j)
		{
			switch (bdrytype(j - 1, i))
			{
			case Free:
				break;
			case Fix:
				if (flg)
					for (int w = 0; w < n; ++w)
						_FixedDisp(Fp, asub(w), j, ptype);

				break;
			case ConstDisp:
				if (flg)
					for (int w = 0; w < n; ++w)
						_ConstDisp(Fp, bdryvalue(j-1, i), asub(w), j);
				break;
			case ConstForc:
				if (flg)
					for (int w = 0; w < n; ++w)
						_ConstForceF(Fp, bdryvalue(j - 1, i), asub(w), j);
				break;
			default:
				printf("Undefined Boundary Type in SetBdryCond() Func. \n");
				return false;
			}
		}
	}
	flg &= ReviseF();
	Fp.output(4);

	return flg;
}

bool Beam_1D::ReviseF(void)
{
	Matrix1<int> asub(nNodeperPart), _nodeID_bdry(2);
	int n;
	bool flg;

	for (int i = 0; i < nBdry; ++i)
	{
		_nodeID_bdry = nodeID_bdry(step(0, 1), i); // i boundary start nodeID and end nodeID
		flg = BiSearchRange(nodeID_seq.v_p, 0, nNodeperPart - 1, _nodeID_bdry.v_p, n, asub.v_p);
		for (int j = 1; j <= dof; ++j)
		{
			switch (bdrytype(j-1,i))
			{
			case Fix:
				if (flg)
					for (int w = 0; w < n; ++w)
						_FixedDispF(Fp, asub(w), i);
				break;
			case ConstDisp:
				if (flg)
					for (int w = 0; w < n; ++w)
						_ConstDispF(Fp, bdryvalue(j - 1, i), asub(w), i);
				break;
			default:
				printf("Undefined Boundary Type in ReviseF() Func. \n");
				return false;
			}
		}
	}
	return flg;
}

void Beam_1D::StaticSolution(void)
{
	if (!Msolver(Ka, Fp, q))
	{
		printf("Static Solving Fail. \n");
		q.setvalue(0);
	}
		

	for (int j = 0; j < dof; ++j)
		for (int i = 0; i < nNodeperPart; ++i)
			p1(i, j) = p0(i, j) + q(dof*i + j);

	p0.output("p0.output", 4);
	p1.output("p1.output", 10);
	q.output("q.output", 4);
	p1.output(4);
}

void Beam_1D::ModeSolution(void)
{
	Matrix1<myTYPE> lam;

	// define search range
	ModePart.wmin = 1;
	ModePart.wmax = 1000;
	// define mode number interested
	ModePart.Nwmax = 40;

	double emin = ModePart.wmin*ModePart.wmin, emax = ModePart.wmax*ModePart.wmax;

	SolveSpEig(Ka, Ma, lam, ModePart.Vn, ModePart.Nwmax, emin, emax, ModePart.Nw);

	if (ModePart.wn.Nv == 0) ModePart.wn.allocate(ModePart.Nw);

	for (int i = 0; i < ModePart.Nw; i++) {
		if (lam(i) < 0) {
			cout << "Warning:  Eig " << "lam(" << i << ") = " << lam(i) << " < 0" << endl;
			ModePart.wn(i) = 0;
			system("pause");
		}
		else {
			ModePart.wn(i) = sqrt(lam(i));
		}
	}

	cout << "Find wn(1:" << ModePart.Nw << ") = " << endl;
	for (int i = 0; i < ModePart.Nw; i++) cout << ModePart.wn(i) << " rad/s = " << ModePart.wn(i) / 2 / PI << " Hz " << endl;
	ModePart.Vn.output("Vn.output", 4);
}