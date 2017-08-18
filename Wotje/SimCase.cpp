#include "stdafx.h"
#include "SimCase.h"

void ProbDefinition::GetProb(void)
{
	// cantilver case
	dof = 2, nbdry = 2;
	frdtype.allocate(dof);
	bdrytype.allocate(dof, nbdry);
	bdryvalue.allocate(dof, nbdry);
	bdryposition.allocate(nbdry, 2);

	probtype = Dynamic, dt = 0.01, Nstep = 1000;
	frdtype(0) = DispPlaneOut;
	frdtype(1) = BendPlaneOut;
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

void Cant::GetProb(void)
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

void Cant::GetModel(void)
{
	myTYPE _rho = 1.0, _E = 1.0, _iflap = 1.0, _A = sqrt(4 * PI);

	length = 1.0, ni = 6;
	ristation.allocate(ni);
	for (int i = 0; i < ni; ++i)
		ristation(i) = length * i / (ni - 1);

	rho_b.allocate(ni);
	rho_b.setvalue(_rho);

	iflap.allocate(ni);
	iflap.setvalue(_iflap*_rho / _A);

	eiflap.allocate(ni);
	eiflap.setvalue(_E*_iflap);
}

void Beam_1D::InitFreeDom(void)
{
	p1.input("p1.output"); // get freedoms as matrix 2
	int k = 0;
	for (int i = 0; i < nNodeperPart; ++i)
		for (int j = 0; j < dof; ++j)
			q0(k++) = p1(i, j); // give q0 as order
}

