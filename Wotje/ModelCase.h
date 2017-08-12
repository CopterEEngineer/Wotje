#pragma once

// header files
#include "MatrixTemplate.h"
#include "BladeCSD.h"


class ModelCase {
private:
	int ni;
	Matrix1<myTYPE> ristation, rho_b, eiflap, eilag, gj, iflap, ilag, ipitch;

public:
	int dof;
	Matrix2<myTYPE> DM, MM;
	Matrix1<FreeDomType> frdtype;

public:
	ModelCase();

	ModelCase(const ModelCase &M);

	~ModelCase();

	inline void GetModel(void);

	inline void CantileverModel(void);

};

inline void GetModel(void)
{
	CantileverModel();
}


ModelCase::ModelCase()
{
	dof = 2;
	frdtype.allocate(dof);
	DM.allocate(dof, dof);
	MM.allocate(dof, dof);
	frdtype.setvalue(Free);

	ni = 5;
	ristation.allocate(ni);
	rho_b.allocate(ni);
	eiflap.allocate(ni);
	eilag.allocate(ni);
	gj.allocate(ni);
	iflap.allocate(ni);
	ilag.allocate(ni);
	ipitich.allocate(ni);
}


ModelCase::ModelCase(const ModelCase &M)
{
	ni = M.ni;
	ristation = M.ristation;
	rho_b = M.rho_b;
	eiflap = M.eiflap;
	eilag = M.eilag;
	gj = M.gj;
	iflap = M.iflap;
	ilag = M.ilag;
	ipitch = M.ipitch;
	dof = M.dof;
	DM = M.DM;
	MM = M.MM;
	frdtype = M.frdtype;
}


ModelCase::~ModelCase()
{
	ni = dof = 0;
	ristation.deallocate();
	rho_b.deallocate();
	eiflap.deallocate();
	eilag.deallocate();
	gj.deallocate();
	iflap.deallocate();
	ilag.deallocate();
	ipitch.deallocate();
	DM.deallocate();
	MM.deallocate();
	frdtype.deallocate();
}


inline void ModelCase::CantileverModel(void)
{
	int ni = 6;

	ristation.allocate(ni);
	for (int i = 0; i < ni; ++i)
		ristation(i) = (double)i / (ni - 1);

	rho_b.allocate(ni);
	rho_b.setvalue(1.0);

	iflap.allocate(ni);
	iflap.setvalue(1.0);

	eiflap.allocate(ni);
	myTYPE E = 1.0;
	eiflap = iflap*E;

	dof = 2;
	frdtype.allocate(dof);
	frdtype(0) = DispPlaneOut;
	frdtype(1) = BendPlaneOut;


}
