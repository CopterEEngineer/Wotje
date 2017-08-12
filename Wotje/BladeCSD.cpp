#include "stdafx.h"
#include "BladeCSD.h"


Blade::Blade()
{
	;
}


Blade::Blade(const Blade &B)
{
	ptype = B.ptype;
	idnode_bdry = B.idnode_bdry;
	bdryvalue = B.bdryvalue;
	bdrytype = B.bdrytype;
	nbd = B.nbd;
	nele = B.nele;
	nnode = B.nnode;
	nnode_all = B.nnode_all;
	dof = B.dof;
	ngs = B.ngs;
	nf = B.nf;
	ni = B.ni;
	omega = B.omega;
	ristation = B.ristation;
	rho_b = B.rho_b;
	eiflap = B.eiflap;
	eilag = B.eilag;
	gj = B.gj;
	iflap = B.iflap;
	ilag = B.ilag;
	ipitch = B.ipitch;
	cent = B.cent;
	nodegrid = B.nodegrid;
	id_node = B.id_node;
	Ka = B.Ka;
	Ma = B.Ma;
	DM = B.DM;
	MM = B.MM;
	Fp = B.Fp;
	q = B.q;
	p0 = B.p0;
	p1 = B.p1;
}


Blade::~Blade()
{
	nele = nnode = nnode_all = 0;
	dof = ngs = nbd = 0;
	nf = ni = 0;
	omega = 0;
	ristation.deallocate();
	rho_b.deallocate();
	eiflap.deallocate();
	eilag.deallocate();
	gj.deallocate();
	iflap.deallocate();
	ilag.deallocate();
	ipitch.deallocate();
	cent.deallocate();

	nodegrid.deallocate();
	id_node.deallocate();
	Ma.deallocate();
	Ka.deallocate();
	DM.deallocate();
	MM.deallocate();
	ptype = Static;
	idnode_bdry.deallocate();
	bdryvalue.deallocate();
	bdrytype.deallocate();
	Fp.deallocate();
	q.deallocate();
	p0.deallocate();
	p1.deallocate();
}



void Blade::SetBdryCond(void)
{
	_setbdrycond();
}


void Blade::EleKMgenrt(Matrix3<myTYPE>& Mele, Matrix3<myTYPE>& Kele)
{
	if (!_eleKMgenrt(Mele, Kele))
	{
		Mele.allocate(2 * dof, 2 * dof, nele);
		Kele.allocate(2 * dof, 2 * dof, nele);
	}

}


void Blade::AssembleKM(SpMtrx<myTYPE>& M, SpMtrx<myTYPE>& K)
{
	Matrix3<myTYPE> mele, kele;
	if (!_eleKMgenrt(mele, kele)) {
		M.allocate(1, 1);
		K.allocate(1, 1);
	}
	else {
		KK2SpM(M, mele, id_node, 0);
		KK2SpM(K, kele, id_node, 0);
	}
}


void Blade::AssembleKM(void)
{
	Matrix3<myTYPE> mele, kele;
	if (!_eleKMgenrt(mele, kele)) {
		Ma.allocate(1, 1);
		Ka.allocate(1, 1);
	}
	else {
		KK2SpM(Ma, mele, id_node, 0);
		KK2SpM(Ka, kele, id_node, 0);
	}
}


void Blade::StaticSolving(void)
{
	if (!Msolver(Ka, Fp, q))
		q.setvalue(0);

	for (int j = 0; j < dof; ++j) 
		for (int i = 0; i < nnode_all; ++i)
			p1(i, j) = p0(i, j) + q(dof*i + j);

	p0.output("p0.output", 4);
	p1.output("p1.output", 4);
	q.output("q.output", 4);
	
	// check F
	/*myTYPE Fx, Fy;
	Fx = Fy = 0;
	for (int i = 0; i < nnode_all; ++i) {
		Fx += Fp(i*dof - 2);
		Fy += Fp(i*dof - 1);
	}*/
}


void Blade::ModeSolving(void)
{
	;
}


void Blade::BladeModel_496(void)
{
	nf = 72;
	ni = 11;
	omega = 54;
	ristation.allocate(ni);
	ristation.input("ristation.txt");

	rho_b.allocate(ni);
	rho_b.setvalue(0.1);

	eiflap.allocate(ni);
	eiflap.setvalue(1.0);

	eilag.allocate(ni);
	eilag.setvalue(1.0);

	gj.allocate(ni);
	gj.setvalue(1.0);

	iflap.allocate(ni);
	iflap.setvalue(1.0);

	ilag.allocate(ni);
	ilag.setvalue(1.0);

	ipitch.allocate(ni);
	ipitch.setvalue(1.0);

	//ristation.output(4);
	Matrix1<int> idnn;
	myTYPE dr = (ristation(ni - 1) - ristation(0)) / (ni - 1);

	cent.allocate(ni);
	for (int i = ni - 1; i >= 0; --i) {
		idnn = step(i, ni - 1);
		cent(i) = omega * omega * (rho_b(idnn)*ristation(idnn)).sum() * dr;
	}

	nele = 10; // number of elements
	nnode = 2;
	dof = 5;
	ngs = 3;
	nbd = 2;
	nnode_all = nnode*nele - (nele - 1);

	Ka.allocate(1, 1);
	Ma.allocate(1, 1);

	nodegrid.allocate(nnode_all); // grid location of element nodes
	id_node.allocate(nele, nnode);

	// start from DIGIT ONE
	for (int i = 0; i < nnode; ++i)
		id_node(0, i) = i + 1;
	for (int i = 1; i < nele; ++i)
	{
		for (int j = 0; j < nnode; ++j)
			if (j > 0)
				id_node(i, j) = id_node(i, j - 1) + 1;
			else
				id_node(i, j) = id_node(i - 1, nnode - 1);
	}

	for (int i = 0; i < nnode_all; ++i)
		nodegrid(i) = i*1.0 / (nnode_all - 1);

	id_node.output();
	nodegrid.output();

	ptype = Dynamic;

	idnode_bdry.allocate(1, nbd);
	bdrytype.allocate(dof, nbd);
	bdryvalue.allocate(dof, nbd);

	Fp.allocate(dof*nnode_all);
	q.allocate(dof*nnode_all);

	p0.allocate(nnode_all, dof);
	p1.allocate(nnode_all, dof);

	frdtype.allocate(dof);
	frdtype(0) = DispPlaneOut;
	frdtype(1) = DispPlaneIn;
	frdtype(2) = BendPlaneIn;
	frdtype(3) = BendPlaneOut;
	frdtype(4) = RotaAxis;

	myTYPE grid;
	DM.allocate(dof, dof, nnode_all);
	MM.allocate(dof, dof, nnode_all);
	for (int i = 0; i < nnode_all; ++i) {
		grid = nodegrid(i);
		MM(0, 0, i) = rho_b.interplinear_fast(ristation, grid);
		MM(1, 1, i) = MM(0, 0, i);
		MM(2, 2, i) = iflap.interplinear_fast(ristation, grid);
		MM(3, 3, i) = ilag.interplinear_fast(ristation, grid);
		MM(4, 4, i) = ipitch.interplinear_fast(ristation, grid);

		DM(0, 0, i) = cent.interplinear_fast(ristation, grid);
		DM(1, 1, i) = DM(0, 0, i);
		DM(2, 2, i) = eiflap.interplinear_fast(ristation, grid);
		DM(3, 3, i) = eilag.interplinear_fast(ristation, grid);
		DM(4, 4, i) = gj.interplinear_fast(ristation, grid);
	}

	DM.output("DM.output", 4);
	MM.output("MM.output", 4);
}


void Blade::CantileverModel(void)
{
	// separate structure, allocate memory
	int ni = 6;

	ristation.allocate(ni);
	for (int i = 0; i < ni; ++i)
		ristation(i) = (double)i / (ni - 1);

	nele = 10; // number of elements
	nnode = 2;
	dof = 2;
	ngs = 2;
	nbd = 2;
	nnode_all = nnode*nele - (nele - 1);

	Ka.allocate(1, 1);
	Ma.allocate(1, 1);

	nodegrid.allocate(nnode_all); // grid location of element nodes
	id_node.allocate(nele, nnode);

	// start from DIGIT ONE
	for (int i = 0; i < nnode; ++i)
		id_node(0, i) = i + 1;
	for (int i = 1; i < nele; ++i)
	{
		for (int j = 0; j < nnode; ++j)
			if (j > 0)
				id_node(i, j) = id_node(i, j - 1) + 1;
			else
				id_node(i, j) = id_node(i - 1, nnode - 1);
	}

	for (int i = 0; i < nnode_all; ++i)
		nodegrid(i) = i*1.0 / (nnode_all - 1);

#ifdef OUTPUT_MODE1
	cout << endl << endl;
	id_node.output();
	nodegrid.output();
	cout << endl << endl;
#endif // OUTPUT_MODE1

	// Set FreeDom
	frdtype.allocate(dof);
	frdtype(0) = DispPlaneOut;
	frdtype(1) = BendPlaneOut;

	q.allocate(dof*nnode_all);
	p0.allocate(nnode_all, dof);
	p1.allocate(nnode_all, dof);

	// Problem definition
	ptype = Static;

	// Boundary definition
	idnode_bdry.allocate(1, nbd);
	bdrytype.allocate(dof, nbd);
	bdryvalue.allocate(dof, nbd);

	idnode_bdry(0, 0) = id_node(0, 0);
	idnode_bdry(0, 1) = id_node(nele - 1, 1);
	for (int j = 0; j < dof; ++j)
	{
		bdrytype(j, 0) = Fix;
		bdrytype(j, 1) = ConstForc;
		bdryvalue(j, 1) = 1;
	}

	// External force
	Fp.allocate(dof*nnode_all);


	// material condition
	rho_b.allocate(ni);
	rho_b.setvalue(1.0);

	iflap.allocate(ni);
	iflap.setvalue(1.0);

	eiflap.allocate(ni);
	myTYPE E = 1.0;
	eiflap = iflap*E;

	DM.allocate(dof, dof, nnode_all);
	MM.allocate(dof, dof, nnode_all);
	myTYPE grid;
	for (int i = nnode_all - 1; i >= 0; --i) {
		grid = nodegrid(i);
		MM(0, 0, i) = rho_b.interplinear_fast(ristation, grid);
		MM(1, 1, i) = iflap.interplinear_fast(ristation, grid);

		//DM(0, 0, i) = eiflap.interplinear_fast(ristation, grid);
		DM(1, 1, i) = eiflap.interplinear_fast(ristation, grid);
	}
#ifdef OUTPUT_MODE1
	cout << endl << endl;
	MM.output();
	DM.output();
	cout << endl << endl;
#endif // OUTPUT_MODE1


}




void Blade::functest(void)
{
	Matrix2<myTYPE> _Ma, _Ka;

	CantileverModel();
	//BladeModel_496();

#ifdef OUTPUT_MODE1
	Matrix3<myTYPE> Mtemp(2 * dof, 2 * dof, nele), Ktemp(2 * dof, 2 * dof, nele);
	EleKMgenrt(Mtemp, Ktemp);
	Mtemp.output("Mtemp.output", 4);
	Ktemp.output("Ktemp.output", 4);
#endif // OUTPUT_MODE1
		
	AssembleKM(Ma, Ka);
#ifdef OUTPUT_MODE1
	SpM2Mtr(_Ma, Ma, dof*nnode_all, dof*nnode_all);
	SpM2Mtr(_Ka, Ka, dof*nnode_all, dof*nnode_all);
	_Ma.output("_Ma1.output", 4);
	_Ka.output("_Ka1.output", 4);

	Ma.Output("SpM1.output");
	Ka.Output("SpK1.output");
#endif // OUTPUT_MODE1

	SetBdryCond();
	
#ifdef OUTPUT_MODE1
	SpM2Mtr(_Ma, Ma, dof*nnode_all, dof*nnode_all);
	SpM2Mtr(_Ka, Ka, dof*nnode_all, dof*nnode_all);
	_Ma.output("_Ma2.output", 4);
	_Ka.output("_Ka2.output", 4);

	Ma.Output("SpM2.output");
	Ka.Output("SpK2.output");
#endif // OUTPUT_MODE1

	StaticSolving();

}