#pragma once
#define BeamCSD_h

// header files
#include "MatrixTemplate.h"
#include "Algorithm.h"
#include "DebugHelper.h"
#include "mkl.h"
#include "mkl_vml.h"
//#include "ModelCase.h"



// macro
#define USE_DOUBLE
#define LARGEVALUE1 1.0e10
#define LARGEVALUE2 1.2e10

#define OUTPUT_MODE1


// enum
enum BDryType {
	Free = 0, Fix = 1, ConstDisp = 2, DistrDisp = 3, ConstPres = 4, DistrPres = 5,
	ConstForc = 6, DistrForc = 7,
};

enum ProbType {
	Static = 0, Dynamic = 1, Mode = 2,
};

enum FreeDomType {
	DispPlaneIn = 0, DispPlaneOut = 1, DispAxis = 2, BendPlaneIn = 3, BendPlaneOut = 4, RotaAxis = 5,
};


class ModeReduced {
public:

	int Nw, Nwmax;
	Matrix2<double> Vn;
	Matrix1<double> wn;
	double wmin, wmax;

	int Nrd;
	Matrix2<double> Krd, Mrd, Crd, KKrd;
	Matrix1<double> qrd, qrd1, dqrd, ddqrd, Frd;
	Matrix1<double> qrd0, dqrd0, ddqrd0;

	~ModeReduced() { ; }

};

class GenArf {
public:
	const double pho = 0.0;
	double af, am, bt, r;
	double ck, c0, c1, c2, c3, c4, c5;

	void InitCoef(double dt) {
		af = pho / (pho + 1.0);
		am = (2.0*pho - 1.0) / (pho + 1.0);
		bt = 4.0 / (1 - am + af) / (1 - am + af);
		r = 0.5 - am + af;

		ck = 1.0 - af;
		c0 = (1 - am)*bt / dt / dt;
		c1 = ck*r*bt / dt;
		c2 = dt*c0;
		c3 = c2*dt / 2.0 - 1.0;
		c4 = ck*r*bt - 1.0;
		c5 = ck*(r / 2.0*bt - 1.0)*dt;
	}

};


class ProbDefinition {
public:
	int dof, nbdry;
	ProbType probtype;
	Matrix1<FreeDomType> frdtype;
	Matrix2<BDryType> bdrytype;
	Matrix2<myTYPE> bdryvalue;
	Matrix2<myTYPE> bdryposition; // (nboudary, start point-end point)

public:
	ProbDefinition() { ; }
	~ProbDefinition() { ; }
	ProbDefinition(const ProbDefinition &P)
	{
		dof = P.dof;
		nbdry = P.nbdry;
		probtype = P.probtype;
		frdtype = P.frdtype;
		bdrytype = P.bdrytype;
		bdryvalue = P.bdryvalue;
		bdryposition = P.bdryposition;
	}

	void GetProb(void);
};


class ModelCase {
public:
	int ni;
	myTYPE length;
	Matrix3<myTYPE> DM, MM;
	Matrix1<myTYPE> ristation, rho_b, eiflap, eilag, gj, iflap, ilag, ipitch;

public:
	ProbDefinition probdef;

public:
	ModelCase() { ; }
	ModelCase(const ModelCase &M)
	{
		ni = M.ni;
		length = M.length;
		ristation = M.ristation;
		rho_b = M.rho_b;
		eiflap = M.eiflap;
		eilag = M.eilag;
		gj = M.gj;
		iflap = M.iflap;
		ilag = M.ilag;
		ipitch = M.ipitch;
		probdef = M.probdef;
		DM = M.DM;
		MM = M.MM;
	}
	~ModelCase() { ; }

	void GetModel(void);

	void GetProb(void);

	void CantileverModel(void);

};




class Beam_1D {
private:
	int dof, nGauss, nEle, nNodeperEle, nNodeperPart, nBdry;
	Matrix2<int> nodeID, nodeID_bdry;
	Matrix1<int> nodeID_seq;
	Matrix1<myTYPE> nodePosition;
	Matrix3<myTYPE> DM, MM, Me, Ke;
	SpMtrx<myTYPE> Ka, Ca, Ma;
	ProbType ptype;
	Matrix2<BDryType> bdrytype;
	Matrix2<myTYPE> bdryvalue;
	Matrix1<FreeDomType> frdtype;
	Matrix1<myTYPE> Fp, q;
	Matrix2<myTYPE> p0, p1;

	ModeReduced ModePart;
	GenArf GAf;
	SpMtrx<myTYPE> *Ka_ptr, *Ma_ptr, *Ca_ptr, *Karf_ptr;
	ModelCase model;
	myTYPE length;
	Matrix1<myTYPE> ristation;

public:
	Beam_1D() { ; }
	~Beam_1D() { ; }
	Beam_1D(const Beam_1D &B)
	{
		dof = B.dof, nGauss = B.nGauss, nEle = B.nEle, nNodeperEle = B.nNodeperEle;
		nNodeperPart = B.nNodeperPart, nBdry = B.nBdry;
		nodeID = B.nodeID, nodeID_bdry = B.nodeID_bdry, nodeID_seq = B.nodeID_seq;
		nodePosition = B.nodePosition;
		DM = B.DM, MM = B.MM, Me = B.Me, Ke = B.Ke;
		Ka = B.Ka, Ca = B.Ca, Ma = B.Ma;
		ptype = B.ptype, bdrytype = B.bdrytype, bdryvalue = B.bdryvalue, frdtype = B.frdtype;
		Fp = B.Fp, q = B.q, p0 = B.p0, p1 = B.p1;
		ModePart = B.ModePart;
		//GAf = B.GAf;
		model = B.model;
		Ka_ptr = B.Ka_ptr, Ma_ptr = B.Ma_ptr, Ca_ptr = B.Ca_ptr, Karf_ptr = B.Karf_ptr;
	}

	void InitBeam1D(ModelCase &M);
	void AssembleKM(void);
	bool ReviseF(void);
	bool SetBdryCond(void);
	void StaticSolution(void);
	void ModeSolution(void);
	void DynamicSolution(void);
	
	template<class _Ty> bool _eleKMgenrt_N2(Matrix3<_Ty> &Me, Matrix3<_Ty> &Ke);
	template<class _Ty> bool _Nmatrix_1d(Matrix3<_Ty> &N, Matrix3<_Ty> &Nd, Matrix1<_Ty> &gauss, Matrix1<FreeDomType> &Dof);
	template<class _Ty> bool _Nmatrix_1d_N2D2(Matrix3<_Ty> &N, Matrix3<_Ty> &Nd, Matrix1<_Ty> &gauss);
	template<class _Ty> bool _Nmatrix_1d_N2D5(Matrix3<_Ty> &N, Matrix3<_Ty> &Nd, Matrix1<_Ty> &gauss);
	template<class _Ty> void _Ninterp_1d(_Ty &n, _Ty &nd, _Ty &ndd, _Ty x, _Ty L, int id);
	template<class _Ty> void _GaussIntegral(Matrix1<_Ty> &w, Matrix1<_Ty> &gp, const int nGs);
	template<class _Ty> bool _eleDMgenrt(Matrix2<_Ty> &M, Matrix2<_Ty> &D, _Ty p);

	template<class _Ty> void _FixedDisp(Matrix1<_Ty> &F, const int nodeid, const int idof, ProbType prob);
	template<class _Ty> void _ConstDisp(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof);
	template<class _Ty> void _FixedDispF(Matrix1<_Ty> &F, const int nodeid, const int idof);
	template<class _Ty> void _ConstDispF(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof);
	template<class _Ty> void _ConstForceF(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof);
	template<class _Ty> void AssembleKM(SpMtrx<_Ty> &M, SpMtrx<_Ty> &K);
};


template<class _Ty> 
bool Beam_1D::_eleKMgenrt_N2(Matrix3<_Ty> &Me, Matrix3<_Ty> &Ke)
{
	Matrix1<_Ty> _w(nGauss), _gp(nGauss);
	Matrix3<_Ty> N(dof, 2 * dof, nGauss), Nd(dof, 2 * dof, nGauss);
	Matrix2<_Ty> temp(dof, 2 * dof), _Ke(2 * dof, 2 * dof), _Me(2 * dof, 2 * dof);
	Matrix2<_Ty> mm(dof, dof), dm(dof, dof);
	_Ty _nodeP;
	int nodeid;
	bool flg = true;

	_GaussIntegral(_w, _gp, nGauss);
	_gp = (_gp*0.5 + 0.5)*length / nEle; // [0, L]
	_w = _w * 0.5;

	flg &= _Nmatrix_1d(N, Nd, _gp, frdtype);

	Me.allocate(2 * dof, 2 * dof, nEle);
	Ke.allocate(2 * dof, 2 * dof, nEle);
	for (int i = 0; i < nEle; ++i)
	{
		nodeid = nodeID(i, 0); // 2 node, select left one
		_nodeP = nodePosition(nodeid - 1);
		_Ke.setvalue(0);
		_Me.setvalue(0);
		for (int g = 0; g < nGauss; ++g)
		{
			flg &= _eleDMgenrt(mm, dm, _nodeP + _gp(g));
			temp = N(2, g);
			_Me += temp.transpose().matrixmultiplyTP(mm).matrixmultiplyTP(temp)*_w(g);
			temp = Nd(2, g);
			_Ke += temp.transpose().matrixmultiplyTP(dm).matrixmultiplyTP(temp)*_w(g);
		}

		for (int k = 0; k < 2 * dof; ++k)
		{
			for (int j = 0; j < 2 * dof; ++j)
			{
				Me(k, j, i) = _Me(k, j)*length / nEle;
				Ke(k, j, i) = _Ke(k, j)*length / nEle;
			}
		}
	}
	
	
	// Me matrix check
	_Ty _sum = 0;
	for (int i = 0; i < nEle; ++i)
		for (int j = 0; j < 2 * dof; j += 2)
			for (int k = 0; k < 2 * dof; k += 2)
				_sum += Me(k, j, i);
	cout << _sum << endl;

	return flg;
}

template<class _Ty> 
void Beam_1D::_GaussIntegral(Matrix1<_Ty> &_w, Matrix1<_Ty> &_gp, const int nGs)
{
	_w.allocate(nGs);
	_gp.allocate(nGs);

	switch (nGs) {
	case 1:
		_gp(0) = 0.0;
		_w(0) = 2.0;
		break;
	case 2:
		_gp(0) = -0.577350269189626;
		_gp(1) = 0.577350269189626;
		_w(0) = 1.0;
		_w(1) = 1.0;
		break;
	case 3:
		_gp(0) = -0.774596669241483;
		_gp(1) = 0.0;
		_gp(2) = 0.774596669241483;
		_w(0) = 0.55555555555555;
		_w(1) = 0.88888888888888;
		_w(2) = 0.55555555555555;
		break;
	case(4):
		_gp(0) = -0.861136311594053;
		_gp(1) = -0.339981043584856;
		_gp(2) = 0.339981043584856;
		_gp(3) = 0.861136311594053;
		_w(0) = 0.347854845137454;
		_w(1) = 0.652145154862546;
		_w(2) = 0.652145154862546;
		_w(3) = 0.347854845137454;
		break;
	case(5):
		_gp(0) = -0.906179845938664;
		_gp(1) = -0.538469310105683;
		_gp(2) = 0.0;
		_gp(3) = 0.538469310105683;
		_gp(4) = 0.906179845938664;
		_w(0) = 0.236926885056189;
		_w(1) = 0.478628670499366;
		_w(2) = 0.568888888888889;
		_w(3) = 0.478628670499366;
		_w(4) = 0.236926885056189;
		break;
	case(6):
		_gp(0) = -0.932469514203152;
		_gp(1) = -0.6612093864662646;
		_gp(2) = -0.2386191860831968;
		_gp(3) = 0.2386191860831968;
		_gp(4) = 0.6612093864662646;
		_gp(5) = 0.932469514203152;
		_w(0) = 0.1713244923791709;
		_w(1) = 0.3607615730481379;
		_w(2) = 0.4679139345726913;
		_w(3) = 0.4679139345726913;
		_w(4) = 0.3607615730481379;
		_w(5) = 0.1713244923791709;
		break;
	case(12):
		_gp(0) = -0.981560634246732;
		_gp(1) = -0.904117256370452;
		_gp(2) = -0.7699026741943177;
		_gp(3) = -0.5873179542866143;
		_gp(4) = -0.3678314989981804;
		_gp(5) = -0.12523340851114688;
		_gp(6) = 0.12523340851114688;
		_gp(7) = 0.3678314989981804;
		_gp(8) = 0.5873179542866143;
		_gp(9) = 0.7699026741943177;
		_gp(10) = 0.904117256370452;
		_gp(11) = 0.981560634246732;
		_w(0) = 0.04717533638647547;
		_w(1) = 0.1069393259953637;
		_w(2) = 0.1600783285433586;
		_w(3) = 0.2031674267230672;
		_w(4) = 0.2334925365383534;
		_w(5) = 0.2491470458134027;
		_w(6) = 0.2491470458134027;
		_w(7) = 0.2334925365383534;
		_w(8) = 0.2031674267230672;
		_w(9) = 0.1600783285433586;
		_w(10) = 0.1069393259953637;
		_w(11) = 0.04717533638647547;
		break;
	default:
		printf("Undefined NGS in _gauss_intg.");
		_gp(1) = 0.0;
		_w(1) = 2.0;
		break;
	}
}

template<class _Ty> 
bool Beam_1D::_Nmatrix_1d(Matrix3<_Ty> &N, Matrix3<_Ty> &Nd, Matrix1<_Ty> &gauss, Matrix1<FreeDomType> &Dof)
{
	_Ty n[6], nd[6], ndd[6], s, _L;
	Matrix1<int> _dof(dof);
	int index;
	_L = length / nEle;
	for (int i = 0; i < dof; ++i)
		_dof(i) = static_cast<int>(Dof(i));

	for (int g = 0; g < nGauss; ++g)
	{
		s = gauss(g);
		for (int id = 0; id < 6; id++)
			_Ninterp_1d(n[id], nd[id], ndd[id], s, _L, id+1);

		for (int i = 0; i < dof; ++i)
		{
			switch (Dof(i))
			{
			case DispPlaneOut:				
				N(i, i, g) = n[0];
				N(i, i + dof, g) = n[2];
				Nd(i, i, g) = nd[0];
				Nd(i, i + dof, g) = nd[2];
				index = BiSearch(_dof.v_p, 0, dof - 1, BendPlaneOut);
				if (index > -1)
				{
					N(i, index, g) = n[1];
					N(i, index + dof, g) = n[3];
					Nd(i, index, g) = nd[1];
					Nd(i, index + dof, g) = nd[3];
				}
				break;
			case DispPlaneIn:				
				N(i, i, g) = n[0];
				N(i, i + dof, g) = n[2];
				Nd(i, i, g) = nd[0];
				Nd(i, i + dof, g) = nd[2];
				index = BiSearch(_dof.v_p, 0, dof - 1, BendPlaneIn);
				if (index > -1)
				{
					N(i, index, g) = -n[1];
					N(i, index + dof, g) = -n[3];
					Nd(i, index, g) = -nd[1];
					Nd(i, index + dof, g) = -nd[3];
				}
				break;
			case BendPlaneOut:
				N(i, i, g) = nd[1];
				N(i, i + dof, g) = nd[3];
				Nd(i, i, g) = ndd[1];
				Nd(i, i + dof, g) = ndd[3];
				index = BiSearch(_dof.v_p, 0, dof - 1, DispPlaneOut);
				if (index > -1)
				{
					N(i, index, g) = nd[0];
					N(i, index + dof, g) = nd[2];
					Nd(i, index, g) = ndd[0];
					Nd(i, index + dof, g) = ndd[2];
				}
				break;
			case BendPlaneIn:
				N(i, i, g) = nd[1];
				N(i, i + dof, g) = nd[3];
				Nd(i, i, g) = ndd[1];
				Nd(i, i + dof, g) = ndd[3];
				index = BiSearch(_dof.v_p, 0, dof - 1, DispPlaneIn);
				if (index > -1)
				{
					N(i, index, g) = -nd[0];
					N(i, index + dof, g) = -nd[2];
					Nd(i, index, g) = -ndd[0];
					Nd(i, index + dof, g) = -ndd[2];
				}
				break;
			case RotaAxis:
				N(i, i, g) = 1 - s / _L;
				N(i, i + dof, g) = s / _L;
				Nd(i, i, g) = -1.0 / _L;
				Nd(i, i + dof, g) = 1.0 / _L;
				break;
			case DispAxis:
				printf("Undefined DOF DispAxis in _Nmatrix_1d(). \n");
				return false;
			}
		}
	}
	return true;
}

template<class _Ty> 
bool Beam_1D::_eleDMgenrt(Matrix2<_Ty> &M, Matrix2<_Ty> &D, _Ty p)
{
	_Ty a;
	int index;
	index = BiSearchFloor(ristation.v_p, 0, ristation.Nv - 1, p);
	if (index > -1 && index < ristation.Nv - 1)
	{
		a = (p - ristation(index)) / (ristation(index + 1) - ristation(index));
		for (int i = 0; i < dof; ++i)
		{
			for (int j = 0; j < dof; ++j)
			{
				M(j, i) = MM(j, i, index)*(1 - a) + MM(j, i, index + 1)*a;
				D(j, i) = DM(j, i, index)*(1 - a) + DM(j, i, index + 1)*a;
			}
		}
		return true;
	}
	else
		return false;

}

template<class _Ty> 
void Beam_1D::_Ninterp_1d(_Ty &n, _Ty &nd, _Ty &ndd, _Ty x, _Ty L, int id)
{
	_Ty ks = x / L, ks2 = ks*ks, ks3 = ks2*ks, L2 = L*L;
	switch (id)
	{
	case 1:
		n = 1 - 3 * ks2 + 2 * ks3;
		nd = 1.0 / L*(-6 * ks + 6 * ks2);
		ndd = 1.0 / L2*(-6 + 12 * ks);
		break;
	case 2:
		n = L*(ks - 2 * ks2 + ks3);
		nd = 1 - 4 * ks + 3 * ks2;
		ndd = 1.0 / L*(-4 + 6 * ks);
		break;
	case 3:
		n = 3 * ks2 - 2 * ks3;
		nd = 1.0 / L*(6 * ks - 6 * ks2);
		ndd = 1.0 / L2*(6 - 12 * ks);
		break;
	case 4:
		n = L*(-ks2 + ks3);
		nd = -2 * ks + 3 * ks2;
		ndd = 1.0 / L*(-2 + 6 * ks);
		break;
	default:
		printf("Undefined id %d of _Niterp_1d() Func. \n", id);
		n = nd = ndd = 0;
		break;
	}
}

template<class _Ty> 
void Beam_1D::AssembleKM(SpMtrx<_Ty> &M, SpMtrx<_Ty> &K)
{
	Matrix3<myTYPE> Me, Ke;
	if (!_eleKMgenrt_N2(Me, Ke))
	{
		M.allocate(1, 1);
		K.allocate(1, 1);
	}
	else
	{
		KK2SpM(M, Me, nodeID, 0);
		KK2SpM(K, Ke, nodeID, 0);
	}
#ifdef OUTPUT_MODE1
	Matrix2<myTYPE> _Ma, _Ka;
	SpM2Mtr(_Ma, M, dof*nNodeperPart, dof*nNodeperPart);
	SpM2Mtr(_Ka, K, dof*nNodeperPart, dof*nNodeperPart);
	_Ma.output("_Ma1.output", 8);
	_Ka.output("_Ka1.output", 8);

	M.Output("SpM1.output");
	K.Output("SpK1.output");
#endif // OUTPUT_MODE1
}

template<class _Ty> 
void Beam_1D::_FixedDisp(Matrix1<_Ty> &F, const int nodeid, const int idof, ProbType prob)
{
	int i, j;
	_Ty kam;

	j = (nodeid - 1)*dof + idof;
	kam = Ka(j, j)*LARGEVALUE1;

	for (i = 1; i <= nNodeperPart*dof; ++i) {
		Ka(i, j) = 0;
		if (prob != Static) Ma(i, j) = 0;
	}

	i = (nodeid - 1)*dof + idof;
	for (j = 1; j <= nNodeperPart*dof; ++j) {
		Ka(i, j) = 0;
		if (prob != Static) Ma(i, j) = 0;
	}

	Ka(i, i) = kam*kam;
	if (prob != Static) Ma(i, i) = 1;
	F(i - 1) = 0;
}

template<class _Ty> 
void Beam_1D::_ConstDisp(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof)
{
	int j;
	j = (nodeid - 1)*dof + idof;
	F(j - 1) = Ka(j, j)*LARGEVALUE2*disp;
	Ka(j, j) *= LARGEVALUE2;
}

template<class _Ty> 
void Beam_1D::_ConstDispF(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof)
{
	int j = dof*(nodeid - 1) + idof;
	//cout << Ka(j, j) << endl;
	F(j - 1) = Ka(j, j)*disp;
}

template<class _Ty> 
void Beam_1D::_FixedDispF(Matrix1<_Ty> &F, const int nodeid, const int idof)
{
	F(dof*(nodeid - 1) + idof - 1) = 0;
}

template<class _Ty> 
void Beam_1D::_ConstForceF(Matrix1<_Ty> &F, _Ty disp, const int nodeid, const int idof)
{
	int j = dof*(nodeid - 1) + idof;
	F(j - 1) += disp;
}