#pragma once
#ifndef BladeCSD_h
#define BladeCSD_h

// header files
#include "MatrixTemplate.h"
#include "DebugHelper.h"
#include "mkl.h"
#include "mkl_vml.h"


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





class Blade {
private:
	int nf, ni;
	int dof, ngs, nele, nnode, nnode_all, nbd;
	myTYPE omega;
	Matrix1<myTYPE> ristation, rho_b, eiflap, eilag, gj, iflap, ilag, ipitch, cent;
	Matrix1<myTYPE> nodegrid;
	Matrix2<int> id_node;
	Matrix3<myTYPE> DM, MM;
	SpMtrx<myTYPE> Ka, Ma;
	ProbType ptype;
	Matrix2<BDryType> bdrytype;
	Matrix2<myTYPE> bdryvalue;
	Matrix2<int> idnode_bdry;
	Matrix1<myTYPE> Fp, q;
	Matrix2<myTYPE> p0, p1;
	Matrix1<FreeDomType> frdtype;

private:
	template <class Type> inline bool _eleDMgenrt(Matrix3<Type> &mele, Matrix3<Type> &dele, const Matrix1<Type> elegrid, int dof);

	template <class Type> inline bool __eleDMgenrt(Matrix2<Type> &mele, Matrix2<Type> &dele, const Type elegrid, int dof);

	template <class Type> inline bool _eleDMgenrt(Matrix2<Type> &mele, Matrix2<Type> &dele, const Type elegrid, const int nodeid);

	template <class Type> inline bool _eleKMgenrt(Matrix3<Type> &Mele, Matrix3<Type> &Kele);
	
	template <class Type> inline bool _ninterp_1d(Type &n, Type &nd, Type &ndd, Type ipt, int id);
	
	template <class Type> inline bool _nmatrix_1d_node2dof1(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);
	
	template <class Type> inline bool _nmatrix_1d_node2dof2(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);
	
	template <class Type> inline bool _nmatrix_1d_node2dof3(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);

	template <class Type> inline bool _nmatrix_1d_node2dof4(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);

	template <class Type> inline bool _nmatrix_1d_node2dof5(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);

	template <class Type> inline bool _nmatrix_1d_node2dof6(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss);

	template <class Type> inline bool _nmatrix_1d(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss, const int dof);

	template <class Type> inline bool _nmatrix_1d(Matrix2<Type> &N, Matrix2<Type> &Nd, const Type &gauss, const int dof);

	template <class Type> inline void _gauss_intg(Matrix1<Type> &_w, Matrix1<Type> &_gp, const int ngs);

	template <class Type> inline void Mat2Arr(Type **A, Matrix2<Type> &B);

	template <class Type> inline void _fixeddisp(Matrix1<Type>& F, const int id_nd, const int idf, ProbType prob);

	template <class Type> inline void _constdisp(Matrix1<Type>& F, Type disp, const int id_nd, const int idf);

	template <class Type> inline void _fixeddispF(Matrix1<Type>&F, const int id_nd, const int idf);

	template <class Type> inline void _constdispF(Matrix1<Type>& F, Type disp, const int id_nd, const int idf);

	inline void _extractbdy(void);

	inline void _fixeddispF(void);

	inline void _setbdrycond(void);

public:
	Blade();

	Blade(const Blade &B);

	~Blade();

	void EleKMgenrt(Matrix3<myTYPE> &Mele, Matrix3<myTYPE> &Kele);

	void AssembleKM(SpMtrx<myTYPE> &M, SpMtrx<myTYPE> &K);

	void AssembleKM(void);

	void SetBdryCond(void);

	void StaticSolving(void);

	void ModeSolving(void);

	void BladeModel_496(void);

	void BladeModel_UH60A(void);

	void CantileverModel(void);

	void functest(void);
};



template<class Type>
inline bool Blade::_eleDMgenrt(Matrix3<Type>& mele, Matrix3<Type>& dele, const Matrix1<Type> elegrid, int dof)
{
	// mele.allocate(dof, dof, elegrid.Nv);
	// dele.allocate(dof, dof, elegrid.Nv);
	Type grid = 0;
	int num = elegrid.Nv;
	if (dof != 5) { 
		printf("Wrong DOF defined.");
		return false;
	}
	for (int i = num - 1; i >= 0; --i) {
		grid = elegrid.v_p[i];
		mele(0, 0, i) = rho_b.interplinear_fast(ristation, grid);
		mele(1, 1, i) = mele(0, 0, i);
		mele(2, 2, i) = iflap.interplinear_fast(ristation, grid);
		mele(3, 3, i) = ilag.interplinear_fast(ristation, grid);
		mele(4, 4, i) = ipitch.interplinear_fast(ristation, grid);

		dele(0, 0, i) = cent.interplinear_fast(ristation, grid);
		dele(1, 1, i) = kele(0, 0, i);
		dele(2, 2, i) = eiflap.interplinear_fast(ristation, grid);
		dele(3, 3, i) = eilag.interplinear_fast(ristation, grid);
		dele(4, 4, i) = gj.interplinear_fast(ristation, grid);
	}
	return true;
}


template <class Type> 
inline bool Blade::__eleDMgenrt(Matrix2<Type> &mele, Matrix2<Type> &dele, const Type elegrid, int dof)
{

	mele.allocate(dof, dof);
	dele.allocate(dof, dof);
	if (dof != 5) { 
		printf("Wrong DOF defined."); 
		return false;
	}

	mele(0, 0) = rho_b.interplinear_fast(ristation, elegrid);
	mele(1, 1) = mele(0, 0);
	mele(2, 2) = iflap.interplinear_fast(ristation, elegrid);
	mele(3, 3) = ilag.interplinear_fast(ristation, elegrid);
	mele(4, 4) = ipitch.interplinear_fast(ristation, elegrid);

	dele(0, 0) = cent.interplinear_fast(ristation, elegrid);
	dele(1, 1) = dele(0, 0);
	dele(2, 2) = eiflap.interplinear_fast(ristation, elegrid);
	dele(3, 3) = eilag.interplinear_fast(ristation, elegrid);
	dele(4, 4) = gj.interplinear_fast(ristation, elegrid);

	return true;
}


template <class Type>
inline bool Blade::_eleDMgenrt(Matrix2<Type> &mele, Matrix2<Type> &dele, const Type elegrid, const int nodeid)
{
	//int indexL = 0, indexR = nnode_all - 1, indexM;
	Type atemp;
	//for (; indexL <= indexR;) {
	//	indexM = indexL + (indexR - indexL >> 1);
	//	if (ristation(indexM) == elegrid)
	//		break;
	//	else {
	//		if (ristation(indexM) > elegrid)
	//			indexR = indexM - 1;
	//		else
	//			indexL = indexM + 1;
	//	}
	//}
	//indexL = Max(Min(indexR, indexL), 0);
	//indexR = Min(nnode_all - 1, indexL + 1);
	atemp = (elegrid - nodegrid(nodeid-1)) / (nodegrid(nodeid) - nodegrid(nodeid-1));
	
	for (int i = 0; i < dof; ++i)
	{
		for (int j = 0; j < dof; ++j) {
			mele(j,i) = MM(j, i, nodeid - 1) * (1 - atemp) + MM(j, i, nodeid) * atemp;
			dele(j,i) = DM(j, i, nodeid - 1) * (1 - atemp) + DM(j, i, nodeid) * atemp;
		}
	}
			
	return true;
}


template <class Type>
inline bool Blade::_eleKMgenrt(Matrix3<Type> &Mele, Matrix3<Type> &Kele)
{
	Matrix1<Type> _weight(ngs), _gaussp(ngs);
	Matrix2<Type> mi(dof, dof), dm(dof, dof);
	Matrix3<Type> N(dof, 2 * dof, ngs), Nd(dof, 2 * dof, ngs);
	Matrix2<Type> temp(2 * dof, dof), temp_m(2 * dof, 2 * dof), temp_k(2 * dof, 2 * dof), temp_N(dof, 2 * dof);
	Type temp_grid;
	int nodeid;
	bool flg = true;

	_gauss_intg(_weight, _gaussp, ngs);
	_gaussp = (_gaussp*0.5 + 0.5) / nele;

	_weight = _weight * 0.5;

	flg &= _nmatrix_1d(N, Nd, _gaussp, dof);
	N.output("N.output", 4);
	Nd.output("Nd.output", 4);

	Mele.allocate(2 * dof, 2 * dof, nele);
	Kele.allocate(2 * dof, 2 * dof, nele);
	for (int i = nele - 1; i >= 0; --i) {
		nodeid = id_node(i, 0);
		temp_grid = nodegrid(nodeid-1);
		temp_m.setvalue(0.0);
		temp_k.setvalue(0.0);
		for (int g = ngs - 1; g >= 0; --g) {
			//flg &= _eleDMgenrt(mi, dm, temp_grid + _gaussp(g), dof);
			cout << "Gauss Point Position: " << temp_grid + _gaussp(g) << endl;
			flg &= _eleDMgenrt(mi, dm, temp_grid + _gaussp(g), nodeid);
			//dm.output(4);

			temp_N = N(2, g);
			//temp_N.output(4);
			//temp_N.output("N.output", 4);

			temp = temp_N.transpose().matrixmultiplyTP(mi);
			temp_m += temp.matrixmultiplyTP(temp_N) * _weight(g);				
			//temp.output();
			//temp_m.output();

			temp_N = Nd(2, g);
			//temp_N.output(4);
			//dm.output(4);
			//temp_N.output("Nd.output", 4);

			temp = temp_N.transpose().matrixmultiplyTP(dm);
			temp_k += temp.matrixmultiplyTP(temp_N) * _weight(g);
			
			//temp_k.output();

			cout << "*****************************************" << endl;
		}
		temp_m.output();
		temp_k.output();

		/*Mele(:, : , i) = temp_m;
		Kele(:, : , i) = temp_k;*/
		for (int j = 2 * dof - 1; j >= 0; --j) {
			for (int k = 2 * dof - 1; k >= 0; --k) {
				Mele(k, j, i) = temp_m(k, j);
				Kele(k, j, i) = temp_k(k, j);
			}
		}
	}
	return flg;
}


template <class Type> 
inline bool Blade::_ninterp_1d(Type &n, Type &nd, Type &ndd, Type ipt, int id)
{
	switch (id) {
	case 1:
		//n = 1 - 3.0 * ipt*ipt + 2.0*ipt*ipt*ipt;
		n = 1 - ipt*ipt*(3 - 2 * ipt);
		//nd = -6 * ipt + 6 * ipt*ipt;
		nd = -6 * ipt * (1 - ipt);
		ndd = -6 + 12 * ipt;
		break;
	case 2:
		//n = ipt - 2 * ipt*ipt + ipt*ipt*ipt;
		n = ipt * (1 - ipt * (2 - ipt));
		//nd = 1 - 4 * ipt + 3 * ipt*ipt;
		nd = 1 - ipt*(4 - 3 * ipt);
		ndd = -4 + 6 * ipt;
		break;
	case 3:
		//n = 3 * ipt*ipt - 2 * ipt*ipt*ipt;
		n = ipt*ipt*(3 - 2 * ipt);
		nd = 6 * ipt *(1 - ipt);
		ndd = 6 - 12 * ipt;
		break;
	case 4:
		n = -ipt*ipt*(1 - ipt);
		//nd = -2 * ipt + 3 * ipt*ipt;
		nd = -ipt*(2 - 3 * ipt);
		ndd = -2 + 6 * ipt;
		break;
	default:
		printf("Wrong id gotten in nmatrix_1d.");
		return false;
	}
	return true;
}


template <class Type>
inline bool Blade::_nmatrix_1d_node2dof5(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss)
{
	N.allocate(dof, 2 * dof, gauss.Nv);
	Nd.allocate(dof, 2 * dof, gauss.Nv);
	Type temp, dtemp, ddtemp, node;
	bool flg = false;
	temp = dtemp = ddtemp = node = 0;

	for (int i = 0; i < dof; ++i)
		flg |= (frdtype(i) == DispAxis);

	if (!flg)
	{
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);

			for (int j = 1; j < 3; ++j) {
				_ninterp_1d(temp, dtemp, ddtemp, node, 2 * j - 1);
				N(0, 5 * j - 5, i) = temp;
				N(1, 5 * j - 4, i) = temp;
				N(2, 5 * j - 4, i) = dtemp;
				N(3, 5 * j - 5, i) = -dtemp;
				Nd(0, 5 * j - 5, i) = dtemp;
				Nd(1, 5 * j - 4, i) = dtemp;
				Nd(2, 5 * j - 4, i) = ddtemp;
				Nd(3, 5 * j - 5, i) = -ddtemp;

				_ninterp_1d(temp, dtemp, ddtemp, node, 2 * j);
				N(0, 5 * j - 2, i) = -temp;
				N(1, 5 * j - 3, i) = temp;
				N(2, 5 * j - 3, i) = dtemp;
				N(3, 5 * j - 2, i) = dtemp;
				Nd(0, 5 * j - 2, i) = -dtemp;
				Nd(1, 5 * j - 3, i) = dtemp;
				Nd(2, 5 * j - 3, i) = ddtemp;
				Nd(3, 5 * j - 2, i) = ddtemp;
			}
			N(4, 4, i) = 1 - node;
			N(4, 9, i) = node;
			Nd(4, 4, i) = -1;
			Nd(4, 9, i) = 1;
		}
	}
	return (!flg);
}


template <class Type>
inline bool Blade::_nmatrix_1d_node2dof2(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss)
{
	N.allocate(dof, 2 * dof, gauss.Nv);
	Nd.allocate(dof, 2 * dof, gauss.Nv);
	Type temp, dtemp, ddtemp, node;
	temp = dtemp = ddtemp = node = 0;


	switch (frdtype(0))
	{
	case DispPlaneOut:
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);

			_ninterp_1d(temp, dtemp, ddtemp, node, 1);
			N(0, 0, i) = temp;
			N(1, 0, i) = dtemp;
			Nd(0, 0, i) = dtemp;
			Nd(1, 0, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 2);
			N(0, 1, i) = temp;
			N(1, 1, i) = dtemp;
			Nd(0, 1, i) = dtemp;
			Nd(1, 1, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 3);
			N(0, 2, i) = temp;
			N(1, 2, i) = dtemp;
			Nd(0, 2, i) = dtemp;
			Nd(1, 2, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 4);
			N(0, 3, i) = temp;
			N(1, 3, i) = dtemp;
			Nd(0, 3, i) = dtemp;
			Nd(1, 3, i) = ddtemp;
		}
		break;
	
	case DispPlaneIn:
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);
			for (int j = 1; j < 3; ++j) {
				node = gauss(i);
				_ninterp_1d(temp, dtemp, ddtemp, node, 1);
				N(0, 0, i) = temp;
				N(1, 0, i) = -dtemp;
				Nd(0, 0, i) = dtemp;
				Nd(1, 0, i) = -ddtemp;

				_ninterp_1d(temp, dtemp, ddtemp, node, 2);
				N(0, 1, i) = -temp;
				N(1, 1, i) = dtemp;
				Nd(0, 1, i) = -dtemp;
				Nd(1, 1, i) = ddtemp;

				_ninterp_1d(temp, dtemp, ddtemp, node, 3);
				N(0, 2, i) = temp;
				N(1, 2, i) = -dtemp;
				Nd(0, 2, i) = dtemp;
				Nd(1, 2, i) = -ddtemp;

				_ninterp_1d(temp, dtemp, ddtemp, node, 4);
				N(0, 3, i) = -temp;
				N(1, 3, i) = dtemp;
				Nd(0, 3, i) = -dtemp;
				Nd(1, 3, i) = ddtemp;
			}
		}
		break;

	case BendPlaneOut:
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);

			_ninterp_1d(temp, dtemp, ddtemp, node, 1);
			N(1, 0, i) = temp;
			N(0, 0, i) = dtemp;
			Nd(1, 0, i) = dtemp;
			Nd(0, 0, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 2);
			N(1, 1, i) = temp;
			N(0, 1, i) = dtemp;
			Nd(1, 1, i) = dtemp;
			Nd(0, 1, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 3);
			N(1, 2, i) = temp;
			N(0, 2, i) = dtemp;
			Nd(1, 2, i) = dtemp;
			Nd(0, 2, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 4);
			N(1, 3, i) = temp;
			N(0, 3, i) = dtemp;
			Nd(1, 3, i) = dtemp;
			Nd(0, 3, i) = ddtemp;
		}
		break;

	case BendPlaneIn:
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);

			_ninterp_1d(temp, dtemp, ddtemp, node, 1);
			N(1, 0, i) = temp;
			N(0, 0, i) = -dtemp;
			Nd(1, 0, i) = dtemp;
			Nd(0, 0, i) = -ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 2);
			N(1, 1, i) = -temp;
			N(0, 1, i) = dtemp;
			Nd(1, 1, i) = -dtemp;
			Nd(0, 1, i) = ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 3);
			N(1, 2, i) = temp;
			N(0, 2, i) = -dtemp;
			Nd(1, 2, i) = dtemp;
			Nd(0, 2, i) = -ddtemp;

			_ninterp_1d(temp, dtemp, ddtemp, node, 4);
			N(1, 3, i) = -temp;
			N(0, 3, i) = dtemp;
			Nd(1, 3, i) = -dtemp;
			Nd(0, 3, i) = ddtemp;
		}
		break;

	default:
		printf("Unfined Dof Type in Func _node2dof2. \n");
		return false;
	}
	return true;
}


template <class Type>
inline bool Blade::_nmatrix_1d(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss, const int dof)
{
	switch (dof)
	{
	case 2:
		_nmatrix_1d_node2dof2(N, Nd, gauss);
		break;
	case 5:
		_nmatrix_1d_node2dof5(N, Nd, gauss);
		break;
	default:
		printf("Undefined Dof %s in func _nmatrix_1d(). \n");
		return false;
	}
	return true;
}


template <class Type> 
inline void Blade::_gauss_intg(Matrix1<Type> &_w, Matrix1<Type> &_gp, const int ngs) 
{
	_w.allocate(ngs);
	_gp.allocate(ngs);

	switch (ngs) {
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


template <class Type> 
inline void Blade::Mat2Arr(Type **A, Matrix2<Type> &B) 
{
	// B is col major
	for (int i = 0; i < B.NI; ++i) {
		*A = &B(i, 0) + i*B.NJ;
		for (int j = 0; j < B.NJ; ++j) {
			*(*A + j) = B(i, j);
		}
	}
}


template<class Type>
inline void Blade::_fixeddisp(Matrix1<Type>& F, const int id_nd, const int idf, ProbType prob)
{
	int i, j;
	Type kam;

	j = (id_nd-1)*dof + idf;
	kam = Ka(j, j)*LARGEVALUE1;

	for (i = 1; i <= nnode_all*dof; ++i) {
		Ka(i, j) = 0;
		if (prob != Static) Ma(i, j) = 0;
	}

	i = (id_nd - 1)*dof + idf;
	for (j = 1; j <= nnode_all*dof; ++j) {
		Ka(i, j) = 0;
		if (prob != Static) Ma(i, j) = 0;
	}

	Ka(i, i) = kam*kam;
	if (prob != Static) Ma(i, i) = 1;
	F(i-1) = 0;
}


template<class Type>
inline void Blade::_constdisp(Matrix1<Type>& F, Type disp, const int id_nd, const int idf)
{
	int j;
	j = (id_nd - 1)*dof + idf;
	F(j-1) = Ka(j, j)*LARGEVALUE2*disp;
	Ka(j, j) *= LARGEVALUE2;
}


template<class Type>
inline void Blade::_fixeddispF(Matrix1<Type>& F, const int id_nd, const int idf)
{
	F(dof*(id_nd - 1) + idf-1) = 0;
}


template<class Type>
inline void Blade::_constdispF(Matrix1<Type>& F, Type disp, const int id_nd, const int idf)
{
	int j = dof*(id_nd - 1) + idf;
	//cout << Ka(j, j) << endl;
	F(j-1) = Ka(j, j)*disp;
}


inline void Blade::_fixeddispF(void)
{
	int j;
	for (int ib = 0; ib < nbd; ++ib) {
		for (int idof = 1; idof <= dof; ++idof) {
			switch (bdrytype(idof - 1, ib))
			{
			case Fix:
				_fixeddispF(Fp, idnode_bdry(0, ib), idof);
				//j = dof*(idnode_bdry(0, ib) - 1) + idof;
				//Fp(j) = 0;
				break;
			case ConstDisp:
				_constdispF(Fp, bdryvalue(idof - 1, ib), idnode_bdry(0, ib), idof);
				//j = dof*(idnode_bdry(0, ib) - 1) + idof;
				//Fp(j) = Ka(j, j)*bdryvalue(idof - 1, ib);
				break;
			}
		}
	}
}


inline void Blade::_extractbdy(void)
{
	idnode_bdry(0, 0) = id_node(0, 0);
	idnode_bdry(0, 1) = id_node(nele - 1, nnode - 1);
	for (int i = 0; i < dof; ++i)
		bdrytype(i, 0) = Fix;
	for (int i = 0; i < dof; ++i)
		bdrytype(i, 1) = Free;
	
	bdryvalue.setvalue(0);
}


inline void Blade::_setbdrycond(void)
{
	//_extractbdy();
	Fp.setvalue(0);

	for (int ib = 0; ib < nbd; ++ib) {
		for (int idof = 1; idof <= dof; ++idof) {
			switch (bdrytype(idof-1, ib)) 
			{
			case Free:
				break;
			case Fix:
				_fixeddisp(Fp, idnode_bdry(0, ib), idof, ptype);
				break;
			case ConstDisp:
				_constdisp(Fp, bdryvalue(idof - 1, ib), idnode_bdry(0, ib), idof);
				break;
			case ConstForc:
				_constdispF(Fp, bdryvalue(idof - 1, ib), idnode_bdry(0, ib), idof);
				break;
			default:
				break;
			}
		}
	}

	Fp.output("Fp.output", 4);
	//_fixeddispF();
}



#endif // !BladeCSD_h