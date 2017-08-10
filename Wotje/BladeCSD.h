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


class Blade {
private:
	int nf, ni;
	myTYPE omega;
	Matrix1<myTYPE> ristation, rho_b, eiflap, eilag, gj, iflap, ilag, ipitch, cent;

private:
	template <class Type> void _eleDMgenrt(Matrix3<Type> &mele, Matrix3<Type> &dele, const Matrix1<Type> elegrid, int dof);


	template <class Type> void _eleDMgenrt(Matrix2<Type> &mele, Matrix2<Type> &dele, const Type elegrid, int dof) {

		mele.allocate(dof, dof);
		dele.allocate(dof, dof);
		if (dof != 5) { print_wrong_msg("Wrong DOF defined."); }

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
	}

	
	template <class Type> void _ninterp_1d(Type &n, Type &nd, Type &ndd, Type ipt, int id) {
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
			print_wrong_msg("Wrong id gotten in nmatrix_1d.");

		}
			

	}
	
	
	template <class Type> void _nmatrix_1d(Matrix3<Type> &N, Matrix3<Type> &Nd, Matrix1<Type> &gauss, const int dof) {
		N.allocate(dof, 2 * dof, gauss.Nv);
		Nd.allocate(dof, 2 * dof, gauss.Nv);
		Type temp, dtemp, ddtemp, node;
		temp = dtemp = ddtemp = node = 0;
		
		if (dof != 5) { print_wrong_msg("Undefined DOF for N and Nd matrix."); }
		for (int i = gauss.Nv - 1; i >= 0; --i) {
			node = gauss(i);

			for (int j = 1; j < 3; ++j) {
				_ninterp_1d(temp, dtemp, ddtemp, node, 2*j-1);
				N(0, 5*j-5, i) = temp;
				N(1, 5*j-4, i) = temp;
				N(2, 5*j-4, i) = dtemp;
				N(3, 5*j-5, i) = -dtemp;
				Nd(0, 5*j-5, i) = dtemp;
				Nd(1, 5*j-4, i) = dtemp;
				Nd(2, 5*j-4, i) = ddtemp;
				Nd(3, 5*j-5, i) = -ddtemp;

				_ninterp_1d(temp, dtemp, ddtemp, node, 2*j);
				N(0, 5*j-2, i) = -temp;
				N(1, 5*j-3, i) = temp;
				N(2, 5*j-3, i) = dtemp;
				N(3, 5*j-2, i) = dtemp;
				Nd(0, 5*j-2, i) = -dtemp;
				Nd(1, 5*j-3, i) = dtemp;
				Nd(2, 5*j-3, i) = ddtemp;
				Nd(3, 5*j-2, i) = ddtemp;
			}
			N(4, 4, i) = 1 - node;
			N(4, 9, i) = node;
			Nd(4, 4, i) = -1;
			Nd(4, 9, i) = 1;
		}
	}


	/*template <class Type> void _nmatrix_1d(Matrix2<Type> &N, Matrix2<Type> &Nd, const Type &gauss, const int dof) {
		;
	}
*/

	template <class Type> void _gauss_intg(Matrix1<Type> &_w, Matrix1<Type> &_gp, const int ngs) {
		_w.allocate(ngs);
		_gp.allocate(ngs);

		switch (ngs) {
		case 1:
			_gp(0) = 0.0;
			_w(0) = 2.0;
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
			print_wrong_msg("Undefined NGS in _gauss_intg.");
			_gp(1) = 0.0;
			_w(1) = 2.0;
		}
	}


	template <class Type> void Mat2Arr(Type **A, Matrix2<Type> &B) {
		// B is col major
		for (int i = 0; i < B.NI; ++i) {
			*A = &B(i, 0) + i*B.NJ; 
			for (int j = 0; j < B.NJ; ++j) {
				*(*A + j) = B(i, j);
			}
		}
	}


public:
	Blade() {
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
		myTYPE dr = (ristation(ni-1) - ristation(0)) / (ni - 1);

		cent.allocate(ni);
		for (int i = ni - 1; i >= 0; --i) {
			idnn = step(i, ni - 1);
			cent(i) = omega * omega * (rho_b(idnn)*ristation(idnn)).sum() * dr;
		}
	}

	Blade(const Blade &B) {
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
	}

	~Blade() {
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
	}


	template <class Type> void EleKMgenrt(Matrix3<Type> &Mele, Matrix3<Type> &Kele, Matrix1<Type> elegrid, int dof) {
		int eleNum = elegrid.Nv - 1;
		int ngs = 3;
		Matrix1<Type> _weight(ngs), _gaussp(ngs);
		Matrix2<Type> mi(dof, dof), dm(dof, dof);
		Matrix3<Type> N(dof, 2 * dof, ngs), Nd(dof, 2 * dof, ngs);
		Matrix2<Type> temp(2 * dof, dof), temp_m(2 * dof, 2 * dof), temp_k(2 * dof, 2 * dof), temp_N(dof, 2 * dof);
		Type temp_grid;

		_gauss_intg(_weight, _gaussp, ngs);
		_gaussp.output(4);
		_gaussp = (_gaussp*0.5 + 0.5) / eleNum;
		_gaussp.output(4);

		_weight *= _weight * 0.5;
		
		_nmatrix_1d(N, Nd, _gaussp, dof);
		
		Mele.allocate(2 * dof, 2 * dof, eleNum);
		Kele.allocate(2 * dof, 2 * dof, eleNum);
		for (int i = eleNum-1; i >= 0; --i) {
			temp_grid = elegrid(i);
			temp_m.setvalue(0.0);
			temp_k.setvalue(0.0);
			for (int g = ngs - 1; g >= 0; --g) {
				_eleDMgenrt(mi, dm, temp_grid + _gaussp(g), dof);
				dm.output(4);
#ifdef USE_DOUBLE
				temp_N = N(2, g);
				//temp_N.output(4);
				//temp_N.output("N.output", 4);

				//Type **temptemp;
				//Mat2Arr(temptemp, temp_N);

				temp = temp_N.transpose().matrixmultiplyTP(mi);
				temp_m = temp.matrixmultiplyTP(temp_N);

				//cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*dof, dof, dof, 1.0, *temp_N_temp, 2*dof, *temp_mk, dof, 0.0, *temp_temp, dof);
				//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2*dof, 2*dof, dof, 1.0, *temp_temp, dof, *temp_N_temp, 2*dof, 0.0, *temp_mk_temp, 2*dof);
				
				//temp_m.output(4);

				temp_N = Nd(2, g);
				//temp_N.output(4);
				//temp_N.output("Nd.output", 4);

				temp = temp_N.transpose().matrixmultiplyTP(dm);
				temp_k = temp.matrixmultiplyTP(temp_N);

				//cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2*dof, dof, dof, 1.0, *temp_N_temp, 2*dof, *temp_mk, dof, 0.0, *temp_temp, dof);
				//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2*dof, 2*dof, dof, 1.0, *temp_temp, dof, *temp_N_temp, 2*dof, 0.0, *temp_mk_temp, 2*dof);


#else
				cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, dof, dof, 2 * dof, 1.0, N.v_p, 2 * dof, mi.v_p, 2 * dof, 0.0, temp.v_p, dof);
				cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2 * dof, dof, dof, 1.0, temp.v_p, 2 * dof, N.v_p, dof, 0.0, temp_m.v_p, 2 * dof);

				cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, dof, dof, 2 * dof, 1.0, Nd.v_p, 2 * dof, dm.v_p, 2 * dof, 0.0, temp.v_p, dof);
				cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2 * dof, dof, dof, 1.0, temp.v_p, 2 * dof, Nd.v_p, dof, 0.0, temp_k.v_p, 2 * dof);
#endif // USE_DOUBLE

				
				temp_m += temp_m*_weight(g);
				temp_k += temp_k*_weight(g);
				//temp_m.output(4);

				cout << "*****************************************" << endl;
			}

			/*Mele(:, : , i) = temp_m;
			Kele(:, : , i) = temp_k;*/
			for (int j = 2*dof - 1; j >= 0; --j) {
				for (int k = 2*dof - 1; k >= 0; --k) {
					Mele(k, j, i) = temp_m(k, j);
					Kele(k, j, i) = temp_k(k, j);
				}
			}
		}
	}



	void functest(void) {
		int dof = 5;
		int num_gauss = 1;
		Matrix1<myTYPE> elegrid;
		Matrix3<myTYPE> mele, kele;
		//Matrix3<myTYPE> N(dof, 2 * dof, num_gauss), Nd(dof, 2 * dof, num_gauss);
		//Matrix1<myTYPE> gauss(num_gauss), weight(num_gauss);

		elegrid.input("blade_grid.txt");
		//elegrid.input("ristation.txt");
		mele.allocate(2 * dof, 2 * dof, elegrid.Nv - 1);
		kele.allocate(2 * dof, 2 * dof, elegrid.Nv - 1);
		
		EleKMgenrt(mele, kele, elegrid, dof);
		mele.output("Mele.output", 4);
		kele.output("Kele.output", 4);
		
		/*_gauss_intg(gauss, weight, num_gauss);

		_nmatrix_1d(N, Nd, gauss*0.5+0.5, dof);*/

		/*N.output("N.output", 4);
		Nd.output("Nd.output", 4);
		N.output(3);
		Nd.output(3);*/
		/*mele.allocate(dof, dof, elegrid.Nv);
		kele.allocate(dof, dof, elegrid.Nv);

		_eleDMgenrt(mele, kele, elegrid, dof);
		mele.output("mele.output", 4);
		kele.output("kele.output", 4);*/
	}
};


template <class Type> void _eleDMgenrt(Matrix3<Type> &mele, Matrix3<Type> &dele, const Matrix1<Type> elegrid, int dof) {
	// mele.allocate(dof, dof, elegrid.Nv);
	// dele.allocate(dof, dof, elegrid.Nv);
	Type grid = 0;
	int num = elegrid.Nv;
	if (dof != 5) { print_wrong_msg("Wrong DOF defined."); }
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

}
#endif // !BladeCSD_h

