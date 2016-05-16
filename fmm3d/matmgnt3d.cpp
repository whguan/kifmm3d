/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#include "matmgnt3d.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

// ---------------------------------------------------------------------- 
double MatMgnt3d::_wsbuf[16384];

// ---------------------------------------------------------------------- 
MatMgnt3d::MatMgnt3d(): _np(6)
{
}
// ---------------------------------------------------------------------- 
MatMgnt3d::~MatMgnt3d()
{
}

double MatMgnt3d::alt()
{
  return pow(0.1, _np+1);
}
// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
int MatMgnt3d::setup()
{
  //--------------------------------------------------------
  _hom = _knl.homogeneous();
  if(_hom==true) {
	 _knl.homogeneousDeg(_degVec); iA(_degVec.size()==srcDOF());
  }
  iC( samPosCal(_np,   1.0, _samPos[UE]) );
  iC( samPosCal(_np+2, 3.0, _samPos[UC]) );
  iC( samPosCal(_np,   3.0, _samPos[DE]) );
  iC( samPosCal(_np,   1.0, _samPos[DC]) );
  
  iC( regPosCal(_np,   1.0, _regPos    ) ); //only one regPos
  //--------------------------------------------------------
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::report()
{
  cerr << "matrix map size"<<endl;
  cerr << _upwEqu2UpwChk.size() << " ";
  cerr << _upwChk2UpwEqu.size() << " ";
  cerr << _dwnChk2DwnEqu.size() << " ";
  cerr << _dwnEqu2DwnChk.size() << " ";
  cerr << _upwEqu2DwnChk.size() << endl;
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::plnDatSze(int tp)
{
  if(tp==UE || tp==DE)
	 return _samPos[tp].n()*srcDOF();
  else
	 return _samPos[tp].n()*trgDOF();
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effDatSze(int tp)
{
  //int et = eq().et();
  int effNum = (2*_np+2)*(2*_np)*(2*_np);
  if(tp==UE || tp==DE)
	 return effNum*srcDOF();
  else
	 return effNum*trgDOF();
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwChk2UpwEqu_dgemv(int l, const DblNumVec& chk, DblNumVec& den)
{
  DblNumMat& _UC2UE = (_hom==true) ? _upwChk2UpwEqu[0] : _upwChk2UpwEqu[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0,l);
  //---------compute matrix
  if(_UC2UE.m()==0) {	 //cerr<<"UpwChk2UpwEqu compute"<<endl;
	 //set matrix
	 DblNumMat ud2c(plnDatSze(UC), plnDatSze(UE));
	 DblNumMat chkPos(dim(),samPos(UC).n());	 clear(chkPos);	 iC( daxpy(R, samPos(UC), chkPos) ); //scale
	 DblNumMat denPos(dim(),samPos(UE).n());	 clear(denPos);	 iC( daxpy(R, samPos(UE), denPos) ); //scale
	 
	 iC( _knl.kernel(denPos, denPos, chkPos, ud2c) );
	 _UC2UE.resize(plnDatSze(UE), plnDatSze(UC));
	 iC( pinv(ud2c, alt(), _UC2UE) );
  }
  //---------matvec
  if(_hom==true) {
	 //matvec
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(UE).n(), false, _wsbuf);	 clear(tmpDen);
	 iC( dgemv(1.0, _UC2UE, chk, 1.0, tmpDen) );
	 //scale
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, - l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i < samPos(UE).n(); i++)
		for(int s=0; s < srcDOF; s++) {
		  den(cnt) = den(cnt) + tmpDen(cnt) * sclvec[s];
		  cnt++;
		}
  } else {
	 iC( dgemv(1.0, _UC2UE, chk, 1.0, den) );
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwEqu2UpwChk_dgemv(int l, Index3 idx, const DblNumVec& den, DblNumVec& chk)
{
  NumTns<DblNumMat>& _UE2UC = (_hom==true) ? _upwEqu2UpwChk[0] : _upwEqu2UpwChk[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0, l);  
  if(_UE2UC.m()==0)	 _UE2UC.resize(2,2,2);
  DblNumMat& _UE2UCii = _UE2UC(idx(0), idx(1), idx(2));
  //---------compute matrix
  if(_UE2UCii.m()==0) {	 //cerr<<"UpwEqu2UpwChk compute"<<endl;
	 _UE2UCii.resize(plnDatSze(UC), plnDatSze(UE)); //_memused[1] += plnnum(UC)*dof()*plnnum(UE)*dof()*sizeof(double);
	 DblNumMat chkPos(dim(),samPos(UC).n());	 clear(chkPos);	 iC( daxpy(2.0*R, samPos(UC), chkPos) ); //scale
	 DblNumMat denPos(dim(),samPos(UE).n());	 clear(denPos);	 iC( daxpy(R, samPos(UE), denPos) ); //scale
	 for(int i=0; i<dim(); i++) for(int j=0; j<samPos(UE).n(); j++)	denPos(i,j) = denPos(i,j) + (2*idx(i)-1)*R;//shift
	 
	 iC( _knl.kernel(denPos, denPos, chkPos, _UE2UCii) );
  }
  //---------matvec
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(UE).n(), false, _wsbuf);	 clear(tmpDen);
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(UE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = den(cnt) * sclvec[s];
		  cnt++;
		}
	 iC( dgemv(1.0, _UE2UCii, tmpDen, 1.0, chk) );
  } else {
	 iC( dgemv(1.0, _UE2UCii, den, 1.0, chk) );
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::DwnChk2DwnEqu_dgemv(int l, const DblNumVec& chk, DblNumVec& den)
{
  DblNumMat& _DC2DE = (_hom==true) ? _dwnChk2DwnEqu[0]: _dwnChk2DwnEqu[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0,l);
  //---------compute matrix
  if(_DC2DE.m()==0) {	 //cerr<<"DwnChk2DwnEqu compute"<<endl;
	 DblNumMat dd2c(plnDatSze(DC), plnDatSze(DE));
	 DblNumMat chkPos(dim(),samPos(DC).n());		clear(chkPos);	 iC( daxpy(R, samPos(DC), chkPos) ); //scale
	 DblNumMat denPos(dim(),samPos(DE).n());		clear(denPos);	 iC( daxpy(R, samPos(DE), denPos) ); //scale
	 
	 iC( _knl.kernel(denPos, denPos, chkPos, dd2c) );//matrix
	 _DC2DE.resize(plnDatSze(DE), plnDatSze(DC)); //_memused[2] += plndnenum()*dof()*plndncnum()*dof()*sizeof(double);
	 iC( pinv(dd2c, alt(), _DC2DE) );
  }
  //---------matvec
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(DE).n(), false, _wsbuf);	 clear(tmpDen);
	 iC( dgemv(1.0, _DC2DE, chk, 1.0, tmpDen) );
	 //scale
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, - l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(DE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  den(cnt) = den(cnt) + tmpDen(cnt) * sclvec[s];
		  cnt++;
		}
  } else {
	 iC( dgemv(1.0, _DC2DE, chk, 1.0, den) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::DwnEqu2DwnChk_dgemv(int l, Index3 idx, const DblNumVec& den, DblNumVec& chk)
{
  NumTns<DblNumMat>& _DE2DC = (_hom==true) ? _dwnEqu2DwnChk[0] : _dwnEqu2DwnChk[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0, l);  //OffTns<DblNumMat>& _DwnEqu2DwnChk = _de2DwnChk[l];  double R       = 1.0/pow(2.0, l);
  if(_DE2DC.m()==0)	 _DE2DC.resize(2,2,2);
  DblNumMat& _DE2DCii = _DE2DC(idx[0], idx[1], idx[2]);
  
  //---------compute matrix
  if(_DE2DCii.m()==0) {	 //cerr<<"DwnEqu2DwnChk compute"<<endl;
	 _DE2DCii.resize(plnDatSze(DC), plnDatSze(DE)); //_memused[3] += plndncnum()*dof()*plndnenum()*dof()*sizeof(double);
	 DblNumMat denPos(dim(),samPos(DE).n());		  clear(denPos);	 iC( daxpy(R, samPos(DE), denPos) ); //scale
	 DblNumMat chkPos(dim(),samPos(DC).n());		  clear(chkPos);	 iC( daxpy(0.5*R, samPos(DC), chkPos) ); //scale
	 for(int i=0; i<dim(); i++) for(int j=0; j<samPos(DC).n(); j++) chkPos(i,j) = chkPos(i,j) + (double(idx(i))-0.5)*R;
	 
	 iC( _knl.kernel(denPos, denPos, chkPos, _DE2DCii) );
  }
  //---------matvec
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(DE).n(), false, _wsbuf);	 clear(tmpDen);
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(DE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = den(cnt) * sclvec[s];
		  cnt++;
		}
	 iC( dgemv(1.0, _DE2DCii, tmpDen, 1.0, chk) );
  } else {
	 iC( dgemv(1.0, _DE2DCii, den, 1.0, chk) );
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::plnDen2EffDen(int l, const DblNumVec& plnDen, DblNumVec& effDen)
{
  DblNumVec regDen(regPos().n()*srcDOF()); clear(regDen);
  //iC( samDen2RegularDen(plnDen, regDen) );
  //rfftwnd_real_to_complex(_forplan, srcDOF(), regDen.data(), srcDOF(), 1, (fftw_complex*)(effDen.data()), srcDOF(), 1);
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(UE).n(), false, _wsbuf);	 clear(tmpDen);
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(UE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = plnDen(cnt) * sclvec[s];
		  cnt++;
		}
	 iC( samDen2RegDen(tmpDen, regDen) );
  } else {
	 iC( samDen2RegDen(plnDen, regDen) );
  }
  
  int nnn[3];  nnn[0] = 2*_np;  nnn[1] = 2*_np;  nnn[2] = 2*_np;
  fftw_plan forplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF(), regDen.data(),NULL, srcDOF(),1, (fftw_complex*)(effDen.data()),NULL, srcDOF(),1, FFTW_ESTIMATE);
  fftw_execute(forplan);
  fftw_destroy_plan(forplan);
  //rfftwnd_real_to_complex(_forplan, srcDOF(), regDen.data(), srcDOF(), 1, (fftw_complex*)(effDen.data()), srcDOF(), 1);
  
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::samDen2RegDen(const DblNumVec& samDen, DblNumVec& regDen)
{
  int np = _np;
  int rgnum = 2*np;
  int srcDOF = this->srcDOF();
  int cnt=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<srcDOF; f++) {
				regDen(srcDOF*rgoff + f) += samDen(srcDOF*cnt + f);
			 }
			 cnt++;
		  }
		}  //iC( PetscLogFlops(np*np*np*dof) );
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effVal2PlnVal(int level, const DblNumVec& effVal, DblNumVec& plnVal)
{
  DblNumVec regVal(regPos().n()*trgDOF());
  int nnn[3];  nnn[0] = 2*_np;  nnn[1] = 2*_np;  nnn[2] = 2*_np;
  fftw_plan invplan = fftw_plan_many_dft_c2r(3,nnn,trgDOF(), (fftw_complex*)(effVal.data()),NULL, trgDOF(),1, regVal.data(),NULL, trgDOF(),1, FFTW_ESTIMATE);
  fftw_execute(invplan);
  fftw_destroy_plan(invplan);  //rfftwnd_complex_to_real(_invplan, trgDOF(), (fftw_complex*)(effVal.data()), trgDOF(), 1, regVal.data(), trgDOF(), 1);
  
  iC( regVal2SamVal(regVal, plnVal) );
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::regVal2SamVal(const DblNumVec& regVal, DblNumVec& samVal)
{
  int np = _np;
  int rgnum = 2*np;
  int trgDOF = this->trgDOF();
  int cnt=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<trgDOF; f++) {
				samVal(trgDOF*cnt + f) += regVal(trgDOF*rgoff + f);
			 }
			 cnt++;
		  }
		}
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwEqu2DwnChk_dgemv(int l, Index3 idx, const DblNumVec& effDen, DblNumVec& effVal)
{
  OffTns<DblNumMat>& _UpwEqu2DwnChk = (_hom==true) ? _upwEqu2DwnChk[0] : _upwEqu2DwnChk[l];
  double R = (_hom==true) ? 1.0 : 1.0/pow(2.0, l); //OffTns< DblNumMat >& _UpwEqu2DwnChk = _upwEqu2DwnChk[l];  double R       = 1.0/pow(2.0, l);
  if(_UpwEqu2DwnChk.m()==0)	 _UpwEqu2DwnChk.resize(7,7,7,-3,-3,-3);
  //DblNumMat* _UpwEqu2DwnChkii = new DblNumMat;
  DblNumMat& _UpwEqu2DwnChkii = _UpwEqu2DwnChk(idx[0], idx[1], idx[2]);
  
  int effNum = (2*_np+2)*(2*_np)*(2*_np);
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  
  //---------compute matrix
  if(_UpwEqu2DwnChkii.m()==0) { //compute it if necessary
	 //-----------------------	 //cerr<<"UpwEqu2DwnChk(FFT) compute"<<endl;	 //COMPUTE FFT	 //Index3 idx = ii.idx();	 
	 iA( idx.linfty()>1 );
	 DblNumMat denPos(dim(),1);	 for(int i=0; i<dim(); i++)		denPos(i,0) = double(idx(i))*2.0*R; //shift
	 DblNumMat chkPos(dim(),regPos().n());	 clear(chkPos);	 iC( daxpy(R, regPos(), chkPos) );
	 DblNumMat tt(regPos().n()*trgDOF, srcDOF);
	 iC( _knl.kernel(denPos, denPos, chkPos, tt) );
	 // move data to tmp
	 DblNumMat tmp(trgDOF,regPos().n()*srcDOF);
	 for(int k=0; k<regPos().n();k++) {
		for(int i=0; i<trgDOF; i++)
		  for(int j=0; j<srcDOF; j++) {
			 tmp(i,j+k*srcDOF) = tt(i+k*trgDOF,j);
		  }
	 }
	 _UpwEqu2DwnChkii.resize(trgDOF*srcDOF, effNum); //_memused[4] += dof*dof*effNum(UE)*sizeof(double);
	 //forward FFT from tmp to _UpwEqu2DwnChkii;
	 
	 int nnn[3];  nnn[0] = 2*_np;  nnn[1] = 2*_np;  nnn[2] = 2*_np;
	 fftw_plan forplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF*trgDOF, tmp.data(),NULL, srcDOF*trgDOF,1, (fftw_complex*)(_UpwEqu2DwnChkii.data()),NULL, srcDOF*trgDOF,1, FFTW_ESTIMATE);
	 fftw_execute(forplan);
	 fftw_destroy_plan(forplan);

	 //rfftwnd_real_to_complex(_forplan, srcDOF*trgDOF, tmp.data(), srcDOF*trgDOF, 1, (fftw_complex*)(_UpwEqu2DwnChkii.data()), srcDOF*trgDOF, 1);
  }
  //---------matvec
  /*
  //scale
  DblNumVec tmpDen(effDen.m(), false, _wsbuf);  clear(tmpDen);
  fftw_complex* denptr = NULL;// = (fftw_complex*)(tmpDen.data());
  if(_hom==true) {
	 vector<double> sclvec(srcDOF);	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(2.0, l*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<effDen.m()/srcDOF; i++)
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = effDen(cnt) * sclvec[s];		cnt++;
		}
	 denptr = (fftw_complex*)(tmpDen.data());
  } else {
	 denptr = (fftw_complex*)(effDen.data());
  }
  */
  //fft mult
  double nrmfc = 1.0/double(regPos().n());
  fftw_complex* matPtr  = (fftw_complex*)(_UpwEqu2DwnChkii.data());
  fftw_complex* denPtr = (fftw_complex*)(effDen.data());
  fftw_complex* chkPtr   = (fftw_complex*)(effVal.data());
  int matStp  = srcDOF*trgDOF;
  int denStp = srcDOF;
  int chkStp   = trgDOF;
  
  double newalpha = nrmfc;
  for(int i=0; i<trgDOF; i++)
	 for(int j=0; j<srcDOF; j++) {
		int matOff = j*trgDOF + i;
		int denOff = j;
		int chkOff = i;
		iC( cptwvv(effNum/2, newalpha, matPtr+matOff, matStp, denPtr+denOff, denStp, chkPtr+chkOff, chkStp) );
	 }
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::localPos(int tp, Point3 center, double radius, DblNumMat& positions)
{
  const DblNumMat& bas = samPos(tp);
  positions.resize(dim(), bas.n());
  for(int i=0; i<dim(); i++)
	 for(int j=0; j<positions.n(); j++)
		positions(i,j) = center(i) + radius * bas(i,j);
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::samPosCal(int np, double R, DblNumMat& positions)
{
  int n = np*np*np - (np-2)*(np-2)*(np-2);
  positions.resize(dim(),n);
  double step = 2.0/(np-1);
  double init = -1.0;
  int cnt = 0;
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 double x = init + i*step;
			 double y = init + j*step;
			 double z = init + k*step;
			 positions(0,cnt) = R*x;
			 positions(1,cnt) = R*y;
			 positions(2,cnt) = R*z;
			 cnt++;
		  }
		}
  iA(cnt==n);
  return 0;
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::regPosCal(int np, double R, DblNumMat& positions)
{
  int n = 2*np*2*np*2*np;
  positions.resize(dim(), n);
  double step = 2.0/(np-1);
  int cnt = 0;
  for(int k=0; k<2*np; k++)
	 for(int j=0; j<2*np; j++)
		for(int i=0; i<2*np; i++) {
		  int gi = (i<np) ? i : i-2*np;
		  int gj = (j<np) ? j : j-2*np;
		  int gk = (k<np) ? k : k-2*np;
		  positions(0, cnt) = R * gi*step;
		  positions(1, cnt) = R * gj*step;
		  positions(2, cnt) = R * gk*step;
		  cnt ++;
		}
  iA(cnt==n);
  return 0;
}

// ---------------------------------------------------------------------- 
inline int MatMgnt3d::cptwvv(int n, double alpha, fftw_complex* x, int incrementX, fftw_complex* y, int incrementY, fftw_complex* z, int incrementZ)
{
  for(int i=0; i<n; i++) {
	 (*z)[0] += alpha * ( (*x)[0] * (*y)[0] - (*x)[1] * (*y)[1]);
	 (*z)[1] += alpha * ( (*x)[0] * (*y)[1] + (*x)[1] * (*y)[0]);
	 x = x + incrementX;
	 y = y + incrementY;
	 z = z + incrementZ;
  }  //iC( PetscLogFlops( 10*n ) );
  return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vector<MatMgnt3d> MatMgnt3d::_mmvec;

MatMgnt3d* MatMgnt3d::getmmptr(Kernel3d knl, int np)
{
  for(int i=0; i<_mmvec.size(); i++)
	 if(_mmvec[i].knl()==knl && _mmvec[i].np()==np)
		return &(_mmvec[i]);
  
  _mmvec.push_back( MatMgnt3d() );
  int last = _mmvec.size()-1;
  MatMgnt3d* tmp = &(_mmvec[last]); //get the last one
  tmp->knl() = knl;  tmp->np() = np;
  tmp->setup();
  return tmp;
}
