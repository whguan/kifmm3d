#ifndef _MEXAUX_HPP_
#define _MEXAUX_HPP_

#include "mex.h"
#include "matrix.h"

#include "fmm3d.hpp"

using std::cerr;
using std::endl;

//int
inline void mex2cpp(const mxArray*& md, int& cd);
inline void cpp2mex(const int& cd, mxArray*& md);
//bool
inline void mex2cpp(const mxArray*& md, bool& cd);
inline void cpp2mex(const bool& cd, mxArray*& md);
//double
inline void mex2cpp(const mxArray*& md, double& cd);
inline void cpp2mex(const double& cd, mxArray*& md);
//dblnummat
inline void mex2cpp(const mxArray*& md, DblNumMat& cd);
inline void cpp2mex(const DblNumMat& cd, mxArray*& md);
//dblnumvec
inline void mex2cpp(const mxArray*& md, DblNumVec& cd);
inline void cpp2mex(const DblNumVec& cd, mxArray*& md);
//point3
inline void mex2cpp(const mxArray*& md, Point3& cd);
inline void cpp2mex(const Point3& cd, mxArray*& md);
//Kernel3d
inline void mex2cpp(const mxArray*& md, Kernel3d& cd);
inline void cpp2mex(const Kernel3d& cd, mxArray*& md);
//vector<double>
inline void mex2cpp(const mxArray*& md, vector<double>& cd);
inline void cpp2mex(const vector<double>& cd, mxArray*& md);
//vector<T>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd);
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md);
//NumTns<T>
template <class T> inline void mex2cpp(const mxArray*& md, NumTns<T>& cd);
template <class T> inline void cpp2mex(const NumTns<T>& cd, mxArray*& md);
//OffTns<T>
template <class T> inline void mex2cpp(const mxArray*& md, OffTns<T>& cd);
template <class T> inline void cpp2mex(const OffTns<T>& cd, mxArray*& md);
//map<S,T>
template <class S, class T> inline void mex2cpp(const mxArray*& md, map<S,T>& cd);
template <class S, class T> inline void cpp2mex(const map<S,T>& cd, mxArray*& md);
//matmgnt
inline void mex2cpp(const mxArray*& md, MatMgnt3d& cd);
inline void cpp2mex(const MatMgnt3d& cd, mxArray*& md);


//----------------------int
inline void mex2cpp(const mxArray*& md, int& cd)
{
  cd = int(mxGetScalar(md));
  return;
}
inline void cpp2mex(const int& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}
//----------------------int
inline void mex2cpp(const mxArray*& md, bool& cd)
{
  cd = bool(mxGetScalar(md));
  return;
}
inline void cpp2mex(const bool& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}
//----------------------double
inline void mex2cpp(const mxArray*& md, double& cd)
{
  cd = mxGetScalar(md);
  return;
}
inline void cpp2mex(const double& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}
//----------------------DblNumMat
inline void mex2cpp(const mxArray*& md, DblNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  cd.resize(m,n);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		cd(i,j) = xr[cnt];
		cnt++;
	 }
  return;
}
inline void cpp2mex(const DblNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		xr[cnt] = cd(i,j);
		cnt++;
	 }
  return;
}
//----------------------dblnumvec
inline void mex2cpp(const mxArray*& md, DblNumVec& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md); iA(n==1);
  double* xr = mxGetPr(md);
  cd.resize(m);
  for(int i=0; i<m; i++)
	 cd(i) = xr[i];
  return;
}
inline void cpp2mex(const DblNumVec& cd, mxArray*& md)
{
  int m = cd.m();
  md = mxCreateDoubleMatrix(m, 1, mxREAL);
  double* xr = mxGetPr(md);
  for(int i=0; i<m; i++)
	 xr[i] = cd(i);
  return;
}
//----------------------Point3
inline void mex2cpp(const mxArray*& md, Point3& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);  iA(m==3 && n==1);
  double* xr = mxGetPr(md);
  cd = Point3(xr[0], xr[1], xr[2]);
  return;
}
inline void cpp2mex(const Point3& cd, mxArray*& md)
{
  md = mxCreateDoubleMatrix(3, 1, mxREAL);
  double* xr = mxGetPr(md);
  xr[0] = cd(0);  xr[1] = cd(1);  xr[2] = cd(2);
  return;
}
//----------------------Kernel3d
inline void mex2cpp(const mxArray*& md, Kernel3d& cd)
{
  int m = mxGetM(md); iA(m==2);
  int n = mxGetN(md); iA(n==1);  //cerr<<m<<" "<<n<<endl;
  { const mxArray* tt = mxGetCell(md, 0);  mex2cpp(tt, cd.kt()); }
  { const mxArray* tt = mxGetCell(md, 1);  mex2cpp(tt, cd.coefs()); }
  return;
}
inline void cpp2mex(const Kernel3d& cd, mxArray*& md)
{
  int m = 2;
  int n = 1;
  md = mxCreateCellMatrix(m,n);
  mxArray* ss;
  cpp2mex(cd.kt(), ss);  mxSetCell(md, 0, ss);
  cpp2mex(cd.coefs(), ss);  mxSetCell(md, 1, ss);
  return;
}
//----------------------vector<double>
inline void mex2cpp(const mxArray*& md, vector<double>& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md); iA(n==1);  //cerr<<m<<"   "<<n<<endl;
  double* xr = mxGetPr(md);
  cd.resize(m*n);
  for(int i=0; i<m*n; i++)
	 cd[i] = xr[i];
  return;
}
inline void cpp2mex(const vector<double>& cd, mxArray*& md)
{
  int m = cd.size();
  int n = 1;
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  for(int i=0; i<m*n; i++)
	 xr[i] = cd[i];
  return;
}
//----------------------vector<...>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);  assert(n==1);
  cd.resize(m*n);
  for(int ci=0; ci<m*n; ci++) {
	 const mxArray*tt = mxGetCell(md, ci);
	 mex2cpp(tt, cd[ci]);
  }
  return;
}
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(n, 1);
  for(int ci=0; ci<n; ci++) {
	 mxArray* ss;	 cpp2mex(cd[ci], ss);
	 mxSetCell(md, ci, ss);
  }
  return;
}
//----------------------numtns<...>
template <class T> inline void mex2cpp(const mxArray*& md, NumTns<T>& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];  int n = dims[1];  int p = dims[2];
  cd.resize(m,n,p);
  int cnt = 0;
  for(int k=0; k<p; k++)
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  const mxArray* tt = mxGetCell(md, cnt);		  mex2cpp(tt, cd(i,j,k));
		  cnt++;
		}
  return;
}
template <class T> inline void cpp2mex(const NumTns<T>& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int p = cd.p();
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateCellArray(3,dims);
  int cnt = 0;
  for(int k=0; k<p; k++)
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  mxArray* ss;		cpp2mex(cd(i,j,k), ss);		mxSetCell(md, cnt, ss);
		  cnt++;
		}
  return;
}
//----------------------offtns<...>
template <class T> inline void mex2cpp(const mxArray*& md, OffTns<T>& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];  int n = dims[1];  int p = dims[2];
  int s = -m/2;  int t = -n/2;  int u = -p/2;
  cd.resize(m,n,p,s,t,u);
  int cnt = 0;
  for(int k=u; k<u+p; k++)
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  const mxArray* tt = mxGetCell(md, cnt);		  mex2cpp(tt, cd(i,j,k));
		  cnt++;
		}
  return;
}
template <class T> inline void cpp2mex(const OffTns<T>& cd, mxArray*& md)
{
  int m = cd.m();  int n = cd.n();  int p = cd.p();
  int s = -m/2;  int t = -n/2;  int u = -p/2;
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateCellArray(3,dims);
  int cnt = 0;
  for(int k=u; k<u+p; k++)
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  mxArray* ss;		cpp2mex(cd(i,j,k), ss);		mxSetCell(md, cnt, ss);
		  cnt++;
		}
  return;
}
//----------------------map<S,T>
template <class S, class T> inline void mex2cpp(const mxArray*& md, map<S,T>& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);  assert(n==1);
  for(int ci=0; ci<m*n; ci++) {
	 const mxArray* now = mxGetCell(md, ci);
	 const mxArray* tmp;
	 S ss;	 tmp = mxGetCell(now,0); mex2cpp(tmp, ss);
	 T tt;	 tmp = mxGetCell(now,1); mex2cpp(tmp, tt);
	 cd[ss] = tt;
  }
  return;
}
template <class S, class T> inline void cpp2mex(const map<S,T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(n,1);
  int cnt = 0;
  for(map<S,T>::const_iterator mi=cd.begin(); mi!=cd.end(); mi++) {
	 mxArray* now = mxCreateCellMatrix(2,1);
	 mxArray* ss;	 cpp2mex((*mi).first, ss);	 mxSetCell(now, 0, ss);
	 mxArray* tt;	 cpp2mex((*mi).second, tt);	 mxSetCell(now, 1, tt);
	 mxSetCell(md, cnt, now);
	 cnt++;
  }
  assert(cnt==n);
  return;
}
//----------------------matmgnt3d
inline void mex2cpp(const mxArray*& md, MatMgnt3d& cd)
{
  const mxArray* tt;
  tt = mxGetCell(md,0);  mex2cpp(tt, cd._knl);
  tt = mxGetCell(md,1);  mex2cpp(tt, cd._np);
  tt = mxGetCell(md,2);  mex2cpp(tt, cd._hom);
  tt = mxGetCell(md,3);  mex2cpp(tt, cd._degvec);
  tt = mxGetCell(md,4);  mex2cpp(tt, cd._uc2ue);
  tt = mxGetCell(md,5);  mex2cpp(tt, cd._ue2uc);
  tt = mxGetCell(md,6);  mex2cpp(tt, cd._dc2de);
  tt = mxGetCell(md,7);  mex2cpp(tt, cd._de2dc);
  tt = mxGetCell(md,8);  mex2cpp(tt, cd._ue2dc);
  tt = mxGetCell(md,9);  mex2cpp(tt, cd._splpos[0]);
  tt = mxGetCell(md,10); mex2cpp(tt, cd._splpos[1]);
  tt = mxGetCell(md,11); mex2cpp(tt, cd._splpos[2]);
  tt = mxGetCell(md,12); mex2cpp(tt, cd._splpos[3]);
  tt = mxGetCell(md,13); mex2cpp(tt, cd._regpos);
  return;
}
inline void cpp2mex(const MatMgnt3d& cd, mxArray*& md)
{
  md = mxCreateCellMatrix(1,14);
  
  mxArray* ss;
  cpp2mex(cd._knl, ss);  mxSetCell(md, 0, ss);
  cpp2mex(cd._np, ss);  mxSetCell(md, 1, ss);
  cpp2mex(cd._hom, ss);  mxSetCell(md, 2, ss);
  cpp2mex(cd._degvec, ss);  mxSetCell(md, 3, ss);
  cpp2mex(cd._uc2ue, ss);  mxSetCell(md, 4, ss);
  cpp2mex(cd._ue2uc, ss);  mxSetCell(md, 5, ss);
  cpp2mex(cd._dc2de, ss);  mxSetCell(md, 6, ss);
  cpp2mex(cd._de2dc, ss);  mxSetCell(md, 7, ss);
  cpp2mex(cd._ue2dc, ss);  mxSetCell(md, 8, ss);
  cpp2mex(cd._splpos[0], ss);  mxSetCell(md, 9, ss);
  cpp2mex(cd._splpos[1], ss);  mxSetCell(md, 10, ss);
  cpp2mex(cd._splpos[2], ss);  mxSetCell(md, 11, ss);
  cpp2mex(cd._splpos[2], ss);  mxSetCell(md, 12, ss);
  cpp2mex(cd._regpos, ss);  mxSetCell(md, 13, ss);
  
  return;
}


#endif
