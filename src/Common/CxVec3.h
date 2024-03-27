/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#ifndef CX_VEC3_H
#define CX_VEC3_H

#include <cstddef>
#include <cmath>

#define TVec3 TVector3<T>

namespace ct {

// don't ask.
inline float cx_sqrt(float z) { using std::sqrt; return sqrt(z); }
inline double cx_sqrt(double z) { using std::sqrt; return sqrt(z); }

template<class T>
inline T sqr(T x) { return x*x; }

template<class scalar_t>
inline scalar_t _ThreshAlmostZero(); // { return 0; };

// default thresholds for AlmostEqualQ and assuming zero lengths in Normalize()/Normalized()
template<> inline float _ThreshAlmostZero<float>() { return 3.81e-6f; }
template<> inline double _ThreshAlmostZero<double>() { return 1.08e-12; }
#ifdef FLOAT_128_NATIVE_EXT
   template<> inline __float128 _ThreshAlmostZero<__float128>() { return 3.07e-26q; }
#endif
// ^-- Note: default thresholds computed based on premise that numerical calculation may cause loss of 25% of manissa bits of precision.
// >>> 2**-(3*24/4) # for binary32 float ("single precision")
// 3.814697265625e-06
// >>> 2**-(3*53/4) # for binary64 float ("double precision")
// 1.0815775704056441e-12
// >>> 2**-(3*113/4) # for binary128 float ("quadruple precision")
// 3.074028343251155e-26
template <class scalar_t> inline scalar_t _ThreshAlmostZero(scalar_t const &) {
   return _ThreshAlmostZero<scalar_t>();
}


// that is only here to avoid external dependencies...
template<class T>
struct TVector3
{
   typedef T value_type;
   T m[3];

   enum {
      static_size = 3
   };

   TVector3() {};
   TVector3(T x_, T y_, T z_) {m[0] = x_; m[1] = y_; m[2] = z_;}
   explicit TVector3(T const *xyz) {m[0] = xyz[0]; m[1] = xyz[1]; m[2] = xyz[2];}
//   explicit TVector3(T xyz) {m[0] = xyz; m[1] = xyz; m[2] = xyz;}
   explicit TVector3(int xyz) {m[0] = T(xyz); m[1] = T(xyz); m[2] = T(xyz);}
   // ^- this one is mainly meant for TVec3<...> v(0);
   template<class S> explicit TVector3(TVector3<S> const &other) {m[0] = static_cast<T>(other[0]); m[1] = static_cast<T>(other[1]); m[2] = static_cast<T>(other[2]);}


   void operator += (TVec3 const &other) {this->m[0] += other.m[0]; this->m[1] += other.m[1]; this->m[2] += other.m[2];}
   void operator -= (TVec3 const &other) {this->m[0] -= other.m[0]; this->m[1] -= other.m[1]; this->m[2] -= other.m[2];}
   void operator *= (T f) {this->m[0] *= f; this->m[1] *= f; this->m[2] *= f;}
   void operator /= (T f) {*this *= T(1)/f;}

//   T &operator[] (unsigned i) { return this->m[i]; }
//   T const &operator[] (unsigned i) const { return this->m[i]; }
//   T &operator[] (int i) { return this->m[i]; }
//   T const &operator[] (int i) const { return this->m[i]; }
//   T &operator[] (unsigned long i) { return this->m[i]; }
//   T const &operator[] (unsigned long i) const { return this->m[i]; }
//   T &operator[] (long i) { return this->m[i]; }
//   T const &operator[] (long i) const { return this->m[i]; }
// ^- in VC long is *also* 32bit, which means that there is no way of
// defining these functions in such a way that int, long, and size_t will work at the same time.
// (because size_t may also be typedef'd to one of those)
   T &operator[] (size_t i) { return this->m[i]; }
   T const &operator[] (size_t i) const { return this->m[i]; }

//    operator T* () {return &this->m[0];}
//    operator T const* () const {return &this->m[0];}
// ^- these ones can also lead to ambiguities... in operator [](int): T* vs actual [].

   T LengthSq() const { return m[0]*m[0] + m[1]*m[1] + m[2]*m[2]; }
   T Length() const { return cx_sqrt(this->LengthSq()); }
   // normalize *this
   void Normalize(T thresh = _ThreshAlmostZero<T>()) {
      T lensq = this->LengthSq();
      if (thresh != 0 && lensq < sqr(thresh))
         *this = Zero();
      else
         (*this) *= value_type(1)/cx_sqrt(lensq);
   }
   // return normalized copy of *this
   TVec3 Normalized(T thresh = _ThreshAlmostZero<T>()) const { TVec3 cp(*this); cp.Normalize(thresh); return cp; }

   T const &x () const { return m[0]; }
   T const &y () const { return m[1]; }
   T const &z () const { return m[2]; }
   T &x () { return m[0]; }
   T &y () { return m[1]; }
   T &z () { return m[2]; }

   size_t size() const { return 3; }
   T *begin() { return &m[0]; }
   T const *begin() const { return &m[0]; }
   T *end() { return &m[size()]; }
   T const *end() const { return &m[size()]; }

   static TVector3 Zero() { return TVector3(T(0),T(0),T(0)); }
   static TVector3 One() { return TVector3(T(1),T(1),T(1)); }
};

template<class T> inline bool operator == (TVec3 const &a, TVec3 const &b) { return a[0]==b[0] && a[1]==b[1] && a[2]==b[2]; }
template<class T> inline bool operator != (TVec3 const &a, TVec3 const &b) { return a[0]!=b[0] || a[1]!=b[1] || a[2]!=b[2]; }
template<class T> inline TVec3 operator + (TVec3 const &a, TVec3 const &b) { return TVec3(a[0]+b[0], a[1]+b[1], a[2]+b[2]); }
template<class T> inline TVec3 operator - (TVec3 const &a, TVec3 const &b) { return TVec3(a[0]-b[0], a[1]-b[1], a[2]-b[2]); }
template<class T> inline TVec3 operator - (TVec3 const &a) { return TVec3(-a[0], -a[1], -a[2]); }
template<class T> inline TVec3 operator * (T f, TVec3 const &b) { return TVec3(f*b[0], f*b[1], f*b[2]); }
template<class T> inline TVec3 operator * (TVec3 const &b, T f) { return TVec3(f*b[0], f*b[1], f*b[2]); }
template<class T> inline T Dot(TVec3 const &a, TVec3 const &b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<class T> inline T LengthSq(TVec3 const &a) { return Dot(a,a); }
template<class T> inline T Length(TVec3 const &a) { return cx_sqrt(Dot(a,a)); }
template<class T> inline T DistSq(TVec3 const &a, TVec3 const &b) { return LengthSq(a-b); }
template<class T> inline T Dist(TVec3 const &a, TVec3 const &b) { return cx_sqrt(DistSq(a,b)); }
template<class T> inline TVec3 Normalized(TVec3 const &a, T thresh = _ThreshAlmostZero<T>()) { return a.Normalized(thresh); }


template<class T> inline T DistSq(T const a[3], T const b[3]) { return sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]); }
template<class T> inline T DistSq3(T const *a, T const *b) { return sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]); }
template<class T> inline T DistSq2(T const *a, T const *b) { return sqr(a[0]-b[0]) + sqr(a[1]-b[1]); }

template<class T>
static void Cross(TVector3<T> &Out, TVector3<T> const &a, TVector3<T> const &b)
{
   Out[0] = a[1]*b[2] - b[1] * a[2];
   Out[1] = a[2]*b[0] - b[2] * a[0];
   Out[2] = a[0]*b[1] - b[0] * a[1];
}

template<class T>
static TVector3<T> Cross(TVector3<T> const &a, TVector3<T> const &b)
{
   TVector3<T> Out;
   Cross(Out, a, b);
   return Out;
}


template<class scalar_t>
bool AlmostEqualQ(ct::TVector3<scalar_t> const &va, ct::TVector3<scalar_t> const &vb, scalar_t thresh = _ThreshAlmostZero<scalar_t>()) {
   return DistSq(va, vb) < ct::sqr(thresh);
}


} // namespace ct

#undef TVec3

#endif // CX_VEC3_H
