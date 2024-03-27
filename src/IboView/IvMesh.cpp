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

#include "Iv.h"
#include <algorithm>
#include <GL/glew.h>
#include "IvMesh.h"
#include "CxPodArray.h"
#include "CxVec3.h"


template<class FVertex>
TIndexedTriangleList<FVertex>::TIndexedTriangleList(FVec3d *Pos, FVec3d *Norm, size_t nVertices, unsigned *Indices, size_t nIndices, FVertex const &RefVertex)
{
   // translate to target location and convert to GL format
   Vertices.reserve(nVertices);
   for (unsigned iVert = 0; iVert < nVertices; ++ iVert) {
      FVertex
         vOut(RefVertex);
      FVec3d const
         &vIn = Pos[iVert],
         &vNorm = Norm[iVert];
      vOut.vNorm = vec3f(vNorm[0], vNorm[1], vNorm[2]);
      vOut.vPos  = vec3f(vIn[0], vIn[1], vIn[2]);
      Vertices.push_back(vOut);
   }

   assert(nIndices % 3 == 0);
   Triangles.reserve(nIndices);
   for (unsigned iTri = 0; iTri < nIndices/3; ++ iTri) {
      Triangles.push_back(vec3ui(Indices[3*iTri+0], Indices[3*iTri+1], Indices[3*iTri+2]));
   }
}


template<class FVertex>
void TIndexedTriangleList<FVertex>::Append(TIndexedTriangleList const &Other)
{
   unsigned
      iVertexOffset = Vertices.size(),
      iTriangle0 = Triangles.size();
   Vertices.insert(Vertices.end(), Other.Vertices.begin(), Other.Vertices.end());
   Triangles.insert(Triangles.end(), Other.Triangles.begin(), Other.Triangles.end());
   // fix up triangle indices of new triangles.
   for (unsigned i = iTriangle0; i < Triangles.size(); ++ i)
      for (unsigned iVx = 0; iVx < 3; ++ iVx)
         Triangles[i][iVx] += iVertexOffset;
}


namespace SubdivSphere {
   using ct::TArray;
   typedef TArray<FTriangle1>
      FTriangleArray;

   // see Icosahedron.py
   static const double fIcosahedronCoordinates[12][3] = {
      {0.85065080835203999, 0.0, 0.52573111211913359},
      {0.85065080835203999, 0.0, -0.52573111211913359},
      {-0.85065080835203999, 0.0, 0.52573111211913359},
      {-0.85065080835203999, 0.0, -0.52573111211913359},
      {0.0, 0.52573111211913359, 0.85065080835203999},
      {0.0, -0.52573111211913359, 0.85065080835203999},
      {0.0, 0.52573111211913359, -0.85065080835203999},
      {0.0, -0.52573111211913359, -0.85065080835203999},
      {0.52573111211913359, 0.85065080835203999, 0.0},
      {-0.52573111211913359, 0.85065080835203999, 0.0},
      {0.52573111211913359, -0.85065080835203999, 0.0},
      {-0.52573111211913359, -0.85065080835203999, 0.0}
   };
   static const unsigned iIcosahedronTriangles[20][3] = {{0,1,8},{0,4,5},{0,5,10},{0,8,4},{0,10,1},
      {1,6,8},{1,7,6},{1,10,7},{2,3,11},{2,4,9},{2,5,4},{2,9,3},{2,11,5},{3,6,7},{3,7,11},{3,9,6},
      {4,8,9},{5,11,10},{6,9,8},{7,10,11}};


   // make a non-indexed triangle list out of a given list of positions and indices.
   void MakeInitialMesh(FTriangleArray &Out, double const (*vCoords)[3], unsigned const (*iTriangles)[3], unsigned nTriangles) {
      Out.reserve(nTriangles);
      for (unsigned iTri = 0; iTri < nTriangles; ++ iTri) {
         unsigned const
            *ii = &iTriangles[iTri][0];
         FTriangle1 t;
         t.v[0] = FVec3d(vCoords[ii[0]]);
         t.v[1] = FVec3d(vCoords[ii[1]]);
         t.v[2] = FVec3d(vCoords[ii[2]]);
         Out.push_back(t);
      }
   }

   // subdivide the given triangles by bisecting their edges and projecting
   // the new points onto the unit sphere.
   //
   // (note: there are some valid subdivisions which cannot be accessed with this:
   //  technically we could divide the edges into any number of equal length segments,
   //  and project those back.)
   void BinarilySubdivideUnitSphere(FTriangleArray &Out, FTriangleArray const &In)
   {
      Out.clear();
      Out.reserve(In.size() * 4);
      for (unsigned iTri = 0; iTri < In.size(); ++ iTri) {
         //            0
         //
         //
         //      5          3
         //
         //
         //  2        4        1
         FVec3d v[6];
         for (unsigned i = 0; i < 3; ++i)
            v[i] = In[iTri].v[i];
         for (unsigned i = 3; i < 6; ++i) {
            v[i] = .5 * (v[i-3] + v[(i-2)%3]); // e.g., new vertex 3 goes between old vertices 0 and 1.
            // project new coordinate to unit sphere.
            v[i] /= ct::Length(v[i]);
         }
         Out.push_back(FTriangle1(v[5], v[0], v[3]));
         Out.push_back(FTriangle1(v[1], v[4], v[3]));
         Out.push_back(FTriangle1(v[3], v[4], v[5]));
         Out.push_back(FTriangle1(v[5], v[4], v[2]));
      }
   }

   // make non-indexed triangle list representing a unit sphere by recursively
   // sub-dividing an icosahedron.
   void MakeSubdivUnitSphere(FTriangleArray &Out, unsigned nSubdivLevel)
   {
      // take an icosahedron and divide it along the edges, each time normalizing
      // the triangles back onto the sphere.
      FTriangleArray
         T0, T1;
      MakeInitialMesh(T0, fIcosahedronCoordinates, iIcosahedronTriangles, 20);
      for (unsigned iSubdivLevel = 0; iSubdivLevel < nSubdivLevel; ++ iSubdivLevel) {
         // subdivide T0 into T1.
         BinarilySubdivideUnitSphere(T1, T0);
         // exchange them.
         T1.swap(T0);
      }
      Out.swap(T0);
   }

   // define some auxiliary structures for joining equal vertices on different triangles.
   // (to turn a direct triangle list into an indexed triangle list)
   static const FVec3d RandomPlane(1.2344591,0.0823435367,-.323451231);
   static const double fEpsilon = 1e-8;

   struct FVertexInfo {
      FVec3d vPos;
      unsigned iOrig;
      double fSortPred;
      FVertexInfo(FVec3d vPos_, unsigned iOrig_)
         : vPos(vPos_), iOrig(iOrig_), fSortPred(ct::Dot(RandomPlane, vPos_))
      {}

      // induce a weak ordering to simplify searching for equal vertices.
      bool operator < (FVertexInfo const &other) const {
         // sort by distance to a random plane.
         return this->fSortPred < other.fSortPred;
      }

      bool operator == (FVertexInfo const &other) const {
         // are both are sufficiently close to be considered equal?
         return DistSq(this->vPos, other.vPos) < fEpsilon*fEpsilon;
         // note: this test is not compatible with operator <.
      }
   };


   // find unique positions in PosIn and return them in PosOut.
   // Indices is an array of length PosIn.size() such that PosOut[Indices[i]] == PosIn[i].
   void FindUniquePositions(TArray<FVec3d> &PosOut, TArray<unsigned> &Indices, TArray<FVec3d> const &PosIn)
   {
      TArray<FVertexInfo>
         vs;
      vs.reserve(PosIn.size());
      for (unsigned i = 0; i < PosIn.size(); ++ i)
         vs.push_back(FVertexInfo(PosIn[i], i));
      std::sort(vs.begin(), vs.end());

      PosOut.clear();
      PosOut.reserve(PosIn.size());
      Indices.resize(PosIn.size());
      unsigned
         iNext;
      for (unsigned i = 0; i < vs.size(); i = iNext) {
         iNext = i + 1;
         while (iNext != vs.size() && vs[i] == vs[iNext])
            iNext += 1;
         // ^- note: that is not quite right. Vertices which are too far apart
         // according to our sorting criterion cannot be equal, but vertices
         // which are equal need non necessarily end up rigt next to each other
         // in the metric. we should search through *all* vertices which lie close
         // to the current one according to the sorting metric.
         for (unsigned k = i; k != iNext; ++ k)
            Indices[vs[k].iOrig] = PosOut.size();
         PosOut.push_back(vs[i].vPos);
      }
   }


   void TranslateAndScale(ct::TArray<FVec3d> &Pos, double fScale, FVec3d vDisplacement) {
      for (unsigned i = 0; i < Pos.size(); ++ i)
         Pos[i] = fScale * Pos[i] + vDisplacement;
   }

} // namespace SubdivSphere


void MakeIndexedTriangles(ct::TArray<FVec3d> &PosOut, ct::TArray<unsigned> &Indices, ct::TArray<FTriangle1> const &TrianglesIn)
{
   using namespace SubdivSphere;

   // flatten the triangle array to turn it into a list of position vectors.
   TArray<FVec3d>
      Triangles;
   Triangles.reserve(3*TrianglesIn.size());
   for (unsigned iTri = 0; iTri < TrianglesIn.size(); ++ iTri)
      for (unsigned iv = 0; iv != 3; ++ iv)
         Triangles.push_back(TrianglesIn[iTri].v[iv]);
   // make indexed triangles.
   TArray<unsigned>
      VertexIndices;
   FindUniquePositions(PosOut, VertexIndices, Triangles);
   Indices.clear();
   Indices.reserve(3*TrianglesIn.size());
   for (unsigned iTri = 0; iTri < TrianglesIn.size(); ++ iTri)
      for (unsigned iv = 0; iv != 3; ++ iv)
         Indices.push_back(VertexIndices[3*iTri + iv]);
}


FIndexedTriangleListPtr MakeSubdivSphere(vec3f vPos, float fRadius, uint nSubdivLevel, FBaseVertex const &RefVertex)
{
   using namespace SubdivSphere;

   // make vertex positions and triangle indices of unit sphere
   ct::TArray<FVec3d> Pos;
   ct::TArray<unsigned> Indices;
   {
      FTriangleArray
         T0;
      MakeSubdivUnitSphere(T0, nSubdivLevel);
      MakeIndexedTriangles(Pos, Indices, T0);
   }

   // translate to target location and convert to GL format
   TranslateAndScale(Pos, (double)fRadius, FVec3d(vPos[0], vPos[1], vPos[2]));
   return FIndexedTriangleListPtr(new FIndexedTriangleList(&Pos[0], &Pos[0], Pos.size(), &Indices[0], Indices.size(), RefVertex));
}

FIndexedTriangleListPtr MakeCylinder(float fBaseRadius, float fTopRadius, float fHeight, uint nSlices, FBaseVertex const &RefVertex)
{
   using namespace SubdivSphere;

   // make vertex positions and triangle indices of cylinder
   ct::TArray<FVec3d> Pos, Norm;
   ct::TArray<unsigned> Indices;

   unsigned n = nSlices;
   // base and top circle
   Pos.resize(2*n);
   Norm.resize(2*n);
   // a circle with nSlices has nSlices edges, and each of those leads to
   // one face composed of two triangles.
   Indices.resize(3*2*n);

   //       \         t
   //     ny \     /<-->|
   //         \   /     |
   //       <->\ /      |
   //        nx /       |
   //          /        |
   //         /<------->|
   //              b

   double
      nx = fHeight,
      ny = -(fTopRadius - fBaseRadius),
      nf = 1./std::sqrt(nx*nx + ny*ny);
   nx *= nf;
   ny *= nf;

   for (unsigned iAngle = 0; iAngle < n; ++ iAngle) {
      double
         fAngle = 2 * M_PI * iAngle/float(n),
         cx = std::cos(fAngle),
         cy = std::sin(fAngle);
      Pos[iAngle] = FVec3d(fBaseRadius*cx, fBaseRadius*cy, 0.);
      Pos[iAngle+n] = FVec3d(fTopRadius*cx, fTopRadius*cy, fHeight);
      Norm[iAngle] = FVec3d(nx*cx, nx*cy, ny);
      Norm[iAngle+n] = FVec3d(nx*cx, nx*cy, ny);
//       Norm[iAngle+n] = FVec3d(cx, cy, 0.); // <- to give a bent look?

      {
         //    iAng-----iAng-1
         //    |    /
         //    |   /
         //    iAng+n   iAng+n-1
         Indices[6*iAngle+0] = iAngle;
         Indices[6*iAngle+1] = (iAngle + 1)%n;
         Indices[6*iAngle+2] = iAngle + n;
         Indices[6*iAngle+3] = iAngle + n;
         Indices[6*iAngle+4] = (iAngle + 1)%n;
         Indices[6*iAngle+5] = (iAngle + 1)%n + n;
      }
   }
   return FIndexedTriangleListPtr(new FIndexedTriangleList(&Pos[0], &Norm[0], Pos.size(), &Indices[0], Indices.size(), RefVertex));
}

FIndexedTriangleListPtr MakeIcosahedron(float fRadius, FBaseVertex const &RefVertex)
{
   using namespace SubdivSphere;
   // we don't bother with indexing it up here. we just create three vertices
   // and three indices per face.
   ct::TArray<FVec3d> Pos, Norm;
   ct::TArray<unsigned> Indices;

   unsigned nFaces = 20;
   Pos.resize(3*nFaces);
   Norm.resize(3*nFaces);
   Indices.reserve(3*nFaces);

   for (size_t iFace = 0; iFace < 20; ++ iFace) {
      // we don't bother with indexing it up. we just create three vertices and three indices
      // per face.
      FVec3d
         &v0 = Pos[3*iFace],
         &v1 = Pos[3*iFace+1],
         &v2 = Pos[3*iFace+2];

      unsigned const
         *pFaceIndices = iIcosahedronTriangles[iFace];
      for (size_t ixyz = 0; ixyz < 3; ++ ixyz) {
         v0[ixyz] = fRadius * fIcosahedronCoordinates[pFaceIndices[0]][ixyz];
         v1[ixyz] = fRadius * fIcosahedronCoordinates[pFaceIndices[1]][ixyz];
         v2[ixyz] = fRadius * fIcosahedronCoordinates[pFaceIndices[2]][ixyz];
      }
      FVec3d
         vFaceNormal = Normalized(Cross(v1-v0,v2-v0));
      for (size_t iVert = 0; iVert < 3; ++ iVert) {
         Norm[3*iFace + iVert] = vFaceNormal;
         Indices.push_back(3*iFace + iVert);
      }
   }
   assert(Pos.size() == 60 && Norm.size() == 60 && Indices.size() == 60);
   return FIndexedTriangleListPtr(new FIndexedTriangleList(&Pos[0], &Norm[0], Pos.size(), &Indices[0], Indices.size(), RefVertex));
}



template struct TIndexedTriangleList<FBaseVertex>;
