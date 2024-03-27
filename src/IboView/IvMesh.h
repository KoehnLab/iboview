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

#ifndef GL_MESH_H
#define GL_MESH_H

#include "Iv.h"

typedef unsigned int
   FIndexType; // <- type we'll use for our index buffers.

// a basic vertex: position, normal, and color.
struct FBaseVertex {
   vec3f vPos;
   vec3f vNorm;
   uint32_t dwColor;
   uint32_t iRole; // used to distinguish sub-surfaces. don't ask .oO (maybe turn this into a float to store scalar value?).
   static void AssignVertexAttributes();
};


// direct, non-indexed triangle.
struct FTriangle1 {
   FVec3d v[3];
   FTriangle1() {};
   FTriangle1(FVec3d v0, FVec3d v1, FVec3d v2) { v[0] = v0; v[1] = v1; v[2] = v2; };
};

// convert list of explicit triangles to indexed triangles.
void MakeIndexedTriangles(ct::TArray<FVec3d> &PosOut, ct::TArray<unsigned> &Indices, ct::TArray<FTriangle1> const &TrianglesIn);


template<class FVertex>
struct TIndexedTriangleList : public ct::FIntrusivePtrDest {
   ct::TArray<FVertex>
      Vertices;
   ct::TArray<vec3ui>
      Triangles;
   TIndexedTriangleList() {};

   // construct from list of vertex positions, vertex normals, and indices
   TIndexedTriangleList(FVec3d *Pos, FVec3d *Norm, size_t nVertices, unsigned *Indices, size_t nIndices, FVertex const &RefVertex);

   // append triangles of other list to this one.
   void Append(TIndexedTriangleList const &Other);
};
typedef TIndexedTriangleList<FBaseVertex>
   FIndexedTriangleList;

typedef ct::TIntrusivePtr<FIndexedTriangleList>
   FIndexedTriangleListPtr;


// make a sphere at vPos, with radius fRadius, by subdividing an icosahedron.
// Number of triangles is 20 * (4^nSubDivLevel).
// Basic vertex attributes (i.e., everything except for position and normal) are copied from RefVertex.
FIndexedTriangleListPtr MakeSubdivSphere(vec3f vPos, float fRadius, uint nSubDivLevel, FBaseVertex const &RefVertex);
// make a cylinder oriented along the Z axis, with radius fBaseRadius at z=0 and fTopRadius at z=Height.
FIndexedTriangleListPtr MakeCylinder(float fBaseRadius, float fTopRadius, float fHeight, uint nSlices, FBaseVertex const &RefVertex);
// similar to MakeSubdivSphere with nSubDivLevel = 0, but creates hard vertex normals.
FIndexedTriangleListPtr MakeIcosahedron(float fRadius, FBaseVertex const &RefVertex);






#endif // GL_MESH_H
