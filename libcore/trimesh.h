//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the software accompanying the following publication
//
// @INPROCEEDINGS{Balzeretal2014,
// author = {J. Balzer and D. Acevedo-Feliz and S. Soatto and S. H\"ofer and M. Hadwiger and J. Beyerer},
// title = {Cavlectometry: Towards Holistic Reconstruction of Large Mirror Objects},
// booktitle = {International Conference on 3D Vision (3DV)},
// year = {2014}}
//
// This file contains free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This source code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this file. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef TRIMESH_H_
#define TRIMESH_H_

#include <set>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Adaptive/Composite/CompositeT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::Subdivider::Adaptive::CompositeTraits> TriangleMesh;

#include "linalg.h"
#include "cam.h"
#include "bbox.h"

/*! \brief triangle mesh
 *
 *
 *
 */
class CTriangleMesh:public TriangleMesh {

public:

	//! Standard constructor.
	CTriangleMesh();

    //! Area.
    float FaceArea(FaceHandle fh);

     //! Computes the integration weight for a vertex.
    float VertexQuadratureWeight(VertexHandle vh);

    //! Computes the Voronoi area around a vertex.
    float VoronoiArea(VertexHandle vh);

	//! Returns barycenter of a mesh face.
    vec3f Barycenter(FaceHandle fh);

    //! Returns the barycenter of the entire mesh.
    vec3f Barycenter() const;

	//! Exterior normal at a boundary vertex.
    CTriangleMesh::Normal ExteriorNormal(VertexHandle vh);

	//! Refines the mesh uniformly with Loop's scheme.
    void UniformMeshRefinement(size_t n);

    //! Collapses edges with aspect ratio to high to do FE-analysis.
    void EdgeCollapse(double threshold);

    /*! Smoothes the boundary of mesh.
     *
     * \details Only works for surfaces with topology of a disk.
     *
     */
    bool SmoothBoundary(size_t n);

    /*! Finds a vertex handle that is on the boundary of the mesh.
     *
     * \details If the boundary is multiply-connected, the first hit will be returned.
     */
    OpenMesh::VertexHandle FindBoundaryVertex();

    //! Deletes a given set of vertices.
    void DeleteVertices(std::set<VertexHandle> vertices);

    //! Deletes a given set of faces.
    void DeleteFaces(std::set<FaceHandle> faces);

    //! Find isolated faces.
    std::set<FaceHandle> FindIsolatedFaces();

    //! Random normal displacement.
    void Disturb(float mu, float sigma);

    //! Simple mesh smoothing by vertex averaging.
    void SimpleSmooth(u_int n, bool boundary);

    //! Casts point into R4R data structure.
    vec3f Point(VertexHandle vh) const;

    //! Casts normal into R4R data structure.
    vec3f Normal(VertexHandle vh) const;

    //! Casts normal into R4R data structure.
    vec3f Normal(FaceHandle vh) const;

    //! Computes the bounding box of the mesh.
    CBoundingBox<float> BoundingBox() const;

    //! Average edge length.
    float MeanEdgeLength();

    //! Compute dual half edge.
    TriangleMesh::Normal DualEdgeVector(HalfedgeHandle heh);

    //! Gets handle of the half edge opposite to a vertex on a face.
    OpenMesh::HalfedgeHandle OppositeHalfedgeHandle(FaceHandle fh, VertexHandle vh);

    //! Computes the discrete gradient operator for vertex-based functions.
    CCSCMatrix<double,int> ComputeGradientOperator();

    /*! \brief Computes the discrete Laplace-Beltrami operator.
     *
     * \details This is much somehow much faster than squaring the
     * discrete gradient operator.
     *
     */
    CCSCMatrix<double,int> ComputeLaplaceBeltramiOperator();

    //! Cotangent of the angle opposite to a half edge.
    float CotanOppositeAngle(HalfedgeHandle heh);

    //! Save mesh to disk.
    void SaveToFile(const char* filename);

private:

    //! Previous half edge on boundary from vertex.
    OpenMesh::HalfedgeHandle PrevBoundaryHalfedge(VertexHandle vh);

    //! Next half edge on boundary from vertex.
    OpenMesh::HalfedgeHandle NextBoundaryHalfedge(VertexHandle vh);

    //! Next vertex on boundary from vertex.
    OpenMesh::VertexHandle PreviousBoundaryVertex(VertexHandle vh);

    //! Next vertex on boundary from vertex.
    OpenMesh::VertexHandle NextBoundaryVertex(VertexHandle vh);

    //! Checks for an obtuse triangle.
    bool IsTriangleObtuse(FaceHandle fh);

    //! Saves as VTK file.
    void SaveAsVTK(const char* filename);

};


#endif // TRIMESH_H_
