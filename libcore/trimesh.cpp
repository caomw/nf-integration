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

#include <assert.h>
#include <limits>
#include <algorithm>

#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>

#include "trimesh.h"

using namespace std;
using namespace OpenMesh;

#define ZERO_TOL 1e-10
#define COTAN_EPSILON 1e-10

CTriangleMesh::CTriangleMesh():
        TriangleMesh() {

}

float CTriangleMesh::CotanOppositeAngle(HalfedgeHandle heh) {

    float sp, cotanAngle;

	cotanAngle = 0;

    VertexHandle vh1, vh2, vh3;

	vh1 = from_vertex_handle(heh);
	vh2 = to_vertex_handle(heh);
	vh3 = opposite_vh(heh);

    if (vh3!=TriangleMesh::InvalidVertexHandle) {

        TriangleMesh::Point v1, v2, v3;
		v1 = point(vh1);
		v2 = point(vh2);
		v3 = point(vh3);

        TriangleMesh::Normal v3v1, v3v2;
		v3v1 = v1 - v3;
		v3v2 = v2 - v3;

		v3v1 = v3v1.normalize();
		v3v2 = v3v2.normalize();
		sp = OpenMesh::dot(v3v1,v3v2);

        cotanAngle = sp / sqrt(1.0f - sp * sp);

        if(std::isnan(cotanAngle))
            cotanAngle = sp / sqrt(1 - sp * sp + COTAN_EPSILON);

	}

	return cotanAngle;

}

float CTriangleMesh::FaceArea(FaceHandle fh) {

    return calc_sector_area(fh_iter(fh));

}

float CTriangleMesh::VertexQuadratureWeight(VertexHandle vh) {

    VertexFaceIter vf_it;

    float w = 0;

    for (vf_it = this->vf_iter(vh); vf_it; ++vf_it)
        w += this->FaceArea(vf_it.handle());

    return w/3.0;

}

bool CTriangleMesh::IsTriangleObtuse(FaceHandle fh) {

    bool obtuse = false;

    FaceHalfedgeIter foh_it;

    for (foh_it=this->fh_iter(fh); foh_it; ++foh_it) {

        if(this->calc_sector_angle(foh_it)>0.5*M_PI) {

            obtuse = true;
            break;

        }

    }

    return obtuse;

}

float CTriangleMesh::VoronoiArea(VertexHandle vh) {

    VertexOHalfedgeIter voh_it;

    float voronoiArea = 0;
    float denom, angle, triArea, previousLength, currentLength, previousWeight, currentWeight;

    TriangleMesh::Normal previousEdge, currentEdge;

    for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(voh_it.handle()) != TriangleMesh::InvalidFaceHandle) {

            calc_edge_vector(voh_it, currentEdge);
            calc_edge_vector(cw_rotated_halfedge_handle(voh_it), previousEdge);

            denom = currentEdge.norm()*previousEdge.norm();

             // degenerate triangle, no area gain
            if(denom > 0) {

                angle = acos(OpenMesh::dot(currentEdge, previousEdge)/denom);

                    // obtuse triangle: area
                    if(IsTriangleObtuse(this->face_handle(voh_it.handle()))) {
                    //if(false) {

                        triArea = this->FaceArea(this->face_handle(voh_it.handle()));

                        // obtuse angle at center vertex
                        if(angle>0.5*M_PI)
                            voronoiArea = voronoiArea + 0.5*triArea;
                        else
                            voronoiArea = voronoiArea + 0.25*triArea;

                    }
                    else {

                        previousLength = this->calc_edge_sqr_length(this->cw_rotated_halfedge_handle(voh_it));
                        currentLength = this->calc_edge_sqr_length(voh_it);
                        previousWeight = this->CotanOppositeAngle(voh_it.handle());
                        currentWeight = this->CotanOppositeAngle(this->opposite_halfedge_handle(voh_it));

                        triArea = (currentLength*currentWeight + previousLength*previousWeight)/8.0;
                        voronoiArea = voronoiArea + triArea;

                    }

            }

        }

    }

    return voronoiArea;
}

vec3f CTriangleMesh::Barycenter(FaceHandle fh) {

    vec3f cog;

    float valence = 0;

    TriangleMesh::FaceVertexIter fv_it;

	for (fv_it = fv_iter(fh); fv_it; ++fv_it) {

        cog += this->Point(fv_it);
    	++valence;

	}

	cog = cog/valence;

	return cog;

}

vec3f CTriangleMesh::Barycenter() const {

    vec3f barycenter;
    size_t counter = 0;

    TriangleMesh::VertexIter v_it;

    for (v_it=vertices_begin(); v_it!=vertices_end(); ++v_it) {

        barycenter = barycenter + Point(v_it);
        counter++;

    }

    return (1.0/float(counter))*barycenter;

}

TriangleMesh::Normal CTriangleMesh::DualEdgeVector(HalfedgeHandle heh) {

    assert(face_handle(heh)!=TriangleMesh::InvalidFaceHandle);

    HalfedgeHandle a, b;
	a = prev_halfedge_handle(heh);
	b = prev_halfedge_handle(a);

    TriangleMesh::Normal aVec, bVec;
	calc_edge_vector(a, aVec);
	calc_edge_vector(b, bVec);
	bVec = -bVec;

    TriangleMesh::Normal result = -(aVec*CotanOppositeAngle(a) + bVec*CotanOppositeAngle(b));

	return result;

}

HalfedgeHandle CTriangleMesh::OppositeHalfedgeHandle(FaceHandle fh, VertexHandle vh) {

    HalfedgeHandle result = TriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::FaceHalfedgeIter fh_it;

	for (fh_it = fh_iter(fh); fh_it; ++fh_it) {

		if(opposite_vh(fh_it) == vh) {

			result = fh_it;
			break;

		}

	}

	return result;

}


HalfedgeHandle CTriangleMesh::PrevBoundaryHalfedge(VertexHandle vh) {

    HalfedgeHandle result;

	if(!is_boundary(vh))
        return TriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::VertexOHalfedgeIter voh_it;

	for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(voh_it.handle())==TriangleMesh::InvalidFaceHandle) {

			result = voh_it.handle();
			break;

		}

	}

	return opposite_halfedge_handle(result);

}

HalfedgeHandle CTriangleMesh::NextBoundaryHalfedge(VertexHandle vh) {

    HalfedgeHandle result;

	if(!is_boundary(vh))
        return CTriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::VertexOHalfedgeIter voh_it;

	for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(opposite_halfedge_handle(voh_it.handle()))==TriangleMesh::InvalidFaceHandle) {

			result = voh_it.handle();

			break;

		}

	}

	return result;

}

VertexHandle CTriangleMesh::PreviousBoundaryVertex(VertexHandle vh) {

    HalfedgeHandle res = PrevBoundaryHalfedge(vh);

    if(res != TriangleMesh::InvalidHalfedgeHandle)
		return from_vertex_handle(res);
	else
        return TriangleMesh::InvalidVertexHandle;

}

VertexHandle CTriangleMesh::NextBoundaryVertex(VertexHandle vh) {

    HalfedgeHandle res = NextBoundaryHalfedge(vh);

    if(res != TriangleMesh::InvalidHalfedgeHandle)
		return to_vertex_handle(res);
	else
        return TriangleMesh::InvalidVertexHandle;

}

TriangleMesh::Normal CTriangleMesh::ExteriorNormal(VertexHandle vh) {

    HalfedgeHandle hehNext = NextBoundaryHalfedge(vh);
    HalfedgeHandle hehPrev = PrevBoundaryHalfedge(vh);

	return -((DualEdgeVector(hehNext).normalize() + DualEdgeVector(hehPrev).normalize())*0.5).normalize();

}

void CTriangleMesh::UniformMeshRefinement(size_t n) {

    OpenMesh::Subdivider::Uniform::LoopT<TriangleMesh> sd;

    sd.attach(*this);
    sd(n);
    sd.detach();

    cout << "No of vertices: " << n_vertices() << endl;

}

void CTriangleMesh::EdgeCollapse(double threshold) {

	size_t before = n_vertices();

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::EdgeIter e_it;
    vector<EdgeHandle> toDel;

	double weight;

	for (e_it = edges_begin(); e_it != edges_end(); ++e_it) {

		weight = CotanOppositeAngle(halfedge_handle(e_it,0)) + CotanOppositeAngle(halfedge_handle(e_it,1));

        if(weight>threshold || std::isnan(weight))
			toDel.push_back(e_it.handle());

	}

	for(size_t i=0; i<toDel.size(); i++) {

		if(is_collapse_ok(halfedge_handle(toDel[i],0)))
			collapse(halfedge_handle(toDel[i],0));

		if(is_collapse_ok(halfedge_handle(toDel[i],1)))
			collapse(halfedge_handle(toDel[i],1));

	}


	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();

	size_t after = n_vertices();

	cout << "No of deleted vertices: " << before-after << std::endl;

}

bool CTriangleMesh::SmoothBoundary(size_t n) {

    VertexHandle init, actual, prev, next;
    TriangleMesh::Point newPt;

	init = FindBoundaryVertex();

	for(size_t k=0; k<n; k++) {

		cout << "Smoothing step k=" << k << endl;

		actual = init;

		do {

            if (actual == TriangleMesh::InvalidVertexHandle) {

				cout << "WARNING: Mesh has no boundary..." << endl;
				return 1;

			}

			prev = PreviousBoundaryVertex(actual);
			next = NextBoundaryVertex(actual);

			newPt = (point(prev) + point(next))*0.5;

			point(actual) = newPt;

			actual = next; //NextBoundaryVertex(actual);

		}
		while (actual!=init);

	}

	return 0;

}

VertexHandle CTriangleMesh::FindBoundaryVertex() {

    VertexHandle result = CTriangleMesh::InvalidVertexHandle;

    TriangleMesh::VertexIter v_it, v_end(vertices_end());

	for (v_it = vertices_begin(); v_it != v_end; ++v_it) {

		if(is_boundary(v_it)) {

			result = v_it.handle();
			break;

		}

	}

	return result;

}

void CTriangleMesh::DeleteVertices(set<VertexHandle> vertices) {

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::VertexIter v_it, v_end(vertices_end());

	for (v_it = vertices_begin(); v_it != v_end; ++v_it) {

		//if(property(toDelete, v_it))
		if(vertices.find(v_it) != vertices.end())
			delete_vertex(v_it, true);

	}

	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();
	//remove_property(toDelete);

}

set<FaceHandle> CTriangleMesh::FindIsolatedFaces() {

    set<FaceHandle> handles;

    TriangleMesh::ConstFaceIter f_it, f_end(faces_end());


    for (f_it = faces_begin(); f_it != f_end; ++f_it) {

        if(!this->is_valid_handle(f_it))
            handles.insert(f_it.handle());
    }

    return handles;

}

void CTriangleMesh::DeleteFaces(set<FaceHandle> faces) {

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::FaceIter f_it, f_end(faces_end());

	for (f_it = faces_begin(); f_it != f_end; ++f_it) {

		if(faces.find(f_it) != faces.end())
			delete_face(f_it, true);

	}

	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();


}

void CTriangleMesh::Disturb(float mu, float sigma) {

    // generate random normal displacement
    vecf noise(this->n_vertices());
    noise.RandN(mu,sigma);

    TriangleMesh::VertexIter v_it;

    this->update_normals();

    for (v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
        this->point(v_it) = this->point(v_it) + this->normal(v_it)*noise.Get(v_it.handle().idx());

    this->update_normals();

}

void CTriangleMesh::SimpleSmooth(u_int n, bool boundary) {

    OpenMesh::VPropHandleT<TriangleMesh::Point> newPoints;
    this->add_property(newPoints);

    TriangleMesh::VertexIter v_it, v_end(this->vertices_end());
    TriangleMesh::VertexVertexIter vv_it;

    TriangleMesh::Point sumPoints;

    for (int i=0; i<n; ++i) {

        for (v_it = this->vertices_begin(); v_it != v_end; ++v_it) {

            sumPoints = TriangleMesh::Point(0.0,0.0,0.0);

            for (vv_it = this->vv_iter(v_it); vv_it; ++vv_it)
                sumPoints = sumPoints + this->point(vv_it);

            this->property(newPoints, v_it) = sumPoints/float(this->valence(v_it));

        }

        for (v_it = this->vertices_begin(); v_it != v_end; ++v_it) {

            if(!boundary) {

                if(!this->is_boundary(v_it))  // shrinking?
                    this->point(v_it) = this->property(newPoints, v_it);

            }
            else
                this->point(v_it) = this->property(newPoints, v_it);

        }

    }

    this->remove_property(newPoints);

}

CCSCMatrix<double, int> CTriangleMesh::ComputeGradientOperator() {

    vector<CCSCTriple<double,int> > entries;
    entries.reserve(10*this->n_faces());
    // row pointer
    size_t row = 0;

    TriangleMesh::FaceIter f_it;
    TriangleMesh::FaceVertexIter fv_it;

    for(f_it=this->faces_begin(); f_it!=this->faces_end(); ++f_it) {

        // face area
        float A =  this->calc_sector_area(this->fh_iter(f_it));
        A = sqrt(A);

        if(A>0) {

            for (fv_it = this->fv_iter(f_it); fv_it; ++fv_it) {

                TriangleMesh::Normal grad = this->DualEdgeVector(this->OppositeHalfedgeHandle(f_it,fv_it));

                // gradient is (0.5/A) times sqrt(A) for integration = (0.5/A)*sqrt(A) = 1/2*sqrt(A)
                grad *= (0.5/A);

                // find col index
                size_t col = fv_it.handle().idx();

                for(uint i=0; i<3; i++) {

                    if(abs(grad[i])>ZERO_TOL)
                        entries.push_back(CCSCTriple<double,int>(row+i,col,grad[i]));

                }

            }

        }

        row += 3;

    }

    return CCSCMatrix<double,int>(3*this->n_faces(),this->n_vertices(),entries);

}

CCSCMatrix<double, int> CTriangleMesh::ComputeLaplaceBeltramiOperator() {

    vector<CCSCTriple<double,int> > entries;
    entries.reserve(10*this->n_vertices());

    TriangleMesh::VertexIter v_it, v_end(this->vertices_end());
    TriangleMesh::VertexOHalfedgeIter voh_it;

    size_t  rowptr = 0;

    for (v_it = this->vertices_begin(); v_it != v_end; ++v_it) {

         size_t j0 = v_it.handle().idx();

         double w0 = 0;

         for (voh_it = this->voh_iter(v_it); voh_it; ++voh_it) {

             size_t j1 = this->to_vertex_handle(voh_it.handle()).idx();

             double weight = this->CotanOppositeAngle(voh_it.handle())
                             + this->CotanOppositeAngle(this->opposite_halfedge_handle(voh_it.handle()));

             if(weight>ZERO_TOL) {
                w0 += 0.5*weight;
                entries.push_back(CCSCTriple<double,int>(rowptr,j1,-0.5*weight));
             }

        }

        entries.push_back(CCSCTriple<double,int>(rowptr,j0,w0));

        rowptr++;

    }

    return CCSCMatrix<double,int>(this->n_vertices(),this->n_vertices(),entries);

}


vec3f CTriangleMesh::Point(VertexHandle vh) const {

    vec3f point = { this->point(vh)[0], this->point(vh)[1], this->point(vh)[2] };

    return point;

}

vec3f CTriangleMesh::Normal(VertexHandle vh) const {

    vec3f normal = { this->normal(vh)[0], this->normal(vh)[1], this->normal(vh)[2] };

    return normal;

}

vec3f CTriangleMesh::Normal(FaceHandle fh) const {

    vec3f normal = { this->normal(fh)[0], this->normal(fh)[1], this->normal(fh)[2] };

    return normal;

}

CBoundingBox<float> CTriangleMesh::BoundingBox() const {

    vec3f lower, upper;

    for(u_int i=0; i<3; i++) {

        lower(i) = std::numeric_limits<float>::max();
        upper(i) = -std::numeric_limits<float>::max();

    }

    TriangleMesh::VertexIter v_it;

    for (v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {

        vec3f point = this->Point(v_it);

        for(u_int i=0; i<3; i++) {

            if(point.Get(i)>=upper.Get(i))
                upper(i) = point.Get(i);

            if(point.Get(i)<=lower.Get(i))
                lower(i) = point.Get(i);

        }

    }

    return CBoundingBox<float>(lower,upper);

}

float CTriangleMesh::MeanEdgeLength() {

    TriangleMesh::EdgeIter e_it;

    float length = 0;

    for (e_it = this->edges_begin(); e_it != this->edges_end(); ++e_it)
        length += this->calc_edge_sqr_length(e_it);

    return length/float(this->n_edges());

}

void CTriangleMesh::SaveToFile(const char* filename) {

    string fn(filename);

    string ext = fn.substr(fn.size()-3,3);

    if(ext==string("vtk")) {
        this->SaveAsVTK(filename);
    }
    else if(ext==string("off")) {

        OpenMesh::IO::Options wopt;
        wopt += OpenMesh::IO::Options::FaceColor;

        if(!OpenMesh::IO::write_mesh(*this,filename,wopt))
            throw std::runtime_error("Could not write file.");

    }
    else if(ext==string("ply") || ext==string("stl") || ext==string("obj")) {

        if(!OpenMesh::IO::write_mesh(*this,filename))
            throw std::runtime_error("Could not write file.");

    }

}

void CTriangleMesh::SaveAsVTK(const char* filename) {

    ofstream out(filename);

    out << "# vtk DataFile Version 2.0" << endl;
    out << "created by blz" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << "POINTS " << this->n_vertices() << " double" << endl;

    TriangleMesh::VertexIter v_it;

    for(v_it=this->vertices_begin(); v_it!=this->vertices_end(); ++v_it)
        out << this->point(v_it)[0] << " " << this->point(v_it)[1] << " " << this->point(v_it)[2] << endl;

    out << "CELLS " << this->n_faces() << " " << 4*this->n_faces() << endl;

    TriangleMesh::ConstFaceIter f_it;
    TriangleMesh::ConstFaceVertexIter fv_it;
    for (f_it=this->faces_begin(); f_it!=this->faces_end(); ++f_it) {

        out << "3 ";

        for (fv_it = this->cfv_iter(f_it.handle()); fv_it; ++fv_it)
            out << fv_it.handle().idx() << " ";

        out << endl;

    }

    out << "CELL_TYPES " << this->n_faces() << endl;

    for(size_t k=0; k<this->n_faces(); k++)
        out << "5" << endl;

    OpenMesh::FPropHandleT<vec3f> nt;

    if(this->get_property_handle(nt,"nt")) {

        out << "CELL_DATA " << this->n_faces() << endl;
        out << "VECTORS \"nt\" double " << endl;

        for (f_it=this->faces_begin(); f_it!=this->faces_end(); ++f_it) {
            vec3f n = this->property(nt,f_it);
            out << n.Get(0) << " " << n.Get(1) << " " << n.Get(2) << endl;
        }

    }

    out.close();

}
