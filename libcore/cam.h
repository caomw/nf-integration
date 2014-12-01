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

#ifndef CAM_H_
#define CAM_H_

#include <stdlib.h>

#include "trafo.h"
#include "darray.h"

/*! \brief linear/projective camera model
 *
 * \details CAVEAT: In the orthographic model #m_f is not really a focal lengths
 * but rather a factor which converts metric units into pixels accounting
 * for the aspect ratio of the image. In this case, the units of #m_f are pixels
 * per mm or m (or whatever the unit for lengt in 3d may be). Otherwise the units
 * of #m_f are just pixels.
 */
template<typename T=double>
class CCamera {

public:

    //! Constructor.
    CCamera();

    //! Constructor pinhole projection.
    CCamera(size_t w, size_t h, T fu, T fv, T cu, T cv, bool orthographic = false);

    /*! \brief Projects a point into the image plane.
     *
     * \param[in] x point in camera coordinates
     * \returns point in pixel coordinates
     *
     */
    CVector<T,2> Project(const CVector<T,3>& x) const;

    //! Writes the camera parameters to a stream.
    template<typename U> friend std::ostream& operator << (std::ostream& os, const CCamera<U>& x);

    /*! \brief Projection matrix for use in OpenGL context.
     *
     * OpenGL matrices are also column-major, so the result can be directly sent
     * to the graphics card via a pointer to the data.
     *
     */
    CDenseArray<T> GetOpenGLProjectionMatrix(T znear, T zfar) const;

    //! Access to image size.
    CVector<size_t,2> GetSize() { return { m_size[0], m_size[1] }; }

private:

    bool m_orthographic;        //!< orthographic projection flag
    size_t m_size[2];           //!< size
    T m_f[2];   				//!< focal length
    T m_c[2];       			//!< principle point

};

/*! \brief viewpoint
 *
 *	\details Maintains a rigid-body transformation and its inverse to reduce computational
 *  complexity and facilitate viewpoint changes.
 *
 */
template<typename T=float>
class CViewPoint {

public:

    //! Constructor.
    CViewPoint();

    //! Constructor.
    CViewPoint(const CRigidMotion<T,3>& F, u_int index = -1);

    //! Writes the viewpoint parameters to a stream.
    template<typename U> friend std::ostream& operator << (std::ostream& os, const CViewPoint<U>& x);

    //! Access to the transformation.
    const CRigidMotion<T,3>& GetTransformation() const { return m_F; }

    //! Access to the inverse transformation.
    const CRigidMotion<T,3>& GetInverseTransformation() const { return m_Finv; }

    //! Set transformation.
    void SetTransformation(const CRigidMotion<T,3>& F) { m_F = F; m_Finv = F; m_Finv.Invert(); }

    //! Set transformation.
    void SetInverseTransformation(const CRigidMotion<T,3>& Finv) { m_Finv = Finv; m_F = Finv; m_F.Invert(); }

    //! Gets the location of the projection center in world coordinates.
    CVector<T,3> GetLocation() const;

    //! Gets the principal axis of the projection device in world coordinates.
    CVector<T,3> GetPrincipalAxis() const;

    //! Return index of the view.
    int GetIndex() const { return m_index; }

    /*! Translates the origin of the view.
     *
     * \param[in] t translation vector in world coordinates
     *
     */
    void Translate(const CVector<T,3>& t);

    /*! Translates the origin of the view.
     *
     * \param[in] t translation vector in camera coordinates
     *
     */
    void DifferentialTranslate(const CVector<T,3>& t);

    /*! Orbits the camera around a point.
     *
     * \param[in] center center of orbit in world coordinates
     * \param[in] axis rotation axis in world coordinates
     */
    void Orbit(const CVector<T,3>& center, const CVector<T,3>& axis);

protected:

    CRigidMotion<T,3> m_F;             //!< extrinsic parameters
    CRigidMotion<T,3> m_Finv;          //!< inverse extrinsic parameters
    int m_index;                       //!< view counter

};

#endif
