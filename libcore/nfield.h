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

#ifndef NFIELD_H
#define NFIELD_H

#include <memory>
#include <stdlib.h>

#include "vecn.h"
#include "cam.h"

// forward declarations
template<typename T> class CNormalField;
template<typename T> class CFlowVisualization;

/*! \brief storage of 3-vector fields over a regular Cartesian lattice
 *
 * \details We cannot use CDenseArray here because of the memory layout of OpenEXR
 * files.
 */
template <typename T>
class CRawData {

    friend class CNormalField<T>;
    friend class CFlowVisualization<T>;

public:

    //! Standard constructor.
    CRawData():m_nrows(0),m_ncols(0),m_data(nullptr) {}

    //! Parametrized constructor.
    CRawData(size_t nrows, size_t ncols);

    //! Get value at grid points.
    CVector<T,3> Get(size_t i, size_t j) const;

    //! Get value at intermediate grid points using bilinear interpolation.
    CVector<T,3> Get(const CVector<T,2>& x) const;

private:

    size_t m_nrows;
    size_t m_ncols;
    std::shared_ptr<T> m_data;

};

/*! \brief normal field plus vantage point and visibility information
 *
 */
template<typename T>
class CNormalField {

    friend class CFlowVisualization<T>;

public:

    //! Standard constructor.
    CNormalField():m_cam(),m_viewpoint(),m_raw_data(), m_mask() {}

    /*! \brief Read normal field from an OpenEXR file.
     *
     * \details Depending on the projection matrix read from file, a orthographic
     * or perspective camera model is used, the latter when the lower-right entry of
     * \f$K_33=0\f$, the former if \f$K_33=1\f$.
     *
     *
     *
     */
    void ReadFromFile(const char* filename);

    //! Project to image and compute value using bilinear interpolation.
    CVector<T,3> Get(const CVector<T,3>& x) const;

    //! Get mask.
    uint GetMask(const CVector<T,3>& x) const;

    //! Get scene point and compute deflectometric normal.
    CVector<T,3> GetDeflectometricNormal(const CVector<T,3>& x) const;

    //! Access to camera.
    const CCamera<T>& GetCam() const { return m_cam; }

    //! Access to vantage point.
    const CViewPoint<T>& GetViewpoint() const { return m_viewpoint; }

protected:

    CCamera<T> m_cam;                   //!< intrinsic camera parameters
    CViewPoint<T> m_viewpoint;          //!< extrinsic camera parameters
    CRawData<T> m_raw_data;             //!< raw normal data
    CDenseArray<uint> m_mask;           //!< optional mask

};

#endif // NFIELD_H
