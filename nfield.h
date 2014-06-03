#ifndef NFIELD_H
#define NFIELD_H

#include <memory>
#include <stdlib.h>

#include "vecn.h"
#include "cam.h"

// forward declaration
template<typename T> class CNormalField;

/*! \brief storage of 3-vector fields over a regular Cartesian lattice
 *
 * \details We cannot use CDenseArray here because of the memory layout of OpenEXR
 * files.
 */
template <typename T>
class CRawData {

    friend class CNormalField<T>;

public:

    //! Standard constructor.
    CRawData():m_nrows(0),m_ncols(0),m_data(nullptr) {}

    //! Parametrized constructor.
    CRawData(size_t nrows, size_t ncols);

    //! Get value at grid points.
    CVector<T,3> Get(size_t i, size_t j);

    //! Get value at intermediate grid points using bilinear interpolation.
    CVector<T,3> Get(const CVector<T,2>& x);

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
    virtual void ReadFromFile(const char* filename);

    //! Project to image and compute value using bilinear interpolation.
    virtual CVector<T,3> Get(const CVector<T,3>& x);

    //! Access to camera.
    const CCamera<T>& GetCam() { return m_cam; }

    //! Access to vantage point.
    const CViewPoint<T>& GetViewpoint() { return m_viewpoint; }

protected:

    CCamera<T> m_cam;
    CViewPoint<T> m_viewpoint;
    CRawData<T> m_raw_data;
    CDenseArray<uint> m_mask;

};

template<typename T>
class CDeflectometricNormalField:public CNormalField<T> {

public:

    //! Standard constructor.
    CDeflectometricNormalField():CNormalField<T>::CNormalField() {}

    //! Project point to image plane and compute normal according to law and reflection.
    CVector<T,3> Get(const CVector<T,3>& x);


};


#endif // NFIELD_H
