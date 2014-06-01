#ifndef CAM_H_
#define CAM_H_

#include <stdlib.h>

#include "trafo.h"
#include "darray.h"

/*! \brief pinhole camera
 *
 *
 *
 */
template<typename T=double>
class CPinholeCam {

public:

    //! Constructor.
    CPinholeCam();

    //! Constructor.
    CPinholeCam(size_t w, size_t h);

    //! Constructor.
    CPinholeCam(T fu, T fv, T cu, T cv);

    //! Constructor.
    CPinholeCam(size_t w, size_t h, T fu, T fv, T cu, T cv);

    /*! \brief Projects a point into the image plane.
     *
     * \param[in] x point in camera coordinates
     * \returns point in pixel coordinates
     *
     */
    CVector<T,2> Project(const CVector<T,3>& x) const;

    /*! \brief Computes the projection of a point into the image plane and its Jacobian.
     *
     * \param[in] x point in camera coordinates
     * \param[out] u point in pixel coordinates
     * \param[out] J Jacobian of the projection mapping in u
     *
     */
    void Project(const CVector<T,3>& x, CVector<T,2>& u, CDenseArray<T>& J) const;

    /*! \brief Converts a pixel into a viewing direction.
     *
     * \param[in] u location w.r.t. the pixel coordinate system
     * \returns direction vector, normalized s.t. \f$z\f$-component equals \f$1\f$
     *
     */
    CVector<T,3> Normalize(const CVector<T,2>& u) const;

    /*! \brief Projects differential motion to optical flow vector.
     *
     * \param[in] x point in camera coordinates
     * \param[in] dx differential motion vector in 3d
     * \returns direction vectors
     *
     */
    CVector<T,2> Flow(const CVector<T,3>& x, const CVector<T,3>& dx) const;

    //! Writes the camera parameters to a stream.
    template<typename U> friend std::ostream& operator << (std::ostream& os, const CPinholeCam<U>& x);

    //! Reads the camera parameters from a stream.
    template<typename U> friend std::istream& operator >> (std::istream& is, CPinholeCam<U>& x);

    //! Checks if two cameras are the same.
    bool operator==(const CPinholeCam<T>& cam);

    //! Access to image size.
    CVector<size_t,2> GetSize() const {  return { m_size[0], m_size[1] }; }

    //! Projection matrix.
    CDenseArray<T> GetProjectionMatrix() const;

    /*! \brief Projection matrix for use in OpenGL context.
     *
     * OpenGL matrices are also column-major, so the result can be directly sent
     * to the graphics card via a pointer to the data.
     *
     */
    CDenseArray<T> GetOpenGLProjectionMatrix(T znear, T zfar) const;

private:

    size_t m_size[2];			//!< pixel size
    T m_f[2];   				//!< focal length
    T m_c[2];       			//!< principle point
    T m_alpha;          		//!< skew coefficient
    T m_k[5];               	//!< distortion coefficients

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
