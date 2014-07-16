#include <math.h>
#include <fstream>
#include <string>
#include <typeinfo>

#include "cam.h"

using namespace std;

template<typename T>
CCamera<T>::CCamera():m_orthographic(false) {

    fill_n(m_size,2,0);
    fill_n(m_f,2,1.0);
    fill_n(m_c,2,0.0);

}

template<typename T>
CCamera<T>::CCamera(size_t w, size_t h, T fu, T fv, T cu, T cv, bool orthographic):m_orthographic(orthographic) {

    m_size[0] = w;
    m_size[1] = h;
    m_f[0] = fu;
    m_f[1] = fv;
    m_c[0] = cu;
    m_c[1] = cv;

}

template<typename T>
CVector<T,2> CCamera<T>::Project(const CVector<T,3>& x) const {

    CVector<T,2> xp;

    if(!m_orthographic) {

        if(x.Get(2)==0)
            throw runtime_error("CCamera::Project: Division by zero.");

        xp(0) = m_f[0]*(x.Get(0)/x.Get(2)) + m_c[0];
        xp(1) = m_f[1]*(x.Get(1)/x.Get(2)) + m_c[1];

    }
    else {

        xp(0) = m_f[0]*(x.Get(0)) + m_c[0];
        xp(1) = m_f[1]*(x.Get(1)) + m_c[1];

    }

    return xp;

}

template<typename U>
ostream& operator<< (ostream& os, const CCamera<U>& x) {

    if(x.m_orthographic)
        os << "orthographic camera: " << endl << endl;
    else
        os << "perspective camera: " << endl << endl;

    os << "# size" << endl;
    os << x.m_size[0] << " " << x.m_size[1] << endl;
    os << "# focal lengths" << endl;
    os << x.m_f[0] << " " << x.m_f[1] << endl;
    os << "# principle point" << endl;
    os << x.m_c[0] << " " << x.m_c[1] << endl;

    return os;

}

template ostream& operator <<(ostream&, const CCamera<float>& x);
template ostream& operator <<(ostream&, const CCamera<double>& x);

template<typename T>
CDenseArray<T> CCamera<T>::GetOpenGLProjectionMatrix(T znear, T zfar) const {

    CDenseArray<T> P(4,4);

    if(!m_orthographic) {

        P(0,0) = 2.0*m_f[0] / T(m_size[0]);
        P(0,2) = 2.0*(m_c[0]/ T(m_size[0])) - 1.0;
        P(1,1) = 2.0*m_f[1] /  T(m_size[1]);
        P(1,2) = 2.0*(m_c[1]/  T(m_size[1])) - 1.0;
        P(2,2) = (-zfar - znear)/(zfar - znear);
        P(2,3) = -2*zfar*znear/(zfar - znear);
        P(3,2) = -1.0;

    }
    else {

        P(0,0) = 2.0*m_f[0] / T(m_size[0]);
        P(0,3) = -2.0*(m_c[0]/ T(m_size[0])) + 1.0;
        P(1,1) = 2.0*m_f[1] /  T(m_size[1]);
        P(1,2) = -2.0*(m_c[1]/  T(m_size[1])) + 1.0;
        P(2,3) = (-zfar - znear)/(zfar - znear);
        P(2,2) = -2/(zfar - znear);
        P(3,3) = 1.0;

    }

    return P;

}

template class CCamera<float>;
template class CCamera<double>;

template<typename T>
CViewPoint<T>::CViewPoint():
    m_F(),
    m_Finv(),
    m_index(-1) {

}

template<typename T>
CViewPoint<T>::CViewPoint(const CRigidMotion<T,3>& F, u_int index):
    m_F(F),
    m_Finv(F),
    m_index(index) {

    m_Finv.Invert();

}

template<typename U>
ostream& operator << (ostream& os, const CViewPoint<U>& x) {

    // guarantees polymorphic behavior
    os << x.m_F << endl;
    os << x.m_Finv;

    return os;

}

template <typename T>
CVector<T,3> CViewPoint<T>::GetLocation() const {

    // F is always world->cam, but the origin is translation of cam->world Finv
    return m_Finv.GetTranslation();

}

template <typename T>
CVector<T,3> CViewPoint<T>::GetPrincipalAxis() const {

    // F is always world->cam, but the origin is translation of cam->world Finv
    return m_Finv.GetPrincipalAxis();

}

template <typename T>
void CViewPoint<T>::Translate(const CVector<T,3>& t) {

    for(u_int i=0; i<3; i++)
        m_Finv(i,3) += t.Get(i);

    m_F = m_Finv;
    m_F.Invert();

}

template <typename T>
void CViewPoint<T>::DifferentialTranslate(const CVector<T,3>& t) {

    // first convert into world coordinates
    CVector<T,3> tw = m_Finv.DifferentialTransform(t);

    // now modify frames
    for(u_int i=0; i<3; i++)
        m_Finv(i,3) += tw.Get(i);

    m_F = m_Finv;
    m_F.Invert();

}

template <typename T>
void CViewPoint<T>::Orbit(const CVector<T,3>& center, const CVector<T,3>& axis) {

    // first get location in world
    CVector<T,3> origin = this->GetLocation();

    // lever
    CVector<T,3> d = origin - center;

    // create rotation
    CRigidMotion<T,3> R(0,0,0,axis.Get(0),axis.Get(1),axis.Get(2));

    // rotate lever and convert it back into point
    d = R.Transform(d);
    origin = center + d;

    // rotate axes
    CTransformation<T,3> Finv = R*m_Finv;
    m_Finv = reinterpret_cast<CRigidMotion<T,3>& >(Finv);

    // set new origin
    for(u_int i=0; i<3; i++)
        m_Finv(i,3) = origin.Get(i);

    // also set inverse
    m_F = m_Finv;
    m_F.Invert();

}


template class CViewPoint<float>;
template class CViewPoint<double>;

template ostream& operator << (ostream& os, const CViewPoint<float>& x);
template ostream& operator << (ostream& os, const CViewPoint<double>& x);
