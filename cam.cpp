#include <math.h>
#include <fstream>
#include <string>
#include <typeinfo>

#include "cam.h"

using namespace std;

template<typename T>
CCamera<T>::CCamera() {

    fill_n(m_size,2,0);
    fill_n(m_f,2,0);
    fill_n(m_c,2,0);
    fill_n(m_k,5,0);
    m_alpha = 0;

}

template<typename T>
CCamera<T>::CCamera(T fu, T fv, T cu, T cv) {

    fill_n(m_size,2,0);
    m_f[0] = fu;
    m_f[1] = fv;
    m_c[0] = cu;
    m_c[1] = cv;
    fill_n(m_k,5,0);
    m_alpha = 0;

}

template<typename T>
CCamera<T>::CCamera(size_t w, size_t h, T fu, T fv, T cu, T cv) {

    m_size[0] = w;
    m_size[1] = h;
    m_f[0] = fu;
    m_f[1] = fv;
    m_c[0] = cu;
    m_c[1] = cv;
    fill_n(m_k,5,0);
    m_alpha = 0;

}

template<typename T>
CCamera<T>::CCamera(size_t w, size_t h) {

    m_size[0] = w;
    m_size[1] = h;
    m_f[0] = min(m_size[0],m_size[1]);
    m_f[1] = m_f[0];
    m_c[0] = 0.5*w;
    m_c[1] = 0.5*h;

    for(size_t i=0;i<5;i++)
        m_k[i]=0;

    m_alpha = 0;

}

template<typename T>
CVector<T,2> CCamera<T>::Project(const CVector<T,3>& x) const {

    CVector<T,2> xn = { x.Get(0)/x.Get(2), x.Get(1)/x.Get(2) };

    // radial distortion
    T r = xn.Norm2();

    CVector<T,2> dx;
    dx(0) = 2*m_k[2]*xn.Get(0)*xn.Get(1) + m_k[3]*(r*r + 2*xn.Get(0)*xn.Get(0));
    dx(1) = m_k[2]*(r*r + 2*xn.Get(1)*xn.Get(1)) + 2*m_k[3]*xn.Get(0)*xn.Get(1);

    double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;
    CVector<T,2> xd = xn*fac + dx;

    // transform to pixel coordinates
    CVector<T,2> xp;
    xp(0) = m_f[0]*(xd.Get(0) + m_alpha*xd.Get(1)) + m_c[0];
    xp(1) = m_f[1]*xd.Get(1) + m_c[1];

    return xp;

}

template<typename T>
void CCamera<T>::Project(const CVector<T,3>& x, CVector<T,2>& u, CDenseArray<T>& J) const {

    u = Project(x);

    J(0,0) = m_f[0]/x.Get(2);
    J(0,2) = -(u.Get(0)-m_c[0])/x.Get(2);
    J(1,1) = m_f[1]/x.Get(2);
    J(1,2) = -(u.Get(1)-m_c[1])/x.Get(2);

}

template<typename T>
CVector<T,2> CCamera<T>::Flow(const CVector<T,3>& x, const CVector<T,3>& dx) const {

    CVector<T,2> u = Project(x);

    T J00 = m_f[0]/x.Get(2);
    T J02 = -(u.Get(0)-m_c[0])/x.Get(2);
    T J11 = m_f[1]/x.Get(2);
    T J12 = -(u(1)-m_c[1])/x.Get(2);

    CVector<T,2> du = { J00*dx.Get(0) + J02*dx.Get(2), J11*dx.Get(1) + J12*dx.Get(2) };

    return du;

}

template<typename T>
CVector<T,3> CCamera<T>::Normalize(const CVector<T,2>& u) const {

    // FIXME: add undistortion
    CVector<T,3> x;
    x(0) = (1.0/m_f[0])*(u.Get(0)-m_c[0]);
    x(1) = (1.0/m_f[1])*(u.Get(1)-m_c[1]);
    x(2) = 1.0;

    return x;

}

template<typename U>
ostream& operator<< (ostream& os, const CCamera<U>& x) {

    os << "# dims" << endl;
    os << x.m_size[0] << " " << x.m_size[1] << endl;
    os << "# focal lengths" << endl;
    os << x.m_f[0] << " " << x.m_f[1] << endl;
    os << "# principle point" << endl;
    os << x.m_c[0] << " " << x.m_c[1] << endl;
    os << "# radial distortion coefficients" << endl;
    os << x.m_k[0] << " " << x.m_k[1] << " " << x.m_k[2] << " " << x.m_k[3] << " " << x.m_k[4] << endl;
    os << "# skew coefficient" << endl;
    os << x.m_alpha;

    return os;

}

template ostream& operator <<(ostream&, const CCamera<float>& x);
template ostream& operator <<(ostream&, const CCamera<double>& x);

template<typename U>
istream& operator >> (istream& is, CCamera<U>& x) {

    string linebuffer;

    getline(is,linebuffer);

    is >> x.m_size[0];
    is >> x.m_size[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_f[0];
    is >> x.m_f[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_c[0];
    is >> x.m_c[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_k[0];
    is >> x.m_k[1];
    is >> x.m_k[2];
    is >> x.m_k[3];
    is >> x.m_k[4];
    is.get();

    getline(is,linebuffer);

    is >> x.m_alpha;

    return is;

}

template<typename T>
bool CCamera<T>::operator==(const CCamera<T>& cam) {

    bool result = (m_size[0]==cam.m_size[0] &&
                   m_size[1]==cam.m_size[1] &&
                   m_f[0]==cam.m_f[0] &&
                   m_f[1]==cam.m_f[1] &&
                   m_c[0]==cam.m_c[0] &&
                   m_c[1]==cam.m_c[1] &&
                   m_k[0]==cam.m_k[0] &&
                   m_k[1]==cam.m_k[1] &&
                   m_k[2]==cam.m_k[2] &&
                   m_k[3]==cam.m_k[3] &&
                   m_k[4]==cam.m_k[4]);

    return result;

}

template<typename T>
CDenseArray<T> CCamera<T>::GetProjectionMatrix() const {

    CDenseArray<T> P(3,3);

    P(0,0) = m_f[0];
    P(1,1) = m_f[1];
    P(0,2) = m_c[0];
    P(1,2) = m_c[1];
    P(2,2) = 1;

    return P;

}

template<typename T>
CDenseArray<T> CCamera<T>::GetOpenGLProjectionMatrix(T znear, T zfar) const {

    CDenseArray<T> P(4,4);

    P(0,0) = 2.0*m_f[0] / T(m_size[0]);
    P(0,2) = 2.0*(m_c[0]/ T(m_size[0])) - 1.0;
    P(1,1) = 2.0*m_f[1] /  T(m_size[1]);
    P(1,2) = 2.0*(m_c[1]/  T(m_size[1])) - 1.0;
    P(2,2) = (-zfar - znear)/(zfar - znear);
    P(2,3) = -2*zfar*znear/(zfar - znear);
    P(3,2) = -1.0;

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
