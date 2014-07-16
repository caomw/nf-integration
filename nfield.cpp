#include <stdexcept>

using namespace std;

#ifdef HAVE_EXR
#include <ImfArray.h>
#include <ImfStandardAttributes.h>
#include <ImfMatrixAttribute.h>
#include <ImfAttribute.h>
#include <ImfInputFile.h>
#include <ImathBox.h>
#include <ImfChannelList.h>

using namespace Imf;
using namespace Imath;
#endif

#include "nfield.h"

template<typename T>
CRawData<T>::CRawData(size_t nrows, size_t ncols):
    m_nrows(nrows),
    m_ncols(ncols),
    m_data(new T[3*nrows*ncols]) {

}

template<typename T>
CVector<T,3> CRawData<T>::Get(size_t i, size_t j) const {

    T* pd = m_data.get();
    size_t nelems = m_nrows*m_ncols;
    return { pd[i*m_ncols+j], pd[nelems+i*m_ncols+j], pd[2*nelems+i*m_ncols+j] };

}

template<typename T>
CVector<T,3> CRawData<T>::Get(const CVector<T,2>& x) const {

    T u, v;
    u = x.Get(0);
    v = x.Get(1);

    int i = (int)floor(v);
    int j = (int)floor(u);

    if(i<0 || i>static_cast<int>(m_nrows)-2 || j<0 || j>static_cast<int>(m_ncols)-2)
        return CVector<T,3>();

    T vd = v - i;
    T ud = u - j;

    CVector<T,3> A00, A01, A10, A11, A0, A1;
    A00 = this->Get(i,j);
    A01 = this->Get(i,j+1);
    A10 = this->Get(i+1,j);
    A11 = this->Get(i+1,j+1);

    A0 = A00*(1-ud) + A01*ud;
    A1 = A10*(1-ud) + A11*ud;

    return A0*(1-vd) + A1*vd;

}

template class CRawData<float>;
template class CRawData<double>;

template<>
void CNormalField<float>::ReadFromFile(const char* filename) {

#ifdef HAVE_EXR
    InputFile file (filename);

    // see if the file contains the right channels
    const ChannelList& channels = file.header().channels();
    if(channels.findChannel("x")==nullptr ||
       channels.findChannel("y")==nullptr ||
       channels.findChannel("z")==nullptr)
        throw runtime_error("OpenEXR file is missing one of the channels x,y,z.");

    Box2i dw = file.header().dataWindow();
    size_t width = dw.max.x - dw.min.x + 1;
    size_t height = dw.max.y - dw.min.y + 1;

    // make sure there is enough memory
    if(width!=m_raw_data.m_ncols || height!=m_raw_data.m_nrows)
        m_raw_data = CRawData<float>(height,width);


    FrameBuffer frameBuffer;
    frameBuffer.insert ("x",Slice(FLOAT,
                                  (char*)m_raw_data.m_data.get(),
                                  sizeof(float)*1,
                                  sizeof(float)*width,                          // stride
                                  1, 1,                                         // subsampling
                                  0.0));                                        // fill
    frameBuffer.insert ("y",Slice(FLOAT,
                                  (char*)(m_raw_data.m_data.get()+m_raw_data.m_nrows*m_raw_data.m_ncols),
                                  sizeof(float)*1,
                                  sizeof(float)*width,
                                  1, 1,
                                  0.0));
    frameBuffer.insert ("z",Slice(FLOAT,
                                  (char*)(m_raw_data.m_data.get()+2*m_raw_data.m_nrows*m_raw_data.m_ncols),
                                  sizeof(float)*1,
                                  sizeof(float)*width,
                                  1, 1,
                                  0.0));

    if(channels.findChannel("mask")!=nullptr) {

        m_mask = CDenseArray<uint>(height,width);
        frameBuffer.insert ("mask",Slice(UINT,
                                         (char*)(m_mask.Data().get()),
                                          sizeof(uint)*1,
                                          sizeof(uint)*width,
                                          1, 1,
                                          0.0));

    }

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y,dw.max.y);

    // read vantage points
    const M44fAttribute* Fa = file.header().findTypedAttribute <M44fAttribute> ("F");
    const M33fAttribute* Ka = file.header().findTypedAttribute <M33fAttribute> ("K");

    Matrix44<float> F;
    Matrix33<float> K;

    // check if we find them in the file
    if(Fa!=nullptr && Ka!=nullptr) {

        F = Fa->value();
        K = Ka->value();

        CRigidMotion<float,3> trafo;
        for(uint i=0; i<3; i++) {
            for (uint j=0; j<4; j++)
                trafo(i,j) = F[i][j];
        }

        m_viewpoint = CViewPoint<float>(trafo);
        m_cam = CCamera<float>(width,height,K[0][0],K[1][1],K[0][2],K[1][2],!K[2][2]>0);

    }
    else
        throw runtime_error("Could not find camera parameters in OpenEXR file.");

#else
#error OpenEXR is required to compile this code...
#endif

}

template<typename T>
CVector<T,3> CNormalField<T>::Get(const CVector<T,3>& x) const {

    CVector<T,3> xc = m_viewpoint.GetTransformation().Transform(x);
    CVector<T,2> xp = m_cam.Project(xc);
    return m_raw_data.Get(xp);

}

template<typename T>
uint CNormalField<T>::GetMask(const CVector<T,3>& x) const {

    CVector<T,3> xc = m_viewpoint.GetTransformation().Transform(x);
    CVector<T,2> xp = m_cam.Project(xc);

    int u = int(xp.Get(0)+0.5);
    int v = int(xp.Get(1)+0.5);

    if(v>=0 && v<m_mask.NRows() && u>=0 && u<m_mask.NCols())
        return m_mask.Get(v,u);
    else
        return 0;

}


template<typename T>
CVector<T,3> CNormalField<T>::GetDeflectometricNormal(const CVector<T,3>& x) const {

    CVector<T,3> o = m_viewpoint.GetLocation();

    // transform to local coordinates
    CVector<T,3> xl = m_viewpoint.GetTransformation().Transform(x);

    // project
    CVector<T,2> xp = m_cam.Project(xl);

    // check whether the pixel is valid
    if(xp.Get(0)>=0 && xp.Get(1)>=0 && xp.Get(0)<m_mask.NCols() && xp.Get(1)<m_mask.NRows()
       && m_mask.Get(size_t(xp.Get(1)),size_t(xp.Get(0)))>=200) {

        CVector<T,3> l = m_raw_data.Get(xp);

        CVector<T,3> s = x - o;
        CVector<T,3> r = x - l;

        s.Normalize();
        r.Normalize();
        CVector<T,3> n = s + r;

        if(InnerProduct(n,s)<0.1)
            n = 0.0*n;
        else
            n = (-1.0)*n;

        return n;
    }
     else
        return CVector<T,3>();

}

template class CNormalField<float>;
template class CNormalField<double>;


