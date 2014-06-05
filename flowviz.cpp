#include "flowviz.h"

#include <assert.h>
#include <math.h>
#include <limits>

using namespace std;

template<typename T>
CFlowVisualization<T>::CFlowVisualization(T sat_threshold, int sat_length, int cycle_length):
    m_sat_threshold(sat_threshold),
    m_sat_length(sat_length),
    m_cycle_length(cycle_length) {

    this->AssignColorMap();

}

template<typename T>
double CFlowVisualization<T>::CapSaturation(T f_x) const {

    if (f_x>m_sat_threshold )
        return (m_sat_threshold*(1.0+log(f_x/m_sat_threshold)));
    else if (f_x<-m_sat_threshold)
        return ((-1.0)*m_sat_threshold*(1.0+log((-1.0)*f_x/m_sat_threshold)));
    else
        return f_x;

}

template<typename T>
void CFlowVisualization<T>::GetColor(uchar& f_r_r, uchar& f_g_r, uchar& f_b_r, T f_index) const {

    // get f_index into range [0..1]
    while (f_index < 0.0) f_index += 1.0;
    while (f_index > 1.0) f_index -= 1.0;

    T adr_f = f_index*static_cast<T>(N_COLORS - 2);

    int adr_i = static_cast<int>(adr_f);

    // check for overlow overflow or underflow
    if(N_COLORS - 2<adr_i || 0>adr_i)
        throw runtime_error("CFlowVisualization::GetColor: Index overflow...");

    // interpolate colour values linearly
    T w = adr_f - static_cast<T>(adr_i);

    // interpolate linearly; round by using +0.5 instead of cutoff during casting
    f_r_r = static_cast<uchar>(0.5 + (1.0 - w)*static_cast<T>(m_cmap_red[ adr_i ])
                               + w*static_cast<T>(m_cmap_red[adr_i + 1]));

    f_g_r = static_cast<uchar>(0.5 + (1.0 - w)*static_cast<T>(m_cmap_green[adr_i])
                               + w*static_cast<T>(m_cmap_green[adr_i + 1]));

    f_b_r = static_cast<uchar>(0.5 + (1.0 - w)*static_cast<T>(m_cmap_blue[adr_i])
                               + w*static_cast<T>(m_cmap_blue[adr_i + 1]));

}

template<typename T>
void CFlowVisualization<T>::AssignColorMap()
{

    int offsetColor = static_cast<int>( static_cast<T>(N_COLORS)*COLOR_OFFSET);

    for (int i = 0; i<N_COLORS; i++) {

        int pos = (i + offsetColor)%N_COLORS;

        T h = pos/static_cast<T>(N_COLORS);
        T s = static_cast<T>(1.0);
        T v = static_cast<T>(1.0);

        T r, g, b;

        this->HSV2RGB(h,s,v,r,g,b);

        //round values into unsigned char for colour images
        m_cmap_red[i] = static_cast<uchar>(255.0*r + 0.5);
        m_cmap_green[i] = static_cast<uchar>(255.0*g + 0.5);
        m_cmap_blue[i] = static_cast<uchar>(255.0*b + 0.5);

    }

}

template<typename T>
void  CFlowVisualization<T>::HSV2RGB(T f_h, T f_s, T f_v, T &f_r, T &f_g, T &f_b) {

    if(0.0>f_s || 1.0<f_s ||  0.0>f_v || 1.0<f_v)
        throw runtime_error("CFlowVisualization::HSV2RGB: Range error...");

    // get f_h into range [0..1]
    while (f_h<0.0) f_h += 1.0;
    while (f_h>1.0) f_h -= 1.0;

    if(fabs(f_v)<numeric_limits<T>::epsilon())
        f_r = f_g = f_b = static_cast<T>(0.0);
    else if(fabs(f_s)<numeric_limits<T>::epsilon() )
        f_r = f_g = f_b = f_v;
    else {

        const T hf  = 6.0f*f_h;
        const int i = static_cast<int>(floor(hf));
        const T f   = hf - i;
        const T vs  = f_v*f_s;

        const T pv  = f_v - vs;
        const T qv  = f_v - vs * f;
        const T tv  = f_v + pv - qv;

        switch(i) {

        case 0:
            f_r = f_v;
            f_g = tv;
            f_b = pv;
            break;
        case 1:
            f_r = qv;
            f_g = f_v;
            f_b = pv;
            break;
        case 2:
            f_r = pv;
            f_g = f_v;
            f_b = tv;
            break;
        case 3:
            f_r = pv;
            f_g = qv;
            f_b = f_v;
            break;
        case 4:
            f_r = tv;
            f_g = pv;
            f_b = f_v;
            break;
        case 5:
            f_r = f_v;
            f_g = pv;
            f_b = qv;
            break;
        case 6:
            f_r = f_v;
            f_g = tv;
            f_b = pv;
            break;
        case -1:
            f_r = f_v;
            f_g = pv;
            f_b = qv;
            break;
        default:
            break;
        }
    }
}

template<typename T>
QImage CFlowVisualization<T>::CalcDirectionEncoding(const CNormalField<T>& nf) const {

    QImage img(nf.m_raw_data.m_ncols,nf.m_raw_data.m_nrows,QImage::Format_RGB888);

    for (size_t i = 0; i <nf.m_raw_data.m_nrows; i++) {

        for (size_t j = 0; j<nf.m_raw_data.m_ncols; j++) {

            // get input value
            CVector<T,3> val = nf.m_raw_data.Get(i,j);

            // so that this also works with deflectometric images
            if(val.Norm2()>0)
                val.Normalize();

            bool o = nf.m_mask.Get(i,j)<128;

            // initialize new values as black
            uchar r, g, b;

            // valid and non-zeros flow vectors
            if (!o) {

                // determine angle of motion vector
                T radian = atan2(val.Get(1),val.Get(o))/M_PI; // in [-1,1[

                // determine colour
                this->GetColor(r,g,b,0.5*(radian+1.0));  // in [0, 1[

                // determine saturation
                T saturation = 1.0-acos(val.Get(2))/M_PI; // in [0, 1[

                r = static_cast<uchar>(255.5 - saturation*(255.0 - static_cast<T>(r)));
                g = static_cast<uchar>(255.5 - saturation*(255.0 - static_cast<T>(g)));
                b = static_cast<uchar>(255.5 - saturation*(255.0 - static_cast<T>(b)));

            }
            else {

                r = 0;
                g = 0;
                b = 0;

            }

            img.setPixel(j,i,qRgb(r,g,b));

        }

    }

    return img;

}

/*void CFlowVisualization::calcDirectionEncoding(QImage& img, Array2D<Rgba>& flowfield, size_t w, size_t h) const {

    const double f_maxSaturation = capSaturation(m_fullSaturationLength);

    // flow smaller than this in both components is displayed as zero
    const double EPS = numeric_limits<double>::epsilon();

    for (size_t y = 0; y < h; y++)
    {
        for (size_t x = 0; x < w; x++ )
        {

            // get input value
            Rgba val = flowfield[y][x];
            float u = (float)val.r;
            float v = (float)val.g;
            bool o = (bool)val.b;

            // initialize new values as black
            unsigned char r = 0;
            unsigned char g = 0;
            unsigned char b = 0;

            // valid and non-zeros flow vectors
            if (( !o ) && ((EPS<fabs(u)) || (EPS < fabs(v)))) {

                // determine angle of motion vector
                double radian = inversPi*atan2(static_cast<double>(v),static_cast<double>(u)); // in [-1,1[

                // determine colour
                getColor(r,g,b,0.5*(radian+1.0)); // in [0, 1[

                // determine saturation
                double saturation = capSaturation(sqrt(u*u+v*v));

                // threshold to [0,1]
                saturation = min(1.0,saturation/f_maxSaturation);

                r = static_cast<unsigned char>( 255.5 - ( saturation * (255.0 - static_cast<double>( r) ) ) );
                g = static_cast<unsigned char>( 255.5 - ( saturation * (255.0 - static_cast<double>( g) ) ) );
                b = static_cast<unsigned char>( 255.5 - ( saturation * (255.0 - static_cast<double>( b) ) ) );
            }
            else if (!o) // flow which is too small is white
            {
                r = 255;
                g = 255;
                b = 255;

            }

            img.setPixel(x,y,qRgb(r,g,b));

        }

    }

}*/

template class CFlowVisualization<float>;
template class CFlowVisualization<double>;
