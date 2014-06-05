#ifndef FLOWVIZ_H
#define FLOWVIZ_H

#include <QImage>

#include "nfield.h"

#define N_COLORS 128

template<typename T>
class CFlowVisualization {

    const T COLOR_OFFSET = static_cast<T>(7.0/16.0);

public:

    //! Constructor.
    CFlowVisualization(T sat_threshold = 2.0, int sat_length = 20, int cycle_length = 10);


    //! Getters/settes.
    void SetSaturationThreshold(T sat_threshold) { m_sat_threshold = sat_threshold; }
    double GetSaturationThreshold() { return m_sat_threshold; }
    void SetFullSaturationLength(int sat_length) { m_sat_length = sat_length; }
    int GetFullSaturationLength() { return m_sat_length; }
    void SetCycleLength(int cycle_length) { m_cycle_length = cycle_length; }
    int GetCycleLength() { return m_cycle_length; }

    //! Compute color map for a normal vector field.
    QImage CalcDirectionEncoding(const CNormalField<T>& nf) const;

private:

    //! Initialize the colormap.
    void AssignColorMap();

    //! Map value to color.
    void GetColor(uchar& f_r_r, uchar& f_g_r, uchar& f_b_r, T f_index) const;

    //! Convert HSV to RGB value.
    void HSV2RGB(T f_h, T f_s, T  f_v, T& f_r, T& f_g, T& f_b);

    //! Compress values greater than threshold.
    double CapSaturation(T f_x) const;

    T m_sat_threshold;
    int m_sat_length;
    int m_cycle_length;
    uchar m_cmap_red[N_COLORS];   //!< look-up tables for color values
    uchar m_cmap_green[N_COLORS];
    uchar m_cmap_blue[N_COLORS];

};

#endif // FLOWVIZ_H
