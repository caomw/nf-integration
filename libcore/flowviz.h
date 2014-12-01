///////////////////////////////////////////////////////////////////////////////
//
// This file is an adaption of a file part of the HCI-Correspondence Estimation
// Benchmark.
//
// More information on this benchmark can be found under:
//
//       http://hci.iwr.uni-heidelberg.de/Benchmarks/
//
// Copyright (C) 2011  <Sellent, Lauer>
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

    //! Getters/setters.
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
