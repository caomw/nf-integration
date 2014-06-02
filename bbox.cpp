#include "bbox.h"

using namespace std;

template<typename T>
CBoundingBox<T>::CBoundingBox(const CVector<T,3>& lower, const CVector<T,3>& upper):
    m_lower(lower),
    m_upper(upper) {}


template<typename T>
vector<CVector<T,3> > CBoundingBox<T>::Corners() const {

    vector<CVector<T,3> > result(8);

    result[0] = m_lower;
    result[1] = { m_upper.Get(0), m_lower.Get(1), m_lower.Get(2) };
    result[2] = { m_lower.Get(0), m_upper.Get(1), m_lower.Get(2) };;
    result[3] = { m_upper.Get(0), m_upper.Get(1), m_lower.Get(2) };;
    result[4] = m_upper;
    result[5] = { m_upper.Get(0), m_lower.Get(1), m_upper.Get(2) };
    result[6] = { m_lower.Get(0), m_upper.Get(1), m_upper.Get(2) };;
    result[7] = { m_upper.Get(0), m_upper.Get(1), m_upper.Get(2) };;

    return result;

}

template class CBoundingBox<double>;
template class CBoundingBox<float>;
