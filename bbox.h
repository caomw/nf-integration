#ifndef BBOX_H
#define BBOX_H

#include <vector>

#include "vecn.h"

template<typename T>
class CBoundingBox {

public:

    //! Constructor.
    CBoundingBox():m_lower(), m_upper() {}

    //! Constructor.
    CBoundingBox(const CVector<T,3>& lower, const CVector<T,3>& upper);

    //! Get corners.
    std::vector<CVector<T,3> > Corners() const;

    //! Access to the lower corner.
    const CVector<T,3>& Lower() { return m_lower; }

    //! Access to the upper corner.
    const CVector<T,3>& Upper() { return m_upper; }

    //! Barycenter.
    CVector<T,3>  Barycenter() { return (m_upper + m_lower)*0.5; }

private:

    CVector<T,3> m_lower;          //!< lower corner
    CVector<T,3> m_upper;          //!< upper corner

};



#endif // BBOX_H
