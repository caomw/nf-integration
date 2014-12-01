#ifndef LIBCORE_GLOBAL_H
#define LIBCORE_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(LIBCORE_LIBRARY)
#  define LIBCORESHARED_EXPORT Q_DECL_EXPORT
#else
#  define LIBCORESHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // LIBCORE_GLOBAL_H
