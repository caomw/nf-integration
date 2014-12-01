QT       += core gui opengl

TARGET = libcore
TEMPLATE = lib

SOURCES += \
    bbox.cpp \
    cam.cpp \
    darray.cpp \
    flowviz.cpp \
    libcore.cpp \
    linalg.cpp \
    nfield.cpp \
    trafo.cpp \
    trimesh.cpp \
    types.cpp \
    vecn.cpp \
    viewer.cpp

HEADERS += \
    bbox.h \
    cam.h \
    darray.h \
    flowviz.h \
    libcore.h \
    libcore_global.h \
    linalg.h \
    nfield.h \
    trafo.h \
    trimesh.h \
    types.h \
    vecn.h \
    viewer.h

QMAKE_CXXFLAGS += -std=c++0x -O3


unix:!symbian|win32 {

    packagesExist(OpenEXR) {

        CONFIG += link_pkgconfig
        PKGCONFIG += OpenEXR
        DEFINES += HAVE_EXR

    }
    else {
        error("Could not resolve mandatory dependence on OpenEXR...")
    }

}
