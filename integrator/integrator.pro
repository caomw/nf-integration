QT       += core gui opengl

QMAKE_CXXFLAGS += -std=c++0x -O3

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = integrator
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

QT += testlib

unix:!symbian|win32 {

    LIBS += -L$$OUT_PWD/../libcore/ -llibcore


    # find OpenEXR
    packagesExist(OpenEXR) {

        CONFIG += link_pkgconfig
        PKGCONFIG += OpenEXR
        DEFINES += HAVE_EXR

    }
    else {
        error("Could not resolve mandatory dependence on OpenEXR...")
    }


    # find OpenMesh library
    OM = $$system(find /usr -name libOpenMeshCore* 2>/dev/null)
    isEmpty(OM) {
        error("Could not resolve mandatory dependency on OpenMesh...")
    } else {

        OM = $$first(OM)
        OMLIBPATH = $$dirname(OM)
        DEFINES += HAVE_OM

        LIBS += -L$$OMLIBPATH \
                -lOpenMeshCore \
                -lOpenMeshTools

        # add this so the binary knows the location of OM in non-standard path
        QMAKE_LFLAGS += -Wl,-rpath=$$OMLIBPATH

    }

    LIBS += -lcholmod

}

RESOURCES += \
    icon.qrc

INCLUDEPATH += $$PWD/../libcore
DEPENDPATH += $$PWD/../libcore
