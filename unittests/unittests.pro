QT += testlib
QT -= gui

QMAKE_CXXFLAGS += -std=c++0x -O3 -msse4

CONFIG   += console
CONFIG   -= app_bundle

TARGET = unittests
TEMPLATE = app

SOURCES += main.cpp \
           unittests.cpp

HEADERS += unittests.h

unix:!symbian|win32 {

    LIBS += -L$$OUT_PWD/../libcore/ -llibcore \
            -lcholmod

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


}

INCLUDEPATH += $$PWD/../libcore
DEPENDPATH += $$PWD/../libcore

