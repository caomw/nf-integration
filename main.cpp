#include <QApplication>

#include "unittests.h"
#include "mainwindow.h"
#include "linalg.h"


using namespace std;


int main(int argc, char *argv[]) {

#ifndef HAVE_UNITTESTS
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
#else
    CLinearAlgebraTest test0;
    QTest::qExec(&test0);
    CDenseArrayTest test1;
    QTest::qExec(&test1);

#endif

}
