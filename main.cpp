#include <QApplication>

#include "mainwindow.h"
#include "linalg.h"


using namespace std;

// Visualization of field in main window -> flow class from server!, save image!!
//
//
//
//
//



int main(int argc, char *argv[]) {

    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();


//    vector<CCSCTriple<double,int> > entries;
//    entries.push_back(CCSCTriple<double,int>(0,0,3.0));
//    entries.push_back(CCSCTriple<double,int>(1,1,2.0));
//    entries.push_back(CCSCTriple<double,int>(2,1,4.0));   // this should not show in unsymmetric version
//    entries.push_back(CCSCTriple<double,int>(2,2,1.0));


//    CCSCMatrix<double,int> A(3,3,entries);
//    cout << A << endl;

//    CDenseVector<double> b(3);
//    b(0) = 1;
//    b(1) = 2;
//    b(2) = 3;

//    CCSCMatrix<double,int> AtA = CCSCMatrix<double,int>::Square(A);
//    cout << AtA << endl;

//    // expected solution: 0.11111111,  -2.5       ,  13.

//    CCholeskySolver<double> solver(AtA);

//    CDenseArray<double> x=solver.Solve(AtA,b);
//    cout << x << endl;

}
