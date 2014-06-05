#ifndef UNITTESTS_H
#define UNITTESTS_H


#include <QtTest/QtTest>
#include "linalg.h"

class CLinearAlgebraTest:public QObject {

  Q_OBJECT

public:

  explicit CLinearAlgebraTest(QObject* parent = nullptr);

private:
    CCSCMatrix<double,int> m_A;
    CDenseVector<double> m_b;

private slots:

  //! Init tests.
  void init();

  //! Test number of non-zeros.
  void testNNz();

  //! Test multiplication with transposed sparse matrix.
  void testMultiplication();

  //! Test basic solver.
  void testSolver();

  //! Test normalization.
  void testNormalization();

  //! Clean-up scaffolding.
  void cleanup();

};


class CDenseArrayTest:public QObject {

  Q_OBJECT

public:

  explicit CDenseArrayTest(QObject* parent = nullptr);

private:
    CDenseArray<double> m_double_array;

private slots:

  //! Init tests.
  void init();

  //! Test access with and without transposition.
  void testAccess();

  //! Test file I/O.
  void testFileIO();

  //! Clean-up scaffolding.
  void cleanup();

};
#endif // DARRAYTEST_H
