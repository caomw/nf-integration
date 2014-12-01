//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the software accompanying the following publication
//
// @INPROCEEDINGS{Balzeretal2014,
// author = {J. Balzer and D. Acevedo-Feliz and S. Soatto and S. H\"ofer and M. Hadwiger and J. Beyerer},
// title = {Cavlectometry: Towards Holistic Reconstruction of Large Mirror Objects},
// booktitle = {International Conference on 3D Vision (3DV)},
// year = {2014}}
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
