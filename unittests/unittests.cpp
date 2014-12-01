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

#include <vector>

#include "unittests.h"

using namespace std;

CLinearAlgebraTest::CLinearAlgebraTest(QObject* parent):
    QObject(parent),
    m_A(),
    m_b(3) {

}

void CLinearAlgebraTest::init() {

    vector<CCSCTriple<double,int> > entries;
    entries.push_back(CCSCTriple<double,int>(0,0,2.0));
    entries.push_back(CCSCTriple<double,int>(1,1,1.0));
    entries.push_back(CCSCTriple<double,int>(2,1,4.0));
    entries.push_back(CCSCTriple<double,int>(2,2,2.0));

    m_A = CCSCMatrix<double,int>(3,3,entries);

    m_b(0) = 2;
    m_b(1) = 6;
    m_b(2) = 2;

}

void CLinearAlgebraTest::testNNz() {

    QCOMPARE(m_A.NNz(),4ul);

}

void CLinearAlgebraTest::testSolver() {

    CCSCMatrix<double,int> AtA = CCSCMatrix<double,int>::Square(m_A);

    CCholeskySolver<double> solver(AtA);
    CDenseArray<double> x = solver.Solve(AtA,m_b);

    QCOMPARE(x.Get(0,0),0.5f);
    QCOMPARE(x.Get(1,0),2.0f);
    QCOMPARE(x.Get(2,0),-3.5f);

}

void CLinearAlgebraTest::testMultiplication() {

    CDenseVector<double> Atb(m_A.NCols()+1);
    m_A.Multiply(Atb,m_b);

    QCOMPARE(Atb.Get(0,0),4.0f);
    QCOMPARE(Atb.Get(1,0),14.0f);
    QCOMPARE(Atb.Get(2,0),4.0f);
    QCOMPARE(Atb.Get(3,0),0.0f);

}

void CLinearAlgebraTest::testNormalization() {

    // square matrix, resize, and add normalization
    CCSCMatrix<double,int> AtA = CCSCMatrix<double,int>::Square(m_A);
    AtA.Resize(4,4);
    size_t nnz_old = AtA.NNz();

    // Atb
    CDenseVector<double> Atb(m_A.NCols()+1);
    m_A.Multiply(Atb,m_b);

    vector<int>* rowptr = AtA.GetRowIndices().get();
    vector<int>* colptr = AtA.GetColumnPointer().get();
    vector<double>* valptr = AtA.GetValues().get();

    for(size_t k=0; k<3; k++) {
        rowptr->push_back(k);
        valptr->push_back(1.0);
    }
    colptr->at(4) = nnz_old + 3;

    CCholeskySolver<double> solver(AtA);
    CDenseArray<double> x = solver.Solve(AtA,Atb);

    QCOMPARE(x.Get(0,0),1.66666667f);
    QCOMPARE(x.Get(1,0),3.33333333f);
    QCOMPARE(x.Get(2,0),-5.0f);
    QCOMPARE(x.Get(3,0),-2.66666667f);

}

void CLinearAlgebraTest::cleanup(){


}

CDenseArrayTest::CDenseArrayTest(QObject* parent):
    QObject(parent),
    m_double_array(2,3) {



}

void CDenseArrayTest::init() {

    // this is called before each test!!!
    m_double_array(0,0) = 1.0;
    m_double_array(0,1) = 2.0;
    m_double_array(0,2) = 3.0;
    m_double_array(1,0) = 4.0;
    m_double_array(1,1) = 5.0;
    m_double_array(1,2) = 6.0;

}

void CDenseArrayTest::testAccess() {

    QCOMPARE(m_double_array.NRows(),2ul);
    QCOMPARE(m_double_array.NCols(),3ul);
    QCOMPARE(m_double_array.Get(0,1),2.0);
    QCOMPARE(m_double_array.Get(0,2),3.0);
    QCOMPARE(m_double_array.Get(1,0),4.0);
    QCOMPARE(m_double_array(0,0),1.0);
    QCOMPARE(m_double_array(1,1),5.0);
    QCOMPARE(m_double_array(1,2),6.0);

    m_double_array.Transpose();

    QCOMPARE(m_double_array.NRows(),3ul);
    QCOMPARE(m_double_array.NCols(),2ul);
    QCOMPARE(m_double_array.Get(2,0),3.0);
    QCOMPARE(m_double_array.Get(0,0),1.0);
    QCOMPARE(m_double_array.Get(2,1),6.0);
    QCOMPARE(m_double_array(1,0),2.0);
    QCOMPARE(m_double_array(1,1),5.0);
    QCOMPARE(m_double_array(0,1),4.0);

    m_double_array.Transpose();

}

void CDenseArrayTest::testFileIO() {

    m_double_array.WriteToFile("test.r4r");

    CDenseArray<double> in;
    in.ReadFromFile("test.r4r");

    QCOMPARE(in.Get(0,0),m_double_array.Get(0,0));
    QCOMPARE(in.Get(0,1),m_double_array.Get(0,1));
    QCOMPARE(in.Get(0,2),m_double_array.Get(0,2));
    QCOMPARE(in.Get(1,0),m_double_array.Get(1,0));
    QCOMPARE(in.Get(1,1),m_double_array.Get(1,1));
    QCOMPARE(in.Get(1,2),m_double_array.Get(1,2));

    system("rm test.r4r");

}

void CDenseArrayTest::cleanup() {



}
