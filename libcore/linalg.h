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

#ifndef LINALG_H_
#define LINALG_H_

#include <map>
#include <vector>
#include <stdlib.h>
#include <memory>
#include <suitesparse/cholmod.h>

#include "darray.h"

/*!
 * \brief matrix-market CSC triple
 */
template<typename T,typename U = size_t>
class CCSCTriple {

public:

    //! Constructor.
    CCSCTriple(U i, U j, T v):m_i(i),m_j(j),m_v(v){}

    //! Lexicographic ordering.
    bool operator<(const CCSCTriple& x) const;

    //! Tests for equal location.
    bool operator==(const CCSCTriple& x) const { return m_i==x.m_i && m_j==x.m_j; }

    //! Tests for equal location.
    bool operator!=(const CCSCTriple& x) const { return m_i!=x.m_i || m_j!=x.m_j; }

    //! Access to row index.
    const U& i() const { return m_i; }

    //! Access to column index.
    const U& j() const { return m_j; }

    //! Access to value.
    const T& v() const { return m_v; }

private:

    U m_i;
    U m_j;
    T m_v;

};

/*!
* matrix in compressed-column storage format
*
*  TODO: add a symmetric flag (similar to CHOLMOD structure, for printing/multiplication)
*
*/
template<typename T,typename U = size_t>
class CCSCMatrix {

public:

    //! Standard constructor.
    CCSCMatrix();

    //! Constructor.
    CCSCMatrix(size_t m, size_t n);

    //! Resize.
    void Resize(size_t m, size_t n);

    /*! \brief Constructor which takes MatrixMarket triples as input.
     *
     * The triples can contain duplicate entries. A heap sort and implicit summation
     * is performed internally in the constructor.
     *
     */
    CCSCMatrix(size_t m, size_t n, std::vector<CCSCTriple<T,U> >& data);

    /*! \brief Constructor for external assembly.
     */
    CCSCMatrix(std::shared_ptr<std::vector<U> >& colptr, std::shared_ptr<std::vector<U> >& rows, std::shared_ptr<std::vector<T> >& vals);

    //! Erases the matrix and replaces it with the identity.
    void Eye();

    //! Access number of cols.
    size_t NRows() const {  return m_nrows; }

    //! Access number of cols.
    size_t NCols() const { return m_ncols; }

    //! In-place scalar multiplication.
    void Scale(T scalar);

    //! Counts the number of non-zero entries.
    size_t NNz() const { return m_vals->size(); }

    //! Writes matrix to a stream.
    template<typename V,typename W> friend std::ostream& operator << (std::ostream& os, const CCSCMatrix<V,W>& x);

    //! Multiplies the transpose (!) of the current object with an array from the right.
    template<class Matrix> Matrix operator*(const Matrix& array) const;

    //! Multiplies the transpose (!) of the current object with an array from the right.
    template<class Matrix> void Multiply(Matrix& out, const Matrix& in) const;

    //! Form the square of a matrix in an efficient way.
    static CCSCMatrix<T,U> Square(const CCSCMatrix<T,U>& A, T  lambda = 0);

    //! Low-level access to the column pointer.
    std::shared_ptr<std::vector<U> > GetColumnPointer() const { return m_colptr; }

    //! Low-level access to the row index array.
    std::shared_ptr<std::vector<U> > GetRowIndices() const { return m_rows; }

    //! Low-level access to the value array.
    std::shared_ptr<std::vector<T> > GetValues() const { return m_vals; }

    //! Write data to file.
    void SaveToFile(const char* filename);

private:

    size_t m_nrows;                                         //!< number of rows
    size_t m_ncols;                                         //!< number of cols
    std::shared_ptr<std::vector<U> > m_colptr;              //!< indicates the beginning of columns in #m_val and #m_rows
    std::shared_ptr<std::vector<U> > m_rows;                //!< row index
    std::shared_ptr<std::vector<T> > m_vals;                //!< value

};


/*!
 * wrapper for the sparse Cholesky decomposition implemented in the
 * CHOLMOD library
 */
template<typename T>
class CCholeskySolver {

public:

    //! Constructor.
    CCholeskySolver(CCSCMatrix<T,int>& AtA);

    //! Destructor.
    ~CCholeskySolver();

    /*! Perform Cholesky factorization then solve by back-substitution.
     *
     * TODO: Separate both steps?
     */
    CDenseArray<T> Solve(CCSCMatrix<T,int>& AtA, CDenseArray<T>& AtB);

    //! Create a sparse CHOLMOD matrix header.
    static cholmod_sparse CreateSparseMatrixHeader(CCSCMatrix<T,int>& M);

    //! Create a dense CHOLMOD matrix header.
    static cholmod_dense CreateDenseMatrixHeader(CDenseArray<T>& M);

private:

    cholmod_common m_c;
    cholmod_factor* m_pattern;

};

#endif
