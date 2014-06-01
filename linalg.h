#ifndef LINALG_H_
#define LINALG_H_

#include <map>
#include <vector>
#include <stdlib.h>
#include <memory>
#include <suitesparse/cholmod.h>

#include "darray.h"

/* Cholmod wrapper class
 * member is structure
 * update of numeric values and right-hand side
 * MEMBERS: are headers!!!!!
 * setter:
 * SetMatrix(const CSymmetricCSCArray& )
 * SetRightHandSide(CDenseVector)
 *
 * Was ist mit reassembly? Lohnt sich ws nicht, da nur die Allokation gespart wird.
 * NORMALIZATION: Add columns of ONES to symmetric CSC matrix!!!
 * DenseVector -> Multiplication -> create a bigger vector for output!
 */

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
    static CCSCMatrix<T,U> Square(const CCSCMatrix<T,U>& A);

    //! Low-level access to the column pointer.
    U* GetColumnPointer() const { return &(*m_colptr.get())[0]; }

    //! Low-level access to the row index array.
    U* GetRowIndices() const { return &(*m_rows.get())[0]; }

    //! Low-level access to the value array.
    T* GetValues() const { return &(*m_vals.get())[0]; }

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
