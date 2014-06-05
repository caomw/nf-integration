#include <string.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "linalg.h"

using namespace std;

template<typename T,typename U>
bool CCSCTriple<T,U>::operator<(const CCSCTriple<T,U>& x) const {

    bool result;

    m_j>x.m_j ? result=true : result = (m_j==x.m_j && m_i>=x.m_i);

    return result;

}

template class CCSCTriple<float,size_t>;
template class CCSCTriple<double,size_t>;
template class CCSCTriple<float,int>;
template class CCSCTriple<double,int>;

template<typename T, typename U>
CCSCMatrix<T,U>::CCSCMatrix():
    m_nrows(),
    m_ncols(),
    m_colptr(new vector<U>(1)),
    m_rows(new vector<U>()),
    m_vals(new vector<T>()) {

    m_colptr->at(0) = 0;

}

template<typename T, typename U>
CCSCMatrix<T,U>::CCSCMatrix(size_t m, size_t n):
    m_nrows(m),
    m_ncols(n),
    m_colptr(new vector<U>(n+1)),
    m_rows(new vector<U>()),
    m_vals(new vector<T>()) {

    fill_n(m_colptr->begin(),n+1,0);
    m_rows->reserve(10*(m+n));
    m_vals->reserve(10*(m+n));
}

template<typename T, typename U>
CCSCMatrix<T,U>::CCSCMatrix(std::shared_ptr<std::vector<U> >& colptr, std::shared_ptr<std::vector<U> >& rows, std::shared_ptr<std::vector<T> >& vals):
    m_nrows(*max_element(rows->begin(),rows->end())),
    m_ncols(colptr->size()-1),
    m_colptr(colptr),
    m_rows(rows),
    m_vals(vals) {}

template<typename T, typename U>
CCSCMatrix<T,U>::CCSCMatrix(size_t m, size_t n, vector<CCSCTriple<T,U> >& data):
    m_nrows(m),
    m_ncols(n),
    m_colptr(new vector<U>()),
    m_rows(new vector<U>()),
    m_vals(new vector<T>()) {

    // make sure we do not have to resize the containers
    m_colptr->reserve(m_ncols+1);
    m_rows->reserve(data.size());
    m_vals->reserve(data.size());

    // make the vector a min heap (max heap by definition of comparison operator)
    std::make_heap(data.begin(),data.end());

    m_colptr->push_back(0);

    U nnz = 0;

    CCSCTriple<T,U> lastentry = data.front();

    // add trailing zeros if the first nonzero entry comes in cols > 0
    for(U k=0; k<lastentry.j(); k++)
        m_colptr->push_back(nnz);

    while(!data.empty()) {

        // col change?
        if(data.front().j()>lastentry.j())
            m_colptr->push_back(nnz);

        // check whether we have to add to the last element or insert a new one
        if(data.front()!=lastentry || m_rows->empty()) {

            m_vals->push_back(data.front().v());
            m_rows->push_back(data.front().i());
            nnz++;

        } else
            m_vals->back() += data.front().v();

        // store last element
        lastentry = data.front();

        // pop element and restore heap property
        std::pop_heap(data.begin(),data.end());
        data.pop_back();

    }

    // need to make sure that rowptr has nrows+1 elements
    while(m_colptr->size()!=m_ncols+1)
        m_colptr->push_back(nnz);

}

template<typename T, typename U>
void CCSCMatrix<T,U>::Eye() {

    // clear old contents
    m_colptr.reset(new vector<U>(m_ncols+1));
    fill_n(m_colptr->begin(),m_ncols+1,0);
    m_rows->clear();
    m_vals->clear();

    // make sure we do not have to resize the containers
    size_t ld = max(m_nrows,m_ncols);
    m_rows->reserve(ld);
    m_vals->reserve(ld);

    for(size_t i=1; i<=ld; i++) {

        m_colptr->at(i) = i;
        m_rows->push_back(U(i-1));
        m_vals->push_back(T(1.0));

    }

}

template<typename T, typename U>
void CCSCMatrix<T,U>::Scale(T scalar) {

    typename vector<T>::iterator it;

    for(it=m_vals->begin(); it!=m_vals->end(); ++it)
        *it *= scalar;

}

template<typename T, typename U>
void CCSCMatrix<T,U>::Resize(size_t m, size_t n) {

    if(m<=m_nrows || n<=m_ncols)
        return;

    m_nrows = m;

    U old_col_ptr = m_colptr->back();
    for(int k=0; k<n-m_ncols; k++)
        m_colptr->push_back(old_col_ptr);

    m_ncols = n;

}


template<typename V,typename W>
ostream& operator << (ostream& os, const CCSCMatrix<V,W>& x) {

    os << "CSC Matrix of size " << x.m_nrows << "x" << x.m_ncols << ":" << endl << endl;

    for(size_t i=0; i<x.m_ncols; i++) {

        for(size_t j=x.m_colptr->at(i); j<x.m_colptr->at(i+1); j++)
            os << "(" << x.m_rows->at(j) << "," << i << "," << x.m_vals->at(j) << ")" << endl;

        os << endl;

    }

    return os;

}

template ostream& operator<<(ostream& os, const CCSCMatrix<float,size_t>& x);
template ostream& operator<<(ostream& os, const CCSCMatrix<double,size_t>& x);
template ostream& operator<<(ostream& os, const CCSCMatrix<float,int>& x);
template ostream& operator<<(ostream& os, const CCSCMatrix<double,int>& x);

template<typename T, typename U>
template<class Matrix>
Matrix CCSCMatrix<T,U>::operator*(const Matrix& array) const {

    if(this->NRows()!=array.NRows())
        throw runtime_error("CSCMatrix::operator*: Dimension mismatch.");

    // zero matrix
    Matrix result(this->NCols(),array.NCols());

    // columns of the input matrix
    for(size_t k=0; k<array.NCols(); k++) {

        for(size_t i=0; i<m_ncols; i++) {

            for(size_t j=m_colptr->at(i); j<m_colptr->at(i+1); j++) {
                result(i,k) += m_vals->at(j)*array.Get(m_rows->at(j),k);

            }

        }

    }

    return result;

}

template CDenseArray<double> CCSCMatrix<double,size_t>::operator*(const CDenseArray<double>& array) const;
template CDenseVector<double> CCSCMatrix<double,size_t>::operator*(const CDenseVector<double>& array) const;
template CDenseArray<float> CCSCMatrix<float,size_t>::operator*(const CDenseArray<float>& array) const;
template CDenseVector<float> CCSCMatrix<float,size_t>::operator*(const CDenseVector<float>& array) const;
template CDenseArray<double> CCSCMatrix<double,int>::operator*(const CDenseArray<double>& array) const;
template CDenseVector<double> CCSCMatrix<double,int>::operator*(const CDenseVector<double>& array) const;
template CDenseArray<float> CCSCMatrix<float,int>::operator*(const CDenseArray<float>& array) const;
template CDenseVector<float> CCSCMatrix<float,int>::operator*(const CDenseVector<float>& array) const;

template<typename T, typename U>
template<class Matrix>
void CCSCMatrix<T,U>::Multiply(Matrix& out, const Matrix& in) const {

    if(this->NRows()!=in.NRows())
        throw runtime_error("CSCMatrix::Multiply: Dimension mismatch.");

    // columns of the input matrix
    for(size_t k=0; k<in.NCols(); k++) {

        for(size_t i=0; i<m_ncols; i++) {

            for(size_t j=m_colptr->at(i); j<m_colptr->at(i+1); j++) {
                out(i,k) += m_vals->at(j)*in.Get(m_rows->at(j),k);

            }

        }

    }

}

template void CCSCMatrix<double,size_t>::Multiply(CDenseArray<double>& out, const CDenseArray<double>& in) const;
template void CCSCMatrix<double,size_t>::Multiply(CDenseVector<double>& out, const CDenseVector<double>& array) const;
template void CCSCMatrix<float,size_t>::Multiply(CDenseArray<float>& out, const CDenseArray<float>& array) const;
template void CCSCMatrix<float,size_t>::Multiply(CDenseVector<float>& out, const CDenseVector<float>& array) const;
template void CCSCMatrix<double,int>::Multiply(CDenseArray<double>& out, const CDenseArray<double>& array) const;
template void CCSCMatrix<double,int>::Multiply(CDenseVector<double>& out, const CDenseVector<double>& array) const;
template void CCSCMatrix<float,int>::Multiply(CDenseArray<float>& out, const CDenseArray<float>& array) const;
template void CCSCMatrix<float,int>::Multiply(CDenseVector<float>& out, const CDenseVector<float>& array) const;


template<typename T, typename U>
CCSCMatrix<T,U> CCSCMatrix<T,U>::Square(const CCSCMatrix<T,U>& A, T lambda) {

    CCSCMatrix<T,U> AtA(A.NCols(),A.NCols());

    // non-zero element counter
    size_t counter = 0;
    T prod;

    // multiply col i of x with all other cols j>=i
    for(size_t i=0; i<A.NCols(); i++) {

        // this is one row
        for(size_t j=0; j<=i; j++) {

            prod = 0;
            size_t k = A.m_colptr->at(i);
            size_t l = A.m_colptr->at(j);

            // x.m_colptr->at(i) -> x.m_colptr->at(i+1) picks range in x.m_vals
            // x.m_colptr->at(j) -> x.m_colptr->at(j+1) picks range in x.m_vals
            // for cols i respectively j.
            // but overlap is determined by x.m_row

            if(A.m_rows->at(k)<=A.m_rows->at(A.m_colptr->at(j+1)-1) && A.m_rows->at(l)<=A.m_rows->at(A.m_colptr->at(i+1)-1)) {

                while(true) {

                    if(A.m_rows->at(k)==A.m_rows->at(l)) {
                        prod += A.m_vals->at(k)*A.m_vals->at(l);
                        k++;
                        l++;
                    }
                    else if(A.m_rows->at(k)<A.m_rows->at(l))
                        k++;
                    else
                        l++;

                    if (k==A.m_colptr->at(i+1) || l==A.m_colptr->at(j+1))
                        break;

                }

                if(prod!=0) {

                    AtA.m_rows->push_back(j);

                    if(i==j)
                        AtA.m_vals->push_back(prod+lambda);
                    else
                        AtA.m_vals->push_back(prod);

                    counter++;

                }



            }

        }

        // after row done, save current number of nnz
        // we have already initialized
        //m_rowptr->push_back(counter);
        //AtA.m_colptr->push_back(counter);
        AtA.m_colptr->at(i+1) = counter;

    }

    return AtA;

}


template<typename T, typename U>
void CCSCMatrix<T,U>::SaveToFile(const char* filename) {

    ofstream out(filename);

    out << m_nrows << " " << m_ncols << endl;

    // write m_ncols + 1 values
    for(typename vector<U>::const_iterator it=m_colptr->begin(); it!=m_colptr->end(); ++it)
        out << *it << endl;

    for(size_t k=0; k<m_vals->size(); k++)
        out << m_rows->at(k) << " " << m_vals->at(k) << endl;

    out.close();

}


template class CCSCMatrix<float,size_t>;
template class CCSCMatrix<double,size_t>;
template class CCSCMatrix<float,int>;
template class CCSCMatrix<double,int>;

template<typename T>
CCholeskySolver<T>::CCholeskySolver(CCSCMatrix<T, int>& AtA):
    m_c(),
    m_pattern(nullptr) {

    cholmod_start(&m_c);
    m_c.print = 5;
    //m_c.default_nesdis = 1;
    m_c.supernodal = CHOLMOD_SIMPLICIAL;

    // pre-analyze
    cholmod_sparse header = CCholeskySolver<T>::CreateSparseMatrixHeader(AtA);
    cout << "CHOLMOD:Analyzing..." << endl;
    m_pattern = cholmod_analyze(&header,&m_c);

}

template<typename T>
CCholeskySolver<T>::~CCholeskySolver() {

    cholmod_free_factor(&m_pattern,&m_c);
    cholmod_finish(&m_c);

}

template<>
cholmod_sparse CCholeskySolver<float>::CreateSparseMatrixHeader(CCSCMatrix<float,int>& M) {

    cholmod_sparse header;
    header.nrow = M.NRows();
    header.ncol = M.NCols();
    header.nzmax = M.NNz();
    header.p = M.GetColumnPointer()->data();
    header.i = M.GetRowIndices()->data();
    header.x = M.GetValues()->data();
    header.stype = 1;                           // use upper triangular part only
    header.itype = CHOLMOD_INT;
    header.xtype = CHOLMOD_REAL;
    header.dtype = CHOLMOD_SINGLE;
    header.packed = 1;
    header.sorted = 1;

    return header;

}

template<>
cholmod_sparse CCholeskySolver<double>::CreateSparseMatrixHeader(CCSCMatrix<double,int>& M) {

    cholmod_sparse header;
    header.nrow = M.NRows();
    header.ncol = M.NCols();
    header.nzmax = M.NNz();
    header.p = M.GetColumnPointer()->data();
    header.i = M.GetRowIndices()->data();
    header.x = M.GetValues()->data();
    header.stype = 1;                           // use upper triangular part only
    header.itype = CHOLMOD_INT;
    header.xtype = CHOLMOD_REAL;
    header.dtype = CHOLMOD_DOUBLE;
    header.packed = 1;
    header.sorted = 1;

    return header;

}

template<>
cholmod_dense CCholeskySolver<float>::CreateDenseMatrixHeader(CDenseArray<float>& M) {

    cholmod_dense header;
    header.nrow = M.NRows();
    header.ncol = M.NCols();
    header.d = M.NRows();                               // leading dimension = nrows
    header.nzmax = M.NElems();
    header.dtype = CHOLMOD_SINGLE;
    header.xtype = CHOLMOD_REAL;
    header.x = M.Data().get();

    return header;

}

template<>
cholmod_dense CCholeskySolver<double>::CreateDenseMatrixHeader(CDenseArray<double>& M) {

    cholmod_dense header;
    header.nrow = M.NRows();
    header.ncol = M.NCols();
    header.d = M.NRows();
    header.nzmax = M.NElems();
    header.dtype = CHOLMOD_DOUBLE;
    header.xtype = CHOLMOD_REAL;
    header.x = M.Data().get();

    return header;

}

template<typename T>
CDenseArray<T> CCholeskySolver<T>::Solve(CCSCMatrix<T,int>& AtA, CDenseArray<T>& AtB) {

    cholmod_sparse AtAhdr = CCholeskySolver<T>::CreateSparseMatrixHeader(AtA);
    cholmod_dense AtBhdr = CCholeskySolver<T>::CreateDenseMatrixHeader(AtB);

    cout << "CHOLMOD:Factorizing..." << endl;
    cholmod_factorize(&AtAhdr,m_pattern,&m_c);

    cout << "CHOLMOD:Solving by backsubstitution..." << endl;
    cholmod_dense* sol = cholmod_solve(CHOLMOD_A,m_pattern,&AtBhdr,&m_c);

    shared_ptr<T> pd(static_cast<T*>(sol->x));

    return CDenseArray<T>(AtB.NRows(),AtB.NCols(),pd);

}

template class CCholeskySolver<float>;
template class CCholeskySolver<double>;

