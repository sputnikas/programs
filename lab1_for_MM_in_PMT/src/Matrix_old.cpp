#include "Matrix.h"

#ifndef STRINGBUFFER
#define STRINGBUFFER
std::string strbuf;
#endif


Matrix::Matrix(int n_row, int n_column, float *elements, float scale = 1.0f):n_row(n_row),n_column(n_column){
    lupped = 0;
    p = new int[1];
    if (IsMatrix()){
        a = new float[n_row*n_column];
        if (elements!=NULL_MASSIVE){
            for (int i = 0; i<n_row; i++){
                for (int j = 0; j<n_column; j++){
                    a[i*n_column+j] = elements[i*n_column+j]*scale;
                }
            }
        }
    } else {
        a = new float[1];
    }
};

Matrix::Matrix(int n_row, int n_column):n_row(n_row),n_column(n_column){
    lupped = 0;
    p = new int[1];
    a = new float[n_row*n_column];
    for (int i = 0; i<n_row; i++){
        for (int j = 0; j<n_column; j++)
            a[i*n_column+j] = 0;
    }
};

Matrix::~Matrix(){
    delete [] p;
    delete [] a;
    n_row = 0;
    n_column = 0;
};

Matrix::Matrix(){
    a = new float[1];
    p = new int[1];
    n_row = 0;
    n_column = 0;
};

const Matrix& Matrix::operator=(const Matrix &m){
    if (this != &m){
        delete [] a;
        delete [] p;
        lupped = m.lupped;
        n_row = m.n_row;
        n_column = m.n_column;
        a = new float[n_row*n_column];
        (lupped==0) ? p = new int[n_row] : p = new int[1];
        for (int i = 0; i<n_row; i++){
            if (lupped){
                p[i] = m.p[i];
            }
            for (int j = 0; j<n_column; j++)
                a[i*n_column+j] = m.a[i*n_column+j];
        }
    }
    return *this;
};

Matrix Matrix::operator-(){
    if (lupped) FromLUP();
    return Matrix(n_row, n_column, a, -1.0f);;
};

Matrix Matrix::operator*(float scale){
    if (lupped) FromLUP();
    return Matrix(n_row, n_column, a, scale);
};

Matrix Matrix::operator/(float scale){
    if (lupped) FromLUP();
    if (scale!=0)
        return Matrix(n_row, n_column, a, 1/scale);
    else return Matrix();
};

Matrix Matrix::operator*(const Matrix &m){
    if (lupped) FromLUP();
    float *b;
    if ((n_column==m.n_row)&&(IsMatrix())&&(m.IsMatrix())){
        b = new float[n_row * m.n_column];
        for (int i=0; i<n_row; i++){
            for (int j=0; j<m.n_column; j++){
                b[i*m.n_column+j] = 0;
                for (int k=0; k<n_column; k++){
                    b[i*m.n_column+j] += a[i*n_column+k] * m.a[k*m.n_column+j];
                }
            }
        }
        return Matrix(n_row, m.n_column, b);
    } else return Matrix();
};

Matrix Matrix::operator+(const Matrix &m){
    if (lupped) FromLUP();
    float *b;
    if ((n_column==m.n_column)&&(n_row==m.n_row)&&(IsMatrix())){
        b = new float[n_row*n_column];
        for (int i=0; i<n_row; i++){
            for (int j=0; j<n_column; j++){
                b[i*n_column+j] = m.a[i*n_column+j]+a[i*n_column+j];
            }
        }
        return Matrix(n_row, m.n_column, b);
    } else return Matrix();
};

Matrix Matrix::operator-(const Matrix &m){
    if (lupped) FromLUP();
    float *b;
    if ((n_column==m.n_column)&&(n_row==m.n_row)&&(IsMatrix())){
        b = new float[n_row*n_column];
        for (int i=0; i<n_row; i++){
            for (int j=0; j<n_column; j++){
                b[i*n_row+j] = a[i*n_row+j]-m.a[i*n_row+j];
            }
        }
        return Matrix(n_row, m.n_column, b);
    } else return Matrix();
};

bool Matrix::operator==(const Matrix &m){
    if (lupped) FromLUP();
    bool result = (n_row == m.n_row)&&(n_column == m.n_column);
    if (result)
        for (int i = 0; i<n_row; i++)
            for (int j = 0; j<n_column; j++)
                result = result&&(a[i*n_column+j] == m.a[i*n_column+j]);
    return result;
};

bool Matrix::IsMatrix()const{
    return !((n_row==0)||(n_column==0));
};

bool Matrix::IsQuad()const{
    return ((n_row==n_column)&&(IsMatrix()));
};

void Matrix::Transpose(){
    if (lupped) FromLUP();
    float *b;
    int n;
    if (IsMatrix()){
        b = a;
        n = n_row; n_row = n_column; n_column = n;
        a = new float[n_row*n_column];
        for (int i=0; i<n_row; i++){
            for (int j=0; j<n_column; j++){
                a[i*n_column+j] = b[j*n_row+i];
            }
        }
        delete [] b;
    }
};

std::string Matrix::ToString(const char *ch)const{
    std::ostringstream str_stream;
    for (int i=0; i<n_row; i++){
        for (int j=0; j<n_column; j++){
            str_stream << a[i*n_column+j] << " " ;
        }
        str_stream <<  "\n";
    }
    str_stream << ch;
    return str_stream.str();
};

void Matrix::Identity(int n){
    delete [] a;
    n_row = n;
    n_column = n;
    a = new float[n*n];
    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            (i==j) ? a[i*n+j] = 1 : a[i*n+j] = 0;
};

void Matrix::Diagonal(int n, float *elements){
    n_row = n;
    n_column = n;
    a = new float[n*n];
    for (int i = 0; i<n; i++)
        a[i*n+i] = elements[i];
};

void Matrix::SwappingOfRows(int i, int j){
    float f;
    if ((i<n_row)&&(j<n_row)&&(i!=j)){
        for (int k = 0; k<n_column; k++){
            f = a[i*n_column+k];
            a[i*n_column+k] = a[j*n_column+k];
            a[j*n_column+k] = f;
        }
    }
};

void Matrix::SwappingOfColumns(int i, int j){
    float f;
    if ((i<n_column)&&(j<n_column)&&(i!=j)){
        for (int k = 0; k<n_row; k++){
            f = a[k*n_column+i];
            a[k*n_column+i] = a[k*n_column+j];
            a[k*n_column+j] = f;
        }
    }
};

Matrix Matrix::LowerTriangular()const{
    int max_column = (n_column>n_row) ? n_row : n_column;
    Matrix result = Matrix(n_row, max_column);
    for (int i = 0; i<n_row; i++){
        for (int j = 0; j<i; j++){
            result.a[i*max_column + j] = a[i*n_column + j];
        }
    }
    return result;
};

Matrix Matrix::UpperTriangular()const{
    int max_row = (n_column>n_row) ? n_row : n_column;
    Matrix result = Matrix(max_row, n_column);
    for (int i = 0; i<max_row; i++){
        for (int j = i; j<n_column; j++){
            result.a[i*n_column + j] = a[i*n_column + j];
        }
    }
    return result;
};

float Matrix::A(int i, int j){
    return a[i*n_column +j];
};

void Matrix::ToLUP(){
    if (!lupped){
        lupped = 1;
        float f;
        delete [] p;
        p = new int[n_row];
        for (int i = 0; i<n_row; i++) p[i] = i;

        int max_size = (n_column>n_row) ? n_row : n_column;
        for (int j = 0; j<max_size; j++){
            f = A(j,j);
            for (int i = j; i<n_row; i++){
                if (fabs(A(i,j))>fabs(f)) {
                    f = A(i,j);
                    p[j] = i;
                }
            }
            if (f != 0.0){
                SwappingOfRows(j, p[j]);
                for (int i = j; i<n_column; i++){
                    for (int k = j+1; k<n_row; k++){
                        if (i == j) a[k*n_column+i] = A(k,i)/f;
                        if (i>j)   a[k*n_column+i] = A(k,i) - A(j,i)*A(k,j);
                    }
                }
            }
        }
    }
}

void Matrix::FromLUP(){
    if (lupped){
        *this = P()*L()*U();
    }
};

Matrix Matrix::L(){
    if (!lupped) ToLUP();
    Matrix m;
    int max_size = (n_column>n_row) ? n_row : n_column;
    m = LowerTriangular();
    for (int i = 0; i<max_size; i++){
        m.a[i*max_size+i] = 1;
    }
    return m;
};

Matrix Matrix::U(){
    if (!lupped) ToLUP();
    return UpperTriangular();
};

Matrix Matrix::P(){
    if (!lupped) ToLUP();
    Matrix m;
    m.Identity(n_row);
    for (int i = 0; i<n_row; i++){
        m.SwappingOfRows(i, p[i]);
    }
    return m;
};

float Matrix::Det(){
    if (!lupped) ToLUP();
    float det = 1.0;
    int pp = 0;
    if (IsQuad()){
        for (int i = 0; i<n_column; i++){
            det *= A(i,i);
            if (p[i]!=i) pp++;
        }
        return (pp % 2 == 0) ? det : -det;
    } else return 0;
};

////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////

Matrix operator*(const float scale, const Matrix &m){
    return Matrix(m.n_row, m.n_column, m.a, scale);
};


void TestMatrix(){
    float a[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    float b[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix m1 = Matrix(3, 3, a); std::cout << m1.ToString() << std::endl;
    Matrix m2 = Matrix(3, 3, b); std::cout << m2.ToString() << std::endl;
    Matrix m3;
    //for (int i = 0; i<10000; i++){
    m3 = m1 + m2; std::cout << m3.ToString() << std::endl;
    m3 = m2 - m1; std::cout << m3.ToString() << std::endl;
    m3 = m1*2;    std::cout << m3.ToString() << std::endl;
    m3 = 3*m1;    std::cout << m3.ToString() << std::endl;
    m3 = m2*m1;    std::cout << m3.ToString() << std::endl;
    m3 = m2/3;    std::cout << m3.ToString() << std::endl;
    m3.Transpose(); std::cout << m3.ToString() << std::endl;
    // for (int i = 0; i<10000; i++){
    m3 = m1 + m2; //std::cout << m3.ToString() << std::endl;
    m3 = m2 - m1; //std::cout << m3.ToString() << std::endl;
    m3 = m1*2;    //std::cout << m3.ToString() << std::endl;
    m3 = 3*m1;    //std::cout << m3.ToString() << std::endl;
    m3 = m1*m2;    std::cout << m3.ToString() << std::endl;
    m3 = m2/3;    //std::cout << m3.ToString() << std::endl;
    m3.Transpose(); //std::cout << m3.ToString() << std::endl;
    m3 = m2;
    m3.SwappingOfColumns(0,2); std::cout << m3.ToString() << std::endl;
    m3 = m2;
    m3.SwappingOfRows(0,2); std::cout << m3.ToString() << std::endl;
    m3 = m2.LowerTriangular(); std::cout << m3.ToString() << std::endl;
    m3 = m2.UpperTriangular(); std::cout << m3.ToString() << std::endl;
    m3.Identity(4);std::cout << m3.ToString() << std::endl;
};

void TestLUP(Matrix &m){
    std::cout << m.ToString() << std::endl;
    m.ToLUP();
    std::cout << m.ToString() << std::endl;
    m.FromLUP();
    std::cout << m.ToString() << std::endl;
};

void TestLUP(){
    float a[9] = {2, 7, -6, 8, 2, 1, 7, 4, 2};
    float b[9] = {0, 8, -6, 0, 2, 1, 0, 4, 2};
    float c[12] = {2, 7, -6, 8, 2, 1, 7, 4, 2, 1, 2, 3};
    Matrix m1 = Matrix(3, 3, a);
    Matrix m2 = Matrix(3, 3, b);
    Matrix m3 = Matrix(3, 4, c);
    Matrix m4 = Matrix(4, 3, c);
    Matrix m5; m5.Identity(3);
    TestLUP(m1);
    TestLUP(m2);
    TestLUP(m3);
    TestLUP(m4);
    TestLUP(m5);
};

