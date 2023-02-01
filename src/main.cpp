/*
Log:

Currently test the Vector class, Matrix class and Conjugate gradient method
Spec Test is as follows:
============================================
Test score: 5/100

Always Fail: 90/100
Heat tests failed to compile: Correct signature for the heat class required.
============================================

TODO:
Finite difference discretization of the heat equation

*/



// header files
// ===================================================================================================
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>
// ===================================================================================================



// Vector template class
// ===================================================================================================
template <typename T>
class Vector
{
public:
    // Default constructor
    Vector() : length(0), data(nullptr) {}

    // Copy constructor
    Vector(const Vector& other) : length(other.length), data(new T[length])
    {
        std::copy(other.data, other.data + length, data);
    }

    // Move constructor
    Vector(Vector&& other) : length(other.length), data(other.data)
    {
        other.length = 0;
        other.data = nullptr;
    }

    // Constructor that takes a length
    Vector(int length) : length(length), data(new T[length]) { initialize(); }

    // Constructor that takes an initializer list
    Vector(std::initializer_list<T> list) : length((int)list.size()), data(new T[length])
    {
        std::copy(list.begin(), list.end(), data);
    }

    // initialize
    void initialize()
    {
        for (int i = 0; i < length; ++i)data[i] = 0;
    }

    // Destructor
    ~Vector()
    {
        length = 0;
        delete[] data;
        data = nullptr;
    }


    // operators =======================================


    // copy assignment operator
    Vector& operator=(const Vector& other)
    {
        if (this != &other)
        {
            delete[] data;
            length = other.length;
            data = new T[length];
            for (int i = 0; i < other.length; ++i)
                data[i] = other.data[i];
        }
        return *this;
    }

    // move assignment operator
    Vector& operator=(Vector&& other)
    {
        if (this != &other)
        {
            delete[] data;
            length = other.length;
            data = other.data;
            other.length = 0;
            other.data = nullptr;
        }
        return *this;
    }

    // operator []
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; } // for const objects

    // operator+
    template<typename U>
    Vector<typename std::common_type<T, U>::type> operator+(const Vector<U>& rhs) const
    {
        if (length != rhs.len())
        {
            throw std::invalid_argument("Different vector lengths");
        }

        using result_type = typename std::common_type<T, U>::type; // necessary
        Vector<result_type> result(length);

        for (int i = 0; i < length; ++i) { result[i] = data[i] + rhs[i]; }

        return result;
    }

    // operator-
    template<typename U>
    Vector<typename std::common_type<T, U>::type> operator-(const Vector<U>& rhs) const
    {
        if (length != rhs.len())
        {
            throw std::invalid_argument("Different vector lengths");
        }

        using result_type = typename std::common_type<T, U>::type; // necessary
        Vector<result_type> result(length);

        for (int i = 0; i < length; ++i) { result[i] = data[i] - rhs[i]; }

        return result;
    }

    // operator* between a vector and a scalar (w = v * s)
    template<typename U>
    Vector<typename std::common_type<T, U>::type> operator*(const U& scalar) const
    {
        using result_type = typename std::common_type<T, U>::type;
        Vector<result_type> result(length);

        for (int i = 0; i < length; ++i)
        {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    // operator* between a scalar and a Vector, a friend function, not a member function
    template <typename U>
    friend Vector<typename std::common_type<T, U>::type> operator*(const U& scalar, const Vector<T>& vector)
    {
        return vector * scalar;
    }


    // functions =======================================


    int len() { return length; }
    const int len() const { return length; } // for const objects

    // get data
    T* Data() { return data; }
    const T* Data() const { return data; }

    // for cout
    friend std::ostream& operator<<(std::ostream& os, const Vector<T>& vec)
    {
        os << "{";
        for (int i = 0; i < vec.len(); ++i)
            os << vec.Data()[i] << ((i + 1) != vec.len() ? "," : "");
        os << "}";
        return os;
    }

private:
    int length;
    T* data;
};
// ===================================================================================================



// dot function - for inner production of Vectors
// ===================================================================================================
template<typename T, typename U>
typename std::common_type<T, U>::type
dot(const Vector<T>& lhs,
    const Vector<U>& rhs)
{
    if (lhs.len() != rhs.len()) // Throw exception if the vectors have different length
        throw std::invalid_argument("Vectors have different lengths!");

    // Calculate the dot product
    using result_type = typename std::common_type<T, U>::type;
    result_type result = 0;

    for (int i = 0; i < lhs.len(); ++i)
        result += lhs[i] * rhs[i];
    return result;
}
// ===================================================================================================



// Matrix
// ===================================================================================================
template <typename T>
class Matrix
{
public:

    // default constructor
    Matrix(): m_nrows(0), m_ncols(0) {}

    // constructor
    Matrix(int rows, int cols) : m_nrows(rows), m_ncols(cols) {}

    // dextructor
    ~Matrix()
    {
        m_nrows = 0;
        m_ncols = 0;
        m_data.clear(); // clear the map
    }

    // operator []
    T& operator[](const std::pair<int, int>& ij) {
        // boundary check
        if (ij.first > m_nrows || ij.second > m_ncols) throw "Matrix [] operator: index out of range";
        return m_data[ij];
    }

    // operator ()
    const T& operator()(const std::pair<int, int>& ij) const {
        // boundary check
        if (ij.first > m_nrows || ij.second > m_ncols) throw "Matrix () operator: index out of range";

        auto it = m_data.find(ij);
        if (it == m_data.end()) { // not found
            throw std::out_of_range("Matrix entry not present");
        }

        // Returns a reference to the mapped value of the element identified with key k.
        // use at() function of std::map or m_data[ij]?
        return m_data.at(ij); 
    }

    // functions

    // iterators
    const auto begin() const { return m_data.begin(); }
    const auto end() const { return m_data.end(); }

    // get the matrix size (number of rows)
    int size() { return m_nrows; }
    const int size() const { return m_nrows; }

    // get the number of rows
    int nrows() { return m_nrows; }
    const int nrows() const { return m_nrows; }

    // get the number of cols
    int ncols() { return m_ncols; }
    const int ncols() const { return m_ncols; }

    // get the data
    std::map<std::pair<int, int>, T> data() { return m_data; }
    const std::map<std::pair<int, int>, T> data() const { return m_data; }

    // for cout
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
    {
        for (int i = 0; i < mat.nrows(); ++i)
        {
            for (int j = 0; j < mat.ncols(); ++j)
            {
                auto x = mat({i, j});   
                if (std::abs(x) < 1.0e-8) x = 0.0;
                os << x << '\t';
            }
            os << '\n';
        }
        return os;
    }

private:
    std::map<std::pair<int, int>, T> m_data; // non-zero matrix entries
    int m_nrows; // number of rows
    int m_ncols; // number of cols
};
// ===================================================================================================



// operator* for Matrix - Vector multiplication
// ===================================================================================================
template<typename T, typename U>
Vector<typename std::common_type<T, U>::type>
operator*(const Matrix<T>& lhs, const Vector<U>& rhs)
{
    // if dimensions are not compatible
    if (lhs.ncols() != rhs.len()) {
        throw std::invalid_argument("Matrix and Vector dimensions not compatible");
    }

    // if compatible, the result of multiplication should be a vector
    using result_type = typename std::common_type<T, U>::type; // necessary
    Vector<result_type> result(lhs.nrows());

    //initialize
    //result.initialize();

    for (auto it = lhs.begin(); it != lhs.end(); ++it) { // need to define begin() and end() function
        //std::cout << it->first.first << " " << it->first.second << " " << it->second << '\n';
        result[it->first.first] = result[it->first.first] + it->second * rhs[it->first.second];
    }

    /*for (auto const& [key, value] : lhs.data()) {
        result[key.first] += value * rhs[key.second];
    }*/

    return result;
}
// ===================================================================================================



// cg function
// ===================================================================================================
template<typename T>
int cg(const Matrix<T>& A,
    const Vector<T>& b,
    Vector<T>& x,
    T                tol = (T)1e-8,
    int              maxiter = 100)
{
    // initialize
    Vector<T> r = b - A * x; // x is usually initialized as vector {0}
    Vector<T> p = r;

    // iterate
    for (int k = 0; k < maxiter; ++k)
    {
        Vector<T> rold = r; // store previous residual
        Vector<T> Ap = A * p; // matrix A times vector P

        T alpha = dot(r, r) / dot(p, Ap);

        x = x + alpha * p; // next estimate of solution
        r = r - alpha * Ap; // residual 

        if (dot(r, r) < tol * tol) // convergence test
            return k + 1;

        T beta = dot(r, r) / dot(rold, rold);
        p = r + beta * p; // next gradient
    }

    return -1;
}
// ===================================================================================================



// Heat class - need to implement
// ===================================================================================================
template <int n, typename T>
class Heat
{
public:
    Heat(T alpha, int m, T dt) : alpha_(alpha), m_(m), dt_(dt)
    {
        // Create the iteration matrix
        int size = (m - 2) * (m - 2);
        Matrix<T> A(size, size);

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (i == j)
                {
                    A[{i, j}] = 1 + 2 * n * alpha * dt / (pow(m - 1, 2));
                }
                else if (i == j - 1 || i == j + 1 || i == j - (m - 2) || i == j + (m - 2))
                {
                    A[{i, j}] = -alpha * dt / (pow(m - 1, 2));
                }
                else
                {
                    A[{i, j}] = 0;
                }
            }
        }
        M_ = A;
        std::cout << "M_: " << '\n';
        std::cout << M_ << '\n';
    }

    Vector<T> exact(T t) const
    {
        Vector<T> result((m_ - 2) * (m_ - 2));
        for (int i = 0; i < m_ - 2; i++)
        {
            for (int j = 0; j < m_ - 2; j++)
            {
                T x = i * pow(m_ - 1, -1);
                T y = j * pow(m_ - 1, -1);
                T sum = 0;
                for (int k = 1; k <= n; k++)
                {
                    for (int l = 1; l <= n; l++)
                    {
                        //sum += pow(-1, k + l) * pow(M_PI, 2) * pow(k, 2) * pow(l, 2) * sin(M_PI * k * x) * sin(M_PI * l * y) * exp(-pow(M_PI, 2) * pow(k, 2) * pow(l, 2) * alpha_ * t);
                    }
                }
                result[i * (m_ - 2) + j] = sum;
            }
        }
        return result;
    }

    Vector<T> solve(T t) const
    {
        int steps = (int)(t / dt_);
        Vector<T> u = exact(0);
        for (int i = 0; i < steps; i++)
        {
            u = M_ * u;
        }
        return u;
    }

private:
    T alpha_;
    int m_;
    T dt_;
    Matrix<T> M_;
};
// ===================================================================================================




// main function
// ===================================================================================================
int main(int argc, char* argv[])
{
    //// test matrix - vector multiplication
    std::cout << "test matrix - vector multiplication" << '\n';
    Matrix<double> m(4, 4);

    
    m[{0,0}] = 1; m[{0,1}] = 0; m[{0,2}] = 0; m[{0,3}] = 1;
    m[{1,0}] = 0; m[{1,1}] = 2; m[{1,2}] = 3; m[{1,3}] = 0;
    m[{2,0}] = 0; m[{2,1}] = 0; m[{2,2}] = 5; m[{2,3}] = 0;
    m[{3,0}] = 1; m[{3,1}] = 0; m[{3,2}] = 0; m[{3,3}] = 1;

    std::cout << m;
    

   /* m[{0, 0}] = 1; m[{0, 3}] = 1;
    m[{1, 1}] = 2; m[{1, 2}] = 3;
    m[{2, 2}] = 5;
    m[{3, 0}] = 1; m[{3, 3}] = 1;*/

    Vector<double> n(4);
    n[0] = 1.22958; n[1] = -0.196732; n[2] = -0.196732; n[3] = 1.22958;
    auto mv = m * n;
    std::cout << "mv: " << mv << '\n';
    //std::cout << "result type of matrix - vector multiplication: " << typeid(mv[0]).name() << '\n';
    //mv.print();
    //std::cout << "== Matrix test end == " << std::endl;
    //std::cout << '\n';
    //std::cout << '\n';
    //---------------- Matrix test end----------------

    for (int i = 0; i < 50; ++i)
    {
        // ---------------- cg function test start----------------
        std::cout << "== cg test start == " << std::endl;

        //Ax = b
        Matrix<double> cg_A(4, 4); // symmetric positive definite matrix A

        // must define the matrix A explicitly - why?

        /*
        cg_A[{0, 0}] = 1; cg_A[{0, 1}] = 0; cg_A[{0, 2}] = 0; cg_A[{0, 3}] = 1;
        cg_A[{1, 0}] = 0; cg_A[{1, 1}] = 2; cg_A[{1, 2}] = 3; cg_A[{1, 3}] = 0;
        cg_A[{2, 0}] = 0; cg_A[{2, 1}] = 0; cg_A[{2, 2}] = 5; cg_A[{2, 3}] = 0;
        cg_A[{3, 0}] = 1; cg_A[{3, 1}] = 0; cg_A[{3, 2}] = 0; cg_A[{3, 3}] = 1;
        */



        cg_A[{0, 0}] = 1; cg_A[{0, 3}] = 1;
        cg_A[{1, 1}] = 2; cg_A[{1, 2}] = 3;
        cg_A[{2, 2}] = 5;
        cg_A[{3, 0}] = 1; cg_A[{3, 3}] = 1;


        Vector<double> cg_x = { 0, 0, 0, 0 }; // initial guess of cg function

        Vector<double> cg_b(4); // known in advance
        cg_b[0] = 2;
        cg_b[1] = 5;
        cg_b[2] = 5;
        cg_b[3] = 2;

        // A, x, b must have the same type
        // x should be: [1, 1, 1, 1] 
        int cg_1 = cg<double>(cg_A, cg_b, cg_x);
        std::cout << "cg test 1: " << cg_1 << '\n';
        std::cout << "x: " << cg_x << '\n';
        std::cout << '\n';


        // another test
        Matrix<double> A_(2, 2);
        A_[{0, 0}] = 4; A_[{0, 1}] = 1;
        A_[{1, 0}] = 1; A_[{1, 1}] = 3;
        Vector<double> b_ = { 1, 2 };
        Vector<double> x_ = { 0, 0 };

        int cg_2 = cg<double>(A_, b_, x_);
        std::cout << "cg test 2: " << cg_2 << '\n';
        std::cout << "x: " << x_ << '\n';

        std::cout << "== cg test end == " << std::endl;
        std::cout << '\n';
        std::cout << '\n';
        // ---------------- cg function test end----------------

    }


    //Heat<1, double> heat(0.3125, 99, 0.001);

    return 0;
}