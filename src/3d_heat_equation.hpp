#pragma once

//https://github.com/DJMvanBemmelen/3D-Heat-Equation/blob/5d9a7a2c4930f0044c0caf91b330063561c05202/heat_equation.cpp

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

namespace heat_equation_3d
{


#define Print(x)(std::cout<< x << std::endl)


    /*
    Dingen die mogelijk verbeterd kunnen worden:
    1. Recursion
    Bij vector + vector + vector ....  (Geen idee of dit aan de orde is)

    2. Imaginaire getallen mogelijk maken
    Alle vermenigvuldigings dingen veranderen (Geen idee of dit aan de orde is)

    3.
    */

    template <typename T>
    class Vector
    {
        // Your implementation of the Vector class starts here
        // Attributes
        int length = 0;
        T* array;

    public:

        // Constructors and Destructors ////////////////////////
        // Default constructor
        Vector()
            : length(0), array(nullptr)
        { }

        // Copy constructor for up and down casting 
        template<typename U>
        Vector(const Vector<U>& other)
            : Vector(other.len())
        {
            // Coppying other data to new array
            array = new T[other.len()];
            for (auto i = 0; i < other.len(); i++)
                array[i] = other[i];
        }

        // Copy constructor for equal type
        Vector(const Vector& other)
            : Vector(other.len())
        {
            // Coppying other data to new array
            array = new T[other.len()];
            for (auto i = 0; i < other.len(); i++)
                array[i] = other[i];
        }

        // Move constructor
        Vector(Vector<T>&& other)
            : length(other.length), array(other.array)
        {
            // Deleting other array
            other.length = 0;
            other.array = nullptr;
        }

        // Constructor (does not allow uniform initialisation)
        Vector(int length)
            : length(length), array(new T[length]()) // Added the {0} such that the vector is deleted in the for loop new T[length]{0}
        { }

        // Constructor (with initialiser list and delegate constructor)
        Vector(std::initializer_list<T> list)
            : Vector((int)list.size())
        {
            std::uninitialized_copy(list.begin(), list.end(), array);
        }

        // Destructor
        ~Vector()
        {
            length = 0;
            delete[] array;
            array = nullptr;
        }

        // Operators and Functions
        // Function Len that returns the length 
        int len() const
        {
            return length;
        }
        // Operator[](int i) and its const overload
        auto& operator[](int i)
        {
            return array[i];
        }
        // @@@@ The const vector can still be changed ????@@@@@@@@
        auto& operator[](int i) const
        {
            return array[i];
        }

        // Copy assignment operator for up and down casting
        template<typename U>
        Vector& operator=(const Vector<U>& other)
        {
            if (this != (Vector*)&other) {
                delete[] array;
                length = other.len();
                array = new T[other.len()];
                for (auto i = 0; i < other.len(); i++)
                    array[i] = other[i];
            }
            return *this;
        }

        // Copy assignment operator for equal type
        Vector& operator=(const Vector& other)
        {
            if (this != (Vector*)&other) {
                delete[] array;
                length = other.len();
                array = new T[other.len()];
                for (auto i = 0; i < other.len(); i++)
                    array[i] = other[i];
            }
            return *this;
        }

        // Move assignment opertator
        Vector& operator=(Vector<T>&& other)
        {
            if (this != (Vector*)&other) {
                delete[] array;
                length = other.len();
                array = other.array;
                other.length = 0; /// CHANGED
                other.array = nullptr;
            }
            return *this;
        }

        // Arithmetic operator+
        template<typename U>
        Vector<typename std::common_type<T, U>::type> operator+(const Vector<U>& other) const
        {
            // Throw exception if the vectors have different length
            if (length != other.len()) throw "Vectors have different size!";

            Vector<typename std::common_type<T, U>::type> Vsum(length);
            for (auto i = 0; i < length; i++)
                Vsum[i] = array[i] + other[i];
            return Vsum;

        }

        // Arithmetic operator- 
        template<typename U>
        Vector<typename std::common_type<T, U>::type> operator-(const Vector<U>& other) const
        {
            // Throw exception if the vectors have different length
            if (length != other.len()) throw "Vectors have different size!";

            Vector<typename std::common_type<T, U>::type> Vdif(length);
            for (auto i = 0; i < length; i++)
                Vdif[i] = array[i] - other[i];
            return Vdif;

        }

        // Numpy like display of vector
        void display()
        {
            int mult_10 = 1;
            std::cout << "[";
            for (auto i = 0; i < length - 1; i++) {
                if (i == mult_10 * 10) {
                    std::cout << array[i] << ", " << std::endl;
                    mult_10 += 1;
                }
                else {
                    std::cout << array[i] << ", ";
                }
            }
            std::cout << array[length - 1] << "]" << std::endl;
        }


    };

    //// ARTIHMATIC OPERATORS FOR VECTOR * SCALAR ///////////////////////////////////////////////////////////////////////

    // Arithmatic opperator* for scalar * vector 
    template<typename A, typename B,
        typename = typename std::enable_if<std::is_fundamental<A>::value>::type> // Cheching if A is a scalar and not a Vector/Matrix
        Vector<typename std::common_type<A, B>::type> operator*(const A lhs, const Vector<B> rhs)
    {
        Vector<typename std::common_type<A, B>::type> Vprod(rhs.len());
        for (auto i = 0; i < rhs.len(); i++)
            Vprod[i] = lhs * rhs[i];
        return Vprod;
    }

    // Arithmatic opperator* for vector * scalar
    template<typename A, typename B,
        typename = typename std::enable_if<std::is_fundamental<B>::value>::type> // Cheching if B is a scalar and not a Vector/Matrix
        Vector<typename std::common_type<A, B>::type> operator*(const Vector<A> lhs, const B rhs)
    {
        Vector<typename std::common_type<A, B>::type> Vprod(lhs.len());
        for (auto i = 0; i < lhs.len(); i++)
            Vprod[i] = lhs[i] * rhs;
        return Vprod;
    }

    ///// DOT PRODUCT /////////////////////////////////////////////////////////////////////////////////////

    // Dot product function 
    template<typename T, typename U>
    typename std::common_type<T, U>::type
        dot(const Vector<T>& lhs,
            const Vector<U>& rhs)
    {
        if ((int)lhs.len() != (int)rhs.len()) throw "Vectors have different size!";

        typename std::common_type<T, U>::type result = 0;
        for (auto i = 0; i < (int)lhs.len(); i++)
            result += lhs[i] * rhs[i];

        return result;
    }


    //// MATRIX CLASS ///////////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    class Matrix
    {
        // Attributes 
        int n_rows;
        int n_colums;
        std::map< std::pair<int, int >, T > M;

    public:

        // Empty constructor
        Matrix()
            : n_rows(0), n_colums(0)
        { }

        // Constructor
        Matrix(int n_rows, int n_colums)
            : n_rows(n_rows), n_colums(n_colums)
        { }

        // Destructor 
        ~Matrix()
        {
            n_rows = 0;
            n_colums = 0;
            M.clear(); // Erasing the full map
        }
        // Function that returns the number of row elements
        const int nrows() const
        {
            return n_rows;
        }
        // Function that returns the number of row elements
        const int ncolumns() const
        {
            return n_colums;
        }

        //  Function that returns the begin iterator of the map @@@@
        const auto begin() const
        {
            return M.begin();
        }

        // Function that returns the end iterator of the map @@@@
        const auto end() const
        {
            return M.end();
        }

        auto& operator[](const std::pair<int, int>& ij)
        {
            // Make sure that i < nrows en j < ncolumns 
            if (ij.first > n_rows || ij.second > n_colums) throw "Pair outside the range of the Matrix";
            return M[ij];
        }

        const auto& operator()(const std::pair<int, int>& ij) const
        {
            // Make sure that i < nrows en j < ncolumns 
            if (ij.first > n_rows || ij.second > n_colums) throw "Pair outside the range of the Matrix";
            // Throw exempt if the entry is not presant
            if (M.find(ij) == M.end()) throw "No value assinged to key";
            return M.at(ij);
        }

        friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
        {
            for (int i = 0; i < mat.nrows(); ++i)
            {
                for (int j = 0; j < mat.ncolumns(); ++j)
                {
                    try
                    {
                        auto x = mat({ i, j });
                        os << x << '\t';
                    }
                    catch (...)
                    {
                        os << 0.0 << '\t';
                    }
                }
                os << '\n';
            }
            return os;
        }

    };

    template<typename T, typename U>
    Vector<typename std::common_type<T, U>::type>
        operator*(const Matrix<T>& lhs,
            const Vector<U>& rhs)
    {
        // Make sure that n_colmns of Matrix = Vector length 
        if (lhs.ncolumns() != rhs.len()) throw "Matrix and Vector do not have the right dimensions";
        // Making a vector of length Matrix n_rows
        Vector<typename std::common_type<T, U>::type> Vprod(lhs.nrows());

        for (auto it = lhs.begin(); it != lhs.end(); ++it) {
            Vprod[it->first.first] = Vprod[it->first.first] + it->second * rhs[it->first.second];
        }
        return Vprod;
    }

    //BEGIN JORT

    //Finite difference discretization of the heat equation

    template <int n, typename T>
    class Heat
    {
        //Attributes
    public:
        T alpha;
        int m;
        T dt;
        Matrix<T> M;

        //Constructor
        //alpha diffusion coefficient, m number of points per dimension, dt time-step size
        Heat(T alpha, int m, T dt)
            : alpha(alpha), m(m), dt(dt)
        {
            //M is the matrix from the introduction

            M = { int(std::pow(m,n)), int(std::pow(m,n)) };

            auto fact = alpha * dt * (m + 1) * (m + 1);

            for (auto i = 0; i < std::pow(m, n); i++)
            {
                for (auto j = 0; j < std::pow(m, n); j++)
                {
                    if (i == j)
                    {
                        M[{i, j}] = 1;
                    }

                    for (auto k = 0; k < n; k++)
                    {
                        if (i == j)
                        {
                            M[{i, j}] = M[{i, j}] + fact * 2;
                        }

                        if ((j - i == std::pow(m, k)) && ((j % int(std::pow(m, k + 1))) / (int(std::pow(m, k))) != 0))
                        {
                            M[{i, j}] = M[{i, j}] - fact;
                        }

                        if ((i - j == std::pow(m, k)) && ((i % int(std::pow(m, k + 1))) / (int(std::pow(m, k))) != 0))
                        {
                            M[{i, j}] = M[{i, j}] - fact;
                        }
                    }
                    //double temp = M[{i,j}];
                    //Print("i: " << i+1 << ", j: " << j+1 << ", M: " << temp);
                }
            }
        }

        //Returns the exact solutions for the given boundary conditions at time t.
        //u_exact(x,t) = exp(-n*pi^2*alpha*t)*u(x,0)
        Vector<T> exact(T t) const
        {
            Vector<T> sol((int)std::pow(m, n));

            for (auto l = 0; l < std::pow(m, n); l++)
            {
                auto x = (l % m + 1) / (m + 1.0);                   //1.0 to force double division
                auto y = (l % (m * m) + m - (l % m)) / (m * (m + 1.0));
                auto z = (l % (m * m * m) + m * m - (l % (m * m))) / (m * m * (m + 1.0));

                for (auto k = 0; k < n; k++)
                {
                    if (k == 0)
                    {
                        sol[l] = std::exp(-n * M_PI * M_PI * alpha * t) * std::sin(M_PI * x);
                    }
                    if (k == 1)
                    {
                        sol[l] = sol[l] * std::sin(M_PI * y);
                    }
                    if (k == 2)
                    {
                        sol[l] = sol[l] * std::sin(M_PI * z);
                    }
                }
            }
            return sol;
        }

        //Solves the initial boundary value problem given in the introduction until time t,
        //returns the numerical solution at time t.
        Vector<T> solve(T t) const
        {
            Vector<T> uk(std::pow(m, n));

            for (auto l = 0; l < std::pow(m, n); l++)
            {
                // making the initial conditions

                auto x = (l % m + 1) / (m + 1.0);                   //1.0 to force double division
                auto y = (l % (m * m) + m - (l % m)) / (m * (m + 1.0));
                auto z = (l % (m * m * m) + m * m - (l % (m * m))) / (m * m * (m + 1.0));

                for (auto k = 0; k < n; k++)
                {
                    if (k == 0)
                    {
                        uk[l] = std::sin(M_PI * x);
                    }
                    if (k == 1)
                    {
                        uk[l] = uk[l] * std::sin(M_PI * y);
                    }
                    if (k == 2)
                    {
                        uk[l] = uk[l] * std::sin(M_PI * z);
                    }
                }
            }

            Vector<T> uk1 = uk;  // initial guess same as b

            int no = t / dt; //number of time-steps

            // Print(no);

            for (auto p = 0; p < no; p++)
            {
                cg(M, uk, uk1);

                uk = uk1;       // next iteration with solution of previous iteration.
            }

            return uk;
        }
    };


    //Conjugate gradient method
    template<typename T>
    int cg(const Matrix<T>& A,
        const Vector<T>& b,
        Vector<T>& x,
        T          tol = (T)1e-8,
        int        maxiter = 100)
    {
        //translated block of pseudocode.
        Vector<T> r = b - A * x;
        Vector<T> p = r;
        int iter = 0;        // variable to keep track of iterations

        // if (dot(r,r) > tol*tol)
        // {
        //     for (auto i = 0; i < maxiter; i++)
        //     {
        //         iter++;
        //         auto Ap = A*p;
        //         auto alphak = dot(r, r) / dot(Ap, p);
        //         // Print(alphak);
        //         x          = x + alphak*p;
        //         auto dotrk = dot(r, r);     // dummy variable to calculate beta
        //         // Print(dotrk);
        //         r          = r - alphak*(Ap);

        //         if (dot(r, r) < tol*tol)
        //         {
        //             break;
        //         }

        //         auto beta  = dot(r, r) / dotrk;
        //         p          = r + beta*p;

        //         // Print("alphak: " << alphak);
        //         // Print(dotrk);
        //         // Print(beta);
        //     }
        // }



        for (auto i = 0; i < maxiter; i++)
        {
            iter++;
            auto Ap = A * p;
            auto alphak = dot(r, r) / dot(Ap, p);
            // Print(alphak);
            x = x + alphak * p;
            auto dotrk = dot(r, r);     // dummy variable to calculate beta
            // Print(dotrk);
            r = r - alphak * (Ap);

            if (dot(r, r) < tol * tol)
            {
                break;
            }

            auto beta = dot(r, r) / dotrk;
            p = r + beta * p;

            // Print("alphak: " << alphak);
            // Print(dotrk);
            // Print(beta);
        }


        if (iter == maxiter)
        {
            iter = 0;
        }
        // Print(iter);
        return iter - 1;
    }
}