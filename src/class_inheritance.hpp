#pragma once

// Include header file for standard input/output stream library
#include <iostream>

// Change #define to #undef to switch between the two variants
#define USE_STD_FUNCTION

// Include header file to enable use of std::function
#ifdef USE_STD_FUNCTION
#include <functional>
#endif

// Abstract class Quadrature which cannot be instanciated
class Quadrature
{
public:
    Quadrature()
        : weights(nullptr),
        points(nullptr),
        n(0)
    { }

    // Constructor
    Quadrature(int n)
        : weights(new double[n]),
        points(new double[n]),
        n(n)
    { }

    // Destructor
    ~Quadrature()
    {
        delete[] weights;
        delete[] points;
        n = 0;
    }

    // Pure virtual factor (this function must be implemented in any
    // class derived from the abstract class Quadrature)
    virtual double factor(double a, double b) const = 0;

    // Pure virtual mapping (this function must be implemented in any
    // class derived from the abstract class Quadrature)
    virtual double mapping(double xi, double a, double b) const = 0;

#ifdef USE_STD_FUNCTION

    // Virtual integrator using lambda expressions (this function can
    // be overridden by a re-implementation in the class derived from
    // the abstract class Quadrature; if the derived class does not
    // provide a re-implementation then the implementation from the
    // calss Quadrature is used instead)
    virtual double integrate(const std::function<double(double)>& func,
        double a, double b) const
    {
        double integral(0);
        for (auto i = 0; i < n; i++)
            integral += weights[i] * func(mapping(points[i], a, b)); // declaration
        return factor(a, b) * integral;
    }

#else

    // Virtual integrator using function pointers (this function can
    // be overridden by a re-implementation in the class derived from
    // the abstract class Quadrature; if the derived class does not
    // provide a re-implementation then the implementation from the
    // calss Quadrature is used instead)
    virtual double integrate(double (*func)(const double x),
        double a, double b) const
    {
        double integral(0);
        for (auto i = 0; i < n; i++)
            integral += weights[i] * func(mapping(points[i], a, b));
        return factor(a, b) * integral;
    }

#endif

protected:
    // Attributes
    double* weights;
    double* points;
    int n;
};

// Midpoint rule: override the virtual integate functions since this
// method is so simple that it can be easily implemented without
// mapping the reference interval [-1,1] to the physical one [a,b]
class MidpointRule : public Quadrature
{
public:
    // Factor (not used!)
    virtual double factor(double a, double b) const override final
    {
        return 1;
    }

    // Mapping (not used!)
    virtual double mapping(double xi, double a, double b) const override final
    {
        return 0;
    }

#ifdef USE_STD_FUNCTION

    // Integrator: lambda expression
    virtual double integrate(const std::function<double(double)>& func,
        double a, double b) const override final
    {
        double m = 0.5 * (a + b);
        return (b - a) * func(m);
    }

#else

    // Integrator: function pointer
    virtual double integrate(double (*func)(const double x),
        double a, double b) const override final
    {
        double m = 0.5 * (a + b);
        return (b - a) * func(m);
    }

#endif
};

// Simpson rule: override the virtual integate functions since this
// method is so simple that it can be easily implemented without
// mapping the reference interval [-1,1] to the physical one [a,b]
class SimpsonRule : public Quadrature
{
public:
    // Factor (not used!)
    virtual double factor(double a, double b) const override final
    {
        return 1;
    }

    // Mapping (not used!)
    virtual double mapping(double xi, double a, double b) const override final
    {
        return 0;
    }

#ifdef USE_STD_FUNCTION

    // Integrator: lambda expression
    double integrate(const std::function<double(double)>& func,
        double a, double b) const override final
    {
        double m = 0.5 * (a + b);
        return (b - a) / 6.0 * (func(a) + 4.0 * func(m) + func(b));
    }


#else

    // Integrator: function pointer
    double integrate(double (*func)(const double x),
        double a, double b) const override final
    {
        double m = 0.5 * (a + b);
        return (b - a) / 6.0 * (func(a) + 4.0 * func(m) + func(b));
    }

#endif
};

// Rectangle rule: override the virtual integate functions since this
// method is so simple that it can be easily implemented without
// mapping the reference interval [-1,1] to the physical one [a,b]
class RectangleRule : public Quadrature
{
public:
    // Constructor: we do not need the weights and points since the
    // rectangle rule is implemented explicitly in the integrate
    // functions. However, we do need the number of sub-rectangles.
    RectangleRule(int n)
        : Quadrature()
    {
        this->n = n;
    }

    // Factor (not used!)
    virtual double factor(double a, double b) const override final
    {
        return 1;
    }

    // Mapping (not used!)
    virtual double mapping(double xi, double a, double b) const override final
    {
        return 0;
    }

#ifdef USE_STD_FUNCTION

    // Integrator: lambda expression
    double integrate(const std::function<double(double)>& func,
        double a, double b) const override final
    {
        double integral(0);
        double h = (b - a) / (double)n;
        for (auto i = 0; i < n; i++)
            integral += h * func(a + i * h);
        return integral;
    }

#else

    // Integrator: function pointer
    double integrate(double (*func)(const double x),
        double a, double b) const override final
    {
        double integral(0);
        double h = (b - a) / (double)n;
        for (auto i = 0; i < n; i++)
            integral += h * func(a + i * h);
        return integral;
    }

#endif
};

// Gauss rule
class GaussRule : public Quadrature
{
public:
    GaussRule(int n)
        : Quadrature(n)
    {
        switch (n)
        {
        case 1:
            weights[0] = { 2.0 };
            points[0] = { 0.0 };
            break;
        case 2:
            weights[0] = { 1.0 };
            weights[1] = { 1.0 };
            points[0] = { -0.57735026919 };
            points[1] = { 0.57735026919 };
            break;
        case 3:
            weights[0] = { 5.0 / 9.0 };
            weights[1] = { 8.0 / 9.0 };
            weights[2] = { 5.0 / 9.0 };
            points[0] = { -0.774596669241 };
            points[1] = { 0.0 };
            points[2] = { 0.774596669241 };
            break;
        default:
            std::cout << "Invalid parameter" << std::endl;
            exit(1);
        }
    }

    // Factor
    virtual double factor(double a, double b) const override final
    {
        return (b - a) / 2.0;
    }

    // Mapping from reference interval [-1,1] to [a,b]
    virtual double mapping(double xi, double a, double b) const override final
    {
        return 0.5 * (b - a) * xi + 0.5 * (a + b);
    }
};

// Define function to integrate as classical function
double myfunc1(const double x) { return x; }

// Define function to integrate as lambda expression
auto myfunc2 = [](double x) { return x; };

// The global main function that is the designated start of the program