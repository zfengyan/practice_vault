#include <functional>
#include <iostream>

// Your implementation of the abstract base class goes here
class Derivative
{
public: //[DO NOT MODIFY/REMOVE THIS GETTER: IT IS USED IN THE SPECTEST]
    double GetH() const { return h; }
public:
    Derivative():h(1e-8)
    {}

    Derivative(double H):h(H)
    {}

    // pure virtual function
    virtual double differentiate(const std::function<double(double)>& func, double x) const = 0;
protected:
    double h;
};


// Your implementation of the derived class for the central difference scheme goes here
class CentralDifference : public Derivative 
{
public:
    using Derivative::Derivative; // inherit the constructors from the base class
public:
    virtual double differentiate(const std::function<double(double)>& func, double x) const override final
    {
        return ((func(x + h) - func(x - h)) / (2 * h));
    }
};


// Your implementation of the derived class for the forward difference scheme goes here
class ForwardDifference : public Derivative
{
public:
    using Derivative::Derivative; // inherit the constructors from the base class
public:
    virtual double differentiate(const std::function<double(double)>& func, double x) const override final
    {
        return ((func(x + h) - func(x)) / h);
    }
};


double test1(double x) { return x; }
auto test2 = [](double x) { return x; };


int main()
{
    // Your tests go here
    CentralDifference central_diff_default;
    std::cout << central_diff_default.differentiate(test1, 10);
    std::cout << central_diff_default.differentiate(test2, 10);

    ForwardDifference forward_diff_default;
    std::cout << forward_diff_default.differentiate(test1, 10);
    std::cout << forward_diff_default.differentiate(test2, 10);

    CentralDifference central_diff_1(0.1);
    std::cout << central_diff_1.differentiate(test1, 10);
    std::cout << central_diff_1.differentiate(test2, 10);

    ForwardDifference forward_diff_1(0.1);
    std::cout << forward_diff_1.differentiate(test1, 1.0);
    std::cout << forward_diff_1.differentiate(test2, 1.0);

    return 0;
}