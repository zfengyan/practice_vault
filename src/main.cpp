#include <iostream>
#include <typeinfo>
#define PRINT_EXPRESSION(expr) std::cout << #expr << ": " << (expr) \
    << " (type: " << typeid(expr).name() << ")" << std::endl

template <typename T>
T add_simple(const T& a, const T& b)
{
    return (a + b);
}

template <typename T, typename U>
auto add(const T& a, const U& b) -> decltype(a + b)
{
    return (a + b);
}

template <typename T, typename U>
auto multiply(const T& a, const U& b) -> decltype(a + b)
{
    return (a * b);
}

template <typename T, typename U>
auto divide(const T& a, const U& b) -> decltype(a + b)
{
    assert(b != static_cast<U>(0));
    return (a / b);
}

template <typename T>
bool is_int(T t) {
    return false;
}

template <>
bool is_int(int i) {
    return true;
}



template <typename T>
class Number {
public:
    const T value; 

    // constructor that takes a value of type T and initializes the attribute
    Number(const T& value) : value(value) 
    {}

    // operators
    template <typename U>
    auto operator+(const Number<U>& rhs) const{
        return Number<decltype(add(value, rhs.value))>(add(value, rhs.value));
    }

    template <typename U>
    auto operator-(const Number<U>& rhs) const {
        return Number<decltype(add(value, -rhs.value))>(add(value, -rhs.value));
    }

    template <typename U>
    auto operator*(const Number<U>& rhs) const {
        return Number<decltype(multiply(value, rhs.value))>(multiply(value, rhs.value));
    }

    template <typename U>
    auto operator/(const Number<U>& rhs) const {
        return Number<decltype(divide(value, rhs.value))>(divide(value, rhs.value));
    }
};


template <int n>
struct fibonacci {
    static const int value = fibonacci<n - 1>::value + fibonacci<n - 2>::value;
};

// template specialization for the first two Fibonacci numbers
template <>
struct fibonacci<0> {
    static const int value = 0;
};

template <>
struct fibonacci<1> {
    static const int value = 1;
};


int main()
{
    PRINT_EXPRESSION(add_simple(1, 2));
    PRINT_EXPRESSION(add_simple(1.5, 2.2));
    PRINT_EXPRESSION(add(1, 2));
    PRINT_EXPRESSION(add(1.0, 2));
    
    std::cout << (Number<int>(2) * Number<double>(1.2)).value << std::endl;

    // print the first 10 Fibonacci numbers
    std::cout << fibonacci<0>::value << '\n';
    std::cout << fibonacci<1>::value << '\n';
    std::cout << fibonacci<2>::value << '\n';
    std::cout << fibonacci<3>::value << '\n';
    std::cout << fibonacci<4>::value << '\n';
    std::cout << fibonacci<5>::value << '\n';

    std::cout << std::endl;
    return 0;
}