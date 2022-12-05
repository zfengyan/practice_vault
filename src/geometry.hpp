#pragma once

#include <iostream>

template<typename T>
class Point
{
public:

    // constructors
    
    Point():_x(T(0)), _y(T(0))
    {} // default constructor

    Point(T px, T py):_x(px), _y(py)
    {} // constructor with two parameters

    template <typename U> // if the constructor gets different @param types
    Point(U px, U py) :
        _x(static_cast<T>(px)), _y(static_cast<T>(py)) {}


    // operators

    // operator+ : yields a new point
    // a+b => a.operator+(b)
    Point<T> operator+(const Point<T>& rhs)const 
    {        
        return Point(this->x() + rhs.x(), this->y() + rhs.y());
    }

    template <typename U>
    Point<T> operator+(Point<U>& rhs)const
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        return Point(this->x() + rhs_x, this->y() + rhs_y);
    }

    

    // methods
    
    // get the x coordinate
    T& x() { return _x; } // return reference, so that we can both access and modify the x coordinate
    const T& x()const { return _x; } // const version, for const objects

    // get the y coordinate
    T& y() { return _y; } // return reference, so that we can both access and modify the x coordinate
    const T& y()const { return _y; } // const version, for const objects

    // ostream overloading, for printing the info of the Point<T> objects
    friend std::ostream& operator<<(std::ostream& os, const Point& p)
    {
        os << "(" << p.x() << ", " << p.y() << ")";
        return os;
    }
protected:
    T _x;
    T _y;
};
