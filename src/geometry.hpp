#pragma once

#include <iostream>
#include <cmath>

template<typename T>
class Point
{
public:

    // constructors --------------------------------------------------------------------------------
    
    Point():_x(T(0)), _y(T(0))
    {} // default constructor

    Point(T px, T py):_x(px), _y(py)
    {} // constructor with two parameters

    template <typename U> // if the constructor gets different @param types
    Point(U px, U py) :
        _x(static_cast<T>(px)), _y(static_cast<T>(py)) {}


    // operators --------------------------------------------------------------------------------

    // operator+ : yields a new point
    // a+b => a.operator+(b), const qualifier: allow const object to access this method
    Point<T> operator+(const Point<T>& rhs)const 
    {        
        return Point(this->x() + rhs.x(), this->y() + rhs.y());
    }

    // operator+ : yields a new point
    // two different types of points, const qualifier: allow const object to access this method
    template <typename U>
    Point<T> operator+(Point<U>& rhs)const
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        return Point(this->x() + rhs_x, this->y() + rhs_y);
    }

    // operator- : yields a new point
    // a-b => a.operator-(b), const qualifier: allow const object to access this method
    Point<T> operator-(const Point<T>& rhs)const
    {
        return Point(this->x() - rhs.x(), this->y() - rhs.y());
    }

    // operator+ : yields a new point
    // two different types of points, const qualifier: allow const object to access this method
    template <typename U>
    Point<T> operator-(Point<U>& rhs)const
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        return Point(this->x() - rhs_x, this->y() - rhs_y);
    }

    // operator += : yields a new point
    Point<T>& operator+=(const Point<T>& rhs)
    {
        this->x() += rhs.x();
        this->y() += rhs.y();
        return *this;
    }

    // operator += : yields a new point
    // two different types of points
    template <typename U>
    Point<T>& operator+=(const Point<U>& rhs)
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        this->x() += rhs_x;
        this->y() += rhs_y;
        return *this;
    }

    // operator -= : yields a new point
    Point<T>& operator-=(const Point<T>& rhs)
    {
        this->x() -= rhs.x();
        this->y() -= rhs.y();
        return *this;
    }

    // operator -= : yields a new point
    // two different types of points
    template <typename U>
    Point<T>& operator-=(const Point<U>& rhs)
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        this->x() -= rhs_x;
        this->y() -= rhs_y;
        return *this;
    }
    

    // methods --------------------------------------------------------------------------------
    
    // get the x coordinate
    T& x() { return _x; } // return reference, so that we can both access and modify the x coordinate
    const T& x()const { return _x; } // const version, for const objects

    // get the y coordinate
    T& y() { return _y; } // return reference, so that we can both access and modify the x coordinate
    const T& y()const { return _y; } // const version, for const objects

    // get the distance between the current point and the origin((0, 0) by default)
    double distance()const
    {
        message("calculate the distance from the origin(0, 0)");
        Point<T> p0; // origin - (0, 0) by default
        return distance(p0);
    }

    // get the distance between the current point and the point rhs
    double distance(const Point<T>& rhs)const
    {
        return std::sqrt((rhs.x() - this->x()) * (rhs.x() - this->x()) + (rhs.y() - this->y()) * (rhs.y() - this->y()));
    }

    // get the distance if two points are not the same type
    template<typename U>
    double distance(const Point<U>& rhs)const
    {
        T rhs_x = static_cast<T>(rhs.x());
        T rhs_y = static_cast<T>(rhs.y()); // conversion
        Point<T> p_conversion(rhs_x, rhs_y);
        return distance(p_conversion);
    }

    // print a point, const qualifier: allow const object to access this method
    void print()const
    {
        std::cout << "(" << this->x() << ", " << this->y() << ")" << '\n';
    }

    // print a message
    void message(const char* msg)const
    {
        std::cout << "-- " << msg << '\n';
    }

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
