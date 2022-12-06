#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <initializer_list> // for delegate constructor

class Point
{
public:

    // constructors --------------------------------------------------------------------------------
    
    Point():_x(0.0), _y(0.0)
    {} // default constructor

    Point(double px, double py):_x(px), _y(py)
    {} // constructor with two parameters

    template <typename U> // if the constructor gets different @param types
    Point(U px, U py) :
        _x(static_cast<double>(px)), _y(static_cast<double>(py)) 
    {
        message("data type will be converted to double");
    }


    // operators --------------------------------------------------------------------------------

    // operator+ : yields a new point
    // a+b => a.operator+(b), const qualifier: allow const object to access this method
    Point operator+(const Point& rhs)const 
    {        
        return Point(this->x() + rhs.x(), this->y() + rhs.y());
    }

    // operator- : yields a new point
    // a-b => a.operator-(b), const qualifier: allow const object to access this method
    Point operator-(const Point& rhs)const
    {
        return Point(this->x() - rhs.x(), this->y() - rhs.y());
    }

    // operator += : yields a new point
    Point& operator+=(const Point& rhs)
    {
        this->x() += rhs.x();
        this->y() += rhs.y();
        return *this;
    }

    // operator -= : yields a new point
    Point& operator-=(const Point& rhs)
    {
        this->x() -= rhs.x();
        this->y() -= rhs.y();
        return *this;
    }

    // methods --------------------------------------------------------------------------------
    
    // get the x coordinate
    double& x() { return _x; } // return reference, so that we can both access and modify the x coordinate
    const double& x()const { return _x; } // const version, for const objects

    // get the y coordinate
    double& y() { return _y; } // return reference, so that we can both access and modify the x coordinate
    const double& y()const { return _y; } // const version, for const objects

    // get the distance between the current point and the origin((0, 0) by default)
    double distance()const
    {
        message("calculate the distance from the origin(0, 0)");
        Point p0; // origin - (0, 0) by default
        return distance(p0);
    }

    // get the distance between the current point and the point rhs
    double distance(const Point& rhs)const
    {
        return std::sqrt((rhs.x() - this->x()) * (rhs.x() - this->x()) + (rhs.y() - this->y()) * (rhs.y() - this->y()));
    }

    // compute and return a new Point rotated with angle radians in anticlockwise direction around the origin
    // angle: radians
    Point rotated(double angle)
    {
        return rotated(angle, Point());
    }

    // compute and return a new Point rotated with angle radians in anticlockwise direction around the certain point
    // angle: radians
    Point rotated(double angle, const Point& rhs)
    {
        double new_x = (this->x()-rhs.x())*cos(angle) - (this->y()-rhs.y())*sin(angle) + rhs.x();
        double new_y = (this->x()-rhs.x())*sin(angle) + (this->y()-rhs.y())*cos(angle) + rhs.y();
        return Point(new_x, new_y);
    }

    // similar with rotated, return the current instance
    Point& rotate(double angle)
    {
        return rotate(angle, Point());
    }
    Point& rotate(double angle, const Point& rhs) 
    {
        double new_x = (this->x() - rhs.x()) * cos(angle) - (this->y() - rhs.y()) * sin(angle) + rhs.x();
        double new_y = (this->x() - rhs.x()) * sin(angle) + (this->y() - rhs.y()) * cos(angle) + rhs.y();
        this->x() = new_x;
        this->y() = new_y;
        return *this;
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
    double _x;
    double _y;
};
