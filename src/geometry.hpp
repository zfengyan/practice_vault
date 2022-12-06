#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <initializer_list> // for delegate constructor

class Point
{
public:

    // constructors --------------------------------------------------------------------------------
    
    Point():x(0.0), y(0.0)
    {} // default constructor

    Point(double px, double py):x(px), y(py)
    {} // constructor with two parameters


    // operators --------------------------------------------------------------------------------

    // operator+ : yields a new point
    // a+b => a.operator+(b), const qualifier: allow const object to access this method
    Point operator+(const Point& other)const 
    {        
        return Point(x + other.x, y + other.y);
    }

    // operator- : yields a new point
    // a-b => a.operator-(b), const qualifier: allow const object to access this method
    Point operator-(const Point& other)const
    {
        return Point(x - other.x, y - other.y);
    }

    // operator += : yields a new point
    Point& operator+=(const Point& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    // operator -= : yields a new point
    Point& operator-=(const Point& other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // methods --------------------------------------------------------------------------------
    
    // get the x coordinate
    //double& X() { return x; } // return reference, so that we can both access and modify the x coordinate
    //const double& X()const { return x; } // const version, for const objects

    //// get the y coordinate
    //double& y() { return _y; } // return reference, so that we can both access and modify the x coordinate
    //const double& y()const { return _y; } // const version, for const objects

    // get the distance between the current point and the origin((0, 0) by default)
    double distance()const
    {
        message("calculate the distance from the origin(0, 0)");
        Point p0; // origin - (0, 0) by default
        return distance(p0);
    }

    // get the distance between the current point and the point other
    double distance(const Point& other)const
    {
        return std::sqrt((other.x - x) * (other.x - x) + (other.y - y) * (other.y - y));
    }

    // compute and return a new Point rotated with angle radians in anticlockwise direction around the origin
    // angle: radians
    Point rotated(double angle)const
    {
        return rotated(angle, Point());
    }

    // compute and return a new Point rotated with angle radians in anticlockwise direction around the certain point
    // angle: radians
    Point rotated(double angle, const Point& other)const
    {
        double new_x = (x-other.x)*cos(angle) - (y-other.y)*sin(angle) + other.x;
        double new_y = (x-other.x)*sin(angle) + (y-other.y)*cos(angle) + other.y;
        return Point(new_x, new_y);
    }

    // similar with rotated, return the current instance
    Point& rotate(double angle)
    {
        return rotate(angle, Point());
    }
    Point& rotate(double angle, const Point& other) 
    {
        double new_x = (x - other.x) * cos(angle) - (y - other.y) * sin(angle) + other.x;
        double new_y = (x - other.x) * sin(angle) + (y - other.y) * cos(angle) + other.y;
        x = new_x;
        y = new_y;
        return *this;
    }

    // print a point, const qualifier: allow const object to access this method
    void print()const
    {
        std::cout << "(" << this->x << ", " << this->y << ")" << '\n';
    }

    // print a message
    void message(const char* msg)const
    {
        std::cout << "-- " << msg << '\n';
    }

    // ostream overloading, for printing the info of the Point<T> objects
    friend std::ostream& operator<<(std::ostream& os, const Point& p)
    {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }


    double x;
    double y;
};


class Triangle
{
public:
    // constructors
    Triangle(): a(Point()), b(Point()), c(Point())
    {}

    Triangle(const Point& pa, const Point& pb, const Point& pc):
        a(pa), b(pb), c(pc)
    {}

    // methods
    
    // translated - return a new triangle
    Triangle translated(const Point& t)const
    { 
        double ax = a.x + t.x;
        double ay = a.y + t.y;
        double bx = b.x + t.x;
        double by = b.y + t.y;
        double cx = c.x + t.x;
        double cy = c.y + t.y;
        return Triangle(Point(ax, ay), Point(bx, by), Point(cx, cy));
    }

    // translate - return the modified triangle
    // return by reference - return the current instance
    Triangle& translate(const Point& t)
    {
        a.x += t.x;
        a.y += t.y;
        b.x += t.x;
        b.y += t.y;
        c.x += t.x;
        c.y += t.y;
        return *this;
    }

    // rotated - return a new triangle
    Triangle rotated(double angle)const
    {
        return Triangle(a.rotated(angle), b.rotated(angle), c.rotated(angle));
    }
    Triangle rotated(double angle, const Point& other)const
    {
        return Triangle(a.rotated(angle, other), b.rotated(angle, other), c.rotated(angle, other));
    }

    // rotate - return the current instance
    Triangle& rotate(double angle)
    {
        a.rotate(angle);
        b.rotate(angle);
        c.rotate(angle);
        return *this;
    }
    Triangle& rotate(double angle, const Point& other)
    {
        a.rotate(angle, other);
        b.rotate(angle, other);
        c.rotate(angle, other);
        return *this;
    }

    // area
    double area()const
    {
        return std::abs(0.5 * ((a.x*b.y+b.x*c.y+c.x*a.y) - (a.y*b.x+b.y*c.x+c.y*a.x)));
    }

    // print the triangle
    void print()const
    {
        a.print();
        b.print();
        c.print();
    }

public:
    Point a;
    Point b;
    Point c;
};
