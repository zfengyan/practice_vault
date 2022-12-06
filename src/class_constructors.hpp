#pragma once

// Include header file for standard input/output stream library
#include <iostream>

// Include header file for initializer list
#include <initializer_list>

// Include header file for memory
#include <memory>

// Class encapsulating array and its lengths
class Container
{
public:
    // Constructor
    explicit Container(int length)
        : data(new double[length]),
        length(length)
    {
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Unified initialization Constructor
    explicit Container(const std::initializer_list<double>& list)
        : Container((int)list.size())
    {
        std::uninitialized_copy(list.begin(), list.end(), data);
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Copy Constructor
    explicit Container(const Container& other)
        : Container(other.length)
    {
        for (auto i = 0; i < other.length; i++)
            data[i] = other.data[i];
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Move Constructor
    Container(Container&& other)
        : data(other.data), length(other.length)
    {
        other.length = 0;
        other.data = nullptr;
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Destructor
    ~Container()
    {
        delete[] data;
        data = nullptr;
        length = 0;
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Copy assignment
    Container& operator=(const Container& other)
    {
        if (this != &other)
        {
            delete[] data;
            data = nullptr;
            length = other.length;
            data = new double[other.length];
            for (auto i = 0; i < other.length; i++)
                data[i] = other.data[i];
        }
        std::cout << __FUNCTION__ << " called" << std::endl;
        return *this;
    }

    // Move assignment
    Container& operator=(Container&& other)
    {
        if (this != &other)
        {
            delete[] data;
            length = other.length;
            data = other.data;
            other.length = 0;
            other.data = nullptr;
        }
        std::cout << __FUNCTION__ << " called" << std::endl;
        return *this;
    }

    // Print container content
    void print()
    {
        std::cout << "Container: ";
        for (auto i = 0; i < length; i++)
            std::cout << data[i] << " ";
        std::cout << std::endl;
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

    // Print container info
    void info()
    {
        std::cout << "Length of data pointer:  " << length << std::endl;
        std::cout << "Address of data pointer: " << &data << std::endl;
        std::cout << "Data pointer:            " << data << std::endl;
        std::cout << __FUNCTION__ << " called" << std::endl;
    }

private:
    // Attributes
    double* data;
    int length;
};