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
public: //[DO NOT MODIFY/REMOVE THESE GETTERS AND SETTERS: THEY ARE USED IN THE SPECTEST]
    int GetLength() const { return length; }
    double* GetData() const { return data; }
    void SetLength(const int length) { this->length = length; }
    void SetData(double* data) { this->data = data; }

public:

    // default constructor
    explicit Container():
        length(0), data(nullptr)
    {
        print("default constructor called!");
    }


    // conversion constructor
    explicit Container(int length): 
        length(length),
        data(new double[length])
    {
        print("conversion constructor called!");
    }


    // unified initialization constructor
    explicit Container(const std::initializer_list<double>& list): 
        Container((int)list.size())
    {
        std::uninitialized_copy(list.begin(), list.end(), data);
        print("unified initialization constructor called!");
    }


    // copy constructor
    explicit Container(const Container& other): 
        Container(other.length)
    {
        for (int i = 0; i < other.length; ++i)
            data[i] = other.data[i];
        print("copy constructor called!");
    }


    // move constructor
    Container(Container&& other): 
        length(other.length),
        data(other.data) // shallow copy: copy the pointer, pointing to the same resources
    {
        other.length = 0;
        other.data = nullptr; // IMPORTANT step
        print("move constructor called!");
    }


    // destructor
    ~Container()
    {
        delete[] data;
        data = nullptr;
        length = 0;
        print("destructor called!");
    }


    // copy assignment operator
    Container& operator=(const Container& other)
    {
        if (this != &other) // guard: a = a is meaningless
        {
            delete[] data; // IMPORTANT: clean up the resources held by the old object
            data = nullptr; // data is going to be used again later, do not make it a wild pointer, set it to nullptr
            length = other.length;
            data = new double[other.length];
            for (int i = 0; i < other.length; ++i)
                data[i] = other.data[i];
        }
        print("copy assignment operator called!");
        return *this;
    }


    // move assignment operator
    Container& operator=(Container&& other)
    {
        if (this != &other)
        {
            delete[] data; // clean up the resources held by the old object
            length = other.length;
            data = other.data; // shallow copy: copy the pointer, now data also points to the memory of other.data
            other.length = 0;
            other.data = nullptr; // set other.data to nullptr
        }
        print("move assignment operator called!");
        return *this;
    }


    // operator+
    Container operator+(const Container& other)
    {
        Container res;
        res.length = this->length;
        for (int i = 0; i < length; ++i)
        {
            res.data[i] = this->data[i] + other.data[i];
        }
        print("operator+ called!");
        return res;
    }


    // Print function
    void print(const char* info)
    {
        std::cout << info << '\n';
    }

private:
    // Attributes
    double* data;
    int length;
};