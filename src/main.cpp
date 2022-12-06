#include "class_constructors.hpp"

// The global main function that is the designated start of the program
int main() {

    Container a({ 1, 2, 3, 4 });
    a.info();
    a.print();

    // Conversion constructor is no longer available due to the use of
    // the explicit specifier in the constructor defined at line 23
    //Container b = {1, 2, 3, 4};
    //b.info(); b.print();

    Container c(3);
    c.info();
    c.print();

    // Conversion constructor is no longer available due to the use of
    // the explicit specifier in the constructor defined at line 15
    //Container d = 3;
    //d.info();
    //d.print();

    // User-defined copy constructor
    Container e(a);
    e.info();
    e.print();

    // User-defined move constructor
    Container f(std::move(a));
    f.info();
    f.print();

    // Container a should be empty now!
    a.info();
    a.print();

    // User-defined copy assignment operator
    Container g(0);
    g = e;
    g.info();
    g.print();

    // User-defined move assignment operator
    Container h(0);
    h = std::move(f);
    h.info();
    h.print();

    // Container f should be empty now!
    f.info();
    f.print();

    // Return code 0 to the operating system (= no error)
    return 0;
}
