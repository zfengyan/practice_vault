#include "geometry.hpp"

int main()
{
	Point<int> p1(2.0, 2.0);
	Point<double> p2(1.0, 1.0);
	std::cout << (p1 + p2);
	return 0;
}