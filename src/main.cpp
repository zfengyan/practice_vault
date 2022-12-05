#include "geometry.hpp"

int main()
{
	Point<int> p1(2.0, 2.0);
	Point<double> p2(1.0, 2.0);
	std::cout << p2.distance() << '\n';
	std::cout << p2.distance(p1) << '\n';
	return 0;
}