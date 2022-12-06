#include "geometry.hpp"

int main()
{
	Point p1(2.0, 2.0);
	int a = 1, b = 2;
	Point p2(a, b);
	Point p = { 1,2 };
	std::cout << p2.distance() << '\n';
	std::cout << p2.distance(p1) << '\n';

	Point p_old(2.0, 0.0);
	Point p_new = p_old.rotated(M_PI / 2);
	p_new.print();

	Point p_self_old(2.0, 0.0);
	p_self_old.rotate(M_PI / 2);
	p_self_old.print();

	return 0;
}