#include "binarysearch.h"
#include <iostream>

int main()
{
	vector<int> sorted_array = { 1,2,3,4,5,6 };
	int target = 1;
	BinarySearch<int> bs(sorted_array, target);
	std::cout << bs.search() << '\n';
	return 0;
}
