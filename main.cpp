#include "binarysearch.h"
#include "Maze.h"
#include <iostream>

/* container_type */
typedef vector<int> container_type;

/* example: a generic traversing function */
std::ostream& operator<<(std::ostream& os, const container_type& c)
{
	for (container_type::const_iterator it = c.begin(); it != c.end(); ++it)
		os << *it << " ";
	return os;
}



int main()
{
	/*vector<int> sorted_array = { 1,2,3,4,5,6 };
	std::cout << "input sorted array is\n";
	std::cout << sorted_array << '\n';
	int target = 1;
	std::cout << "target value is: " << target << '\n';
	BinarySearch<int> bs(sorted_array, target);
	std::cout << "search result (return index if exists): " << bs.search() << '\n';
	return 0;*/

	bool construct_default_maze = true;
	Maze<int> maze(construct_default_maze);  // construct the default maze, type: int

	// alter the parameters here: maze.solve_maze_backtracking(start, end)
	std::list<int> path = maze.solve_maze_backtracking(0, 8);
	// alter done

	std::cout << "the size of path list is: " << path.size() << '\n';
	std::cout << "the path is: " << '\n';
	std::list<int>::iterator iter = path.begin();
	while (iter != path.end())
	{
		std::cout << *iter << " ";
		++iter;
	}
		
	return 0;
}
