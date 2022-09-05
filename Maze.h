/*
use list to define a maze:

simple maze:
start from position 0, end at position 8
|--------|
 0  1  2 |            
|  ------|            
|3  4  5 |     
| |   ---|            
|6| 7  8
|--------|

maze[0] stores 1 and 3 --> indicating 1 and 3 are the viable neighbors of pos 0;
maze[2] stores only 1  --> indicating 1 is the only viable neighbor of pos 2;

*/


#pragma once


#include<list>
#include<unordered_set>
#include<iterator>


template<typename T>
class Maze
{
public:
	Maze() = default;  // default constructor
	Maze(bool flag)
	{
		if (flag)  // if flag is true, construct the default maze
		{
			maze[0].emplace_back(); maze[0].back() = T(1);
			maze[0].emplace_back(); maze[0].back() = T(3);

			maze[1].push_back(0);
			maze[1].push_back(2);

			maze[2].push_back(1);

			maze[3].push_back(0);
			maze[3].push_back(4);
			maze[3].push_back(6);

			maze[4].push_back(3);
			maze[4].push_back(5);
			maze[4].push_back(7);

			maze[5].push_back(4);

			maze[6].push_back(3);

			maze[7].push_back(4);
			maze[7].push_back(8);

			maze[8].push_back(7);

			/*
			* use list to store the maze
			* why use list?
			* @pro:
			*	fast inserting / removing elements:
			*	lists perform generally better in inserting, extracting and moving elements in any position
			*	within the container for which an iterator has already been obtained,
			*	and therefore also in algorithms that make intensive use of these, like sorting algorithms.
			* @con:
			*	(1) they lack direct access to the elements by their position
			*	(2) consume some extra memory to keep the linking information associated to each element
			*		(which may be an important factor for large lists of small-sized elements).
			*
			* path list needs to be dynamically maintained (inserting / removing)
			* and the example maze we use is not a big maze, thus using list seems to be a relatively good choice.
			*/
		}
	}


	std::list<T> get_neighbors(T pos) const  // get the neighbors for one specific position
	{
		return maze[pos];
	}


	/*
	solve the maze using backtrack
	@param:
		maze, start position, end position
	@return:
		a list containing the positions that make up the path (if any)
	 

	simple maze:
	start from position 0, end at position 8
	|--------|
	 0  1  2 |            process of backtracking:
	|  ------|            path: [0 3 4 5]  visited: [0 1 2 3 4 5] 5 -> wall
	|3  4  5 |            backtrack to the "latest" and "can reach other positions" node -> 4
	| |   ---|            path: [0 3 4]    visited: [0 1 2 3 4 5]
	|6| 7  8              
	|--------|
	need two list: (a)path list  (b)visited list
	e.g.
	path: [0 3 4]  visited: [0 1 2 3 4]
	in some search problems, keeping a visited list is for improving the performance
	in this problem, keeping a visited list is meant for avoiding endless loop
	*/
	std::list<T> solve_maze_backtracking(T start, T end)
	{
		std::unordered_set<T> visited;  // visited list
		std::list<T> path;  // sotre the path

		// initialize
		path.emplace_back();
		path.back() = start;  // add the start position to the path
		T currentPos = start;
		visited.insert(currentPos);  // add the current position to the visited list

		// loop for processing: when reaching the end position or the path is empty
		while (path.back() != end && path.empty() == false)
		{
			typename std::list<T>::iterator iter = maze[currentPos].begin();  // iterate all neighbors of currentPos
			bool foundOutlet = false;  // when an outlet is found, mark it as true

			// inner loop: for currentPos, when reaching the last neighbor or an outlet is found
			while (iter != maze[currentPos].end() && (!foundOutlet))
			{
				// check if connection leads to unvisited point
				if (visited.count(*iter) == 0)
				{
					foundOutlet = true;
				}
				else
				{
					iter++;
				}
			}

			// outside the inner loop
			if (foundOutlet)  // if an outlet is found, add it to path list and visited set
			{
				path.emplace_back();
				path.back() = *iter;  // add the outlet to the path list
				visited.insert(*iter);  // add the outlet to the visited set
			}
			else 
			{
				path.pop_back();
			}

			// last step: set the currentPos to the last element of the path list
			currentPos = path.back();
		}
		// end of the outer loop

		// return the path
		return path;

	}


	virtual ~Maze() = default;  // virtual key word for inheritance reason

protected:
	std::list<T> maze[9];  // a list store 9 lists, i.e. maze[0] is a list, stores 1 and 3(its neighbors)
};