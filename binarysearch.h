#pragma once
#include <vector>
#include <cassert>

using std::vector;

template<typename T>
class BinarySearch
{
public:
	BinarySearch() = default;  /* constructor can't be virtual */
	BinarySearch(
		/* sorted array, store all the values */ const vector<T>& p_array, 
		/* target value */ const T& p_target)
	{
		m_array  = p_array;
		m_target = p_target;
	}
	virtual ~BinarySearch() = default;
public:
	/*
	* the search range is [0, n)
	* example:
	* [0, 1, 2, 3, 4, 5]
	* (1) low = 0, high = 6 (NOT 5) mid = 3
	*     [0, 1, 2] 3 [4, 5]
	*/
	int search()
	{
		int N = m_array.size();
		assert(N >= 2);

		int low = 0, high = N;
		while (low < high)  /* search range: [0, N) */
		{
			int mid = (high + low) / 2;
			if (m_target == m_array[mid])return mid;
			else if (m_target < m_array[mid])high = mid;
			else low = mid + 1;
		}
		return INT_MIN; // not found		
	}

protected:
	vector<T> m_array;  /* sorted array, where everything stores */
	T m_target;  /* the target value */
};
