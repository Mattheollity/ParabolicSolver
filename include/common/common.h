#ifndef COMMON_
#define COMMON_

#include <exception>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <list>
#include <deque>
#include <set>

#include <memory>
#include <utility>
#include <cmath>
#include <cstdlib>

#include <omp.h>

using std::cout;
using std::cin;
using std::endl;

using std::string;
using std::to_string;
using std::vector;
using std::list;
using std::pair;
using std::deque;
using std::set;
using std::unique_ptr;
using std::shared_ptr;

typedef vector<double> point;
typedef vector<int> voxel_point;

class simple_error : public std::exception
{
    public:
    simple_error(string message)
    {
        std::cerr<<message<<endl;
        std::exit(1);
    }
};

#endif