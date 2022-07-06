#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <bits/stdc++.h>
#include <bitset>
#include <omp.h>
//------for time------
#include <chrono>
#include <cstdlib>

// Class that keeps info about initial grid map and relationships between nodes
class Grid
{
public:
  short int value;
  unsigned short int id;
  unsigned short int domain;
  std::vector<unsigned short int> dependedGrids;
  std::vector<unsigned short int> biggerThan;
  std::vector<unsigned short int> smallerThan;

  Grid() {}

  Grid(unsigned short int size, unsigned short int val, unsigned short int pos_id)
  {
    this->domain = pow(2, size) - 1;
    this->value = val;
    this->id = pos_id;
  }

  std::vector<int> possibleValues()
  {
    std::vector<int> possiblevalues;
    auto s = std::bitset<12>(this->domain);
    for (int i = 0; i <= s.size(); i++)
      if (s[i])
        possiblevalues.push_back(i + 1);

    return possiblevalues;
  }
};

// Class that keeps info about grid map contains only value and domain of grids
class miniGrid
{
public:
  short int value;
  unsigned short int domain;

  std::vector<int> possibleValues()
  {
    std::vector<int> possiblevalues;
    auto s = std::bitset<12>(this->domain);
    for (int i = 0; i <= s.size(); i++)
      if (s[i])
        possiblevalues.push_back(i + 1);

    return possiblevalues;
  }
};

// GridMapState that keeps empty nodes IDs and gridmap consisted of miniGrids
class GridState
{
public:
  std::vector<miniGrid> gridMap;
  std::vector<unsigned short int> emptyGrids;

  GridState() {}

  GridState(const GridState &onRight)
  {
    this->gridMap = onRight.gridMap;
    this->emptyGrids = onRight.emptyGrids;
  }
};

// Creating gridmap, arrange domain after every iteration and solver functions.
class PuzzleSolver
{
public:
  std::vector<Grid> initialGridMap;
  std::vector<unsigned short int> initialEmptyGrids;
  GridState final;
  bool done = false;
  unsigned short int domainRange;

  // used at every change
  bool arrangeDomainsForward(Grid &focus, GridState &GridMapInUse, short int value)
  {
    bool anyEmptyGridCheck = true;
    unsigned short int dummy;

    // according to the value of a grid, other grids' domains that have comparision relation with that grid
    // should be pruned accordingly.

    for (auto i = focus.biggerThan.begin(); i != focus.biggerThan.end() && anyEmptyGridCheck; ++i)
    {
      dummy = pow(2, domainRange) - pow(2, value - 1);
      GridMapInUse.gridMap[*i].domain = GridMapInUse.gridMap[*i].domain & (USHRT_MAX ^ dummy);
      anyEmptyGridCheck = anyEmptyGridCheck & (GridMapInUse.gridMap[*i].domain > 0);
    }

    for (auto i = focus.smallerThan.begin(); i != focus.smallerThan.end() && anyEmptyGridCheck; ++i)
    {
      dummy = pow(2, value) - 1;
      GridMapInUse.gridMap[*i].domain = GridMapInUse.gridMap[*i].domain & (USHRT_MAX ^ dummy);
      anyEmptyGridCheck = anyEmptyGridCheck & (GridMapInUse.gridMap[*i].domain > 0);
    }

    for (auto i = focus.dependedGrids.begin(); i != focus.dependedGrids.end() && anyEmptyGridCheck; ++i)
    {
      if (std::bitset<12>(GridMapInUse.gridMap[*i].domain)[value - 1] == 1)
      {
        dummy = pow(2, value - 1);
        GridMapInUse.gridMap[*i].domain = GridMapInUse.gridMap[*i].domain & (USHRT_MAX ^ dummy);
        anyEmptyGridCheck = anyEmptyGridCheck & (GridMapInUse.gridMap[*i].domain > 0);
      }
    }
    return anyEmptyGridCheck;
  }

  // used once at the beginning in createGridMap
  bool arrangeDomainsForward(Grid &focus, std::vector<Grid> &initialGridMap)
  {
    // auto start = std::chrono::high_resolution_clock::now();
    // std::cout<<"begin arrange"<<std::endl;
    bool anyEmptyGridCheck = true;
    unsigned short int dummy;
    // according to the value of a grid, other grids' domains that have comparision relation with that grid
    // should be pruned accordingly.

    for (auto i = focus.biggerThan.begin(); i != focus.biggerThan.end() && anyEmptyGridCheck; ++i)
    {
      dummy = pow(2, domainRange) - pow(2, focus.value - 1);
      initialGridMap[*i].domain = initialGridMap[*i].domain & (USHRT_MAX ^ dummy);
      anyEmptyGridCheck = anyEmptyGridCheck & (initialGridMap[*i].domain > 0);
    }

    for (auto i = focus.smallerThan.begin(); i != focus.smallerThan.end() && anyEmptyGridCheck; ++i)
    {
      dummy = pow(2, focus.value) - 1;
      initialGridMap[*i].domain = initialGridMap[*i].domain & (USHRT_MAX ^ dummy);
      anyEmptyGridCheck = anyEmptyGridCheck & (initialGridMap[*i].domain > 0);
    }

    for (auto i = focus.dependedGrids.begin(); i != focus.dependedGrids.end() && anyEmptyGridCheck; ++i)
    {
      if (std::bitset<12>(initialGridMap[*i].domain)[focus.value - 1] == 1)
      {
        dummy = pow(2, focus.value - 1);
        initialGridMap[*i].domain = initialGridMap[*i].domain & (USHRT_MAX ^ dummy);
        anyEmptyGridCheck = anyEmptyGridCheck & (initialGridMap[*i].domain > 0);
      }
    }

    return anyEmptyGridCheck;
  }

  // Creating grid map and dependencies accordingly
  void createGridMap(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints, int size) // optimal
  {
    this->domainRange = size;
    std::vector<Grid> initial_grids;

    // Initilization of grids
    for (auto index = 0; index < domainRange * domainRange; index++)
    {
      unsigned short int row = index / this->domainRange;
      unsigned short int col = index % this->domainRange;
      Grid new_grid(domainRange, matrix[row][col], index);

      if (matrix[row][col] != -1)
        initial_grids.push_back(new_grid);

      else
        initialEmptyGrids.push_back(index);

      initialGridMap.push_back(new_grid);
    }
    // end of initiazation

    // constraints
    for (auto element : constraints)
    {
      unsigned short int big = domainRange * element.first.first + element.first.second;
      unsigned short int small = domainRange * element.second.first + element.second.second;

      this->initialGridMap[big].biggerThan.push_back(small);
      this->initialGridMap[big].domain = this->initialGridMap[big].domain & (USHRT_MAX ^ 1);
      this->initialGridMap[small].smallerThan.push_back(big);
      this->initialGridMap[small].domain = this->initialGridMap[small].domain & (USHRT_MAX ^ int(pow(2, domainRange - 1)));
    }

    // binding nodes to each other
    for (auto id = 0; id < this->domainRange * this->domainRange; id++)
    {
      for (int row = -id % this->domainRange; row < this->domainRange - id % this->domainRange; row++)
        if (row != 0)
          this->initialGridMap[id].dependedGrids.push_back(id + row);

      for (int col = -id / this->domainRange; col < this->domainRange - id / this->domainRange; col++)
        if (col != 0)
          this->initialGridMap[id].dependedGrids.push_back(id + this->domainRange * col);
    }

    std::vector<Grid> dummy;
    for (auto initial : initial_grids)
      arrangeDomainsForward(this->initialGridMap[initial.id], this->initialGridMap);
  }

  void matrixCallback(int **&matrix) // optimal
  {
    for (int row = 0; row < domainRange; row++)
      for (int col = 0; col < domainRange; col++)
        matrix[row][col] = final.gridMap[domainRange * row + col].value;
  }

  void debug(const Grid &focus, const GridState &currentGridMap, int threadID) // NO NEED FOR PARALLELIZATION
  {
    std::cout << "Focus Pos: " << focus.id << ",(" << int(focus.id / domainRange) << "," << int(focus.id % domainRange) << ")" << std::endl;
    std::cout << "Thread ID: " << threadID << std::endl;
    auto s = std::bitset<12>(focus.domain);
    std::cout << "domain " << s << std::endl;
    std::cout << "----GridMap----" << std::endl;

    for (int i = 0; i < domainRange; i++)
    {
      for (int j = 0; j < domainRange; j++)
      {
        std::cout << currentGridMap.gridMap[domainRange * i + j].value << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------"
              << "Empty Grid Count " << currentGridMap.emptyGrids.size() << std::endl;
  }

  // THIS FUNCTION DOES EVERYTHING
  void parallelRecursiveSolver(GridState &GridMapInUse, unsigned int depth)
  {
    if (GridMapInUse.emptyGrids.size() == 0)
    {
      final = GridMapInUse;
#pragma omp atomic write
      done = true;
    }

    else if (!done)
    {
      std::vector<int> indexes;
      std::vector<unsigned short int> minDomain;

      minDomain.push_back(GridMapInUse.emptyGrids[0]);
      indexes.push_back(0);

      for (int index = 1; index < GridMapInUse.emptyGrids.size(); index++)
      {
        if (GridMapInUse.gridMap[GridMapInUse.emptyGrids[index]].domain == GridMapInUse.gridMap[minDomain[0]].domain)
        {
          minDomain.push_back(GridMapInUse.emptyGrids[index]);
          indexes.push_back(index);
        }

        else if (__builtin_popcount(GridMapInUse.gridMap[GridMapInUse.emptyGrids[index]].domain) < __builtin_popcount(GridMapInUse.gridMap[minDomain[0]].domain))
        {
          minDomain.clear();
          indexes.clear();
          minDomain.push_back(GridMapInUse.emptyGrids[index]);
          indexes.push_back(index);
        }
      }

      unsigned int min_power = 0;
      for (const auto &dependedGrid : this->initialGridMap[minDomain[0]].dependedGrids)
        min_power += __builtin_popcount(GridMapInUse.gridMap[dependedGrid].domain);

      unsigned short int min_power_pos = minDomain[0];
      unsigned short int min_power_index = indexes[0];
      unsigned int min_power_found;

      for (int index = 1; index < minDomain.size(); index++)
      {
        min_power_found = 0;

        for (const auto &i : initialGridMap[minDomain[index]].dependedGrids)
          min_power_found += __builtin_popcount(GridMapInUse.gridMap[i].domain);

        if (min_power_found < min_power)
        {
          min_power = min_power_found;
          min_power_pos = minDomain[index];
          min_power_index = indexes[index];
        }
      }

      GridMapInUse.emptyGrids.erase(GridMapInUse.emptyGrids.begin() + min_power_index);

      auto remaining_domain = GridMapInUse.gridMap[min_power_pos].domain;
      for (const auto &i : initialGridMap[min_power_pos].dependedGrids)
        remaining_domain = remaining_domain & (USHRT_MAX ^ GridMapInUse.gridMap[i].domain);

      if (remaining_domain != 0 && __builtin_popcount(GridMapInUse.gridMap[min_power_pos].domain) - __builtin_popcount(remaining_domain) > 0)
        GridMapInUse.gridMap[min_power_pos].domain = remaining_domain;

      std::vector<int> possible_values;
      possible_values = GridMapInUse.gridMap[min_power_pos].possibleValues();
      auto size = possible_values.size();

      //#BRANCHING STARTS HERE!!!!!
      for (auto i = 0; i < possible_values.size() && !done; i++)
      {
        depth++;
// int priorityValue =pow(2,52) - int(omp_get_wtime()) + depth;
#pragma omp task firstprivate(GridMapInUse, depth) // priority(priorityValue) //final(done) //priority(cycle + possible_values[i]) //if(possible_values.size()>1)
        {

          GridMapInUse.gridMap[min_power_pos].value = possible_values[i];
          if (!done && arrangeDomainsForward(initialGridMap[min_power_pos], GridMapInUse, GridMapInUse.gridMap[min_power_pos].value))
            parallelRecursiveSolver(GridMapInUse, depth);
        }
      }
    }
  }
};

// Adaptation of classes into main function
void solve(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints, int size)
{
  GridState initialMiniGridMap;
  PuzzleSolver solver;
  omp_set_nested(1);
  // omp_set_max_active_levels(10);
  solver.createGridMap(matrix, constraints, size);
  int cycle = 0;

  for (auto i = 0; i < solver.initialGridMap.size(); i++)
  {
    miniGrid dummy;
    dummy.value = solver.initialGridMap[i].value;
    dummy.domain = solver.initialGridMap[i].domain;
    initialMiniGridMap.gridMap.push_back(dummy);
  }

  initialMiniGridMap.emptyGrids = solver.initialEmptyGrids;
  unsigned int depth = 0;
// STARTING THE RECURSIVE FUNCTION
#pragma omp parallel
  {
#pragma omp single
    solver.parallelRecursiveSolver(initialMiniGridMap, depth);
  }
  solver.matrixCallback(matrix);
}

// USE THIS FUNCTION WHILE CHECKING FINAL RESULT
// YOU MAY FURTHER OPTIMIZE THIS FUNCTION TO USE FOR CALCULATION
bool solved(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints, int size)
{
  for (int i = 0; i < constraints.size(); i++)
  {
    if (matrix[constraints[i].first.first][constraints[i].first.second] < matrix[constraints[i].second.first][constraints[i].second.second])
      return false;
  }

  std::vector<int> rows;
  std::vector<int> cols;
  for (int rc = 0; rc < size; rc++)
  {
    for (int s = 0; s < size; s++)
    {
      rows.push_back(matrix[rc][s]);
      cols.push_back(matrix[s][rc]);
    }

    std::sort(rows.begin(), rows.end());
    std::sort(cols.begin(), cols.end());

    if ((rows[0] == -1) || (cols[0] == -1))
      return false;

    for (int i = 0; i < size - 1; i++)
    {
      if ((rows[i] == rows[i + 1]) || (cols[i] == cols[i + 1]))
      {
        return false;
      }
    }

    rows.clear();
    cols.clear();
  }

  return true;
}

// DON'T CHANGE THIS FUNCTION
void get_constraints(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &constraints, std::ifstream &file)
{

  int ctr = 1;
  std::string constraint;
  while (!file.eof())
  {
    std::getline(file, constraint);
    if (constraint != "")
    {
      std::cout << "Constraint " << ctr++ << ": " << constraint << std::endl;
      std::stringstream ss(constraint);
      int val1, val2, val3, val4; // Coordinate (val1, val2) > (val3, val4)
      ss >> val1 >> val2 >> val3 >> val4;
      constraints.push_back(std::pair<std::pair<int, int>, std::pair<int, int>>(std::pair<int, int>(val1 - 1, val2 - 1), std::pair<int, int>(val3 - 1, val4 - 1)));
    }
  }
}

// DON'T CHANGE THIS FUNCTION
void read_matrix(int **&matrix, std::ifstream &file, int size)
{
  matrix = new int *[size];

  for (int i = 0; i < size; i++)
  {
    matrix[i] = new int[size];
  }

  int val;
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      file >> val;
      matrix[i][j] = val;
    }
  }
}

int main(int argc, char **argv)
{
  std::string filename(argv[1]);

  std::ifstream file;
  file.open(filename.c_str());
  int size;

  file >> size;
  std::cout << "Size: " << size << std::endl;

  int **matrix;

  read_matrix(matrix, file, size);

  std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints;
  get_constraints(constraints, file);

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  const size_t bigger_than_cachesize = 10 * 1024 * 1024;
  long *p = new long[bigger_than_cachesize];

  double start = omp_get_wtime();
  solve(matrix, constraints, size);
  double stop = omp_get_wtime();

  std::cout << "Solving Time  "
            << ": " << stop - start << " seconds." << std::endl;

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  // Check if the solution is correct
  if (solved(matrix, constraints, size))
    std::cout << "correctly solved" << std::endl;
  else
    std::cout << "WRONG" << std::endl;

  // DELETE//
  for (int i = 0; i < size; i++)
  {
    delete matrix[i];
  }

  delete[] matrix;
  // DELETE//

  return 0;
}