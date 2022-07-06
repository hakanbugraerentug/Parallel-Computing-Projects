#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
//------for time------
#include <chrono>
#include <cstdlib>

class Grid
{
public:
  int value;
  std::pair<int, int> position;
  std::vector<int> domain; // first index is the size of domain
  std::vector<Grid *> biggerThan;
  std::vector<Grid *> smallerThan;
  std::vector<Grid *> column;
  std::vector<Grid *> row;

  Grid(int size, int val, std::pair<int, int> pos)
  {
    domain.push_back(size);
    for (int i = 1; i <= size; i++)
      domain.push_back(1);

    value = val;
    position = pos;
  }

  bool anyEmptyDomain()
  {
    for (int index = 0; index < biggerThan.size(); index++)
    {
      if (biggerThan[index]->domain[0] == 0)
        return true;
    }

    for (int index = 0; index < smallerThan.size(); index++)
    {
      if (smallerThan[index]->domain[0] == 0)
        return true;
    }

    for (int index = 0; index < column.size(); index++)
    {
      if (column[index]->domain[0] == 0)
        return true;
    }

    for (int index = 0; index < row.size(); index++)
    {
      if (row[index]->domain[0] == 0)
        return true;
    }

    return false;
  }

  std::vector<int> possibleValues()
  {
    std::vector<int> possiblevalues;
    for (int i = 1; i <= domain.size() - 1; i++)
      if (domain[i] == 1)
        possiblevalues.push_back(i);

    return possiblevalues;
  }
};

class Puzzle
{
public:
  int domainRange;
  int emptyGridCount;
  std::vector<std::vector<Grid>> gridMap;
  std::vector<std::pair<int, int>> emptyGrids;
  bool anyEmptyDomain;

  void arrangeDomains(Grid focus, std::vector<Grid> &modifiedGrids, bool forward)
  {
    if (forward)
    {
      modifiedGrids.clear();

      // according to the value of a grid, other grids' domains that have comparision relation with that grid
      // should be pruned accordingly.

      // for the smaller grids
      for (int index = 0; index < focus.biggerThan.size(); index++)
      {
        Grid dummy = *focus.biggerThan[index];
        modifiedGrids.push_back(dummy);
        // std::cout<<"focus.biggerThan's pos: "<<focus.biggerThan[index]->position.first<<","<<focus.biggerThan[index]->position.second<<std::endl;
        for (int number = focus.value; number <= domainRange; number++)
        {

          if (dummy.domain[number] != 0)
          {
            focus.biggerThan[index]->domain[number] = 0;
            focus.biggerThan[index]->domain[0]--;

            if (focus.biggerThan[index]->domain[0] == 0)
            {
              anyEmptyDomain = true;
              return;
            }
          }
        }
      }

      // for the larger grids
      for (int index = 0; index < focus.smallerThan.size(); index++)
      {
        Grid dummy = *focus.smallerThan[index];
        modifiedGrids.push_back(dummy);
        // std::cout<<"focus.smalerThan's pos: "<<focus.smallerThan[index]->position.first<<","<<focus.smallerThan[index]->position.second<<std::endl;
        for (int number = 1; number <= focus.value; number++)
        {
          if (dummy.domain[number] != 0)
          {
            focus.smallerThan[index]->domain[number] = 0;
            focus.smallerThan[index]->domain[0]--;

            if (focus.smallerThan[index]->domain[0] == 0)
            {
              anyEmptyDomain = true;
              return;
            }
          }
        }
      }

      // a number cannot be used again in the same row or column
      // std::cout<<"focus pos: "<<focus.position.first<<","<<focus.position.second<<std::endl;
      // std::cout<<"row.size: "<<focus.row.size()<<std::endl;
      for (int index = 0; index < focus.row.size(); index++) // domain pruning on column
      {
        // std::cout<<"focus.row's pos: "<<focus.row[index]->position.first<<","<<focus.row[index]->position.second<<std::endl;
        if (focus.row[index]->domain[focus.value] != 0)
        {
          Grid dummy = *focus.row[index];
          modifiedGrids.push_back(dummy);
          focus.row[index]->domain[focus.value] = 0;
          focus.row[index]->domain[0]--;

          if (focus.row[index]->domain[0] == 0)
          {
            anyEmptyDomain = true;
            return;
          }
        }
        // std::cout<<"focus.column's pos: "<<focus.column[index]->position.first<<","<<focus.column[index]->position.second<<std::endl;
        if (focus.column[index]->domain[focus.value] != 0)
        {
          Grid dummy = *focus.column[index];
          modifiedGrids.push_back(dummy);
          focus.column[index]->domain[focus.value] = 0;
          focus.column[index]->domain[0]--;

          if (focus.column[index]->domain[0] == 0)
          {
            anyEmptyDomain = true;
            return;
          }
        }
      }
    }

    else
    { // for iteration recalls
      // std::cout<<"modifiedGrid size: "<<modifiedGrids.size()<<std::endl;

      anyEmptyDomain = false;
      for (int index = 0; index < modifiedGrids.size(); index++)
      {
        gridMap[modifiedGrids[index].position.first][modifiedGrids[index].position.second] = modifiedGrids[index];
      }
    }
  }

  void create_GridMap(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints, int size)
  {
    domainRange = size;
    emptyGridCount = domainRange * domainRange;
    anyEmptyDomain = false;

    std::vector<std::pair<int, int>> initial_grids;
    std::vector<Grid> dummy;
    // initializing the gridMap
    for (int row = 0; row < domainRange; row++)
    {
      dummy.clear();
      for (int col = 0; col < domainRange; col++)
      {
        Grid new_grid(domainRange, matrix[row][col], std::make_pair(row, col));
        dummy.push_back(new_grid);

        if (matrix[row][col] != -1)
          initial_grids.push_back(std::make_pair(row, col));

        else
          emptyGrids.push_back(std::make_pair(row, col));
      }

      gridMap.push_back(dummy);
    }
    // end of initiazation

    for (int index = 0; index < constraints.size(); index++)
    {
      // adding to the node the address of smaller nodes
      Grid *bigger = &gridMap[constraints[index].first.first][constraints[index].first.second];
      Grid *smaller = &gridMap[constraints[index].second.first][constraints[index].second.second];

      bigger->biggerThan.push_back(smaller);
      if (bigger->domain[1] != 0)
      {
        bigger->domain[1] = 0;
        bigger->domain[0]--;
      }
      // adding to the node the address of the larger nodes
      smaller->smallerThan.push_back(bigger);

      if (smaller->domain[domainRange] != 0)
      {
        smaller->domain[domainRange] = 0;
        smaller->domain[0]--;
      }
    }

    // binding nodes to each other
    for (int r = 0; r < domainRange; r++)
    {
      for (int c = 0; c < domainRange; c++)
      {
        for (int row = 0; row < domainRange; row++)
          if (row != r)
            gridMap[r][c].column.push_back(&gridMap[row][c]);

        for (int col = 0; col < domainRange; col++)
          if (col != c)
            gridMap[r][c].row.push_back(&gridMap[r][col]);
      }
    }

    for (int index = 0; index < initial_grids.size(); index++)
    {
      arrangeDomains(gridMap[initial_grids[index].first][initial_grids[index].second], dummy, true);
      emptyGridCount--;
    }
  }

  void matrixCallback(int **&matrix)
  {
    for (int row = 0; row < domainRange; row++)
      for (int col = 0; col < domainRange; col++)
        matrix[row][col] = gridMap[row][col].value;
  }

  void debug(const Grid focus)
  {
    std::cout << "Focus Pos: " << focus.position.first << "," << focus.position.second << std::endl;

    std::cout << "----GridMap----" << std::endl;

    for (int i = 0; i < focus.domain.size() - 1; i++)
    {
      for (int j = 0; j < focus.domain.size() - 1; j++)
      {
        std::cout << gridMap[i][j].value << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------"
              << "Empty Grid Count " << emptyGridCount << std::endl;
  }

  // pick the grid that has the minumum domain number and biggest affect to its neighbors
  bool minimum_remaining_max_effect()
  {
    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<int> possible_values;
      std::vector<int> indexes;
      std::vector<Grid> modified_grids;
      std::vector<std::pair<int, int>> dummy;

      int deleted_index;

      std::vector<std::pair<int, int>> minDomain;
      minDomain.push_back(gridMap[emptyGrids[0].first][emptyGrids[0].second].position);
      indexes.push_back(0);

      for (int index = 1; index < emptyGrids.size(); index++)
      {
        if (gridMap[emptyGrids[index].first][emptyGrids[index].second].domain[0] == gridMap[minDomain[0].first][minDomain[0].second].domain[0])
        {
          minDomain.push_back(gridMap[emptyGrids[index].first][emptyGrids[index].second].position);
          indexes.push_back(index);
        }

        else if (gridMap[emptyGrids[index].first][emptyGrids[index].second].domain[0] < gridMap[minDomain[0].first][minDomain[0].second].domain[0])
        {
          minDomain.clear();
          indexes.clear();
          minDomain.push_back(gridMap[emptyGrids[index].first][emptyGrids[index].second].position);
          indexes.push_back(index);
        }
      }

      int min_power = 2 * domainRange * domainRange;
      std::pair<int, int> min_power_pos;
      int min_power_index = 0;
      int min_power_find;

      for (int index = 0; index < minDomain.size(); index++)
      {
        Grid *dummy = &gridMap[minDomain[index].first][minDomain[index].second];
        min_power_find = 0;
        for (int i = 0; i < dummy->column.size(); i++)
        {
          min_power_find += dummy->column[i]->domain[0];
          min_power_find += dummy->row[i]->domain[0];
        }

        if (min_power_find < min_power)
        {
          min_power = min_power_find;
          min_power_pos = minDomain[index];
          min_power_index = indexes[index];
        }
      }

      Grid *focus = &gridMap[min_power_pos.first][min_power_pos.second];
      dummy = emptyGrids;

      possible_values = focus->possibleValues();
      emptyGrids.erase(emptyGrids.begin() + min_power_index);

      for (int index = 0; index < possible_values.size(); index++)
      {
        // debug(*focus);

        modified_grids.clear();
        focus->value = possible_values[index];

        arrangeDomains(*focus, modified_grids, true);
        emptyGridCount--;

        // debug(*focus);

        if (!anyEmptyDomain)
          if (minimum_remaining_max_effect())
            return true;

        arrangeDomains(*focus, modified_grids, false);
        focus->value = -1;

        emptyGridCount++;
      }
      emptyGrids = dummy;
      return false;
    }
  }

  // pick the grid that has the minimum domain number and located in the end of the empty grids set
  bool min_remaining_domain()
  {
    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<int> possible_values;
      std::vector<Grid> modified_grids;
      std::vector<std::pair<int, int>> dummy;

      int deleted_index;

      std::pair<int, int> minDomain = gridMap[emptyGrids[0].first][emptyGrids[0].second].position;
      for (int index = 0; index < emptyGrids.size(); index++)
      {
        if (gridMap[emptyGrids[index].first][emptyGrids[index].second].domain[0] <= gridMap[minDomain.first][minDomain.second].domain[0])
        {
          minDomain = gridMap[emptyGrids[index].first][emptyGrids[index].second].position;
          deleted_index = index;
        }
      }

      Grid *focus = &gridMap[minDomain.first][minDomain.second];
      dummy = emptyGrids;

      possible_values = focus->possibleValues();
      emptyGrids.erase(emptyGrids.begin() + deleted_index);

      for (int index = 0; index < possible_values.size(); index++)
      {
        // debug(*focus);

        modified_grids.clear();
        focus->value = possible_values[index];

        arrangeDomains(*focus, modified_grids, true);
        emptyGridCount--;

        // debug(*focus);

        if (!anyEmptyDomain)
          if (min_remaining_domain())
            return true;

        arrangeDomains(*focus, modified_grids, false);
        focus->value = -1;

        emptyGridCount++;
      }
      emptyGrids = dummy;
      return false;
    }
  }

  // pick the grid randomly among the empty grids that have the minumum domain number
  bool min_remaining_domain_random()
  {
    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<int> possible_values;
      std::vector<int> indexes;
      std::vector<Grid> modified_grids;
      std::vector<std::pair<int, int>> dummy;

      int deleted_index;

      std::vector<std::pair<int, int>> minDomain;
      minDomain.push_back(gridMap[emptyGrids[0].first][emptyGrids[0].second].position);
      indexes.push_back(0);

      for (int index = 1; index < emptyGrids.size(); index++)
      {
        if (gridMap[emptyGrids[index].first][emptyGrids[index].second].domain[0] == gridMap[minDomain[0].first][minDomain[0].second].domain[0])
        {
          minDomain.push_back(gridMap[emptyGrids[index].first][emptyGrids[index].second].position);
          indexes.push_back(index);
        }

        else if (gridMap[emptyGrids[index].first][emptyGrids[index].second].domain[0] < gridMap[minDomain[0].first][minDomain[0].second].domain[0])
        {
          minDomain.clear();
          indexes.clear();
          minDomain.push_back(gridMap[emptyGrids[index].first][emptyGrids[index].second].position);
          indexes.push_back(index);
        }
      }

      int rand = std::rand() % minDomain.size();

      Grid *focus = &gridMap[minDomain[rand].first][minDomain[rand].second];
      dummy = emptyGrids;

      possible_values = focus->possibleValues();
      emptyGrids.erase(emptyGrids.begin() + indexes[rand]);

      for (int index = 0; index < possible_values.size(); index++)
      {
        // debug(*focus);

        modified_grids.clear();
        focus->value = possible_values[index];

        arrangeDomains(*focus, modified_grids, true);
        emptyGridCount--;

        // debug(*focus);

        if (!anyEmptyDomain)
          if (min_remaining_domain_random())
            return true;

        arrangeDomains(*focus, modified_grids, false);
        focus->value = -1;

        emptyGridCount++;
      }
      emptyGrids = dummy;
      return false;
    }
  }

  // pick the grid randomy among the empty grids set
  bool empty_grids_random()
  {

    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<std::pair<int, int>> dummy;
      std::vector<int> possible_values;
      std::vector<Grid> modified_grids;

      int rand = std::rand() % emptyGridCount;
      Grid *focus = &gridMap[emptyGrids[rand].first][emptyGrids[rand].second];

      possible_values = focus->possibleValues();
      dummy = emptyGrids;

      emptyGrids.erase(emptyGrids.begin() + rand);
      for (int index = 0; index < possible_values.size(); index++)
      {
        // debug(*focus);

        modified_grids.clear();
        focus->value = possible_values[index];

        arrangeDomains(*focus, modified_grids, true);
        emptyGridCount--;

        // debug(*focus);

        if (!anyEmptyDomain)
          if (empty_grids_random())
            return true;

        arrangeDomains(*focus, modified_grids, false);
        focus->value = -1;

        emptyGridCount++;
      }
      // emptyGrids.push_back(focus->position);
      emptyGrids = dummy;
      return false;
    }
  }

  // traverse through grid map and pick the first available grid
  bool first_available()
  {
    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<int> possible_values;
      std::vector<Grid> modified_grids;
      for (int row = 0; row < domainRange; row++)
      {
        for (int col = 0; col < domainRange; col++)
        {
          Grid *focus = &gridMap[row][col];

          if (focus->value == -1)
          {
            possible_values = focus->possibleValues();

            for (int index = 0; index < possible_values.size(); index++)
            {
              // debug(*focus);

              modified_grids.clear();
              focus->value = possible_values[index];
              arrangeDomains(*focus, modified_grids, true);
              emptyGridCount--;

              // debug(*focus);

              if (!focus->anyEmptyDomain())
                if (first_available())
                  return true;

              arrangeDomains(*focus, modified_grids, false);
              focus->value = -1;
              emptyGridCount++;
            }
            return false;
          }
        }
      }

      return false;
    }
  }

  // pick the first grid in empty grids set
  bool empty_grids_first()
  {
    if (emptyGridCount == 0)
    {
      return true;
    }
    else
    {
      std::vector<int> possible_values;
      std::vector<Grid> modified_grids;
      std::vector<std::pair<int, int>> dummy;

      Grid *focus = &gridMap[emptyGrids[0].first][emptyGrids[0].second];
      dummy = emptyGrids;

      possible_values = focus->possibleValues();
      emptyGrids.erase(emptyGrids.begin());
      for (int index = 0; index < possible_values.size(); index++)
      {
        // debug(*focus);

        modified_grids.clear();
        focus->value = possible_values[index];

        arrangeDomains(*focus, modified_grids, true);
        emptyGridCount--;

        // debug(*focus);

        if (!anyEmptyDomain)
          if (empty_grids_first())
            return true;

        arrangeDomains(*focus, modified_grids, false);
        focus->value = -1;

        emptyGridCount++;
      }
      // emptyGrids.push_back(focus->position);
      emptyGrids = dummy;
      return false;
    }
  }
};

void solve(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> constraints, int size)
{
  Puzzle obj;
  obj.create_GridMap(matrix, constraints, size);

  // *********** Possible Functions **************
  obj.minimum_remaining_max_effect();
  // obj.min_remaining_domain();
  // obj.min_remaining_domain_random();
  // obj.empty_grids_random();
  // obj.first_available();
  // obj.empty_grids_first();
  // *********** Possible Functions ***************

  obj.matrixCallback(matrix);
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

// Very Brute Force
bool solve_very_brute(int **&matrix, const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &constraints, int size)
{
  int counter = 1;

  if (solved(matrix, constraints, size))
    return true;
  else
  {
    for (int row = 0; row < size; row++)
    {
      for (int col = 0; col < size; col++)
      {
        if (matrix[row][col] == -1)
        {
          while (counter <= size)
          {
            matrix[row][col] = counter;
            if (solve_very_brute(matrix, constraints, size))
              return true;

            matrix[row][col] = -1;
            counter++;
          }
          return false;
        }
      }
    }
    return false;
  }
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

  // solve(matrix, constraints, size)
  // Measure only solution function's execution time

  auto start = std::chrono::high_resolution_clock::now();
  // solve(matrix, constraints, size);
  solve_very_brute(matrix, constraints, size);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Solving Time " << duration.count() << " microseconds." << std::endl;

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }

  // DELETE//
  for (int i = 0; i < size; i++)
  {
    delete matrix[i];
  }

  delete[] matrix;
  // DELETE//

  return 0;
}