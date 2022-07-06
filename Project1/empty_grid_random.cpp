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
      for (int index = 0; index < focus.biggerThan.size() && !anyEmptyDomain; index++)
      {
        Grid dummy = *focus.biggerThan[index];
        modifiedGrids.push_back(dummy);
        // std::cout<<"focus.biggerThan's pos: "<<focus.biggerThan[index]->position.first<<","<<focus.biggerThan[index]->position.second<<std::endl;
        for (int number = focus.value; number <= domainRange && !anyEmptyDomain; number++)
        {

          if (dummy.domain[number] != 0)
          {
            focus.biggerThan[index]->domain[number] = 0;
            focus.biggerThan[index]->domain[0]--;

            if (focus.biggerThan[index]->domain[0] == 0)
              anyEmptyDomain = true;
          }
        }
      }

      // for the larger grids
      for (int index = 0; index < focus.smallerThan.size() && !anyEmptyDomain; index++)
      {
        Grid dummy = *focus.smallerThan[index];
        modifiedGrids.push_back(dummy);
        // std::cout<<"focus.smalerThan's pos: "<<focus.smallerThan[index]->position.first<<","<<focus.smallerThan[index]->position.second<<std::endl;
        for (int number = 1; number <= focus.value && !anyEmptyDomain; number++)
        {
          if (dummy.domain[number] != 0)
          {
            focus.smallerThan[index]->domain[number] = 0;
            focus.smallerThan[index]->domain[0]--;

            if (focus.smallerThan[index]->domain[0] == 0)
              anyEmptyDomain = true;
          }
        }
      }

      // a number cannot be used again in the same row or column
      // std::cout<<"focus pos: "<<focus.position.first<<","<<focus.position.second<<std::endl;
      // std::cout<<"row.size: "<<focus.row.size()<<std::endl;
      for (int index = 0; index < focus.row.size() && !anyEmptyDomain; index++) // domain pruning on column
      {
        // std::cout<<"focus.row's pos: "<<focus.row[index]->position.first<<","<<focus.row[index]->position.second<<std::endl;
        if (focus.row[index]->domain[focus.value] != 0)
        {
          Grid dummy = *focus.row[index];
          modifiedGrids.push_back(dummy);
          focus.row[index]->domain[focus.value] = 0;
          focus.row[index]->domain[0]--;

          if (focus.row[index]->domain[0] == 0)
            anyEmptyDomain = true;
        }
        // std::cout<<"focus.column's pos: "<<focus.column[index]->position.first<<","<<focus.column[index]->position.second<<std::endl;
        if (focus.column[index]->domain[focus.value] != 0)
        {
          Grid dummy = *focus.column[index];
          modifiedGrids.push_back(dummy);
          focus.column[index]->domain[focus.value] = 0;
          focus.column[index]->domain[0]--;

          if (focus.column[index]->domain[0] == 0)
            anyEmptyDomain = true;
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
        emptyGrids.push_back(std::make_pair(row, col));

        if (matrix[row][col] != -1)
        {
          initial_grids.push_back(std::make_pair(row, col));
          emptyGrids.erase(emptyGrids.end());
        }
      }

      gridMap.push_back(dummy);
    }
    // end of initiazation

    for (int index = 0; index < constraints.size(); index++)
    {
      // adding to the node the address of smaller nodes
      gridMap[constraints[index].first.first][constraints[index].first.second].biggerThan.push_back(&gridMap[constraints[index].second.first][constraints[index].second.second]);

      if (gridMap[constraints[index].first.first][constraints[index].first.second].domain[1] != 0)
      {
        gridMap[constraints[index].first.first][constraints[index].first.second].domain[1] = 0;
        gridMap[constraints[index].first.first][constraints[index].first.second].domain[0]--;
      }
      // adding to the node the address of the larger nodes
      gridMap[constraints[index].second.first][constraints[index].second.second].smallerThan.push_back(&gridMap[constraints[index].first.first][constraints[index].first.second]);

      if (gridMap[constraints[index].second.first][constraints[index].second.second].domain[domainRange] != 0)
      {
        gridMap[constraints[index].second.first][constraints[index].second.second].domain[domainRange] = 0;
        gridMap[constraints[index].second.first][constraints[index].second.second].domain[0]--;
      }
    }

    // binding nodes to each other

    for (int r = 0; r < domainRange; r++)
    {
      for (int c = 0; c < domainRange; c++)
      {

        for (int row = 0; row < domainRange; row++)
        {
          if (row != r)
          {

            gridMap[r][c].column.push_back(&gridMap[row][c]);
          }
        }

        for (int col = 0; col < domainRange; col++)
        {
          if (col != c)
          {
            // std::cout<<"r - c -col: "<<r<<"-"<<c<<"-"<<col<<std::endl;
            gridMap[r][c].row.push_back(&gridMap[r][col]);
          }
          // std::cout<<"size of row "<<gridMap[r][c].row.size()<<std::endl;
        }
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
    {
      for (int col = 0; col < domainRange; col++)
      {
        matrix[row][col] = gridMap[row][col].value;
      }
    }
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

  bool puzzleSolver()
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
          if (puzzleSolver())
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
  obj.puzzleSolver();
  obj.matrixCallback(matrix);
  // std::cout<<"Empty Grid Count "<<obj.emptyGridCount<<std::endl;
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

  // solve_very_brute(matrix, constraints, size);
  auto start = std::chrono::high_resolution_clock::now();
  solve(matrix, constraints, size);
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