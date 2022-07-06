#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

#define SIZE 5 // Matrix size
#define INPUTSIZE 2306451
#define BLOCKSIZE 1125 // 1125 * 12 = 144000
#define THREADPERBLOCK 128

inline cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess)
  {
    printf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    if (result != cudaSuccess)
    {
      exit(1);
    }
  }
#endif
  return result;
}

__device__ bool constraintCheck(const int *matrix, const int *constraints, int *constraint_size, int *CSI, int *RC, int *index)
{

  for (int i = 0; i < *constraint_size; i++)
  {

    if ((SIZE * (constraints[(*CSI) + 4 * i] - 1) + constraints[(*CSI) + 4 * i + 1] - 1 == *RC) ||
        (SIZE * (constraints[(*CSI) + 4 * i + 2] - 1) + constraints[(*CSI) + 4 * i + 3] - 1 == *RC))
    {

      if (matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i] - 1) + constraints[(*CSI) + 4 * i + 1] - 1] != -1 &&
          matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i + 2] - 1) + constraints[(*CSI) + 4 * i + 3] - 1] != -1 &&
          matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i] - 1) + constraints[(*CSI) + 4 * i + 1] - 1] <
              matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i + 2] - 1) + constraints[(*CSI) + 4 * i + 3] - 1])
      {
        return false;
      }

      if (matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i] - 1) + constraints[(*CSI) + 4 * i + 1] - 1] == 1 ||
          matrix[25 * (*index) + SIZE * (constraints[(*CSI) + 4 * i + 2] - 1) + constraints[(*CSI) + 4 * i + 3] - 1] == SIZE)
      {
        return false;
      }
    }
  }
  return true;
}

__device__ bool checkRowColumn(const int *matrix, int *id, int *value, int *index)
{
  for (int row = -(*id % SIZE); row < SIZE - (*id % SIZE); row++)
  {
    if (row != 0 && *value == matrix[25 * (*index) + *id + row])
    {
      return true;
    }
  }

  for (int col = -(*id / SIZE); col < SIZE - (*id / SIZE); col++)
  {
    if (col != 0 && *value == matrix[25 * (*index) + *id + SIZE * col])
    {
      return true;
    }
  }

  return false;
}

__device__ void iterative(int *matrix, int *constraints, int *constraint_size, int *constraintStartIndex, int *index)
{

  int RC = 0;
  bool keepIteration = true;
  int value;

  int valuesArray[SIZE * SIZE];
  for (int i = 0; i < SIZE * SIZE; i++)
  {
    valuesArray[i] = 0;
    if (matrix[25 * (*index) + i] == -1)
    {
      valuesArray[i] = 1;
    }
  }

  while (RC < SIZE * SIZE)
  {
    if (valuesArray[RC] != 0)
    {
      while (valuesArray[RC] <= SIZE && keepIteration)
      {
        value = valuesArray[RC];
        valuesArray[RC]++;
        if (!checkRowColumn(matrix, &RC, &value, index))
        {
          matrix[25 * (*index) + RC] = value;
          if (constraintCheck(matrix, constraints, constraint_size, constraintStartIndex, &RC, index))
          {
            keepIteration = false;
          }
        }
      }

      if (!keepIteration)
      {
        RC++;
        keepIteration = true;
      }
      else
      {
        for (int i = RC; i < SIZE * SIZE; i++)
        {
          if (valuesArray[i] != 0)
          {
            valuesArray[i] = 1;
          }
        }
        matrix[25 * (*index) + RC] = -1;
        RC--;
        while (valuesArray[RC] == 0)
        {
          RC--;
        }
      }
    }
    else
    {
      RC++;
    }
  }
}

__global__ void solver(int *matrix_in, int *constraints_in, long int *constraintSums_in, int *inputGridNumber)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < *inputGridNumber)
  {
    int sizeOfConstraint = constraintSums_in[index + 1] - constraintSums_in[index];
    int constraintStartIndex = constraintSums_in[index] * 4;

    iterative(matrix_in, constraints_in, &sizeOfConstraint, &constraintStartIndex, &index);
  }
}

int isSolved(int *&matrix, int *constraints, int ConstraintSize)
{

  for (int i = 0; i < ConstraintSize; i++)
  {
    if (matrix[SIZE * (constraints[4 * i] - 1) + constraints[4 * i + 1] - 1] < matrix[SIZE * (constraints[4 * i + 2] - 1) + constraints[4 * i + 3] - 1])
      return 0;
  }

  std::vector<int> rows;
  std::vector<int> cols;
  for (int rc = 0; rc < SIZE; rc++)
  {
    for (int s = 0; s < SIZE; s++)
    {
      rows.push_back(matrix[SIZE * rc + s]);
      cols.push_back(matrix[SIZE * s + rc]);
    }

    std::sort(rows.begin(), rows.end());
    std::sort(cols.begin(), cols.end());

    if ((rows[0] == -1) || (cols[0] == -1))
      return 0;

    for (int i = 0; i < SIZE - 1; i++)
    {
      if ((rows[i] == rows[i + 1]) || (cols[i] == cols[i + 1]))
      {
        return 0;
      }
    }

    rows.clear();
    cols.clear();
  }

  return 1;
}

int main(int argc, char **argv)
{
  printf("in Oguzhan Ilter's futoshiki_hw3.cu \n");
  cudaSetDevice(0);

  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // ------------------------------- Reading Data From txt Start -----------------------------------------
  std::string filename(argv[1]);
  std::ifstream file(filename.c_str());
  std::ifstream scout(filename.c_str());

  int no_grids;
  file >> no_grids;

  int dummy;
  scout >> dummy;

  int **grids = new int *[no_grids];
  int **constraints = new int *[no_grids];

  for (int i = 0; i < no_grids; i++)
  {
    grids[i] = new int[SIZE * SIZE];
  }

  int elem0, elem1, elem2, elem3, elem4;
  int pre_cursor = 0;
  int cursor = 0;
  int csize = 0;

  std::string file_line;
  std::string scout_line;

  int *constraint_sizes = new int[no_grids];
  int maxNumberofConstraint = 0;
  std::getline(scout, scout_line); // These are for spare lines
  std::getline(scout, scout_line);
  for (int i = 0; i < INPUTSIZE; i++)
  {
    std::getline(scout, scout_line);
    if (scout_line == "-------")
    {
      csize = i - pre_cursor - 5;
      constraint_sizes[cursor] = csize;
      if (csize > maxNumberofConstraint)
      {
        maxNumberofConstraint = csize;
      }
      cursor++;
      pre_cursor = i + 1;
    }
  }

  std::ofstream solved("solved.txt");

  solved << no_grids << std::endl;
  solved << "-------" << std::endl;

  for (int i = 0; i < no_grids; i++)
  {
    constraints[i] = new int[constraint_sizes[i] * 4];
  }

  std::getline(file, file_line);
  for (int i = 0; i < no_grids; i++)
  {
    std::getline(file, file_line);
    for (int j = 0; j < SIZE; j++)
    {
      std::getline(file, file_line);
      // std::cout << "i: " << i << " file_line: " << file_line << std::endl;
      std::istringstream iss(file_line);
      iss >> elem0 >> elem1 >> elem2 >> elem3 >> elem4;
      grids[i][SIZE * j + 0] = elem0;
      grids[i][SIZE * j + 1] = elem1;
      grids[i][SIZE * j + 2] = elem2;
      grids[i][SIZE * j + 3] = elem3;
      grids[i][SIZE * j + 4] = elem4;
    }
    for (int c = 0; c < constraint_sizes[i]; c++)
    {
      std::getline(file, file_line);
      // std::cout << "i: " << i << "c line: " << file_line << std::endl;
      std::istringstream iss(file_line);
      iss >> elem0 >> elem1 >> elem2 >> elem3;
      constraints[i][4 * c] = elem0;
      constraints[i][4 * c + 1] = elem1;
      constraints[i][4 * c + 2] = elem2;
      constraints[i][4 * c + 3] = elem3;
    }
  }

  long int *constraintsSums = new long int[no_grids + 1];
  constraintsSums[0] = 0;
  for (int i = 1; i < no_grids + 1; i++)
  {
    constraintsSums[i] += constraintsSums[i - 1] + constraint_sizes[i - 1];
  }

  int *constraintsArray = new int[4 * constraintsSums[no_grids]];
  int *gridsArray = new int[SIZE * SIZE * no_grids];

  for (int y = 0, k = 0, j = 0; y < no_grids; y++)
  {
    for (int z = 0; z < constraint_sizes[y]; z++)
    {
      constraintsArray[k++] = constraints[y][4 * z];
      constraintsArray[k++] = constraints[y][4 * z + 1];
      constraintsArray[k++] = constraints[y][4 * z + 2];
      constraintsArray[k++] = constraints[y][4 * z + 3];
    }

    for (int i = 0; i < SIZE * SIZE; i++)
    {
      gridsArray[j++] = grids[y][i];
    }
  }

  // ------------------------------- Reading Data From txt Done -----------------------------------------

  // ------------------------------- GPU Memcpy and Kernel Call -----------------------------------------

  int *matrix_in, *constraints_in, *numberOfInput;
  long int *constraintSums_in;
  size_t Smatrix_in = SIZE * SIZE * no_grids * sizeof(int);
  size_t Sconstraints_in = 4 * constraintsSums[no_grids] * sizeof(int);
  size_t SconstraintSums_in = (no_grids + 1) * sizeof(long int);

  cudaEventRecord(start, 0);

  checkCuda(cudaMalloc((void **)&matrix_in, Smatrix_in));
  checkCuda(cudaMalloc((void **)&constraints_in, Sconstraints_in));
  checkCuda(cudaMalloc((void **)&constraintSums_in, SconstraintSums_in));
  checkCuda(cudaMalloc((void **)&numberOfInput, sizeof(int)));

  checkCuda(cudaMemcpy(matrix_in, gridsArray, Smatrix_in, cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(constraints_in, constraintsArray, Sconstraints_in, cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(constraintSums_in, constraintsSums, SconstraintSums_in, cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(numberOfInput, &no_grids, sizeof(int), cudaMemcpyHostToDevice));

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf(" Host2Device %lf sec ", elapsedTime / 1000);

  cudaEventRecord(start, 0);
  solver<<<BLOCKSIZE, THREADPERBLOCK>>>(matrix_in, constraints_in, constraintSums_in, numberOfInput);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf(" Kernel %lf sec ", elapsedTime / 1000);
  cudaDeviceSynchronize();
  checkCuda(cudaPeekAtLastError());

  cudaEventRecord(start, 0);
  checkCuda(cudaMemcpy(gridsArray, matrix_in, Smatrix_in, cudaMemcpyDeviceToHost));

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf(" Device2Host %lf sec ", elapsedTime / 1000);
  // ------------------------------- **************************************************** -----------------------------------------

  // ------------------------------- Correctness Check and Solved Matrices txt Write Start-----------------------------------------
  int *matrix = new int[SIZE * SIZE];
  bool isAllCorrect = true;
  for (int index = 0; index < no_grids; index++)
  {
    for (int j = 0; j < 5; j++)
    {

      matrix[5 * j] = gridsArray[25 * index + 5 * j];
      matrix[5 * j + 1] = gridsArray[25 * index + 5 * j + 1];
      matrix[5 * j + 2] = gridsArray[25 * index + 5 * j + 2];
      matrix[5 * j + 3] = gridsArray[25 * index + 5 * j + 3];
      matrix[5 * j + 4] = gridsArray[25 * index + 5 * j + 4];

      solved << matrix[5 * j] << " " << matrix[5 * j + 1] << " " << matrix[5 * j + 2]
             << " " << matrix[5 * j + 3] << " " << matrix[5 * j + 4] << std::endl;
    }

    solved << "-------" << std::endl;

    if (!isSolved(matrix, constraints[index], constraint_sizes[index]))
    {
      printf("Puzzle %i is WRONG!!! \n", index);
      isAllCorrect = false;
    }
  }

  if (isAllCorrect)
  {
    printf(" \n ***ALL DONE !!! ALL CORRECT*** \n");
    printf("Output file name is solved.txt \n");
  }
  // ------------------------------- Correctness Check and Solved Matrices txt Write Done-----------------------------------------

  checkCuda(cudaFree(matrix_in));
  checkCuda(cudaFree(constraints_in));
  checkCuda(cudaFree(constraintSums_in));

  return 0;
}