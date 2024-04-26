/**
 * @file    life.cpp
 *
 * @author Pavel Kratochvil
 *         xkrato61@vutbr.cz
 * 
 * @date   26 April  2024
 * 
 * @brief  The program reads the file passed as the first argument and performs a given number of iterations
 *         passed as the second argument. The number of CPUs is decided by the `test.sh` script to be the
 *         min(PHYSICAL_CORES, ROW_COUNT). This way, all the cores in the system are utilized so there
 *         is no need for oversubscription and software threads. The program works with arbitrary number of
 *         spawned process as long as ROW_COUNT >= NP. The program also works with arbitrarily sized input.
 *        
 *         Bonus feature: I tried to make the assignment to be as efficient as possible. Therefore, I implemented
 *         overlaying of communication and computation. Firstly, the data is distributed by blocks of rows.
 *         The overlay can be seen in the main function, where I first asynchronously send the data on the
 *         borders, then use the otherwise wasted time in communication for center data calculation and then
 *         wait for the communication to finish. Only then I calculate the next state for the parts which are
 *         dependent on the neighboring ranks. All of this is only possible iff the thickness of the block of
 *         rows is greater then 5, otherwise, there can be no overlay and the computation has to made
 *         synchronously in iterations by calculating all of the data in one rank in one go.
 */

#include <thread>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include <mpi.h>

// message tags for P2P communication
constexpr int TO_ABOVE = 0;
constexpr int TO_BELLOW = 1;
constexpr int FROM_ABOVE = 1;
constexpr int FROM_BELLOW = 0;

constexpr int TO_ABOVE_REQUEST = 0;
constexpr int TO_BELLOW_REQUEST = 1;
constexpr int FROM_ABOVE_REQUEST = 2;
constexpr int FROM_BELLOW_REQUEST = 3;

// threshold for utilizing communication and computation overlay
constexpr int COMM_OVERLAY_THRESHOLD = 5;

// macros for neighbor rank calculation (accounting for wrap-around)
#define ABOVE_RANK ((rank - 1 + size) % size)
#define BELLOW_RANK ((rank + 1) % size)

#define MPI_ROOT_RANK 0

/**
 * Calculates the sum of the neighbors of a given cell in a 2D grid.
 * The grid is represented as a 1D vector for efficiency. The function takes into 
 * account the boundary conditions and wraps around the x coordinate if necessary.
 * 
 * @param rs: The 1D vector representing the 2D grid.
 * @param idx: The index of the cell in the 1D vector.
 * @param my_counts: Cell count of the rank.
 * @param row_length: The number of cells in a row of the 2D grid.
 * @return The sum of the neighbors of the cell.
 */
int sum_neighbours(std::vector<int>& rs, int idx, int my_counts, int row_length) {
  int sum = 0;
  int x = idx % row_length;
  int y = idx / row_length;
  int nx, ny;
  int field_length = (my_counts / row_length) + 2;

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      // skip the current cell
      if (!dx && !dy) continue;
      nx = x + dx, ny = y + dy;

      // wrap around the x coordinate, y coordinate does not need wrap is there is always the neighbors' line
      if (nx > row_length - 1 || nx < 0) {
        nx = ((nx % row_length) + row_length) % row_length;
      }
      sum += rs[ny * row_length + nx];
    }
  }
  return sum;
}

/**
 * Performs one iteration of the game. The function calculates the sum of the neighbors
 * of a given range of cells.
 * 
 * 
 * @param rs: The 1D vector representing the current state of the 2D grid.
 * @param ns: The 1D vector representing the neighbor sum of the 2D grid.
 * @param my_counts: Cell count of the rank.
 * @param row_length: The number of cells in a row of the 2D grid.
 * @param start: The index of the first cell to be updated.
 * @param end: The index of the last cell to be updated.
 */
void gol_iteration(std::vector<int>& rs, std::vector<int>& ns, int my_counts, int row_length, int start, int end) {
  for (int i = start; i < end; i++) {
    ns[i] = sum_neighbours(rs, i, my_counts, row_length);
  }
}

/**
 * Determines the next state of the cell given the sum of the neighbors in `ns`.
 * 
 * each live cell with less then 2 neighbors dies
 * each live cell with 2 or 3  neighbors stays alive
 * each live cell with more than 3 neighbors dies
 * each dead with exactly 3 neighbors alives
 * 
 * @param rs: The 1D vector representing the current state of the 2D grid.
 * @param ns: The 1D vector representing the neighbor sum of the 2D grid.
 * @param my_counts: Cell count of the rank.
 * @param row_length: The number of cells in a row of the 2D grid.
 */
void gol_next_state(std::vector<int>& rs, std::vector<int>& ns, int my_counts, int row_length) {
  for (int i = row_length; i < my_counts + row_length; i++) {
    // cell is alive
    if (rs[i]) {
      if (ns[i] < 2) {
        rs[i] = 0;
      } else if (ns[i] > 3) {
        rs[i] = 0;
      }
    } 
    // cell is dead
    else {
      if (ns[i] == 3) rs[i] = 1;
    }
  }
}

/**
 * Reads a file containing a 2D grid for the Game of Life.
 * The grid is represented as a sequence of '0's and '1's, where '0' represents a dead cell
 * and '1' represents a live cell. Each row of the grid is terminated by a newline character ('\n').
 * 
 * The function reads the file character by character and pushes the cell states into a 1D vector.
 * It also calculates the number of rows and the length of each row in the grid.
 * 
 * If the file cannot be opened, or if it contains any character other than '0', '1', or '\n',
 * or if the rows have different lengths, the function returns -1. Otherwise, it returns 0.
 * 
 * @param file_name: The name of the file to be read.
 * @param field: The 1D vector to store the cell states.
 * @param row_count: A pointer to an integer to store the number of rows in the grid.
 * @param row_length: A pointer to an integer to store the length of each row in the grid.
 * @return 0 if the file is successfully read, -1 otherwise.
 */
int read_file(std::string file_name, std::vector<int>& field, int* row_count, int* row_length) {
    std::ifstream file(file_name);
    if (!file) return -1;

    bool first_length = false;
    int length = 0;
    char value;
    while (file.get(value)) {
      switch (value) {
        case '0':
          field.push_back(0);
          length++;
          break;
        case '1':
          field.push_back(1);
          length++;
          break;
        case '\n':
          if (!first_length) {
            *row_length = length;
            first_length = true;
          }
          // check for various length lines
          if (length != *row_length) return -1;
          length = 0;
          break;
        default:
          // any other character is faulty
          return -1;
      }
    }
    file.close();
    // Assuming the file is a grid of integers, calculate row_count and row_length
    *row_count = field.size() / *row_length;
    return 0;
}

/**
 * This function prints a 2D grid for the Game of Life to the console, along with the rank of each row.
 * Each row of the grid is terminated by a newline character ('\n'). The rank of each row is printed
 * at the beginning of the row.
 * 
 * @param field: The 1D vector representing the 2D grid.
 * @param row_length: The number of cells in a row of the 2D grid.
 * @param counts: The vector containing the number of rows assigned to each rank.
 */
void print_field_ranks(std::vector<int>& field, int row_length, std::vector<int>& counts) {
  int currrent_rank = 0;
  std::vector<int> ranks;
  for (int i = 0; i < counts.size(); i++) {
    for (int r = 0; r < counts[i]; r += row_length) {
      ranks.push_back(currrent_rank);
    }
    currrent_rank++;
  }

  for (int i = 0; i < field.size(); i++) {
    if (i % row_length == 0) printf("%d: ", ranks[i / row_length]);
    printf("%d", field[i]);
    if (i % row_length == row_length - 1) printf("\n");
  }
}

/**
 * Parses the command-line arguments.
 * The program expects two arguments: the path to a file containing the initial state of the game,
 * and the number of iterations to run the game.
 * 
 * The function reads the initial state from the file and stores it in a 1D vector. It also calculates
 * the number of rows and the length of each row in the grid. The number of iterations is converted to an integer.
 * 
 * If the file cannot be read, or if the second argument is not a valid integer, or if the number of arguments
 * is not 2, the function returns -1. Otherwise, it returns 0.
 * 
 * @param argc: The number of command-line arguments.
 * @param argv: The array of command-line arguments.
 * @param field: The 1D vector to store the initial state of the game.
 * @param row_count: The number of rows in the grid.
 * @param row_length: The length of each row in the grid.
 * @param iterations: The umber of iterations.
 * @return 0 if the arguments are successfully parsed, -1 otherwise.
 */
int parse_and_read(int &argc, char **argv, std::vector<int>& field, int &row_count, int &row_length, int &iterations) {
  if (argc == 3) {
      std::string file_path = argv[1];
      if (read_file(file_path, field, &row_count, &row_length)) {
          return -1;
      }
      try {
        iterations = std::stoi(argv[2]);
      } catch (std::invalid_argument& e) {
          std::cerr << "Invalid argument: " << argv[2] << " is not a number." << std::endl;
          return -1;
      } catch (std::out_of_range& e) {
          std::cerr << "Invalid argument: " << argv[2] << " is out of range for an integer." << std::endl;
          return -1;
      }
    } else {
        std::cerr << "Invalid number of arguments" << std::endl;
        return -1;
    }
    return 0;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int row_length, row_count, iterations, rank, size, start, end;
  std::vector<int> field;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // root rank parses args and reads the input file
  if (rank == MPI_ROOT_RANK) {
    if (parse_and_read(argc, argv, field, row_count, row_length, iterations)) {
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // broadcast necessary information to all ranks
  MPI_Bcast(&iterations, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);
  MPI_Bcast(&row_length, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);
  MPI_Bcast(&row_count, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);

  std::vector<int> counts(size), displacements(size);
  displacements[0] = 0;

  // compute the counts and displacements for MPI_Scatterv and MPI_Gatherv
  int row_size = row_count / size, row_rema = row_count % size;
  for (int i = 0; i < size; i++) {
    counts[i] = row_size * row_length;
    if (i < row_rema) counts[i] += row_length;
  }
  for (int i = 1; i < size; i++) {
    displacements[i] = displacements[i-1] + counts[i-1];
  }
  
  // allocate space for local rows with space for neighbors' rows
  std::vector<int> local_rows(counts[rank] + (2 * row_length)), next_state(counts[rank] + (2 * row_length));
  MPI_Scatterv(field.data(), counts.data(), displacements.data(), MPI_INT, local_rows.data() + row_length, counts[rank], MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);

  // saving requests for async sending and then waiting
  auto requests = std::make_unique<MPI_Request[]>(4);

  // main loop
  for (int i = 0; i < iterations; i++) {
    // send up
    MPI_Isend(local_rows.data() + row_length, row_length, MPI_INT, ABOVE_RANK, TO_ABOVE, MPI_COMM_WORLD, &requests[TO_ABOVE_REQUEST]);
    // send down
    MPI_Isend(local_rows.data() + counts[rank], row_length, MPI_INT, BELLOW_RANK, TO_BELLOW, MPI_COMM_WORLD, &requests[TO_BELLOW_REQUEST]);

    // receive from above me
    MPI_Irecv(local_rows.data(), row_length, MPI_INT, ABOVE_RANK, FROM_ABOVE, MPI_COMM_WORLD, &requests[FROM_ABOVE_REQUEST]);
    // receive from bellow me
    MPI_Irecv(local_rows.data() + counts[rank] + row_length, row_length, MPI_INT, BELLOW_RANK, FROM_BELLOW, MPI_COMM_WORLD, &requests[FROM_BELLOW_REQUEST]);

    // calculate the centers if thickness is >= 5
    if (counts[rank] >= COMM_OVERLAY_THRESHOLD * row_length) {
      start = row_length * 2;
      end = counts[rank];
      gol_iteration(local_rows, next_state, counts[rank], row_length, start, end);
    }
  
    MPI_Waitall(4, requests.get(), MPI_STATUSES_IGNORE);
    
    if (counts[rank] >= COMM_OVERLAY_THRESHOLD * row_length) {
      // if the tile is thick enough, there are two remaining areas
      start = row_length;
      end = row_length * 2;
      gol_iteration(local_rows, next_state, counts[rank], row_length, start, end);

      start = counts[rank];
      end = counts[rank] + row_length;
      gol_iteration(local_rows, next_state, counts[rank], row_length, start, end);
    } else {
      // if the tile is too small, do the entire area
      start = row_length;
      end = counts[rank] + row_length;
      gol_iteration(local_rows, next_state, counts[rank], row_length, start, end);
    }

    // calculate the next state
    gol_next_state(local_rows, next_state, counts[rank], row_length);
  }

  // at the end gather all the data to root rank
  MPI_Gatherv(local_rows.data() + row_length, counts[rank], MPI_INT, field.data(), counts.data(), displacements.data(), MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);

  // only root rank prints the result
  if (rank == MPI_ROOT_RANK) {
    print_field_ranks(field, row_length, counts);
  }
  MPI_Finalize();
}
