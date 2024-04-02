/**
 * @file      gol.cpp
 *
 * @author    Pavel Kratochvil \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xkrato61@vutbr.cz
 *
 * @brief     PC lab 2 / MPI Broadcast
 *
 * @version   2023
 *
 * @date      23 February  2020, 11:13 (created) \n
 *
 */

#include <thread>
#include <cstdlib>
#include <cstdio>

#include <string>
#include <fstream>

#include <iostream>

#include <mpi.h>

constexpr int TO_ABOVE = 0;
constexpr int TO_BELLOW = 1;
constexpr int FROM_ABOVE = 1;
constexpr int FROM_BELLOW = 0;

constexpr int TO_ABOVE_REQUEST = 0;
constexpr int TO_BELLOW_REQUEST = 1;
constexpr int FROM_ABOVE_REQUEST = 2;
constexpr int FROM_BELLOW_REQUEST = 3;

constexpr int COMM_OVERLAY_THRESHOLD = 5;

#define ABOVE_RANK (rank - 1)
#define BELLOW_RANK (rank + 1)

#define MPI_ROOT_RANK 0

#define DBG 0

int mpiGetCommRank(MPI_Comm comm) {
  int rank{};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int mpiGetCommSize(MPI_Comm comm) {
  int size{};
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

int sum_neighbours(std::vector<int>& rs, int idx, int my_counts, int row_length) {
  int sum = 0;
  int m = idx % row_length;
#if DBG
  if (mpiGetCommRank(MPI_COMM_WORLD) == 0) {
    if (idx == 4) {
      // look above
      if (rs[idx - row_length]) printf("above\n");
      // look bellow
      if (rs[idx + row_length]) printf("bellow\n");
      // look left
      if (m > 0 && rs[idx - 1]) printf("left\n");
      // look above left
      if (m > 0 && rs[idx - row_length - 1]) printf("above left\n");
      // look above right
      if ((m != row_length - 1) && rs[idx - row_length + 1]) printf("above right\n");
      // look right
      if ((m != row_length - 1) && rs[idx + 1]) printf("right\n");
      // look bellow left
      if (m > 0 && rs[idx + row_length - 1]) printf("bellow left\n");
      // look bellow right
      if ((m != row_length - 1) && rs[idx + row_length + 1]) printf("bellow right\n");
    }
  }
#endif
  // N - look above
  if (rs[idx - row_length]) sum++;
  // S - look bellow
  if (rs[idx + row_length]) sum++;
  // W - look left
  if (m > 0 && rs[idx - 1]) sum++;
  // NW - look above left
  if (m > 0 && rs[idx - row_length - 1]) sum++;
  // NE - look above right
  if ((m != row_length - 1) && rs[idx - row_length + 1]) sum++;
  // E - look right
  if ((m != row_length - 1) && rs[idx + 1]) sum++;
  // SW - look bellow left
  if (m > 0 && rs[idx + row_length - 1]) sum++;
  // SE - look bellow right
  if ((m != row_length - 1) && rs[idx + row_length + 1]) sum++;
  
  return sum;
}

void gol_iteration_inner(std::vector<int>& rs, std::vector<int>& ns, int my_counts, int row_length, int start, int end) {
  // rs == localRows
  // count == how many rows do I have
  // row_length == row length
  // the rs contains at least 5 full rows, the top and bottom are being exchanged
  // we can calculated only layers with a top and bottom padding of at least 2
  // layer 0    - exchanged
  // layer 1    - can be influenced by layer 0 in the next step
  // layer 2..n - can be calculated
  // layer n+1  - can be influenced by layer n+2 in the next step
  // layer n+2  - exchanged

  for (int i = start; i < end; i++) {
    ns[i] = sum_neighbours(rs, i, my_counts, row_length);
  }
}

void gol_next_state(std::vector<int>& rs, std::vector<int>& ns, int my_counts, int row_length) {
  for (int i = row_length; i < my_counts + row_length; i++) {
    // cel is alive
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

void print_field(std::vector<int>& field, int row_length) {
  for (int i = 0; i < field.size(); i++) {
    printf("%d", field[i]);
    if (i % row_length == row_length - 1) printf("\n");
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  
  int row_length, row_count, iterations;
  std::vector<int> field;

  int size = mpiGetCommSize(MPI_COMM_WORLD);
  int rank = mpiGetCommRank(MPI_COMM_WORLD);
  
  if (rank == MPI_ROOT_RANK) {
    if (argc == 3) {
      std::string file_path = argv[1];
      if (read_file(file_path, field, &row_count, &row_length)) {
          MPI_Abort(MPI_COMM_WORLD, -1);
      }
      try {
        iterations = std::stoi(argv[2]);
      } catch (std::invalid_argument& e) {
          std::cerr << "Invalid argument: " << argv[2] << " is not a number." << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
      } catch (std::out_of_range& e) {
          std::cerr << "Invalid argument: " << argv[2] << " is out of range for an integer." << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
      }
    } else {
        std::cerr << "Invalid number of arguments" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  MPI_Bcast(&iterations, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);
  MPI_Bcast(&row_length, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);
  MPI_Bcast(&row_count, 1, MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);

  int row_size = row_count / size;
  int row_rema = row_count % size;

#if DBG
  if (rank == 0) printf("row_size: %d row_count: %d row_rema: %d size: %d\n", row_size, row_count, row_rema, size);
#endif

  std::vector<int> counts(size), displacements(size);
  std::fill_n(displacements.data(), size, 0);

  // compute the counts and displacements for MPI_Scatterv and MPI_Gatherv
  for (int i = 0; i < size; i++) {
    counts[i] = row_size * row_length;
    if (i < row_rema) counts[i] += row_length;
  }
  for (int i = 1; i < size; i++) {
    displacements[i] = displacements[i-1] + counts[i-1];
  }
#if DBG
  if (rank == 0) {
    for (int i = 0; i < size; i++) {
      printf("%d ", counts[i]);
    }
    printf("\n");

    for (int i = 0; i < size; i++) {
      printf("%d ", displacements[i]);
    }
    printf("\n");
  }
#endif

  std::vector<int> local_rows(counts[rank] + (2 * row_length)), next_state(counts[rank] + (2 * row_length));

  MPI_Scatterv(field.data(), counts.data(), displacements.data(), MPI_INT, local_rows.data() + row_length, counts[rank], MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);
  
  // saving requests for async sending and then waiting
  auto requests = std::make_unique<MPI_Request[]>(4);

  for (int i = 0; i < iterations; i++) {
    // send up
    if (rank > MPI_ROOT_RANK) {
      MPI_Isend(local_rows.data() + row_length, row_length, MPI_INT, ABOVE_RANK, TO_ABOVE, MPI_COMM_WORLD, &requests[TO_ABOVE_REQUEST]);
    } else {
      requests[TO_ABOVE_REQUEST] = MPI_REQUEST_NULL;
    }

    // send down
    if (rank < size - 1) {
      MPI_Isend(local_rows.data() + counts[rank], row_length, MPI_INT, BELLOW_RANK, TO_BELLOW, MPI_COMM_WORLD, &requests[TO_BELLOW_REQUEST]);
    } else {
      requests[TO_BELLOW_REQUEST] = MPI_REQUEST_NULL;
    }

    // calculate the centers if thickness is >= 5
    if (counts[rank] >= COMM_OVERLAY_THRESHOLD * row_length) {
      int start = row_length * 2;
      int end = counts[rank];
#if DBG
      if (rank == 0) printf("#0 start: %d  end: %d\n", start, end);
#endif
      gol_iteration_inner(local_rows, next_state, counts[rank], row_length, start, end);
    }

    // receive from above me
    if (rank > 0) {
      MPI_Irecv(local_rows.data(), row_length, MPI_INT, ABOVE_RANK, FROM_ABOVE, MPI_COMM_WORLD, &requests[FROM_ABOVE_REQUEST]);
    } else {
      requests[FROM_ABOVE_REQUEST] = MPI_REQUEST_NULL;
    }

    // receive from bellow me
    if (rank < size - 1) {
      MPI_Irecv(local_rows.data() + counts[rank] + row_length, row_length, MPI_INT, BELLOW_RANK, FROM_BELLOW, MPI_COMM_WORLD, &requests[FROM_BELLOW_REQUEST]);
    } else {
      requests[FROM_BELLOW_REQUEST] = MPI_REQUEST_NULL;
    }

  MPI_Waitall(4, requests.get(), MPI_STATUSES_IGNORE);

#if DBG
  if (rank == 0) {
    printf("[%d] local rows:\n", rank);
    for (int n = 0; n < counts[rank] + (2 * row_length); n++) {
      printf("(%d, %d)", local_rows[n], next_state[n]);
      if (n % row_length == row_length - 1) printf("\n");
    }
    printf("\n");
  }
#endif

  int start, end;
  if (counts[rank] >= COMM_OVERLAY_THRESHOLD * row_length) {
    start = row_length;
    end = row_length * 2;
#if DBG
    if (rank == 0) printf("#1 start: %d  end: %d\n", start, end);
#endif
    gol_iteration_inner(local_rows, next_state, counts[rank], row_length, start, end);

    start = counts[rank];
    end = counts[rank] + row_length;
#if DBG
    if (rank == 0) printf("#2 start: %d  end: %d\n", start, end);
#endif
    gol_iteration_inner(local_rows, next_state, counts[rank], row_length, start, end);
  } else {
    start = row_length;
    end = counts[rank] + row_length;
#if DBG
    if (rank == 0) printf("start: %d  end: %d\n", start, end);
#endif
    gol_iteration_inner(local_rows, next_state, counts[rank], row_length, start, end);
  }

#if DBG
  if (rank == 0) {
    printf("[%d] next state:\n", rank);
    for (int n = 0; n < counts[rank] + (2 * row_length); n++) {
      printf("(%d, %d)", local_rows[n], next_state[n]);
      if (n % row_length == row_length - 1) printf("\n");
    }
    printf("\n");
  }
#endif

  gol_next_state(local_rows, next_state, counts[rank], row_length);

  // if (rank == 0) {
  //   printf("[%d] next state:\n", rank);
  //   for (int n = 0; n < counts[rank] + (2 * row_length); n++) {
  //     printf("(%d, %d)", local_rows[n], next_state[n]);
  //     if (n % row_length == row_length - 1) printf("\n");
  //   }
  //   printf("\n");
  // }

  //   if (rank == 0) {
  //   printf("[%d] next state local rows:\n", rank);
  //   for (int n=0; n < counts[rank] + (2 * row_length); n++) {
  //     printf("%d ", local_rows[n]);
  //     if (n % row_length == row_length - 1) printf("\n");
  //   }
  //   printf("\n");
  // }
  }

  MPI_Gatherv(local_rows.data() + row_length, counts[rank], MPI_INT, field.data(), counts.data(), displacements.data(), MPI_INT, MPI_ROOT_RANK, MPI_COMM_WORLD);

  if (rank == MPI_ROOT_RANK) {
    print_field(field, row_length);
  }

  MPI_Finalize();
}
