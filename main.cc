#include "lib.h"

#include <iostream>
#include <random>

template <typename T, typename U>
void print_output(const T& input, const U& output) {
  std::cout << "row strategy" << std::endl;
  for (int i = 0; i < input.rows; ++i) {
    std::cout << output.row_strategy[i] << ", ";
  }
  std::cout << '\n';
  std::cout << "col strategy" << std::endl;
  for (int i = 0; i < input.cols; ++i) {
    std::cout << output.col_strategy[i] << ", ";
  }
  std::cout << '\n';
}

int main(int argc, char **argv) {

  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " rows cols den trials\n";
    return 1; // indicate error
  }
  FastInput input;
  input.rows = std::atoi(argv[1]);
  input.cols = std::atoi(argv[2]);
  input.den = std::atoi(argv[3]);
  const int trials = std::atoi(argv[4]);
  const int entries = input.rows * input.cols;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(0, input.den + 1);


  for (int t = 0; t < trials; ++t) {

    std::vector<int> data{};
    data.resize(entries);
    input.data = data.data();

    for (int i = 0; i < entries; ++i) {
      data[i] = dis(gen);
    }

    FloatOneSumOutput output{};
    std::vector<float> row_strategy{};
    std::vector<float> col_strategy{};
    row_strategy.resize(input.rows + 2);
    col_strategy.resize(input.cols + 2);
    output.row_strategy = row_strategy.data();
    output.col_strategy = col_strategy.data();
    
    solve_fast(&input, &output);

    // print_output(input, output);
  }

  return 0;
}
