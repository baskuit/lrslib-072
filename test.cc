#include "src/lib.h"

#include "assert.h"
#include <algorithm>
#include <iostream>
#include <random>

template <typename T, typename U>
void print_output(const T &input, const U &output) {
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

template <typename T, typename U>
void check_output(const T &input, const U &output) {

  const static float eps = .0001;

  std::vector<float> row_scores{};
  std::vector<float> col_scores{};
  row_scores.resize(input.rows);
  col_scores.resize(input.cols);
  float value = 0;
  int k = 0;
  for (int i = 0; i < input.rows; ++i) {
    for (int j = 0; j < input.cols; ++j) {
      int d = input.data[k];
      row_scores[i] += d * output.col_strategy[j];
      col_scores[j] += d * output.row_strategy[i];
      value += output.row_strategy[i] * output.col_strategy[j] * d;
      ++k;
    }
  }

  value /= input.den;

  float row_max = *std::max_element(row_scores.begin(), row_scores.end());
  float col_min = *std::min_element(col_scores.begin(), col_scores.end());

  float expl = row_max - col_min;

  // std::cout << "row_max: " << row_max << std::endl;
  // std::cout << "col_min: " << col_min << std::endl;

  // std::cout << "expl: " << expl << std::endl;
  // std::cout << "calulated value: " << value << std::endl;
  // std::cout << "output value: " << output.value << std::endl;

  assert(std::abs(value - output.value) < eps);
  assert(std::abs(expl) < eps);
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
    // check_output(input, output);
    // std::cout << std::endl;
  }

  return 0;
}
