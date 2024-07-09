#include "src/lib.h"

#include "assert.h"
#include <algorithm>
#include <iostream>
#include <random>

int main(int argc, char **argv) {

  FastInput input;
  input.rows = 9;
  input.cols = 9;
  input.den = 80;
  const int trials = 100000;
  const int entries = input.rows * input.cols;

//   std::random_device rd{};
  std::mt19937 gen(100000001);
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
  }

  return 0;
}
