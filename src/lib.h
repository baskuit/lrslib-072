struct lrs_dic;
struct lrs_dat;

struct ratnum {
  long num;
  long den;
};

struct FastInput {
  long rows;
  long cols;
  int *data;
  int den;
};

struct FloatOneSumOutput {
  float *row_strategy;
  float *col_strategy;
  float value;
};

void solve_fast(const FastInput *g, FloatOneSumOutput *gg);