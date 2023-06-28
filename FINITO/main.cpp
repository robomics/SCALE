#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <straw.h>
#include <string>
#include <sys/timeb.h>
#include <unistd.h>

#include "scale/scale.hpp"

static void usage(const char *argv0) {
  fprintf(stderr,
          "Usage: %s [-e (oe)][-D (diag)][-p perc][-v verbose][-t tol][-I "
          "max_iterations][-A All_iterations][-T threads] <hicfile> "
          "<chromosome> <resolution> <outfile>\n",
          argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chr>: chromosome\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

struct Config {
  std::string norm{"NONE"};

  std::string ob{"observed"};
  std::string norm_name{"SCALE"};
  double tol = 1.0e-4;
  double del = 5.0e-2;
  double perc = 0.01;
  int low;
  int maxiter = 1000;
  int All_iter = 2000;
  int threads = 1;
  int verb = 1;
  bool diag = false;
  int opt;
};

[[nodiscard]] Config parse_cli_args(int argc, char **argv) {
  Config c{};

  while ((c.opt = getopt(argc, argv, "ep:v:t:Dd:I:A:T:h")) != -1) {
    switch (c.opt) {
    case 'e':
      c.ob = "oe";
      break;
    case 'D':
      c.diag = true;
      break;
    case 'p':
      c.perc = std::stof(optarg);
      break;
    case 'd':
      c.del = std::stof(optarg);
      break;
    case 't':
      c.tol = std::stof(optarg);
      break;
    case 'v':
      c.verb = std::stoi(optarg);
      break;
    case 'A':
      c.All_iter = std::stoi(optarg);
      break;
    case 'I':
      c.maxiter = std::stoi(optarg);
      break;
    case 'T':
      c.threads = std::stoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }

    return c;
  }

  if (argc - optind < 4) {
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {

  const auto c = parse_cli_args(argc, argv);
  time_t t0, t1;

  std::string fname = argv[optind++];
  if (c.verb)
    printf("Reading hic file from %s\n", argv[optind - 1]);
  std::string chrom = argv[optind++];
  int binsize = atoi(argv[optind++]);
  char *out_name = argv[optind++];

  auto interactions = getSingleMatrix(fname, binsize, chrom);

  int width = -1;
  if (c.diag) {
    width = (int)(10000 / binsize);
    if (width < 2)
      width = 2;
  }

  auto &ii = interactions.ii;
  auto &jj = interactions.jj;
  auto &xx = interactions.xx;
  const auto& m = ii.size();

  for (std::size_t p = 0; p < m; p++)
    if (jj[p] == ii[p])
      xx[p] *= 0.5;

  std::vector<double> b((interactions.chrom1.length + binsize - 1) / binsize,
                        std::numeric_limits<double>::quiet_NaN());
  int totalIt = c.All_iter;
  std::vector<double> report(totalIt + 3 + 100);
  std::vector<int> allIters(totalIt + 3 + 100);

  int iter;
  time(&t0);

  int low{};
  double perc = c.perc;
  iter = balance(m, ii.data(), jj.data(), xx.data(), b.data(), report.data(), allIters.data(), c.tol, &low, c.maxiter, c.del,
                 &totalIt, c.threads, b.size(), &perc, width);


  double **space = (double **) malloc(c.threads*sizeof(double *));
  for (std::size_t p=0;p<c.threads; p++) space[p] = (double *) malloc(b.size()*sizeof(double));
  double sum = ppNormVector(m,ii.data(),jj.data(),xx.data(),b.data(), b.size(),c.threads, space);
  for (std::size_t p=0;p<c.threads; p++) free(space[p]);
  free(space);

  for (const auto&w: b) {
    printf("%f\n", w);
  }
  time(&t1);

    printf("iterations took %ld seconds\n", t1 - t0);
      printf("%d: %30.15lf\n", allIters.back(), report.back());
    printf("total %d iterations; final perc = %g; final cutoff = %d\n", totalIt,
           perc, low);
    printf("final error in scaling vector is %g and in row sums is %g\n",
           report[totalIt + 1], report[totalIt + 2]);
/*
  double **space = (double **)malloc(threads * sizeof(double *));
  for (p = 0; p < threads; p++)
    space[p] = (double *)malloc(k * sizeof(double));
  double sum = ppNormVector(m, ii, jj, xx, b, k, threads, space);

  fprintf(fout, "vector %s %s %d BP\n", const_cast<char *>(norm_name.c_str()),
          const_cast<char *>(chr.c_str()), binsize);
  for (int p = 0; p < k; p++) {
    if (!std::isnan(b[p]))
      fprintf(fout, "%20.10f\n", b[p]);
    else
      fprintf(fout, "%20s\n", "NaN");
  }
  fclose(fout);
  return (iter);
  */
}
