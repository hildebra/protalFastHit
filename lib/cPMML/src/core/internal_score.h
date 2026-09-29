
/*******************************************************************************
 * Copyright 2019 AMADEUS. All rights reserved.
 * Author: Paolo Iannino
 *******************************************************************************/

#ifndef CPMML_SCORE_H
#define CPMML_SCORE_H

#include <cerrno>
#include <cstdlib>
#include <string>
#include "core/value.h"
#include "utils/utils.h"

/**
 * @class InternalScore
 *
 * Abstract class internally representing the prediction produced by the scoring
 * process.
 *
 * It contains both double and literal representations of the score, as well as
 * the associated probabilities and all values produced by Output.
 */
class InternalScore {
 public:
  bool empty = true;
  std::string score;
  double double_score = double_min();
  std::unordered_map<std::string, double> probabilities;
  std::unordered_map<std::string, double> num_outputs;
  std::unordered_map<std::string, std::string> str_outputs;

  InternalScore() = default;

  explicit InternalScore(const double &score) : empty(false), score(std::to_string(score)), double_score(score) {}

  explicit InternalScore(const std::string &score) : empty(false), score(score), double_score(parse_score(score)) {}

  // ::to_double (std::stod) without exceptions, same values: internal tree nodes have no score and
  // leaves hold class labels ("TRUE"), so stod threw, and to_double threw again, for every node.
  static double parse_score(const std::string &score) {
    errno = 0;
    char *end = nullptr;
    double const value = std::strtod(score.c_str(), &end);
    if (end == score.c_str() || errno == ERANGE) return double_min();
    return value;
  }

  InternalScore(const double &score, const std::unordered_map<std::string, double> &probabilities)
      : empty(false), score(std::to_string(score)), double_score(score), probabilities(probabilities) {}

  InternalScore(const std::string &score, const std::unordered_map<std::string, double> &probabilities)
      : empty(false), score(score), double_score(parse_score(score)), probabilities(probabilities) {}

  InternalScore(const InternalScore &) = default;

  InternalScore(InternalScore &&) = default;

  InternalScore &operator=(const InternalScore &) = default;

  InternalScore &operator=(InternalScore &&) = default;

  virtual ~InternalScore() = default;
};

#endif
