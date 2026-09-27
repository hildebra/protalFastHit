
/*******************************************************************************
 * Copyright 2019 AMADEUS. All rights reserved.
 * Author: Paolo Iannino
 *******************************************************************************/

#include "cPMML.h"
#include "core/internal_evaluator.h"
#include "core/internal_score.h"
#include "core/modelbuilder.h"
#include "utils/csvreader.h"
#include "utils/utils.h"

namespace cpmml {
Model::Model(const std::string &model_filepath) : evaluator(ModelBuilder::build(model_filepath, false)) {}

Model::Model(const std::string &model_filepath, const bool zipped = false)
    : evaluator(ModelBuilder::build(model_filepath, zipped)) {}

Model Model::from_string(const std::string &pmml) {
  std::vector<char> data(pmml.begin(), pmml.end());
  data.push_back('\0');
  Model model;
  model.evaluator = ModelBuilder::build(std::move(data));
  return model;
}

bool Model::validate(const std::unordered_map<std::string, std::string> &sample) const {
  return evaluator->validate(sample);
}

Prediction Model::score(const std::unordered_map<std::string, std::string> &sample) const {
  return Prediction(evaluator->score(sample));
}

std::string Model::predict(const std::unordered_map<std::string, std::string> &sample) const {
  return evaluator->predict(sample);
}
}  // namespace cpmml
