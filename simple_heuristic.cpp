/* 
 * File:   simple_heuristic.h
 * Author: Shunsuke Iwamoto
 *
 * Created on October 30, 2015, 9:38 PM
 */

#include <iostream>
#include <ostream>

#include <math.h>

#include "simple_heuristic.h"
#include "best_score_calculator.h"

heuristics::SimpleHeuristic::SimpleHeuristic(int variableCount){
  this->variableCount = variableCount;
}

heuristics::SimpleHeuristic::~SimpleHeuristic() {
    
}

void heuristics::SimpleHeuristic::initialize(std::vector<bestscorecalculators::BestScoreCalculator*> &spgs){
  this->spgs = spgs;
}

float heuristics::SimpleHeuristic::h(const varset &variables, bool &complete){
  float h = 0;

  VARSET_NEW(mask, variableCount);
  VARSET_SET_ALL(mask, variableCount);
  VARSET_COPY(VARSET_NOT(variables), remaining);

  // use the mask to make sure "remaining" only has valid variables left
  remaining = VARSET_AND(remaining, mask);
  for(int i = 0; i < variableCount; i++){
    if (VARSET_GET(remaining, i)){
      h += spgs[i]->getScore(mask);
    }
  }

  return h;
}
