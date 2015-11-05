/* 
 * File:   simple_heuristic.h
 * Author: Shunsuke Iwamoto
 *
 * Created on October 30, 2015, 9:38 PM
 */

#ifndef SIMPLE_HEURISTIC_H
#define SIMPLE_HEURISTIC_H

#include <vector>

#include "typedefs.h"
#include "heuristic.h"

namespace heuristics {
class SimpleHeuristic : public Heuristic {
public:
    SimpleHeuristic() {}
    SimpleHeuristic(int variableCount);
    ~SimpleHeuristic();
    
    void initialize(std::vector<bestscorecalculators::BestScoreCalculator*> &spgs);

    int size() {};
    float h(const varset &variables, bool &complete);
    void print() {}
    void printStatistics() {
        printf("Simple heuristic does not collect statistics\n");
    }

private:
    int variableCount;
    std::vector<bestscorecalculators::BestScoreCalculator*> spgs;
    
    std::vector<varset> variableSets;
};

}
#endif	/* SIMPLE_HEURISTIC_H */

