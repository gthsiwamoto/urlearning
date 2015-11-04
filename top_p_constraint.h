/* 
 * File:   top_p_constraint.h
 * Author: bmmalone
 *
 * Created on April 21, 2014, 6:20 PM
 */

#ifndef TOP_P_CONSTRAINT_H
#define	TOP_P_CONSTRAINT_H

#include <vector>

#include "typedefs.h"

namespace bestscorecalculators {
    class BestScoreCalculator;
}

namespace heuristics {
    
    class TopPConstraint {
    
    public:
        
    std::vector<bitset> findParentSets(std::vector<bestscorecalculators::BestScoreCalculator*> &spgs, int maxScc);
    
    };
}


#endif	/* TOP_P_CONSTRAINT_H */

