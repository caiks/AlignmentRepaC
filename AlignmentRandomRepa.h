#ifndef ALIGNMENTRANDOMREPA_H
#define ALIGNMENTRANDOMREPA_H

#include "AlignmentRepa.h"

namespace Alignment
{
    // historyRepasShuffle_u :: HistoryRepa -> Int -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasShuffle_u(const HistoryRepa&, unsigned);
}



#endif