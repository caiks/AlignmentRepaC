#ifndef ALIGNMENTPRACTICABLEIOREPA_H
#define ALIGNMENTPRACTICABLEIOREPA_H

#include "AlignmentRepa.h"

namespace Alignment
{
    // parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u ::
    //   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
    //   [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa-> HistogramRepaRed -> Integer ->
    //   IO(SystemRepa, FudRepa, [(Double, [VariableRepa]])
    std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, const SizeList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&, std::size_t, SystemRepa&);

}



#endif