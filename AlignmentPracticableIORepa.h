#ifndef ALIGNMENTPRACTICABLEIOREPA_H
#define ALIGNMENTPRACTICABLEIOREPA_H

#include "AlignmentPracticableRepa.h"

namespace Alignment
{
    // parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u ::
    //   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
    //   [VariableRepa] -> HistoryRepa -> HistoryRepa-> Integer ->
    //   IO (SystemRepa, FudRepa, [(Double, [VariableRepa]])
    std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, const SizeList&, const HistoryRepa&, const HistoryRepa&, std::size_t, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa ::
    //   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
    //   Integer -> Integer ->
    //   [VariableRepa] -> HistoryRepa ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, const SizeList&, const HistoryRepa&, SystemRepa&);

}



#endif