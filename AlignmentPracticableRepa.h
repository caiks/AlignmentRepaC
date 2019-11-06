#ifndef ALIGNMENTPRACTICABLEREPA_H
#define ALIGNMENTPRACTICABLEREPA_H

#include "AlignmentRepa.h"

namespace Alignment
{
    // parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui ::
    //   Integer -> Integer -> Integer -> Integer -> [VariableRepa] -> Fud ->
    //   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
    //   ([[VariableRepa]],Integer)
    std::tuple<std::unique_ptr<SizeListList>, std::size_t> parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(std::size_t, std::size_t, std::size_t, std::size_t, const SizeList&, const FudRepa&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);

    // parametersSystemsPartitionerMaxRollByMRepa_ui ::
    //   Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistoryRepa -> 
    //   ([[[VariableRepa]]],Integer)
    std::tuple<std::unique_ptr<SizeListListList>, std::size_t> parametersSystemsPartitionerMaxRollByMRepa_ui(std::size_t, std::size_t, std::size_t, const SizeList&, const HistoryRepa&, const HistoryRepa&);


    // parametersRollerMaximumRollExcludedSelfRepa_i ::
    //   [[VariableRepa]] -> -> HistogramRepa -> -> HistogramRepa ->
    //   ([[[Int]]],Integer)
    std::tuple<std::unique_ptr<SizeListListList>, std::size_t> parametersRollerMaximumRollExcludedSelfRepa_i(const SizeListList&, const HistogramRepa&, const HistogramRepa&, double);

    // parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui ::
    //   Integer -> Integer -> [VariableRepa] -> Fud ->
    //   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
    //   ([(Double, [VariableRepa]],Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(std::size_t, std::size_t, const SizeList&, const FudRepa&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);


}



#endif