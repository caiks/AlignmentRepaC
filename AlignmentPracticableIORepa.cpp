#include "AlignmentPracticableIORepa.h"

using namespace Alignment;

// parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa-> HistogramRepaRed -> Integer ->
//   IO(SystemRepa, FudRepa, [(Double, [VariableRepa]])
std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> Alignment::parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, const SizeList& vv, const HistoryRepa& xx, const HistogramRepaRed& xxp, const HistoryRepa& xxrr, const HistogramRepaRed& xxrrp, std::size_t f, SystemRepa& ur)
{
    auto fr = std::make_unique<FudRepa>();
    auto xx1 = std::make_unique<DoubleSizeListPairList>();

    return std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>>(std::move(fr), std::move(xx1));
}



