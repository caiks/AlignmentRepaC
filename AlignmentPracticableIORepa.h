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
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t seed, const SizeList&, const HistoryRepa&, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa ::
    //   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
    //   Integer -> Integer ->
    //   [VariableRepa] -> FudRepa -> HistoryRepa ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsFudRepasHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t seed, const SizeList&, const FudRepa&, const HistoryRepa&, SystemRepa&);


    // parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa ::
    //   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
    //   Integer -> Integer -> Integer ->
    //   [VariableRepa] -> FudRepa -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, double znnmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t smin, std::size_t seed, const SizeList&, const FudRepa&, const HistoryRepa&, int, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u ::
    //  Integer ->
    //   [VariableRepa] -> VariableRepa -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u(std::size_t fmax, const SizeList&, std::size_t, const HistoryRepa&, int, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u_1 ::
    //  Integer ->
    //   [VariableRepa] -> VariableRepa -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u_1(std::size_t fmax, const SizeList&, std::size_t, const HistoryRepa&, int, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_u ::
    //  Integer ->
    //   [VariableRepa] -> VariableRepa -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_u(std::size_t fmax, const SizeList&, std::size_t, const HistoryRepa&, int, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up ::
    //  Integer -> Integer ->
    //   [VariableRepa] -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up(std::size_t fmax, std::size_t tint, const SizeList&, std::size_t, const HistoryRepa&, int, SystemRepa&);

    // parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_1 ::
    //  Integer -> Integer ->
    //   [VariableRepa] -> HistoryRepa -> Integer ->
    //   IO (SystemRepa, ApplicationRepa)
    std::unique_ptr<ApplicationRepa> parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_1(std::size_t fmax, std::size_t tint, const SizeList&, std::size_t, const HistoryRepa&, int, SystemRepa&);


}



#endif