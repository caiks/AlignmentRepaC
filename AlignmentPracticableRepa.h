#ifndef ALIGNMENTPRACTICABLEREPA_H
#define ALIGNMENTPRACTICABLEREPA_H

#include "AlignmentRepa.h"

namespace Alignment
{
    // parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui ::
    //   Integer->Integer->Integer->Integer->System->Set.Set Variable->Fud ->
    //   HistoryRepa->HistogramRepaRed->HistoryRepa->HistogramRepaRed ->
    //   ([Set.Set Variable],Integer)
    std::tuple<std::unique_ptr<SizeListList>, std::size_t> parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(std::size_t, std::size_t, std::size_t, std::size_t, const SizeList&, const FudRepa&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);
}



#endif