#ifndef ALIGNMENTRANDOMREPA_H
#define ALIGNMENTRANDOMREPA_H

#include "AlignmentRepa.h"
#include <random>

namespace Alignment
{
	// historyRepasShuffle_u :: HistoryRepa -> Int -> HistoryRepa
	std::unique_ptr<HistoryRepa> historyRepasShuffle_u(const HistoryRepa&, unsigned);
	
	// thread safe version with std::subtract_with_carry_engine
	std::unique_ptr<HistoryRepa> historyRepasShuffle_us(const HistoryRepa&, std::ranlux48_base&);

	// historySparseArrayShuffle_us :: HistorySparseArray -> Int -> HistorySparseArray
	std::unique_ptr<HistorySparseArray> historySparseArrayShuffle_us(const HistorySparseArray&, std::ranlux48_base&);

}



#endif