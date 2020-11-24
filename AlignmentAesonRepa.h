#ifndef ALIGNMENTAESONREPA_H
#define ALIGNMENTAESONREPA_H

#include "AlignmentAeson.h"
#include "AlignmentRepa.h"

namespace Alignment
{
	// systemRepasPersistent :: SystemRepa -> SystemRepaPersistent
	void systemRepasPersistent(const SystemRepa&, std::ostream&);

	// persistentsSystemRepa :: SystemRepaPersistent -> Maybe SystemRepa
	std::unique_ptr<SystemRepa> persistentsSystemRepa(std::istream&, StrVarPtrMap&);

	// historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
	void historyRepasPersistent(const HistoryRepa&, std::ostream&);

	// persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
	std::unique_ptr<HistoryRepa> persistentsHistoryRepa(std::istream&);

	// historyRepasPersistentInitial :: HistoryRepa -> Int -> HistoryRepaPersistentInitial
	// assumes evient
	void historyRepasPersistentInitial(const HistoryRepa&, std::size_t, std::ostream&);

	// persistentInitialsHistoryRepa :: HistoryRepaPersistentInitial -> Maybe HistoryRepa
	std::unique_ptr<HistoryRepa> persistentInitialsHistoryRepa(std::istream&);

	// historySparsesPersistent :: HistorySparse -> HistorySparsePersistent
	void historySparsesPersistent(const HistorySparse&, std::ostream&);

	// persistentsHistorySparse :: HistorySparsePersistent -> Maybe HistorySparse
	std::unique_ptr<HistorySparse> persistentsHistorySparse(std::istream&);

	// historySparseArraysPersistent :: HistorySparseArray -> HistorySparseArrayPersistent
	void historySparseArraysPersistent(const HistorySparseArray&, std::ostream&);

	// persistentsHistorySparseArray :: HistorySparseArrayPersistent -> Maybe HistorySparseArray
	std::unique_ptr<HistorySparseArray> persistentsHistorySparseArray(std::istream&);

	// historySparseArraysPersistentInitial :: HistorySparseArray -> Int -> HistorySparseArrayPersistentInitial
	void historySparseArraysPersistentInitial(const HistorySparseArray&, std::size_t, std::ostream&);

	// persistentInitialsHistorySparseArray :: HistorySparseArrayPersistentInitial -> Maybe HistorySparseArray
	std::unique_ptr<HistorySparseArray> persistentInitialsHistorySparseArray(std::istream&);

	// transformRepasPersistent :: TransformRepa -> TransformRepaPersistent
	void transformRepasPersistent(const TransformRepa&, std::ostream&);

	// persistentsTransformRepa :: TransformRepaPersistent -> Maybe TransformRepa
	std::unique_ptr<TransformRepa> persistentsTransformRepa(std::istream&);

	// fudRepasPersistent :: FudRepa -> FudRepaPersistent
	void fudRepasPersistent(const FudRepa&, std::ostream&);

	// persistentsFudRepa :: FudRepaPersistent -> Maybe FudRepa
	std::unique_ptr<FudRepa> persistentsFudRepa(std::istream&);

	// decompFudRepasPersistent :: DecompFudRepa -> DecompFudRepaPersistent
	void decompFudRepasPersistent(const DecompFudRepa&, std::ostream&);

	// persistentsDecompFudRepa :: DecompFudRepaPersistent -> Maybe DecompFudRepa
	std::unique_ptr<DecompFudRepa> persistentsDecompFudRepa(std::istream&);

	// applicationRepasPersistent :: ApplicationRepa -> ApplicationRepaPersistent
	void applicationRepasPersistent(const ApplicationRepa&, std::ostream&);

	// persistentsApplicationRepa :: ApplicationRepaPersistent -> Maybe ApplicationRepa
	std::unique_ptr<ApplicationRepa> persistentsApplicationRepa(std::istream&);

	// decompFudSlicedRepasPersistent :: DecompFudSlicedRepa -> DecompFudSlicedRepaPersistent
	void decompFudSlicedRepasPersistent(const DecompFudSlicedRepa&, std::ostream&);

	// persistentsDecompFudSlicedRepa :: DecompFudSlicedRepaPersistent -> Maybe DecompFudSlicedRepa
	std::unique_ptr<DecompFudSlicedRepa> persistentsDecompFudSlicedRepa(std::istream&);
}



#endif