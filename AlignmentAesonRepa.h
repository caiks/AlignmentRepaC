#ifndef ALIGNMENTAESONREPA_H
#define ALIGNMENTAESONREPA_H

#include "AlignmentAeson.h"
#include "AlignmentRepa.h"

namespace Alignment
{
    // historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
    void historyRepasPersistent(const HistoryRepa&, std::ostream&);

    // persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
    std::unique_ptr<HistoryRepa> persistentsHistoryRepa(std::istream&, StrVarPtrMap&);

    // transformRepasPersistent :: TransformRepa -> TransformRepaPersistent
    void transformRepasPersistent(const TransformRepa&, std::ostream&);

    // persistentsTransformRepa :: TransformRepaPersistent -> Maybe TransformRepa
    std::unique_ptr<TransformRepa> persistentsTransformRepa(std::istream&, StrVarPtrMap&);

    // fudRepasPersistent :: FudRepa -> FudRepaPersistent
    void fudRepasPersistent(const FudRepa&, std::ostream&);

    // persistentsFudRepa :: FudRepaPersistent -> Maybe FudRepa
    std::unique_ptr<FudRepa> persistentsFudRepa(std::istream&, StrVarPtrMap&);


}



#endif