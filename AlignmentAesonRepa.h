#ifndef ALIGNMENTAESONREPA_H
#define ALIGNMENTAESONREPA_H

#include "AlignmentRepa.h"

namespace Alignment
{
    // historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
    void historyRepasPersistent(const HistoryRepa&, std::ostream&);

    // persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
    std::unique_ptr<HistoryRepa> persistentsHistoryRepa(std::istream&);

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


}



#endif