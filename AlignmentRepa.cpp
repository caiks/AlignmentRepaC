#include "AlignmentRepa.h"
#include <iostream>

using namespace Alignment;

std::unordered_map<Variable, std::size_t>& Alignment::HistogramRepa::mapVarInt()
{
    if (!_mapVarInt)
    {
	_mapVarInt = std::move(std::unique_ptr<VarSizeUMap>(new VarSizeUMap(vectorVar.size())));
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    _mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}