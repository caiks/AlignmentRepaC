#ifndef ALIGNMENTREPA_H
#define ALIGNMENTREPA_H

#include "Alignment.h"

#include <iostream>
#include <memory>
#include <string>

namespace Alignment
{
    typedef std::unordered_map<Variable, std::size_t> VarSizeUMap;
    typedef std::vector<std::size_t> SizeList;
    typedef std::vector<double> DoubleList;

    // data HistogramRepa = HistogramRepa {
    //  histogramRepasVectorVar :: !(V.Vector Variable),
    //   histogramRepasMapVarInt::Map.Map Variable Int,
    //   histogramRepasArray :: !(Array U VShape Double)

    struct HistogramRepa
    {
	VarList vectorVar;

	VarSizeUMap& mapVarInt();
	std::unique_ptr<VarSizeUMap> _mapVarInt;

	SizeList shape;

	DoubleList arr;
    };

}





#endif