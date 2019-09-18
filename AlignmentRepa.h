#ifndef ALIGNMENTREPA_H
#define ALIGNMENTREPA_H

#include "Alignment.h"

#include <iostream>
#include <memory>
#include <string>

namespace Alignment
{
    typedef std::unordered_map<Variable, std::size_t> VarSizeUMap;
    typedef std::unordered_map<Value, std::size_t> ValSizeUMap;
    typedef std::vector<std::size_t> SizeList;
    typedef std::vector<double> DoubleList;
    typedef std::vector<std::vector<double>> DoubleListList;
}

namespace Alignment
{
    // data HistogramRepa = HistogramRepa {
    //   histogramRepasVectorVar :: !(V.Vector Variable),
    //   histogramRepasMapVarInt::Map.Map Variable Int,
    //   histogramRepasArray :: !(Array U VShape Double)

    struct HistogramRepa
    {
	VarList vectorVar;

	VarSizeUMap& mapVarInt() const;
	std::unique_ptr<VarSizeUMap> _mapVarInt;

	SizeList shape;

	DoubleList arr;
    };

    // systemsHistogramsHistogramRepa_u :: System -> Histogram -> Maybe HistogramRepa
    std::unique_ptr<HistogramRepa> systemsHistogramsHistogramRepa_u(const System&, const Histogram&);

    // systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
    std::unique_ptr<Histogram> systemsHistogramRepasHistogram_u(const System&, const HistogramRepa&);

    // setVarsHistogramRepasReduce_u :: Set.Set Variable -> HistogramRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistogramRepasReduce_u(const VarList&, const HistogramRepa&);
}

namespace Alignment
{
    // data HistogramRepaVec = HistogramRepaVec {
    //   histogramRepaVecsVectorVar :: !(V.Vector Variable),
    //   histogramRepaVecsMapVarInt::Map.Map Variable Int,
    //   histogramRepaVecsSize :: !Double,
    //   histogramRepaVecsShape :: !VShape,
    //   histogramRepaVecsArray :: !(V.Vector(UV.Vector Double))

    struct HistogramRepaVec
    {
	VarList vectorVar;

	VarSizeUMap& mapVarInt() const;
	std::unique_ptr<VarSizeUMap> _mapVarInt;

	SizeList shape;

	DoubleListList arr;
    };
}




#endif