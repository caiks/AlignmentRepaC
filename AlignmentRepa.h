﻿#ifndef ALIGNMENTREPA_H
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

	double size;
	SizeList shape;
	DoubleListList arr;
    };
}

namespace Alignment
{
    // data HistoryRepa = HistoryRepa {
    //   historyRepasVectorVar :: !(V.Vector Variable),
    //   historyRepasMapVarInt::Map.Map Variable Int,
    //   historyRepasShape :: !VShape,
    //   historyRepasArray :: !(Array U DIM2 Int)

    struct HistoryRepa
    {
	VarList vectorVar;

	VarSizeUMap& mapVarInt() const;
	std::unique_ptr<VarSizeUMap> _mapVarInt;

	std::size_t size;
	SizeList shape;
	std::unique_ptr<unsigned char[]> arr;
    };

    // systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
    std::unique_ptr<HistoryRepa> systemsHistoriesHistoryRepa_u(const System&, const History&);

    // systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
    std::unique_ptr<History> systemsHistoryRepasHistory_u(const System&, const HistoryRepa&);

    // eventsHistoryRepasHistoryRepaSelection :: [Int] -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> eventsHistoryRepasHistoryRepaSelection_u(const SizeList&, const HistoryRepa&);

    // historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasHistoryRepasHistoryRepaSelection_u(const HistoryRepa&, const HistoryRepa&);

    // setVarsHistoryRepasHistoryRepaReduced_u :: Set.Set Variable -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> setVarsHistoryRepasHistoryRepaReduced_u(const VarList&, const HistoryRepa&);

    // setVarsHistoryRepasReduce_u :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistoryRepasReduce_u(double, const VarList&, const HistoryRepa&);
}

namespace Alignment
{
    // data TransformRepa = TransformRepa{
    //   transformRepasVectorVar :: !(V.Vector Variable),
    //   transformRepasMapVarInt::Map.Map Variable Int,
    //   transformRepasVarDerived :: !Variable,
    //   transformRepasValency :: !Int,
    //   transformRepasArray :: !(Array U VShape Int) }

    struct TransformRepa
    {
	VarList vectorVar;

	VarSizeUMap& mapVarInt() const;
	std::unique_ptr<VarSizeUMap> _mapVarInt;

	std::unique_ptr<Variable> derived;
	unsigned char valency;

	SizeList shape;
	std::unique_ptr<unsigned char[]> arr;
    };

    // systemsTransformsTransformRepa_u :: System -> Transform -> TransformRepa
    std::unique_ptr<TransformRepa> systemsTransformsTransformRepa_u(const System&, const Transform&);

    // systemsTransformRepasTransform_u :: System -> TransformRepa -> Transform
    std::unique_ptr<Transform> systemsTransformRepasTransform_u(const System&, const TransformRepa&);
}

namespace Alignment
{
    typedef std::shared_ptr<TransformRepa> TransformRepaPtr;
    typedef std::vector<TransformRepaPtr> TransformRepaPtrList;
    typedef std::vector<TransformRepaPtrList> TransformRepaPtrListList;
}

namespace Alignment
{
    struct FudRepa
    {
	TransformRepaPtrListList layers;
    };

    // setVariablesListTransformRepasFudRepa_u :: Set.Set Variable -> V.Vector TransformRepa -> FudRepa
    std::unique_ptr<FudRepa> setVariablesListTransformRepasFudRepa_u(const VarUSet&, const TransformRepaPtrList&);

    // systemsFudsFudRepa_u :: System -> Fud -> FudRepa
    std::unique_ptr<FudRepa> systemsFudsFudRepa_u(const System&, const Fud&);

    // systemsFudRepasFud_u :: System -> FudRepa -> Fud
    std::unique_ptr<Fud> systemsFudRepasFud_u(const System&, const FudRepa&);

    // historyRepasFudRepasMultiply_u :: HistoryRepa -> FudRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasFudRepasMultiply_u(const HistoryRepa&, const FudRepa&);
}

namespace Alignment
{
    struct HistoryRepaPtrFudRepaPtrPair
    {
	HistoryRepaPtrFudRepaPtrPair() {}
	HistoryRepaPtrFudRepaPtrPair(std::shared_ptr<HistoryRepa> ss, std::shared_ptr<FudRepa> ff) : _state(ss), _fud(ff) {}

	std::shared_ptr<HistoryRepa> _state;
	std::shared_ptr<FudRepa> _fud;
    };

    struct DecompFudRepa
    {
	Tree<HistoryRepaPtrFudRepaPtrPair> tree;
    };

    // systemsDecompFudsDecompFudRepa_u :: System -> DecompFud -> DecompFudRepa
    std::unique_ptr<DecompFudRepa> systemsDecompFudsDecompFudRepa_u(const System&, const DecompFud&);

    // systemsDecompFudRepasDecompFud_u :: System -> DecompFudRepa -> DecompFud
    std::unique_ptr<DecompFud> systemsDecompFudRepasDecompFud_u(const System&, const DecompFudRepa&);
}


#endif