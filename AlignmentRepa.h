#ifndef ALIGNMENTREPA_H
#define ALIGNMENTREPA_H

#include "Alignment.h"

#include <iostream>
#include <memory>
#include <string>

namespace Alignment
{
    typedef std::pair<Variable, std::size_t> VarSizePair;
    typedef std::vector<VarSizePair> VarSizePairList;
    typedef std::unordered_map<std::size_t, std::size_t> SizeSizeUMap;
    typedef std::unordered_map<Variable, std::size_t> VarSizeUMap;
    typedef std::unordered_map<Value, std::size_t> ValSizeUMap;
    typedef std::vector<std::size_t> SizeList;
    typedef std::vector<SizeList> SizeListList;
    typedef std::vector<SizeListList> SizeListListList;
    typedef std::pair<std::size_t, std::size_t> SizeSizePair;
    typedef std::vector<SizeSizePair> SizeSizePairList;
    typedef std::pair<double, SizeList> DoubleSizeListPair;
    typedef std::vector<DoubleSizeListPair> DoubleSizeListPairList;
    typedef std::unordered_set<std::size_t> SizeUSet;
    typedef std::set<std::size_t> SizeSet;
    typedef std::map<std::size_t, SizeSet> SizeSizeSetMap;
    typedef std::vector<SizeSizeSetMap> SizeSizeSetMapList;
    typedef std::vector<double> DoubleList;
    typedef std::vector<std::vector<double>> DoubleListList;
    typedef Tree<std::size_t> SizeTree;
    typedef std::pair<std::size_t, std::shared_ptr<SizeTree>> SizeSizeTreePair;
}

namespace Alignment
{
    class SystemRepa
    {
    public: SystemRepa();
    private: SystemRepa(const SystemRepa &);
    public: ~SystemRepa();

    private: SystemRepa& operator=(const SystemRepa &);

    public: VarSizePairList listVarSizePair;

    public: VarSizeUMap& mapVarSize() const;
    private: VarSizeUMap* _mapVarSize;
    };

    // systemsSystemRepa :: System -> SystemRepa
    std::unique_ptr<SystemRepa> systemsSystemRepa(const System&);

    // systemsRepasSystem :: SystemRepa -> System
    void systemsRepasSystem(const SystemRepa&, System&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::SystemRepa&);


namespace Alignment
{
    // data HistogramRepa = HistogramRepa {
    //   histogramRepasVectorVar :: !(V.Vector VariableRepa),
    //   histogramRepasMapVarInt::Map.Map VariableRepa Int,
    //   histogramRepasArray :: !(Array U VShape Double)

    class HistogramRepa
    {
    public: HistogramRepa();
    public: HistogramRepa(double);
    private: HistogramRepa(const HistogramRepa &);
    public: ~HistogramRepa();

    private: HistogramRepa& operator=(const HistogramRepa &);

    public: std::size_t dimension;
    public: std::size_t* vectorVar;

    public: SizeSizeUMap& mapVarInt() const;
    private: SizeSizeUMap* _mapVarInt;

    public: std::size_t* shape;
    public: double* arr;

    public: double size() const;
    public: double facLn() const;

    };

    // systemsHistogramsHistogramRepa_u :: System -> Histogram -> Maybe HistogramRepa
    std::unique_ptr<HistogramRepa> systemsHistogramsHistogramRepa_u(const System&, const SystemRepa&, const Histogram&);

    // systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
    std::unique_ptr<Histogram> systemsHistogramRepasHistogram_u(const System&, const SystemRepa&, const HistogramRepa&);

    // setVarsHistogramRepasReduce_u :: [VariableRepa] -> HistogramRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistogramRepasReduce_u(std::size_t m, const std::size_t* kk, const HistogramRepa&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::HistogramRepa&);

namespace Alignment
{
    // data HistogramRepaRed = HistogramRepaRed{
    //   histogramRepaRedsVectorVar :: !(V.Vector VariableRepa),
    //   histogramRepaRedsMapVarInt::Map.Map VariableRepa Int,
    //   histogramRepaRedsShape :: !VShape,
    //   histogramRepaRedsVectorArray :: !(V.Vector(UV.Vector Double)) }

    class HistogramRepaRed
    {
    public: HistogramRepaRed();
    private: HistogramRepaRed(const HistogramRepaRed &);
    public: ~HistogramRepaRed();

    private: HistogramRepaRed& operator=(const HistogramRepaRed &);

    public: std::size_t dimension;
    public: std::size_t* vectorVar;

    public: SizeSizeUMap& mapVarInt() const;
    private: SizeSizeUMap* _mapVarInt;

    public: std::size_t* shape;
    public: double* arr;
    };

    // histogramRepasRed :: Double -> HistogramRepa -> HistogramRepaRed
    std::unique_ptr<HistogramRepaRed> histogramRepasRed(double, const HistogramRepa&);

    // histogramRepaRedsIndependent :: Double -> HistogramRepaRed -> HistogramRepa
    std::unique_ptr<HistogramRepa> histogramRepaRedsIndependent(double, const HistogramRepaRed&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::HistogramRepaRed&);

namespace Alignment
{
    // data HistoryRepa = HistoryRepa {
    //   historyRepasVectorVar :: !(V.Vector VariableRepa),
    //   historyRepasMapVarInt::Map.Map VariableRepa Int,
    //   historyRepasShape :: !VShape,
    //   historyRepasArray :: !(Array U DIM2 Int)

    class HistoryRepa
    {
    public: HistoryRepa();
    private: HistoryRepa(const HistoryRepa &);
    public: ~HistoryRepa();

    private: HistoryRepa& operator=(const HistoryRepa &);

    public: std::size_t dimension;
    public: std::size_t* vectorVar;

    public: SizeSizeUMap& mapVarInt() const;
    private: SizeSizeUMap* _mapVarInt;

    public: std::size_t* shape;
    public: std::size_t size;

    public: void transpose();
    public: unsigned char evient;
    public: unsigned char* arr;
    };

    typedef std::shared_ptr<HistoryRepa> HistoryRepaPtr;
    typedef std::vector<HistoryRepaPtr> HistoryRepaPtrList;

    // systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
    std::unique_ptr<HistoryRepa> systemsHistoriesHistoryRepa_u(const System&, const SystemRepa&, const History&, unsigned char evient = false);

    // systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
    std::unique_ptr<History> systemsHistoryRepasHistory_u(const System&, const SystemRepa&, const HistoryRepa&);

    // eventsHistoryRepasHistoryRepaSelection :: [Int] -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> eventsHistoryRepasHistoryRepaSelection_u(std::size_t, const std::size_t*, const HistoryRepa&);

    // historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasHistoryRepasHistoryRepaSelection_u(const HistoryRepa&, const HistoryRepa&);

    // setVarsHistoryRepasHistoryRepaReduced_u :: [VariableRepa] -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> setVarsHistoryRepasHistoryRepaReduced_u(std::size_t, const std::size_t*, const HistoryRepa&);

    // setVarsHistoryRepasReduce_u :: Double -> [VariableRepa] -> HistoryRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistoryRepasReduce_u(double, std::size_t, const std::size_t*, const HistoryRepa&);

    // historyRepasRed :: HistoryRepa -> HistogramRepaRed
    std::unique_ptr<HistogramRepaRed> historyRepasRed(const HistoryRepa&);

    // vectorHistoryRepasConcat_u :: V.Vector HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> vectorHistoryRepasConcat_u(const HistoryRepaPtrList&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::HistoryRepa&);


namespace Alignment
{
    // data TransformRepa = TransformRepa{
    //   transformRepasVectorVar :: !(V.Vector VariableRepa),
    //   transformRepasMapVarInt::Map.Map VariableRepa Int,
    //   transformRepasVarDerived :: !VariableRepa,
    //   transformRepasValency :: !Int,
    //   transformRepasArray :: !(Array U VShape Int) }

    class TransformRepa
    {
    public: TransformRepa();
    private: TransformRepa(const TransformRepa &);
    public: ~TransformRepa();

    private: TransformRepa& operator=(const TransformRepa &);

    public: std::size_t dimension;
    public: std::size_t* vectorVar;

    public: SizeSizeUMap& mapVarInt() const;
    private: SizeSizeUMap* _mapVarInt;

    public: std::size_t derived;
    public: std::size_t valency;

    public: std::size_t* shape;
    public: unsigned char* arr;
    };

    // systemsTransformsTransformRepa_u :: System -> Transform -> TransformRepa
    std::unique_ptr<TransformRepa> systemsTransformsTransformRepa_u(const System&, const SystemRepa&, const Transform&);

    // systemsTransformRepasTransform_u :: System -> TransformRepa -> Transform
    std::unique_ptr<Transform> systemsTransformRepasTransform_u(const System&, const SystemRepa&, const TransformRepa&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::TransformRepa&);


namespace Alignment
{
    typedef std::shared_ptr<TransformRepa> TransformRepaPtr;
    typedef std::vector<TransformRepaPtr> TransformRepaPtrList;
    typedef std::vector<TransformRepaPtrList> TransformRepaPtrListList;
    typedef std::map<std::size_t, TransformRepaPtr> SizeTransformRepaPtrMap;
}

namespace Alignment
{
    struct FudRepa
    {
	TransformRepaPtrListList layers;
    };

    // setVariablesListTransformRepasFudRepa_u :: [VariableRepa] -> [TransformRepa] -> FudRepa
    std::unique_ptr<FudRepa> setVariablesListTransformRepasFudRepa_u(const SizeUSet&, const TransformRepaPtrList&);

    // systemsFudsFudRepa_u :: System -> Fud -> FudRepa
    std::unique_ptr<FudRepa> systemsFudsFudRepa_u(const System&, const SystemRepa&, const Fud&);

    // systemsFudRepasFud_u :: System -> FudRepa -> Fud
    std::unique_ptr<Fud> systemsFudRepasFud_u(const System&, const SystemRepa&, const FudRepa&);

    // fudRepasSize :: FudRepa -> Integer
    std::size_t fudRepasSize(const FudRepa&);

    // fudRepasSetVar :: FudRepa -> [VariableRepa]
    std::unique_ptr<SizeUSet> fudRepasSetVar(const FudRepa&);

    // fudRepasDerived :: FudRepa -> [VariableRepa]
    std::unique_ptr<SizeUSet> fudRepasDerived(const FudRepa&);

    // fudRepasUnderlying :: FudRepa -> [VariableRepa]
    std::unique_ptr<SizeUSet> fudRepasUnderlying(const FudRepa&);

    // fudRepasDefinitions :: FudRepa -> Map.Map VariableRepa TransformRepa
    std::unique_ptr<SizeTransformRepaPtrMap> fudRepasDefinitions(const FudRepa&);

    // fudsSetVarsDepends :: FudRepa -> [VariableRepa] -> FudRepa
    std::unique_ptr<TransformRepaPtrList> fudsSetVarsDepends(const FudRepa&, const SizeUSet&);

    // historyRepasFudRepasMultiply_u :: HistoryRepa -> FudRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasFudRepasMultiply_u(const HistoryRepa&, const FudRepa&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::FudRepa&);


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
    std::unique_ptr<DecompFudRepa> systemsDecompFudsDecompFudRepa_u(const System&, const SystemRepa&, const DecompFud&);

    // systemsDecompFudRepasDecompFud_u :: System -> DecompFudRepa -> DecompFud
    std::unique_ptr<DecompFud> systemsDecompFudRepasDecompFud_u(const System&, const SystemRepa&, const DecompFudRepa&);
}

namespace Alignment
{
    struct ApplicationRepa
    {
	SizeTree slices;
	FudRepa fud;
	SizeList substrate;
    };
}

std::ostream& operator<<(std::ostream& out, const Alignment::ApplicationRepa&);


namespace Alignment
{
    // parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->([(Double, [VariableRepa])],Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t, std::size_t, const SizeList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);

    // parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa])],Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t, std::size_t, const SizeList&, const SizeListList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);

    // parametersHistogramRepaVecsSetTuplePartitionTopByM_u ::
    //   Integer -> Integer -> Integer -> HistogramRepa -> HistogramRepa -> Double -> Double -> ([[[VariableRepa]]],Integer)
    std::tuple<std::unique_ptr<SizeListListList>, std::size_t> parametersHistogramRepaVecsSetTuplePartitionTopByM_u(std::size_t, std::size_t, std::size_t, const HistogramRepa&, const HistogramRepa&, double, double);

    // histogramRepaVecsRollMax :: [[VariableRepa]] -> HistogramRepa -> HistogramRepa -> Double -> Double -> ([[VariableRepa]],Integer)
    std::tuple<std::unique_ptr<SizeListList>, std::size_t> histogramRepaVecsRollMax(const SizeListList&, const HistogramRepa&, const HistogramRepa&, double);

    // parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u :: Integer -> Integer -> [(Variable, Variable)] -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa]], Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u(std::size_t, std::size_t, const SizeSizePairList&, const SizeList&, const SizeListList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);

}


#endif
