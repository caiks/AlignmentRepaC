#ifndef ALIGNMENTREPA_H
#define ALIGNMENTREPA_H

#include "Alignment.h"

#include <iostream>
#include <memory>
#include <string>

namespace Alignment
{
    typedef std::pair<Variable, unsigned char> VarUCharPair;
    typedef std::vector<VarUCharPair> VarUCharPairList;
    typedef std::unordered_map<std::size_t, std::size_t> SizeSizeUMap;
    typedef std::unordered_map<Variable, std::size_t> VarSizeUMap;
    typedef std::unordered_map<Value, std::size_t> ValSizeUMap;
    typedef std::vector<std::size_t> SizeList;
    typedef std::vector<SizeList> SizeListList;
    typedef std::pair<double, SizeList> DoubleSizeListPair;
    typedef std::vector<DoubleSizeListPair> DoubleSizeListPairList;
    typedef std::unordered_set<std::size_t> SizeUSet;
    typedef std::vector<double> DoubleList;
    typedef std::vector<std::vector<double>> DoubleListList;
}

namespace Alignment
{
    class SystemRepa
    {
    public: SystemRepa();
    private: SystemRepa(const SystemRepa &);
    public: ~SystemRepa();

    private: SystemRepa& operator=(const SystemRepa &);

    public: VarUCharPairList listVarUCharPair;

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
    //   histogramRepasVectorVar :: !(V.Vector Variable),
    //   histogramRepasMapVarInt::Map.Map Variable Int,
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

    public: unsigned char* shape;
    public: double* arr;
    };

    // systemsHistogramsHistogramRepa_u :: System -> Histogram -> Maybe HistogramRepa
    std::unique_ptr<HistogramRepa> systemsHistogramsHistogramRepa_u(const System&, const SystemRepa&, const Histogram&);

    // systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
    std::unique_ptr<Histogram> systemsHistogramRepasHistogram_u(const System&, const SystemRepa&, const HistogramRepa&);

    // setVarsHistogramRepasReduce_u :: Set.Set Variable -> HistogramRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistogramRepasReduce_u(std::size_t m, std::size_t* kk, const HistogramRepa&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::HistogramRepa&);

namespace Alignment
{
    // data HistogramRepaRed = HistogramRepaRed{
    //   histogramRepaRedsVectorVar :: !(V.Vector Variable),
    //   histogramRepaRedsMapVarInt::Map.Map Variable Int,
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

    public: unsigned char* shape;
    public: double* arr;
    };

    // histogramRepasRed :: Double -> HistogramRepa -> HistogramRepaRed
    std::unique_ptr<HistogramRepaRed> histogramRepasRed(double, const HistogramRepa&);

    // histogramRepaRedsIndependent :: Double -> HistogramRepaRed -> HistogramRepa
    std::unique_ptr<HistogramRepa> histogramRepaRedsIndependent(double, const HistogramRepaRed&);
}

std::ostream& operator<<(std::ostream& out, const Alignment::HistogramRepaRed&);


//namespace Alignment
//{
//    // data HistogramRepaVec = HistogramRepaVec {
//    //   histogramRepaVecsVectorVar :: !(V.Vector Variable),
//    //   histogramRepaVecsMapVarInt::Map.Map Variable Int,
//    //   histogramRepaVecsSize :: !Double,
//    //   histogramRepaVecsShape :: !VShape,
//    //   histogramRepaVecsArray :: !(V.Vector(UV.Vector Double))
//
//    class HistogramRepaVec
//    {
//    public: HistogramRepaVec();
//    private: HistogramRepaVec(const HistogramRepaVec &);
//    public: ~HistogramRepaVec();
//
//    private: HistogramRepaVec& operator=(const HistogramRepaVec &);
//
//    public: VarList vectorVar;
//
//    public: VarSizeUMap& mapVarInt() const;
//    private: VarSizeUMap* _mapVarInt;
//
//    public: double size;
//    public: SizeList shape;
//    public: DoubleListList arr;
//    };
//}

namespace Alignment
{
    // data HistoryRepa = HistoryRepa {
    //   historyRepasVectorVar :: !(V.Vector Variable),
    //   historyRepasMapVarInt::Map.Map Variable Int,
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

    public: unsigned char* shape;
    public: std::size_t size;

    public: void transpose();
    public: unsigned char evient;
    public: unsigned char* arr;
    };

    // systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
    std::unique_ptr<HistoryRepa> systemsHistoriesHistoryRepa_u(const System&, const SystemRepa&, const History&, unsigned char evient = false);

    // systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
    std::unique_ptr<History> systemsHistoryRepasHistory_u(const System&, const SystemRepa&, const HistoryRepa&);

    // eventsHistoryRepasHistoryRepaSelection :: [Int] -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> eventsHistoryRepasHistoryRepaSelection_u(std::size_t, std::size_t*, const HistoryRepa&);

    // historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> historyRepasHistoryRepasHistoryRepaSelection_u(const HistoryRepa&, const HistoryRepa&);

    // setVarsHistoryRepasHistoryRepaReduced_u :: Set.Set Variable -> HistoryRepa -> HistoryRepa
    std::unique_ptr<HistoryRepa> setVarsHistoryRepasHistoryRepaReduced_u(std::size_t, std::size_t*, const HistoryRepa&);

    // setVarsHistoryRepasReduce_u :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa
    std::unique_ptr<HistogramRepa> setVarsHistoryRepasReduce_u(double, std::size_t, std::size_t*, const HistoryRepa&);

    // historyRepasRed :: HistoryRepa -> HistogramRepaRed
    std::unique_ptr<HistogramRepaRed> historyRepasRed(const HistoryRepa&);

}

std::ostream& operator<<(std::ostream& out, const Alignment::HistoryRepa&);


namespace Alignment
{
    // data TransformRepa = TransformRepa{
    //   transformRepasVectorVar :: !(V.Vector Variable),
    //   transformRepasMapVarInt::Map.Map Variable Int,
    //   transformRepasVarDerived :: !Variable,
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
    public: unsigned char valency;

    public: unsigned char* shape;
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
}

namespace Alignment
{
    struct FudRepa
    {
	TransformRepaPtrListList layers;
    };

    // setVariablesListTransformRepasFudRepa_u :: Set.Set Variable -> V.Vector TransformRepa -> FudRepa
    std::unique_ptr<FudRepa> setVariablesListTransformRepasFudRepa_u(const SizeUSet&, const TransformRepaPtrList&);

    // systemsFudsFudRepa_u :: System -> Fud -> FudRepa
    std::unique_ptr<FudRepa> systemsFudsFudRepa_u(const System&, const SystemRepa&, const Fud&);

    // systemsFudRepasFud_u :: System -> FudRepa -> Fud
    std::unique_ptr<Fud> systemsFudRepasFud_u(const System&, const SystemRepa&, const FudRepa&);

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
    std::unique_ptr<DecompFudRepa> systemsDecompFudsDecompFudRepa_u(const System&, const SystemRepa&, const DecompFud&);

    // systemsDecompFudRepasDecompFud_u :: System -> DecompFudRepa -> DecompFud
    std::unique_ptr<DecompFud> systemsDecompFudRepasDecompFud_u(const System&, const SystemRepa&, const DecompFudRepa&);
}

namespace Alignment
{
    // parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Integer -> V.Vector Variable -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> (V.Vector ((Double, V.Vector Variable)),Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t, std::size_t, const SizeList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);

    // parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Set.Set Variable -> V.Vector (Set.Set Variable) -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> (V.Vector ((Double, V.Vector Variable)),Integer)
    std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t, std::size_t, const SizeList&, const SizeListList&, const HistoryRepa&, const HistogramRepaRed&, const HistoryRepa&, const HistogramRepaRed&);
}


#endif