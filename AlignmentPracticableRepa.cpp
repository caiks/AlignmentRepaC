#include "AlignmentPracticableRepa.h"

using namespace Alignment;

// parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui ::
//   Integer -> Integer -> Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([[VariableRepa]],Integer)
std::tuple<std::unique_ptr<SizeListList>, std::size_t> Alignment::parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, const SizeList& vv, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
    auto frvars = fudRepasSetVar;
    auto frder = fudRepasDerived;
    auto cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u;
    auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u;

    SizeList vv1;
    vv1.reserve(vv.size());
    {
	auto n = hh.dimension;
	auto& mvv = hh.mapVarInt();
	auto sh = hh.shape;
	auto pxx1 = hhx.arr;
	std::size_t sz = 0;
	SizeList xx;
	xx.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	{
	    xx.push_back(sz);
	    sz += sh[i];
	}
	for (auto& w : vv)
	{
	    auto p = mvv[w];
	    auto e = sh[p];
	    bool first = false;
	    bool second = false;
	    for (std::size_t i = 0; !second && i < e; i++)
	    {
		bool nz = pxx1[xx[p] + i] > 0.0;
		second = second || (first && nz);
		first = first || nz;
	    }
	    if (second)
		vv1.push_back(w);
	}
    }
    std::size_t s = 0;
    auto xx1 = std::make_unique<SizeListList>();
    if (vv1.size() < 2)
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
    DoubleSizeListPairList xx;
    xx.reserve(omax*10);
    if (!fr.layers.size())
    {
	auto t = cross(xmax, omax, vv1, hh, hhx, hhrr, hhrrx);
	auto& xc = std::get<0>(t);
	s += std::get<1>(t);
	xx.insert(xx.end(), xc->begin(), xc->end());
	xx1->reserve(xc->size());
	for (auto& pp : *xc)
	    xx1->push_back(pp.second);
    }
    else
    {
	auto vvf = frvars(fr);
	vv1.reserve(vv1.size() + vvf->size());
	vv1.insert(vv1.end(), vvf->begin(), vvf->end());
	auto wwf = frder(fr);
	xx1->reserve(wwf->size());
	for (auto& w : *wwf)
	    xx1->push_back(SizeList{ w });
    }
    while (xx1->size())
    {
	auto t = append(xmax, omax, vv1, *xx1, hh, hhx, hhrr, hhrrx);
	auto& xa = std::get<0>(t);
	s += std::get<1>(t);
	xx.insert(xx.end(), xa->begin(), xa->end());
	xx1->clear();
	xx1->reserve(xa->size());
	for (auto& pp : *xa)
	    xx1->push_back(pp.second);
    }
    if (!xx.size())
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
    std::sort(xx.begin(), xx.end());
    xx1->reserve(bmax / mmax);
    std::size_t start = xx.size() - 1;
    std::size_t end = xx.size() > bmax/mmax ? xx.size() - bmax/mmax : 0;
    for (long long i = start; i >= end; i--)
	xx1->push_back(xx[i].second);
    return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
}


// parametersSystemsPartitionerMaxRollByMRepa_ui ::
//   Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistoryRepa -> 
//   ([[[VariableRepa]]],Integer)
std::tuple<std::unique_ptr<SizeListListList>, std::size_t> Alignment::parametersSystemsPartitionerMaxRollByMRepa_ui(std::size_t mmax, std::size_t umax, std::size_t pmax, const SizeList& kk, const HistoryRepa& hh, const HistoryRepa& hhrr)
{
    auto hrred = [](double f, const HistoryRepa& hr, const SizeList& kk)
    {
	return setVarsHistoryRepasReduce_u(f, kk.size(), kk.data(), hr);
    };
    auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;

    auto z = (double)hh.size;
    auto zr = (double)hhrr.size;
    auto aa = hrred(1.0, hh, kk);
    auto aarr = hrred(z/zr, hhrr, kk);
    double y1 = aa->facLn() - aarr->facLn();
    return rrvqqy(mmax, umax, pmax, *aa, *aarr, z, y1);
}


// parametersRollerMaximumRollExcludedSelfRepa_i ::
//   [[VariableRepa]] -> -> HistogramRepa -> -> HistogramRepa ->
//   ([[[Int]]],Integer)
std::tuple<std::unique_ptr<SizeListListList>, std::size_t> Alignment::parametersRollerMaximumRollExcludedSelfRepa_i(const SizeListList& pp, const HistogramRepa& aa, const HistogramRepa& aarr, double z)
{
    auto t = histogramRepaVecsRollMax(pp, aa, aarr, z);
    auto tt = std::move(std::get<0>(t));
    auto s = std::get<1>(t);
    auto tt1 = std::make_unique<SizeListListList>();
    bool any = false;
    for (auto& vv : *tt)
	if (*max_element(vv.begin(), vv.end()) < vv.size() - 1)
	{
	    tt1->push_back(*tt);
	    break;
	}
    return std::tuple<std::unique_ptr<SizeListListList>, std::size_t>(std::move(tt1), s);
}

// parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui ::
//   Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([(Double, [VariableRepa]],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(std::size_t wmax, std::size_t omax, const SizeList& vv, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
    auto frdef = fudRepasDefinitions;
    auto frvars = fudRepasSetVar;
    auto frder = fudRepasDerived;
    auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u;

    std::size_t s = 0;
    SizeListList xx1;
    DoubleSizeListPairList xx;
    xx.reserve(omax * 10);
    auto vvf = frvars(fr);
    for (auto& w : vv)
	vvf->erase(w);
    SizeList yy;
    yy.reserve(vvf->size());
    yy.insert(yy.end(), vvf->begin(), vvf->end());
    SizeSizePairList cc;
    auto dd = frdef(fr);
    auto p = dd->size();
    for (auto& pp : *dd)
    {
	auto& w = pp.first;
	auto& m = pp.second->dimension;
	auto& zz = pp.second->vectorVar;
	SizeUSet uu1(p);
	SizeUSet uu2(m);
	for (std::size_t i = 0; i < m; i++)
	    uu2.insert(zz[i]);
	while (uu2.size())
	{
	    SizeList uu3;
	    uu3.reserve(uu2.size());
	    uu3.insert(uu3.end(), uu2.begin(), uu2.end());
	    for (auto& v : uu3)
	    {
		uu2.erase(v);
		if (uu1.find(v) == uu1.end() || dd->find(v) == dd->end())
		    break;
		uu1.insert(v);
		auto& tt = (*dd)[v];
		auto& m1 = tt->dimension;
		auto& zz1 = tt->vectorVar;
		for (std::size_t i = 0; i < m1; i++)
		    uu2.insert(zz1[i]);
	    }
	}
	for (auto& v : uu1)
	    cc.push_back(SizeSizePair(w, v));
    }
    auto wwf = frder(fr);
    xx1.reserve(wwf->size());
    for (auto& w : *wwf)
	xx1.push_back(SizeList{ w });
    while (xx1.size())
    {
	auto t = append(wmax, omax, cc, yy, xx1, hh, hhx, hhrr, hhrrx);
	auto& xa = std::get<0>(t);
	s += std::get<1>(t);
	xx.insert(xx.end(), xa->begin(), xa->end());
	xx1.clear();
	xx1.reserve(xa->size());
	for (auto& pp : *xa)
	    xx1.push_back(pp.second);
    }
    std::sort(xx.begin(), xx.end());
    auto xx2 = std::make_unique<DoubleSizeListPairList>();
    xx2->reserve(1);
    if (xx.size())
	xx2->push_back(xx.back());
    return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(xx2), s);
}

