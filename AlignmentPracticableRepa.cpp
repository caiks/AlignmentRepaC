#include "AlignmentPracticableRepa.h"

using namespace Alignment;

// parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui ::
//   Integer->Integer->Integer->Integer->System->Set.Set Variable->Fud ->
//   HistoryRepa->HistogramRepaRed->HistoryRepa->HistogramRepaRed ->
//   ([Set.Set Variable],Integer)
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
    /**/
	std::cout << "s = " << s << std::endl;
    	for (auto& pp : *xa)
		std::cout << pp.first << "," << pp.second << std::endl;
    /**/
	xx.insert(xx.end(), xa->begin(), xa->end());
	xx1->clear();
	xx1->reserve(xa->size());
	for (auto& pp : *xa)
	    xx1->push_back(pp.second);
    }
    if (!xx.size())
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
    std::sort(xx.begin(), xx.end());
    /**/
    for (auto& pp : xx)
	std::cout << pp.first << "," << pp.second << std::endl;
    /**/
    xx1->reserve(bmax / mmax);
    std::size_t start = xx.size() - 1;
    std::size_t end = xx.size() > bmax/mmax ? xx.size() - bmax/mmax : 0;
    for (long long i = start; i >= end; i--)
	xx1->push_back(xx[i].second);
    return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
}
