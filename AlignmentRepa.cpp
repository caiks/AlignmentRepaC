﻿#include "AlignmentRepa.h"
#include <iostream>

using namespace Alignment;

VarSizeUMap& Alignment::HistogramRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	std::unique_ptr<VarSizeUMap>& mm = (std::unique_ptr<VarSizeUMap>&)_mapVarInt;
	mm = std::move(std::unique_ptr<VarSizeUMap>(new VarSizeUMap(vectorVar.size())));
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    mm->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// systemsHistogramsHistogramRepa_u :: System -> Histogram -> Maybe HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::systemsHistogramsHistogramRepa_u(const System& uu, const Histogram& aa)
{
    auto ar = std::make_unique<HistogramRepa>();
    auto ww = histogramsSetVar(aa);
    auto n = ww->size();
    auto& vv = ar->vectorVar;
    vv.reserve(n);
    vv.insert(vv.end(), ww->begin(), ww->end());
    auto& sh = ar->shape;
    sh.reserve(n);
    std::vector<ValSizeUMap> mm(n);
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& xx = uu.map_u()[vv[i]];
	auto s = xx.size();
	sh.push_back(s);
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& w : xx)
	    yy.insert_or_assign(w, j++);
    }
    auto& rr = ar->arr;
    rr.resize(sz);
    for (auto& sc : aa.map_u())
    {
	auto& sm = sc.first.map_u();
	std::size_t j = 0;
	for (std::size_t i = 0; i < n; i++)
	    j = j*sh[i] + mm[i][sm.find(vv[i])->second];
	rr[j] = sc.second.getDouble();
    }
    return ar;
}

// systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
std::unique_ptr<Histogram> Alignment::systemsHistogramRepasHistogram_u(const System& uu, const HistogramRepa& ar)
{
    auto scalar = histogramScalar_u;
    auto& vv = ar.vectorVar;
    auto n = vv.size();
    auto& sh = ar.shape;
    auto& rr = ar.arr;
    auto sz = rr.size();
    if (!n)
	return scalar(Rational(sz ? rr[0] : 0));
    std::vector<ValList> mm(n);
    for (std::size_t i = 0; i < n; i++)
    {
	auto& xx = uu.map_u()[vv[i]];
	auto s = xx.size();
	auto& yy = mm[i];
	yy.reserve(s);
	for (auto& w : xx)
	    yy.push_back(w);
    }
    auto aa = std::make_unique<Histogram>();
    auto& am = aa->map_u();
    am.reserve(sz);
    SizeList ii(n);
    for (std::size_t j = 0; j < sz; j++)
    {
	std::vector<VarValPair> ss;
	ss.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	    ss.push_back(VarValPair(vv[i], mm[i][ii[i]]));
	am.insert_or_assign(State(ss), Rational(rr[j]));
	for (std::size_t k = n - 1; k >= 0; k--)
	{
	    std::size_t y = ii[k] + 1;
	    if (y == sh[k])
		ii[k] = 0;
	    else
	    {
		ii[k] = y;
		break;
	    }
	}
    }
    return aa;
}

// setVarsHistogramRepasReduce_u :: Set.Set Variable -> HistogramRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistogramRepasReduce_u(const VarList& kk, const HistogramRepa& ar)
{
    auto& vvv = ar.vectorVar;
    auto n = vvv.size();
    auto& svv = ar.shape;
    auto& mvv = ar.mapVarInt();
    auto& rvv = ar.arr;
    auto v = rvv.size();
    auto br = std::make_unique<HistogramRepa>();
    br->vectorVar = kk;
    auto& vkk = br->vectorVar;
    auto m = vkk.size();
    SizeList pkk(m);
    for (std::size_t i = 0; i < m; i++)
	pkk[i] = mvv[vkk[i]];
    auto& skk = br->shape;
    skk.reserve(m);
    std::size_t w = 1;
    for (std::size_t i = 0; i < m; i++)
    {
	auto s = svv[pkk[i]];
	w *= s;
	skk.push_back(s);
    }
    auto& rkk = br->arr;
    rkk.resize(w);
    SizeList ivv(n);
    for (std::size_t j = 0; j < v; j++)
    {
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++)
	    k = skk[i]*k + ivv[pkk[i]];
	rkk[k] += rvv[j];
	for (std::size_t i = n - 1; i >= 0; i--)
	{
	    std::size_t y = ivv[i] + 1;
	    if (y == svv[i])
		ivv[i] = 0;
	    else
	    {
		ivv[i] = y;
		break;
	    }
	}
    }
    return br;
}


VarSizeUMap& Alignment::HistogramRepaVec::mapVarInt() const
{
    if (!_mapVarInt)
    {
	std::unique_ptr<VarSizeUMap>& mm = (std::unique_ptr<VarSizeUMap>&)_mapVarInt;
	mm = std::move(std::unique_ptr<VarSizeUMap>(new VarSizeUMap(vectorVar.size())));
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    mm->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

VarSizeUMap& Alignment::HistoryRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	std::unique_ptr<VarSizeUMap>& mm = (std::unique_ptr<VarSizeUMap>&)_mapVarInt;
	mm = std::move(std::unique_ptr<VarSizeUMap>(new VarSizeUMap(vectorVar.size())));
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    mm->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::systemsHistoriesHistoryRepa_u(const System& uu, const History& hh)
{
    auto hr = std::make_unique<HistoryRepa>();
    auto ww = historiesSetVar(hh);
    auto n = ww->size();
    auto& vv = hr->vectorVar;
    vv.reserve(n);
    vv.insert(vv.end(), ww->begin(), ww->end());
    auto& sh = hr->shape;
    sh.reserve(n);
    std::vector<ValSizeUMap> mm(n);
    hr->size = hh.map_u().size();
    for (std::size_t i = 0; i < n; i++)
    {
	auto& xx = uu.map_u()[vv[i]];
	auto s = xx.size();
	sh.push_back(s);
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& w : xx)
	    yy.insert_or_assign(w, j++);
    }
    auto& rr = hr->arr;
    rr.reserve(hr->size * n);
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    for (auto& is : hh.map_u())
    {
	auto& sm = is.second.map_u();
	for (std::size_t i = 0; i < n; i++)
	{
	    std::size_t k = mm[i][sm.find(vv[i])->second];
	    if (k >= ucmax)
		throw std::out_of_range::out_of_range("systemsHistoriesHistoryRepa_u");
	    rr.push_back(k);
	}
    }
    return hr;
}

// systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
std::unique_ptr<History> Alignment::systemsHistoryRepasHistory_u(const System& uu, const HistoryRepa& hr)
{
    auto& vv = hr.vectorVar;
    auto n = vv.size();
    auto& sh = hr.shape;
    auto& rr = hr.arr;
    auto z = hr.size;
    std::vector<ValList> mm(n);
    for (std::size_t i = 0; i < n; i++)
    {
	auto& xx = uu.map_u()[vv[i]];
	auto s = xx.size();
	auto& yy = mm[i];
	yy.reserve(s);
	for (auto& w : xx)
	    yy.push_back(w);
    }
    auto hh = std::make_unique<History>();
    auto& hm = hh->map_u();
    hm.reserve(z);
    for (std::size_t j = 0; j < z; j++)
    {
	std::vector<VarValPair> ss;
	ss.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	    ss.push_back(VarValPair(vv[i], mm[i][rr[j*n+i]]));
	hm.insert_or_assign(Id(j+1), State(ss));
    }
    return hh;
}

// eventsHistoryRepasHistoryRepaSelection_u :: [Int] -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::eventsHistoryRepasHistoryRepaSelection_u(const SizeList& ll, const HistoryRepa& hr)
{
    auto hr1 = std::make_unique<HistoryRepa>();
    hr1->vectorVar = hr.vectorVar;
    hr1->size = ll.size();
    hr1->shape = hr.shape;
    auto n = hr.vectorVar.size();
    auto& rr = hr.arr;
    auto& rr1 = hr1->arr;
    rr1.reserve(hr1->size * n);
    for (auto j : ll)
    {
	std::size_t jn = j*n;
	for (std::size_t i = 0; i < n; i++)
	    rr1.push_back(rr[jn+i]);
    }
    return hr1;
}

// setVarsHistoryRepasHistoryRepaReduced_u :: Set.Set Variable -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::setVarsHistoryRepasHistoryRepaReduced_u(const VarList& kk, const HistoryRepa& hr)
{
    auto n = hr.vectorVar.size();
    auto& svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto m = kk.size();
    SizeList pkk(m);
    for (std::size_t i = 0; i < m; i++)
	pkk[i] = mvv[kk[i]];
    auto hr1 = std::make_unique<HistoryRepa>();
    hr1->vectorVar = kk;
    hr1->size = z;
    auto& skk = hr1->shape;
    skk.reserve(m);
    for (std::size_t i = 0; i < m; i++)
	skk.push_back(svv[pkk[i]]);
    auto& rr = hr.arr;
    auto& rr1 = hr1->arr;
    rr1.reserve(z*m);
    for (std::size_t j = 0; j < z; j++)
    {
	std::size_t jn = j*n;
	for (std::size_t i = 0; i < m; i++)
	    rr1.push_back(rr[jn + pkk[i]]);
    }
    return hr1;
}

// setVarsHistoryRepasReduce_u :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistoryRepasReduce_u(double f, const VarList& kk, const HistoryRepa& hr)
{
    auto n = hr.vectorVar.size();
    auto& svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto m = kk.size();
    SizeList pkk(m);
    for (std::size_t i = 0; i < m; i++)
	pkk[i] = mvv[kk[i]];
    auto ar1 = std::make_unique<HistogramRepa>();
    ar1->vectorVar = kk;
    auto& skk = ar1->shape;
    skk.reserve(m);
    std::size_t w = 1;
    for (std::size_t i = 0; i < m; i++)
    {
	auto s = svv[pkk[i]];
	w *= s;
	skk.push_back(s);
    }
    auto& rr = hr.arr;
    auto& rr1 = ar1->arr;
    rr1.resize(w);
    if (m > 0)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn = j*n;
	    std::size_t k = rr[jn+pkk[0]];
	    for (std::size_t i = 1; i < m; i++)
		k = skk[i]*k + rr[jn+pkk[i]];
	    rr1[k] += f;
	}
    else
	rr1[0] = f*z;
    return ar1;
}
