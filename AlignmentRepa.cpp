#include "AlignmentRepa.h"
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
    auto sat = statesVarsValue;

    auto ar = std::make_unique<HistogramRepa>();
    auto ww = histogramsSetVar(aa);
    auto n = ww->size();
    auto& vv = ar->vectorVar;
    vv.reserve(n);
    for (auto& v : *ww)
	vv.push_back(v);
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
    auto zvv = rvv.size();
    auto br = std::make_unique<HistogramRepa>();
    br->vectorVar = kk;
    auto& vkk = br->vectorVar;
    auto m = vkk.size();
    SizeList pkk(m);
    for (std::size_t i = 0; i < m; i++)
	pkk[i] = mvv[vkk[i]];
    auto& skk = br->shape;
    skk.reserve(m);
    std::size_t zkk = 1;
    for (std::size_t i = 0; i < m; i++)
    {
	auto s = svv[pkk[i]];
	zkk *= s;
	skk.push_back(s);
    }
    auto& rkk = br->arr;
    rkk.resize(zkk);
    SizeList ivv(n);
    for (std::size_t j = 0; j < zvv; j++)
    {
	std::size_t i = 0;
	for (std::size_t k = 0; k < m; k++)
	    i = i*skk[k] + ivv[pkk[k]];
	rkk[i] += rvv[j];
	for (std::size_t k = n - 1; k >= 0; k--)
	{
	    std::size_t y = ivv[k] + 1;
	    if (y == svv[k])
		ivv[k] = 0;
	    else
	    {
		ivv[k] = y;
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