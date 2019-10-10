#include "AlignmentRepa.h"
#include <iostream>

using namespace Alignment;

SystemRepa::SystemRepa() : _mapVarSize(0)
{
}

SystemRepa::~SystemRepa()
{
    delete _mapVarSize;
}

VarSizeUMap& Alignment::SystemRepa::mapVarSize() const
{
    if (!_mapVarSize)
	const_cast<SystemRepa*>(this)->_mapVarSize = new VarSizeUMap(listVarSizePair.size());
    if (_mapVarSize->size() < listVarSizePair.size())
    {
	for (std::size_t i = _mapVarSize->size(); i < listVarSizePair.size(); i++)
	    const_cast<SystemRepa*>(this)->_mapVarSize->insert_or_assign(listVarSizePair[i].first, i);
    }
    return *_mapVarSize;
}

// systemsSystemRepa :: System -> SystemRepa
std::unique_ptr<SystemRepa> Alignment::systemsSystemRepa(const System& uu)
{
    auto ur = std::make_unique<SystemRepa>();
    auto& ll = ur->listVarSizePair;
    auto mm = sorted(uu.map_u());
    ll.reserve(mm.size());
    for (auto& vww : mm)
	ll.push_back(VarSizePair(vww.first, vww.second.size()));
    return ur;
}

// systemsRepasSystem :: SystemRepa -> System
void Alignment::systemsRepasSystem(const SystemRepa& ur, System& uu)
{
    for (auto& vs : ur.listVarSizePair)
    {
	auto it = uu.map_u().find(vs.first);
	if (it == uu.map_u().end())
	{
	    ValSet ww;
	    for (int i = 0; i < vs.second; i++)
		ww.insert(Value(i));
	    uu.map_u().insert_or_assign(vs.first, ww);
	}
    }
}

HistogramRepa::HistogramRepa() : _mapVarInt(0)
{
}

HistogramRepa::~HistogramRepa()
{
    delete _mapVarInt;
}

VarSizeUMap& Alignment::HistogramRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistogramRepa*>(this)->_mapVarInt = new VarSizeUMap(vectorVar.size());
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    const_cast<HistogramRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
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
	auto xx = systemsVarsSetValue(uu,vv[i]);
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
    auto rrp = rr.data();
    for (auto& sc : aa.map_u())
    {
	auto& sm = sc.first.map_u();
	std::size_t j = 0;
	for (std::size_t i = 0; i < n; i++)
	    j = sh[i]*j + mm[i][sm.find(vv[i])->second];
	rrp[j] = sc.second.getDouble();
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
	auto xx = systemsVarsSetValue(uu,vv[i]);
	auto s = xx.size();
	auto& yy = mm[i];
	yy.reserve(s);
	for (auto& w : xx)
	    yy.push_back(w);
    }
    auto aa = std::make_unique<Histogram>();
    aa->map_u().reserve(sz);
    SizeList ii(n);
    auto rrp = rr.data();
    for (std::size_t j = 0; j < sz; j++)
    {
	std::vector<VarValPair> ss;
	ss.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	{
	    VarValPair pp(vv[i], mm[i][ii[i]]);
	    ss.push_back(pp);
	}
	aa->map_u().insert_or_assign(State(ss),Rational(rrp[j]));
	for (int k = n - 1; k >= 0; k--)
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
    auto rvvp = rvv.data();
    auto rkkp = rkk.data();
    for (std::size_t j = 0; j < v; j++)
    {
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++)
	    k = skk[i]*k + ivv[pkk[i]];
	rkkp[k] += rvvp[j];
	for (int i = n - 1; i >= 0; i--)
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

HistogramRepaVec::HistogramRepaVec() : _mapVarInt(0)
{
}

HistogramRepaVec::~HistogramRepaVec()
{
    delete _mapVarInt;
}

VarSizeUMap& Alignment::HistogramRepaVec::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistogramRepaVec*>(this)->_mapVarInt = new VarSizeUMap(vectorVar.size());
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    const_cast<HistogramRepaVec*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

HistoryRepa::HistoryRepa() : _mapVarInt(0), arr(0)
{
}

HistoryRepa::~HistoryRepa()
{
    delete _mapVarInt;
    delete[] arr;
}

VarSizeUMap& Alignment::HistoryRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistoryRepa*>(this)->_mapVarInt = new VarSizeUMap(vectorVar.size());
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    const_cast<HistoryRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
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
	auto xx = systemsVarsSetValue(uu,vv[i]);
	auto s = xx.size();
	sh.push_back(s);
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& w : xx)
	    yy.insert_or_assign(w, j++);
    }
    hr->arr = new unsigned char[hr->size * n];
    auto rr = hr->arr;
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    auto hm = sorted(hh.map_u());
    std::size_t j = 0;
    for (auto& is : hm)
    {
	auto& sm = is.second.map_u();
	for (std::size_t i = 0; i < n; i++)
	{
	    std::size_t k = mm[i][sm.find(vv[i])->second];
	    if (k > ucmax)
		throw std::out_of_range("systemsHistoriesHistoryRepa_u");
	    rr[j] = k;
	    j++;
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
    auto rr = hr.arr;
    auto z = hr.size;
    std::vector<ValList> mm(n);
    for (std::size_t i = 0; i < n; i++)
    {
	auto xx = systemsVarsSetValue(uu,vv[i]);
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
    auto rr = hr.arr;
    hr1->arr = new unsigned char[hr1->size * n];
    auto rr1 = hr1->arr;
    std::size_t k = 0;
    for (auto j : ll)
    {
	std::size_t jn = j*n;
	for (std::size_t i = 0; i < n; i++)
	{
	    rr1[k] = rr[jn + i];
	    k++;
	}
    }
    return hr1;
}

// historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historyRepasHistoryRepasHistoryRepaSelection_u(const HistoryRepa& ss, const HistoryRepa& hr)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;

    auto& kk = ss.vectorVar;
    auto n = hr.vectorVar.size();
    auto& svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto m = kk.size();
    auto y = ss.size;
    SizeList pkk;
    for (std::size_t i = 0; i < m; i++)
	pkk.push_back(mvv[kk[i]]);
    auto rr1 = ss.arr;
    auto rr = hr.arr;
    SizeList ll;
    for (std::size_t j = 0; j < z; j++)
    {
	std::size_t jn = j*n;
	bool any = false;
	for (std::size_t k = 0; !any && k < y; k++)
	{
	    std::size_t km = k*m;
	    bool all = true;
	    for (std::size_t i = 0; all && i < m; i++)
		all = rr1[km + i] == rr[jn + pkk[i]];
	    any = all;
	}
	if (any)
	    ll.push_back(j);
    }
    return hrsel(ll,hr);
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
    auto rr = hr.arr;
    hr1->arr = new unsigned char[z*m];
    auto rr1 = hr1->arr;
    std::size_t k = 0;
    for (std::size_t j = 0; j < z; j++)
    {
	std::size_t jn = j*n;
	for (std::size_t i = 0; i < m; i++)
	{
	    rr1[k] = rr[jn + pkk[i]];
	    k++;
	}
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
    auto rr = hr.arr;
    auto& rr1 = ar1->arr;
    rr1.resize(w);
    auto rr1p = rr1.data();
    if (m > 0)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn = j*n;
	    std::size_t k = rr[jn+pkk[0]];
	    for (std::size_t i = 1; i < m; i++)
		k = skk[i]*k + rr[jn+pkk[i]];
	    rr1p[k] += f;
	}
    else
	rr1p[0] = f*z;
    return ar1;
}

TransformRepa::TransformRepa() : _mapVarInt(0), derived(0), arr(0)
{
}

TransformRepa::~TransformRepa()
{
    delete _mapVarInt;
    delete derived;
    delete[] arr;
}

VarSizeUMap& Alignment::TransformRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<TransformRepa*>(this)->_mapVarInt = new VarSizeUMap(vectorVar.size());
	for (std::size_t i = 0; i < vectorVar.size(); i++)
	    const_cast<TransformRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// systemsTransformsTransformRepa_u :: System -> Transform -> TransformRepa
std::unique_ptr<TransformRepa> Alignment::systemsTransformsTransformRepa_u(const System& uu, const Transform& tt)
{
    auto tr = std::make_unique<TransformRepa>();
    auto ww = transformsDerived(tt);
    auto it = ww.begin();
    if (it == ww.end())
	return tr;
    tr->derived = new Variable(*it);
    auto zz = systemsVarsSetValue(uu,*tr->derived);
    auto w = zz.size();
    if (w > std::numeric_limits<unsigned char>::max())
	throw std::out_of_range("systemsTransformsTransformRepa_u");
    tr->valency = w;
    auto qq = transformsUnderlying(tt);
    auto n = qq->size();
    auto& vv = tr->vectorVar;
    vv.reserve(n);
    vv.insert(vv.end(), qq->begin(), qq->end());
    auto& sh = tr->shape;
    sh.reserve(n);
    std::vector<ValSizeUMap> mm(n+1);
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto xx = systemsVarsSetValue(uu,vv[i]);
	auto s = xx.size();
	sh.push_back(s);
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& x : xx)
	    yy.insert_or_assign(x, j++);
    }
    {
	auto xx = systemsVarsSetValue(uu, *tr->derived);
	auto s = xx.size();
	auto& yy = mm[n];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& x : xx)
	    yy.insert_or_assign(x, j++);
    }
    tr->arr = new unsigned char[sz];
    auto rr = tr->arr;
    for (auto& sc : tt.histogram_u().map_u())
    {
	auto& sm = sc.first.map_u();
	std::size_t j = 0;
	for (std::size_t i = 0; i < n; i++)
	    j = sh[i]*j + mm[i][sm.find(vv[i])->second];
	rr[j] = mm[n][sm.find(*tr->derived)->second];
    }
    return tr;
}

// systemsTransformRepasTransform_u :: System -> TransformRepa -> Transform
std::unique_ptr<Transform> Alignment::systemsTransformRepasTransform_u(const System& uu, const TransformRepa& tr)
{
    auto& vv = tr.vectorVar;
    auto n = vv.size();
    auto& sh = tr.shape;
    auto rr = tr.arr;
    if (!tr.derived || !n)
	return std::make_unique<Transform>();
    auto& w = *tr.derived;
    std::size_t sz = 1;
    std::vector<ValList> mm(n+1);
    for (std::size_t i = 0; i < n; i++)
    {
	auto xx = systemsVarsSetValue(uu, vv[i]);
	auto s = xx.size();
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	for (auto& x : xx)
	    yy.push_back(x);
    }
    {
	auto xx = systemsVarsSetValue(uu, w);
	auto s = xx.size();
	auto& yy = mm[n];
	yy.reserve(s);
	for (auto& x : xx)
	    yy.push_back(x);
    }
    auto tt = std::make_unique<Transform>();
    tt->derived_u().insert(w);
    auto& am = tt->histogram_u().map_u();
    am.reserve(sz);
    SizeList ii(n);
    for (std::size_t j = 0; j < sz; j++)
    {
	std::vector<VarValPair> ss;
	ss.reserve(n+1);
	for (std::size_t i = 0; i < n; i++)
	    ss.push_back(VarValPair(vv[i], mm[i][ii[i]]));
	ss.push_back(VarValPair(w, mm[n][rr[j]]));
	am.insert_or_assign(State(ss), Rational(1));
	for (int k = n - 1; k >= 0; k--)
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
    return tt;
}

// setVariablesListTransformRepasFudRepa_u :: Set.Set Variable -> V.Vector TransformRepa -> FudRepa
// cf listVariablesListTransformRepasSort :: V.Vector Variable -> V.Vector TransformRepa -> V.Vector
std::unique_ptr<FudRepa> Alignment::setVariablesListTransformRepasFudRepa_u(const VarUSet& vv, const TransformRepaPtrList& ff)
{
    auto fr = std::make_unique<FudRepa>();
    VarUSet vv1;
    vv1.reserve(vv.size() + ff.size());
    vv1.insert(vv.begin(), vv.end());
    TransformRepaPtrList ff1(ff);
    TransformRepaPtrList ff0;
    TransformRepaPtrList* ffa = &ff0;
    TransformRepaPtrList* ffb = &ff1;
    bool found = true;
    while (found)
    {
	found = false;
	VarUSet vv0;
	vv0.reserve(ffb->size());
	for (auto& tt : *ffb)
	{
	    bool layer = true;
	    for (auto& v : tt->vectorVar)
	    {
		auto it = vv1.find(v);
		layer = it != vv1.end();
		if (!layer)
		    break;
	    }
	    if (layer)
	    {
		vv0.insert(*tt->derived);
		if (!found)
		{
		    fr->layers.push_back(TransformRepaPtrList());
		    found = true;
		}
		fr->layers.back().push_back(tt);
	    }
	    else
		ffa->push_back(tt);
	}
	if (found)
	{
	    TransformRepaPtrList* ffc = ffa;
	    ffa = ffb;
	    ffb = ffc;
	    ffa->clear();
	    vv1.insert(vv0.begin(), vv0.end());
	}
    }
    return fr;
}

// systemsFudsFudRepa_u :: System -> Fud -> FudRepa
std::unique_ptr<FudRepa> Alignment::systemsFudsFudRepa_u(const System& uu, const Fud& ff)
{
    auto fund = fudsUnderlying;
    auto tttr = systemsTransformsTransformRepa_u;
    auto llfr = setVariablesListTransformRepasFudRepa_u;

    auto vv = fund(ff);
    TransformRepaPtrList ll;
    for (auto& tt : ff.list_u())
    {
	auto tr = tttr(uu, *tt);
	ll.push_back(std::move(tr));
    }
    return llfr(*vv, ll);
}

// systemsFudRepasFud_u :: System -> FudRepa -> Fud
std::unique_ptr<Fud> Alignment::systemsFudRepasFud_u(const System& uu, const FudRepa& fr)
{
    auto trtt = systemsTransformRepasTransform_u;

    auto ff = std::make_unique<Fud>();
    auto& mm = ff->list_u();
    for (auto& ll : fr.layers)
    {
	for (auto& tr : ll)
	{
	    auto tt = trtt(uu, *tr);
	    mm.push_back(std::move(tt));
	}
    }
    return ff;
}

// historyRepasFudRepasMultiply_u :: HistoryRepa -> FudRepa -> HistoryRepa
// cf historyRepasListTransformRepasApply_u :: HistoryRepa -> V.Vector TransformRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historyRepasFudRepasMultiply_u(const HistoryRepa& hr, const FudRepa& fr)
{
    auto hr1 = std::make_unique<HistoryRepa>();
    hr1->vectorVar = hr.vectorVar;
    hr1->size = hr.size;
    hr1->shape = hr.shape;
    auto& kk = hr1->vectorVar;
    auto& skk = hr1->shape;
    auto z = hr.size;
    auto n = hr.vectorVar.size();
    auto p = n;
    for (auto& ll : fr.layers)
	for (auto& tr : ll)
	{
	    kk.push_back(*tr->derived);
	    skk.push_back(tr->valency);
	    p++;
	}
    auto& mkk = hr1->mapVarInt();
    auto rr = hr.arr;
    hr1->arr = new unsigned char[z*p];
    auto rr1 = hr1->arr;
    for (std::size_t j = 0; j < z; j++)
    {
	std::size_t jn = j*n;
	std::size_t jp = j*p;
	for (std::size_t i = 0; i < n; i++)
	    rr1[jp+i] = rr[jn+i];
    }
    auto q = n;
    for (auto& ll : fr.layers)
	for (auto& tr : ll)
	{
	    auto& vv = tr->vectorVar;
	    auto& svv = tr->shape;
	    auto m = vv.size();
	    auto ar = tr->arr;
	    SizeList pkk;
	    for (auto& v : vv)
		pkk.push_back(mkk[v]);
	    if (m > 0)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t jp = j*p;
		    std::size_t k = rr1[jp + pkk[0]];
		    for (std::size_t i = 1; i < m; i++)
			k = svv[i] * k + rr1[jp + pkk[i]];
		    rr1[jp+q] = ar[k];
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		    rr1[j*p+q] = 0;
	    q++;
	}
    return hr1;
}

typedef std::shared_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> HistoryRepaPtrFudRepaPtrPairTreePtr;
typedef std::pair<HistoryRepaPtrFudRepaPtrPair, HistoryRepaPtrFudRepaPtrPairTreePtr> HistoryRepaPtrFudRepaPtrPairTreePtrPair;

// systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u :: Tree (State,Fud) -> Tree (HistoryRepa,FudRepa)
std::unique_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u(const System& uu, const Tree<StatePtrFudPtrPair>& rr)
{
    auto sshr = [](const System& uu, const State& ss)
    {
	History hh;
	hh.map_u().insert_or_assign(Id(1),ss);
	return systemsHistoriesHistoryRepa_u(uu, hh);
    };
    auto fffr = systemsFudsFudRepa_u;
    auto zzzr = systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u;

    auto tt = std::make_unique<Tree<HistoryRepaPtrFudRepaPtrPair>>();
    for (auto& pp : rr._list)
    {
	auto hr = sshr(uu, *pp.first._state);
	auto fr = fffr(uu, *pp.first._fud);
	HistoryRepaPtrFudRepaPtrPair mm(std::move(hr), std::move(fr));
	if (pp.second)
	{
	    auto qq = zzzr(uu, *pp.second);
	    tt->_list.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mm, std::move(qq)));
	}
	else
	{
	    auto nn = HistoryRepaPtrFudRepaPtrPairTreePtr(new Tree<HistoryRepaPtrFudRepaPtrPair>());
	    tt->_list.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mm, nn));
	}
    }
    return tt;
}

// systemsDecompFudsDecompFudRepa_u :: System -> DecompFud -> DecompFudRepa
std::unique_ptr<DecompFudRepa> Alignment::systemsDecompFudsDecompFudRepa_u(const System& uu, const DecompFud& df)
{
    auto sshr = [](const System& uu, const State& ss)
    {
	History hh;
	hh.map_u().insert_or_assign(Id(1), ss);
	return systemsHistoriesHistoryRepa_u(uu, hh);
    };
    auto fffr = systemsFudsFudRepa_u;
    auto zzzr = systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u;

    auto dr = std::make_unique<DecompFudRepa>();
    auto& rr = df.tree_u();
    auto& tt = dr->tree;
    for (auto& pp : rr._list)
    {
	auto hr = sshr(uu, *pp.first._state);
	auto fr = fffr(uu, *pp.first._fud);
	HistoryRepaPtrFudRepaPtrPair mm(std::move(hr), std::move(fr));
	if (pp.second)
	{
	    auto qq = zzzr(uu, *pp.second);
	    tt._list.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mm, std::move(qq)));
	}
	else
	{
	    auto nn = HistoryRepaPtrFudRepaPtrPairTreePtr(new Tree<HistoryRepaPtrFudRepaPtrPair>());
	    tt._list.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mm, nn));
	}
    }
    return dr;
}


typedef std::shared_ptr<Tree<StatePtrFudPtrPair>> StatePtrFudPtrPairTreePtr;
typedef std::pair<StatePtrFudPtrPair, StatePtrFudPtrPairTreePtr> StatePtrFudPtrPairTreePtrPair;

// systemsHistoryRepaFudRepaPairTreesStateFudPairTree_u :: Tree (HistoryRepa,FudRepa) -> Tree (State,Fud)
std::unique_ptr<Tree<StatePtrFudPtrPair>> systemsHistoryRepaFudRepaPairTreesStateFudPairTree_u(const System& uu, const Tree<HistoryRepaPtrFudRepaPtrPair>& rr)
{
    auto hrss = [](const System& uu, const HistoryRepa& hr)
    {
	auto hh = systemsHistoryRepasHistory_u(uu,hr);
	auto it = hh->map_u().begin();
	if (it != hh->map_u().end())
	    return std::make_unique<State>(it->second);
	return std::make_unique<State>(State());
    };
    auto frff = systemsFudRepasFud_u;
    auto zrzz = systemsHistoryRepaFudRepaPairTreesStateFudPairTree_u;

    auto tt = std::make_unique<Tree<StatePtrFudPtrPair>>();
    for (auto& pp : rr._list)
    {
	auto ss = hrss(uu, *pp.first._state);
	auto ff = frff(uu, *pp.first._fud);
	StatePtrFudPtrPair mm(std::move(ss), std::move(ff));
	if (pp.second)
	{
	    auto qq = zrzz(uu, *pp.second);
	    tt->_list.push_back(StatePtrFudPtrPairTreePtrPair(mm, std::move(qq)));
	}
	else
	{
	    auto nn = StatePtrFudPtrPairTreePtr(new Tree<StatePtrFudPtrPair>());
	    tt->_list.push_back(StatePtrFudPtrPairTreePtrPair(mm, nn));
	}
    }
    return tt;
}

// systemsDecompFudRepasDecompFud_u :: System -> DecompFudRepa -> DecompFud
std::unique_ptr<DecompFud> Alignment::systemsDecompFudRepasDecompFud_u(const System& uu, const DecompFudRepa& dr)
{
    auto hrss = [](const System& uu, const HistoryRepa& hr)
    {
	auto hh = systemsHistoryRepasHistory_u(uu, hr);
	auto it = hh->map_u().begin();
	if (it != hh->map_u().end())
	    return std::make_unique<State>(it->second);
	return std::make_unique<State>(State());
    };
    auto frff = systemsFudRepasFud_u;
    auto zrzz = systemsHistoryRepaFudRepaPairTreesStateFudPairTree_u;

    auto df = std::make_unique<DecompFud>();
    auto& tt = df->tree_u();
    auto& rr = dr.tree;
    for (auto& pp : rr._list)
    {
	auto ss = hrss(uu, *pp.first._state);
	auto ff = frff(uu, *pp.first._fud);
	StatePtrFudPtrPair mm(std::move(ss), std::move(ff));
	if (pp.second)
	{
	    auto qq = zrzz(uu, *pp.second);
	    tt._list.push_back(StatePtrFudPtrPairTreePtrPair(mm, std::move(qq)));
	}
	else
	{
	    auto nn = StatePtrFudPtrPairTreePtr(new Tree<StatePtrFudPtrPair>());
	    tt._list.push_back(StatePtrFudPtrPairTreePtrPair(mm, nn));
	}
    }
    return df;
}

