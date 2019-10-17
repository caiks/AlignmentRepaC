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

std::ostream& operator<<(std::ostream& out, const SystemRepa& ur)
{
    auto& ll = ur.listVarUCharPair;
    auto n = ll.size();
    out << "[";
    for (std::size_t i = 0; i < n; i++)
    {
	if (i) out << ",";
	out << "(" << ll[i].first << "," << (std::size_t)ll[i].second << ")";
    }
    out << "]";
    return out;
}

VarSizeUMap& Alignment::SystemRepa::mapVarSize() const
{
    auto n = listVarUCharPair.size();
    if (!_mapVarSize || _mapVarSize->size() < n)
    {
	if (_mapVarSize)
	    delete const_cast<SystemRepa*>(this)->_mapVarSize;
	const_cast<SystemRepa*>(this)->_mapVarSize = new VarSizeUMap(n);
	for (std::size_t i = 0; i < n; i++)
	    const_cast<SystemRepa*>(this)->_mapVarSize->insert_or_assign(listVarUCharPair[i].first, i);
    }
    return *_mapVarSize;
}

// systemsSystemRepa :: System -> SystemRepa
std::unique_ptr<SystemRepa> Alignment::systemsSystemRepa(const System& uu)
{
    auto ur = std::make_unique<SystemRepa>();
    auto& ll = ur->listVarUCharPair;
    auto mm = sorted(uu.map_u());
    ll.reserve(mm.size());
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    for (auto& vww : mm)
    {
	auto s = vww.second.size();
	if (s > ucmax)
	    throw std::out_of_range("systemsSystemRepa");
	ll.push_back(VarUCharPair(vww.first,(unsigned char)s));
    }
    return ur;
}

// systemsRepasSystem :: SystemRepa -> System
void Alignment::systemsRepasSystem(const SystemRepa& ur, System& uu)
{
    for (auto& vs : ur.listVarUCharPair)
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

HistogramRepa::HistogramRepa() : _mapVarInt(0), dimension(0), vectorVar(0), shape(0), arr(0)
{
}

HistogramRepa::~HistogramRepa()
{
    delete[] arr;
    delete[] shape;
    delete[] vectorVar;
    delete _mapVarInt;
}

std::ostream& operator<<(std::ostream& out, const HistogramRepa& ar)
{
    auto n = ar.dimension;
    auto vv = ar.vectorVar;
    auto sh = ar.shape;
    auto rr = ar.arr;
    out << "(" << n << ",[";
    for (std::size_t i = 0; i < n; i++)
    {
	if (i) out << ",";
	out << vv[i];
    }
    out << "],[";
    std::size_t v = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto s = sh[i];
	v *= s;
	if (i) out << ",";
	out << (std::size_t)s;
    }
    out << "],[";
    for (std::size_t j = 0; j < v; j++)
    {
	if (j) out << ",";
	out << rr[j];
    }
    out << "])";
    return out;
}


SizeSizeUMap& Alignment::HistogramRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistogramRepa*>(this)->_mapVarInt = new SizeSizeUMap(dimension);
	for (std::size_t i = 0; i < dimension; i++)
	    const_cast<HistogramRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// systemsHistogramsHistogramRepa_u :: System -> Histogram -> Maybe HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::systemsHistogramsHistogramRepa_u(const System& uu, const SystemRepa& ur, const Histogram& aa)
{
    auto& vvi = ur.mapVarSize();
    auto ww = histogramsSetVar(aa);
    VarList ww1(ww->begin(), ww->end());
    auto ar = std::make_unique<HistogramRepa>();
    ar->dimension = ww1.size();
    auto n = ar->dimension;
    ar->vectorVar = new std::size_t[n];
    auto vv = ar->vectorVar;
    ar->shape = new unsigned char[n];
    auto sh = ar->shape;
    std::vector<ValSizeUMap> mm(n);
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = ww1[i];
	vv[i] = vvi[v];
	auto xx = systemsVarsSetValue(uu,v);
	auto s = xx.size();
	sh[i] = s;
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& w : xx)
	    yy.insert_or_assign(w, j++);
    }
    ar->arr = new double[sz];
    auto rr = ar->arr;
    for (std::size_t j = 0; j < sz; j++)
	rr[j] = 0.0;
    for (auto& sc : aa.map_u())
    {
	auto& sm = sc.first.map_u();
	std::size_t j = 0;
	for (std::size_t i = 0; i < n; i++)
	    j = sh[i]*j + mm[i][sm.find(ww1[i])->second];
	rr[j] = sc.second.getDouble();
    }
    return ar;
}

// systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
std::unique_ptr<Histogram> Alignment::systemsHistogramRepasHistogram_u(const System& uu, const SystemRepa& ur, const HistogramRepa& ar)
{
    auto scalar = histogramScalar_u;

    auto& ivv = ur.listVarUCharPair;
    auto n = ar.dimension;
    auto vv = ar.vectorVar;
    auto sh = ar.shape;
    auto rr = ar.arr;
    if (!rr)
	return scalar(Rational(0));
    if (!n)
	return scalar(Rational(rr[0]));
    std::size_t sz = 1;
    std::vector<ValList> mm(n);
    unsigned char* ii = new unsigned char[n];
    VarList ww;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = ivv[vv[i]].first;
	ww.push_back(v);
	auto xx = systemsVarsSetValue(uu,v);
	auto& yy = mm[i];
	yy.reserve(xx.size());
	for (auto& w : xx)
	    yy.push_back(w);
	sz *= sh[i];
	ii[i] = 0;
    }
    auto aa = std::make_unique<Histogram>();
    aa->map_u().reserve(sz);
    for (std::size_t j = 0; j < sz; j++)
    {
	std::vector<VarValPair> ss;
	ss.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	{
	    VarValPair pp(ww[i], mm[i][ii[i]]);
	    ss.push_back(pp);
	}
	aa->map_u().insert_or_assign(State(ss),Rational(rr[j]));
	for (int k = n-1; k >= 0; k--)
	{
	    auto y = ii[k] + 1;
	    if (y == sh[k])
		ii[k] = 0;
	    else
	    {
		ii[k] = y;
		break;
	    }
	}
    }
    delete[] ii;
    return aa;
}

// setVarsHistogramRepasReduce_u :: Set.Set Variable -> HistogramRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistogramRepasReduce_u(std::size_t m, std::size_t* kk, const HistogramRepa& ar)
{
    auto n = ar.dimension;
    auto svv = ar.shape;
    auto& mvv = ar.mapVarInt();
    auto rvv = ar.arr;
    auto br = std::make_unique<HistogramRepa>();
    if (!rvv)
	return br;
    std::size_t v = 1;
    unsigned char* ivv = new unsigned char[n];
    for (std::size_t i = 0; i < n; i++)
    {
	v *= svv[i];
	ivv[i] = 0;
    }
    br->dimension = m;
    br->vectorVar = new std::size_t[m];
    auto vkk = br->vectorVar;
    br->shape = new unsigned char[m];
    auto skk = br->shape;
    std::size_t* pkk = new std::size_t[m];
    std::size_t w = 1;
    for (std::size_t i = 0; i < m; i++)
    {
	auto x = kk[i];
	vkk[i] = x;
	auto p = mvv[x];
	pkk[i] = p;
	auto s = svv[p];
	w *= s;
	skk[i] = s;
    }
    br->arr = new double[w];
    auto rkk = br->arr;
    for (std::size_t j = 0; j < w; j++)
	rkk[j] = 0.0;
    for (std::size_t j = 0; j < v; j++)
    {
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++)
	    k = skk[i]*k + ivv[pkk[i]];
	rkk[k] += rvv[j];
	for (int i = n-1; i >= 0; i--)
	{
	    auto y = ivv[i] + 1;
	    if (y == svv[i])
		ivv[i] = 0;
	    else
	    {
		ivv[i] = y;
		break;
	    }
	}
    }
    delete[] pkk;
    delete[] ivv;
    return br;
}

HistogramRepaRed::HistogramRepaRed() : _mapVarInt(0), dimension(0), vectorVar(0), shape(0), arr(0)
{
}

HistogramRepaRed::~HistogramRepaRed()
{
    delete[] arr;
    delete[] shape;
    delete[] vectorVar;
    delete _mapVarInt;
}

std::ostream& operator<<(std::ostream& out, const HistogramRepaRed& ar)
{
    auto n = ar.dimension;
    auto vv = ar.vectorVar;
    auto sh = ar.shape;
    auto rr = ar.arr;
    out << "(" << n << ",[";
    for (std::size_t i = 0; i < n; i++)
    {
	if (i) out << ",";
	out << vv[i];
    }
    out << "],[";
    std::size_t sz = 0;
    for (std::size_t i = 0; i < n; i++)
    {
	auto s = sh[i];
	sz += s;
	if (i) out << ",";
	out << (std::size_t)s;
    }
    out << "],[";
    for (std::size_t j = 0; j < sz; j++)
    {
	if (j) out << ",";
	out << rr[j];
    }
    out << "])";
    return out;
}


SizeSizeUMap& Alignment::HistogramRepaRed::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistogramRepaRed*>(this)->_mapVarInt = new SizeSizeUMap(dimension);
	for (std::size_t i = 0; i < dimension; i++)
	    const_cast<HistogramRepaRed*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// histogramRepasRed_u :: Double -> HistogramRepa -> HistogramRepaRed
std::unique_ptr<HistogramRepaRed> Alignment::histogramRepasRed_u(double z, const HistogramRepa& ar)
{
    auto n = ar.dimension;
    auto vv = ar.vectorVar;
    auto sh = ar.shape;
    auto rr = ar.arr;
    if (!n || !rr || z <= 0.0)
	return std::make_unique<HistogramRepaRed>();
    double f = 1.0 / z;
    auto pr = std::make_unique<HistogramRepaRed>();
    pr->dimension = n;
    pr->vectorVar = new std::size_t[n];
    auto vv1 = pr->vectorVar;
    pr->shape = new unsigned char[n];
    auto sh1 = pr->shape;
    std::size_t v = 1;
    std::size_t sz = 0;
    unsigned char* ii = new unsigned char[n];
    std::size_t* xx = new std::size_t[n];
    for (std::size_t i = 0; i < n; i++)
    {
	vv1[i] = vv[i];
	auto s = sh[i];
	sh1[i] = s;
	xx[i] = sz;
	sz += s;
	v *= s;
	ii[i] = 0;
    }
    pr->arr = new double[sz];
    auto rr1 = pr->arr;
    for (std::size_t j = 0; j < sz; j++)
	rr1[j] = 0.0;
    for (std::size_t j = 0; j < v; j++)
    {
	auto a = f * rr[j];
	for (std::size_t i = 0; i < n; i++)
	    rr1[xx[i]+ii[i]] += a;
	for (int i = n-1; i >= 0; i--)
	{
	    auto y = ii[i] + 1;
	    if (y == sh[i])
		ii[i] = 0;
	    else
	    {
		ii[i] = y;
		break;
	    }
	}
    }
    delete[] xx;
    delete[] ii;
    return pr;
}


//HistogramRepaVec::HistogramRepaVec() : _mapVarInt(0)
//{
//}
//
//HistogramRepaVec::~HistogramRepaVec()
//{
//    delete _mapVarInt;
//}
//
//VarSizeUMap& Alignment::HistogramRepaVec::mapVarInt() const
//{
//    if (!_mapVarInt)
//    {
//	const_cast<HistogramRepaVec*>(this)->_mapVarInt = new VarSizeUMap(vectorVar.size());
//	for (std::size_t i = 0; i < vectorVar.size(); i++)
//	    const_cast<HistogramRepaVec*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
//    }
//    return *_mapVarInt;
//}

HistoryRepa::HistoryRepa() : _mapVarInt(0), dimension(0), vectorVar(0), shape(0), size(0), arr(0)
{
}

HistoryRepa::~HistoryRepa()
{
    delete[] arr;
    delete[] shape;
    delete[] vectorVar;
    delete _mapVarInt;
}

std::ostream& operator<<(std::ostream& out, const HistoryRepa& hr)
{
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto sh = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    out << "(" << n << ",[";
    for (std::size_t i = 0; i < n; i++)
    {
	if (i) out << ",";
	out << vv[i];
    }
    out << "],[";
    for (std::size_t i = 0; i < n; i++)
    {
	auto s = sh[i];
	if (i) out << ",";
	out << (std::size_t)s;
    }
    out << "]," << z << ",[";
    for (std::size_t j = 0; j < n*z; j++)
    {
	if (j) out << ",";
	out << (std::size_t)rr[j];
    }
    out << "])";
    return out;
}


SizeSizeUMap& Alignment::HistoryRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<HistoryRepa*>(this)->_mapVarInt = new SizeSizeUMap(dimension);
	for (std::size_t i = 0; i < dimension; i++)
	    const_cast<HistoryRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}

// systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::systemsHistoriesHistoryRepa_u(const System& uu, const SystemRepa& ur, const History& hh)
{
    auto& vvi = ur.mapVarSize();
    auto ww = historiesSetVar(hh);
    VarList ww1(ww->begin(), ww->end());
    auto hr = std::make_unique<HistoryRepa>();
    hr->dimension = ww1.size();
    auto n = hr->dimension;
    hr->vectorVar = new std::size_t[n];
    auto vv = hr->vectorVar;
    hr->shape = new unsigned char[n];
    auto sh = hr->shape;
    hr->size = hh.map_u().size();
    auto z = hr->size;
    std::vector<ValSizeUMap> mm(n);
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = ww1[i];
	vv[i] = vvi[v];
	auto xx = systemsVarsSetValue(uu, v);
	auto s = xx.size();
	sh[i] = s;
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& w : xx)
	    yy.insert_or_assign(w, j++);
    }
    hr->arr = new unsigned char[z*n];
    auto rr = hr->arr;
    auto hm = sorted(hh.map_u());
    std::size_t j = 0;
    for (auto& is : hm)
    {
	auto& sm = is.second.map_u();
	for (std::size_t i = 0; i < n; i++)
	{
	    rr[j] = mm[i][sm.find(ww1[i])->second];
	    j++;
	}
    }
    return hr;
}

// systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
std::unique_ptr<History> Alignment::systemsHistoryRepasHistory_u(const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
{
    auto& ivv = ur.listVarUCharPair;
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto sh = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    std::vector<ValList> mm(n);
    VarList ww;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = ivv[vv[i]].first;
	ww.push_back(v);
	auto xx = systemsVarsSetValue(uu,v);
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
	    ss.push_back(VarValPair(ww[i], mm[i][rr[j*n+i]]));
	hm.insert_or_assign(Id(j+1), State(ss));
    }
    return hh;
}

// eventsHistoryRepasHistoryRepaSelection_u :: [Int] -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::eventsHistoryRepasHistoryRepaSelection_u(std::size_t z1, std::size_t* ll, const HistoryRepa& hr)
{
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto sh = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    auto hr1 = std::make_unique<HistoryRepa>();
    if (!rr)
	return hr1;
    hr1->dimension = n;
    hr1->vectorVar = new std::size_t[n];
    auto vv1 = hr1->vectorVar;
    hr1->shape = new unsigned char[n];
    auto sh1 = hr1->shape;
    for (std::size_t i = 0; i < n; i++)
    {
	vv1[i] = vv[i];
	sh1[i] = sh[i];
    }
    hr1->size = z1;
    hr1->arr = new unsigned char[z1*n];
    auto rr1 = hr1->arr;
    std::size_t k = 0;
    for (std::size_t p = 0; p < z1; p++)
    {
	std::size_t jn = ll[p]*n;
	for (std::size_t i = 0; i < n; i++)
	{
	    rr1[k] = rr[jn+i];
	    k++;
	}
    }
    return hr1;
}

// historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historyRepasHistoryRepasHistoryRepaSelection_u(const HistoryRepa& ss, const HistoryRepa& hr)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;

    auto n = hr.dimension;
    auto svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto rr = hr.arr;
    if (!rr)
	return std::make_unique<HistoryRepa>();
    auto m = ss.dimension;
    auto kk = ss.vectorVar;
    auto y = ss.size;
    auto rr1 = ss.arr;
    std::size_t* pkk = new std::size_t[m];
    for (std::size_t i = 0; i < m; i++)
	pkk[i] = mvv[kk[i]];
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
    delete[] pkk;
    return hrsel(ll.size(),ll.data(),hr);
}

// setVarsHistoryRepasHistoryRepaReduced_u :: Set.Set Variable -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::setVarsHistoryRepasHistoryRepaReduced_u(std::size_t m, std::size_t* kk, const HistoryRepa& hr)
{
    auto n = hr.dimension;
    auto svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto rr = hr.arr;
    auto hr1 = std::make_unique<HistoryRepa>();
    if (!rr)
	return hr1;
    hr1->dimension = m;
    hr1->vectorVar = new std::size_t[m];
    auto vkk = hr1->vectorVar;
    hr1->shape = new unsigned char[m];
    auto skk = hr1->shape;
    std::size_t* pkk = new std::size_t[m];
    for (std::size_t i = 0; i < m; i++)
    {
	auto x = kk[i];
	vkk[i] = x;
	auto p = mvv[x];
	pkk[i] = p;
	skk[i] = svv[p];
    }
    hr1->size = z;
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
    delete[] pkk;
    return hr1;
}

// setVarsHistoryRepasReduce_u :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistoryRepasReduce_u(double f, std::size_t m, std::size_t* kk, const HistoryRepa& hr)
{
    auto n = hr.dimension;
    auto svv = hr.shape;
    auto& mvv = hr.mapVarInt();
    auto z = hr.size;
    auto rr = hr.arr;
    auto ar = std::make_unique<HistogramRepa>();
    if (!rr)
	return ar;
    ar->dimension = m;
    ar->vectorVar = new std::size_t[m];
    auto vkk = ar->vectorVar;
    ar->shape = new unsigned char[m];
    auto skk = ar->shape;
    std::size_t w = 1;
    std::size_t* pkk = new std::size_t[m];
    for (std::size_t i = 0; i < m; i++)
    {
	auto x = kk[i];
	vkk[i] = x;
	auto p = mvv[x];
	pkk[i] = p;
	auto s = svv[p];
	w *= s;
	skk[i] = s;
    }
    ar->arr = new double[w];
    auto rr1 = ar->arr;
    for (std::size_t j = 0; j < w; j++)
	rr1[j] = 0.0;
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
    delete[] pkk;
    return ar;
}

TransformRepa::TransformRepa() : _mapVarInt(0), dimension(0), vectorVar(0), derived(0), valency(0), shape(0), arr(0)
{
}

TransformRepa::~TransformRepa()
{
    delete[] arr;
    delete[] shape;
    delete[] vectorVar;
    delete _mapVarInt;
}

std::ostream& operator<<(std::ostream& out, const TransformRepa& tr)
{
    auto n = tr.dimension;
    auto vv = tr.vectorVar;
    auto w = tr.derived;
    auto u = tr.valency;
    auto sh = tr.shape;
    auto rr = tr.arr;
    out << "(" << n << ",[";
    for (std::size_t i = 0; i < n; i++)
    {
	if (i) out << ",";
	out << vv[i];
    }
    out << "]," << w << ",[";
    std::size_t v = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto s = sh[i];
	v *= s;
	if (i) out << ",";
	out << (std::size_t)s;
    }
    out << "]," << (std::size_t)u << ",[";
    for (std::size_t j = 0; j < v; j++)
    {
	if (j) out << ",";
	out << (std::size_t)rr[j];
    }
    out << "])";
    return out;
}

SizeSizeUMap& Alignment::TransformRepa::mapVarInt() const
{
    if (!_mapVarInt)
    {
	const_cast<TransformRepa*>(this)->_mapVarInt = new SizeSizeUMap(dimension);
	for (std::size_t i = 0; i < dimension; i++)
	    const_cast<TransformRepa*>(this)->_mapVarInt->insert_or_assign(vectorVar[i], i);
    }
    return *_mapVarInt;
}


// systemsTransformsTransformRepa_u :: System -> Transform -> TransformRepa
std::unique_ptr<TransformRepa> Alignment::systemsTransformsTransformRepa_u(const System& uu, const SystemRepa& ur, const Transform& tt)
{
    auto& vvi = ur.mapVarSize();
    auto ww = transformsDerived(tt);
    auto wit = ww.begin();
    if (wit == ww.end())
	return std::make_unique<TransformRepa>();
    auto qq = transformsUnderlying(tt);
    VarList qq1(qq->begin(), qq->end());
    auto tr = std::make_unique<TransformRepa>();
    tr->dimension = qq1.size();
    auto n = tr->dimension;
    tr->vectorVar = new std::size_t[n];
    auto vv = tr->vectorVar;
    tr->shape = new unsigned char[n];
    auto sh = tr->shape;
    std::vector<ValSizeUMap> mm(n+1);
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = qq1[i];
	vv[i] = vvi[v];
	auto xx = systemsVarsSetValue(uu,v);
	auto s = xx.size();
	sh[i] = s;
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	std::size_t j = 0;
	for (auto& x : xx)
	    yy.insert_or_assign(x, j++);
    }
    {
	auto& v = *wit;
	tr->derived = vvi[v];
	auto xx = systemsVarsSetValue(uu, v);
	auto s = xx.size();
	tr->valency = s;
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
	    j = sh[i]*j + mm[i][sm.find(qq1[i])->second];
	rr[j] = mm[n][sm.find(*wit)->second];
    }
    return tr;
}

// systemsTransformRepasTransform_u :: System -> TransformRepa -> Transform
std::unique_ptr<Transform> Alignment::systemsTransformRepasTransform_u(const System& uu, const SystemRepa& ur, const TransformRepa& tr)
{
    auto& ivv = ur.listVarUCharPair;
    auto n = tr.dimension;
    auto vv = tr.vectorVar;
    auto sh = tr.shape;
    auto rr = tr.arr;
    if (!tr.arr || !n)
	return std::make_unique<Transform>();
    std::size_t sz = 1;
    std::vector<ValList> mm(n+1);
    VarList ww;
    for (std::size_t i = 0; i < n; i++)
    {
	auto& v = ivv[vv[i]].first;
	ww.push_back(v);
	auto xx = systemsVarsSetValue(uu, v);
	auto s = xx.size();
	sz *= s;
	auto& yy = mm[i];
	yy.reserve(s);
	for (auto& x : xx)
	    yy.push_back(x);
    }
    auto& w = ivv[tr.derived].first;
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
	    ss.push_back(VarValPair(ww[i], mm[i][ii[i]]));
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
std::unique_ptr<FudRepa> Alignment::setVariablesListTransformRepasFudRepa_u(const SizeUSet& vv, const TransformRepaPtrList& ff)
{
    auto fr = std::make_unique<FudRepa>();
    SizeUSet vv1;
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
	SizeUSet vv0;
	vv0.reserve(ffb->size());
	for (auto& tt : *ffb)
	{
	    bool layer = true;
	    auto n = tt->dimension;
	    auto ww = tt->vectorVar;
	    for (std::size_t i = 0; i < n; i++)
	    {
		auto it = vv1.find(ww[i]);
		layer = it != vv1.end();
		if (!layer)
		    break;
	    }
	    if (layer)
	    {
		vv0.insert(tt->derived);
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
std::unique_ptr<FudRepa> Alignment::systemsFudsFudRepa_u(const System& uu, const SystemRepa& ur, const Fud& ff)
{
    auto fund = fudsUnderlying;
    auto tttr = systemsTransformsTransformRepa_u;
    auto llfr = setVariablesListTransformRepasFudRepa_u;

    auto& vvi = ur.mapVarSize();
    auto vv = fund(ff);
    SizeUSet vv1;
    for (auto& v : *vv)
	vv1.insert(vvi[v]);
    TransformRepaPtrList ll;
    for (auto& tt : ff.list_u())
    {
	auto tr = tttr(uu, ur, *tt);
	ll.push_back(std::move(tr));
    }
    return llfr(vv1, ll);
}

// systemsFudRepasFud_u :: System -> FudRepa -> Fud
std::unique_ptr<Fud> Alignment::systemsFudRepasFud_u(const System& uu, const SystemRepa& ur, const FudRepa& fr)
{
    auto trtt = systemsTransformRepasTransform_u;

    auto ff = std::make_unique<Fud>();
    auto& mm = ff->list_u();
    for (auto& ll : fr.layers)
    {
	for (auto& tr : ll)
	{
	    auto tt = trtt(uu, ur, *tr);
	    mm.push_back(std::move(tt));
	}
    }
    return ff;
}

// historyRepasFudRepasMultiply_u :: HistoryRepa -> FudRepa -> HistoryRepa
// cf historyRepasListTransformRepasApply_u :: HistoryRepa -> V.Vector TransformRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historyRepasFudRepasMultiply_u(const HistoryRepa& hr, const FudRepa& fr)
{
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto svv = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    auto hr1 = std::make_unique<HistoryRepa>();
    auto n1 = n;
    for (auto& ll : fr.layers)
	n1 += ll.size();
    hr1->dimension = n1;
    hr1->vectorVar = new std::size_t[n1];
    auto kk = hr1->vectorVar;
    hr1->shape = new unsigned char[n1];
    auto skk = hr1->shape;
    for (std::size_t i = 0; i < n; i++)
    {
	kk[i] = vv[i];
	skk[i] = svv[i];
    }
    std::size_t mmax = 1;
    auto p = n;
    for (auto& ll : fr.layers)
	for (auto& tr : ll)
	{
	    kk[p] = tr->derived;
	    skk[p] = tr->valency;
	    p++;
	    if (mmax < tr->dimension)
		mmax = tr->dimension;
	}
    hr1->size = z;
    auto& mkk = hr1->mapVarInt();
    hr1->arr = new unsigned char[z*p];
    auto rr1 = hr1->arr;
    for (std::size_t j = 0; j < z; j++)
    {
	std::size_t jn = j*n;
	std::size_t jp = j*p;
	for (std::size_t i = 0; i < n; i++)
	    rr1[jp+i] = rr[jn+i];
    }
    std::size_t* pkk = new std::size_t[mmax];
    auto q = n;
    for (auto& ll : fr.layers)
	for (auto& tr : ll)
	{
	    auto m = tr->dimension;
	    auto ww = tr->vectorVar;
	    auto sh = tr->shape;
	    auto ar = tr->arr;
	    for (std::size_t i = 0; i < m; i++)
		pkk[i] = mkk[ww[i]];
	    if (m > 0)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t jp = j*p;
		    std::size_t k = rr1[jp + pkk[0]];
		    for (std::size_t i = 1; i < m; i++)
			k = sh[i] * k + rr1[jp + pkk[i]];
		    rr1[jp+q] = ar[k];
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		    rr1[j*p+q] = 0;
	    q++;
	}
    delete[] pkk;
    return hr1;
}

typedef std::shared_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> HistoryRepaPtrFudRepaPtrPairTreePtr;
typedef std::pair<HistoryRepaPtrFudRepaPtrPair, HistoryRepaPtrFudRepaPtrPairTreePtr> HistoryRepaPtrFudRepaPtrPairTreePtrPair;

// systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u :: Tree (State,Fud) -> Tree (HistoryRepa,FudRepa)
std::unique_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u(const System& uu, const SystemRepa& ur, const Tree<StatePtrFudPtrPair>& rr)
{
    auto sshr = [](const System& uu, const SystemRepa& ur, const State& ss)
    {
	History hh;
	hh.map_u().insert_or_assign(Id(1),ss);
	return systemsHistoriesHistoryRepa_u(uu, ur, hh);
    };
    auto fffr = systemsFudsFudRepa_u;
    auto zzzr = systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u;

    auto tt = std::make_unique<Tree<HistoryRepaPtrFudRepaPtrPair>>();
    for (auto& pp : rr._list)
    {
	auto hr = sshr(uu, ur, *pp.first._state);
	auto fr = fffr(uu, ur, *pp.first._fud);
	HistoryRepaPtrFudRepaPtrPair mm(std::move(hr), std::move(fr));
	if (pp.second)
	{
	    auto qq = zzzr(uu, ur, *pp.second);
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
std::unique_ptr<DecompFudRepa> Alignment::systemsDecompFudsDecompFudRepa_u(const System& uu, const SystemRepa& ur, const DecompFud& df)
{
    auto sshr = [](const System& uu, const SystemRepa& ur, const State& ss)
    {
	History hh;
	hh.map_u().insert_or_assign(Id(1), ss);
	return systemsHistoriesHistoryRepa_u(uu, ur, hh);
    };
    auto fffr = systemsFudsFudRepa_u;
    auto zzzr = systemsStateFudPairTreesHistoryRepaFudRepaPairTree_u;

    auto dr = std::make_unique<DecompFudRepa>();
    auto& rr = df.tree_u();
    auto& tt = dr->tree;
    for (auto& pp : rr._list)
    {
	auto hr = sshr(uu, ur, *pp.first._state);
	auto fr = fffr(uu, ur, *pp.first._fud);
	HistoryRepaPtrFudRepaPtrPair mm(std::move(hr), std::move(fr));
	if (pp.second)
	{
	    auto qq = zzzr(uu, ur, *pp.second);
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
std::unique_ptr<Tree<StatePtrFudPtrPair>> systemsHistoryRepaFudRepaPairTreesStateFudPairTree_u(const System& uu, const SystemRepa& ur, const Tree<HistoryRepaPtrFudRepaPtrPair>& rr)
{
    auto hrss = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
    {
	auto hh = systemsHistoryRepasHistory_u(uu,ur,hr);
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
	auto ss = hrss(uu, ur, *pp.first._state);
	auto ff = frff(uu, ur, *pp.first._fud);
	StatePtrFudPtrPair mm(std::move(ss), std::move(ff));
	if (pp.second)
	{
	    auto qq = zrzz(uu, ur, *pp.second);
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
std::unique_ptr<DecompFud> Alignment::systemsDecompFudRepasDecompFud_u(const System& uu, const SystemRepa& ur, const DecompFudRepa& dr)
{
    auto hrss = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
    {
	auto hh = systemsHistoryRepasHistory_u(uu,ur,hr);
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
	auto ss = hrss(uu, ur, *pp.first._state);
	auto ff = frff(uu, ur, *pp.first._fud);
	StatePtrFudPtrPair mm(std::move(ss), std::move(ff));
	if (pp.second)
	{
	    auto qq = zrzz(uu, ur, *pp.second);
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

