#include "AlignmentApprox.h"
#include "AlignmentRepa.h"

#include <thread>
#include <cstring>
#include <iostream>

using namespace Alignment;

#define ECHO(x) std::cout << #x << std::endl; x
#define EVAL(x) std::cout << #x << ": " << (x) << std::endl
#define EVALL(x) std::cout << #x << ": " << std::endl << (x) << std::endl
#define TRUTH(x) std::cout << #x << ": " << ((x) ? "true" : "false") << std::endl

SystemRepa::SystemRepa() : _mapVarSize(0)
{
}

SystemRepa::~SystemRepa()
{
	delete _mapVarSize;
}

std::ostream& operator<<(std::ostream& out, const SystemRepa& ur)
{
	auto& ll = ur.listVarSizePair;
	auto n = ll.size();
	out << "[";
	for (std::size_t i = 0; i < n; i++)
	{
		if (i) out << ",";
		out << "(" << *ll[i].first << "," << ll[i].second << ")";
	}
	out << "]";
	return out;
}

VarSizeUMap& Alignment::SystemRepa::mapVarSize() const
{
	auto n = listVarSizePair.size();
	if (!_mapVarSize || _mapVarSize->size() < n)
	{
		if (_mapVarSize)
			delete const_cast<SystemRepa*>(this)->_mapVarSize;
		const_cast<SystemRepa*>(this)->_mapVarSize = new VarSizeUMap(n);
		for (std::size_t i = 0; i < n; i++)
			const_cast<SystemRepa*>(this)->_mapVarSize->insert_or_assign(*listVarSizePair[i].first, i);
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
		ll.push_back(VarSizePair(std::make_shared<Variable>(vww.first), vww.second.size()));
	return ur;
}

// systemsRepasSystem :: SystemRepa -> System
void Alignment::systemsRepasSystem(const SystemRepa& ur, System& uu)
{
	for (auto& vs : ur.listVarSizePair)
	{
		auto it = uu.map_u().find(*vs.first);
		if (it == uu.map_u().end())
		{
			ValSet ww;
			for (int i = 0; i < vs.second; i++)
				ww.insert(Value(i));
			uu.map_u().insert_or_assign(*vs.first, ww);
		}
	}
}

HistogramRepa::HistogramRepa() : _mapVarInt(0), dimension(0), vectorVar(0), shape(0), arr(0)
{
}

HistogramRepa::HistogramRepa(double z) : _mapVarInt(0), dimension(0), vectorVar(0), shape(0)
{
	arr = new double[1];
	arr[0] = z;
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
		out << s;
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

double Alignment::HistogramRepa::size() const
{
	auto n = dimension;
	auto sh = shape;
	auto rr = arr;
	std::size_t v = 1;
	for (std::size_t i = 0; i < n; i++)
		v *= sh[i];
	double a = 0.0;
	for (std::size_t j = 0; j < v; j++)
		a += rr[j];
	return a;
}

double Alignment::HistogramRepa::facLn() const
{
	auto n = dimension;
	auto sh = shape;
	auto rr = arr;
	std::size_t v = 1;
	for (std::size_t i = 0; i < n; i++)
		v *= sh[i];
	double a = 0.0;
	for (std::size_t j = 0; j < v; j++)
		a += alngam(rr[j] + 1.0);
	return a;
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
	ar->shape = new std::size_t[n];
	auto sh = ar->shape;
	std::vector<ValSizeUMap> mm(n);
	std::size_t sz = 1;
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = ww1[i];
		vv[i] = vvi[v];
		auto xx = systemsVarsSetValue(uu, v);
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
			j = sh[i] * j + mm[i][sm.find(ww1[i])->second];
		rr[j] = sc.second.getDouble();
	}
	return ar;
}

// systemsHistogramRepasHistogram_u :: System -> HistogramRepa -> Maybe Histogram
std::unique_ptr<Histogram> Alignment::systemsHistogramRepasHistogram_u(const System& uu, const SystemRepa& ur, const HistogramRepa& ar)
{
	auto scalar = histogramScalar_u;

	auto& ivv = ur.listVarSizePair;
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
	std::size_t* ii = new std::size_t[n];
	VarList ww;
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = *ivv[vv[i]].first;
		ww.push_back(v);
		auto xx = systemsVarsSetValue(uu, v);
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
		aa->map_u().insert_or_assign(State(ss), Rational(rr[j]));
		for (long long k = n - 1; k >= 0; k--)
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

// setVarsHistogramRepasReduce_u :: [VariableRepa] -> HistogramRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistogramRepasReduce_u(std::size_t m, const std::size_t* kk, const HistogramRepa& ar)
{
	auto n = ar.dimension;
	auto svv = ar.shape;
	auto& mvv = ar.mapVarInt();
	auto rvv = ar.arr;
	auto br = std::make_unique<HistogramRepa>();
	if (!rvv)
		return br;
	std::size_t v = 1;
	std::size_t* ivv = new std::size_t[n];
	for (std::size_t i = 0; i < n; i++)
	{
		v *= svv[i];
		ivv[i] = 0;
	}
	br->dimension = m;
	br->vectorVar = new std::size_t[m];
	auto vkk = br->vectorVar;
	br->shape = new std::size_t[m];
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
			k = skk[i] * k + ivv[pkk[i]];
		rkk[k] += rvv[j];
		for (long long i = n - 1; i >= 0; i--)
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

// histogramRepasEntropy :: HistogramRepa -> Double
double Alignment::histogramRepasEntropy(const HistogramRepa& ar)
{
	auto n = ar.dimension;
	auto sh = ar.shape;
	auto rr = ar.arr;
	std::size_t v = 1;
	for (std::size_t i = 0; i < n; i++)
		v *= sh[i];
	double z = 0.0;
	for (std::size_t j = 0; j < v; j++)
		z += rr[j];
	double e = 0.0;
	for (std::size_t j = 0; j < v; j++)
	{
		auto a = rr[j] / z;
		if (a > 0.0)
			e -= a * log(a);
	}
	return e;
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
		out << s;
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

// histogramRepasRed :: Double -> HistogramRepa -> HistogramRepaRed
std::unique_ptr<HistogramRepaRed> Alignment::histogramRepasRed(double z, const HistogramRepa& ar)
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
	pr->shape = new std::size_t[n];
	auto sh1 = pr->shape;
	std::size_t v = 1;
	std::size_t sz = 0;
	std::size_t* ii = new std::size_t[n];
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
			rr1[xx[i] + ii[i]] += a;
		for (long long i = n - 1; i >= 0; i--)
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

// histogramRepaRedsIndependent :: Double -> HistogramRepaRed -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::histogramRepaRedsIndependent(double z, const HistogramRepaRed& pr)
{
	if (z < 0.0)
		return std::make_unique<HistogramRepa>(0.0);
	auto n = pr.dimension;
	auto vv = pr.vectorVar;
	auto sh = pr.shape;
	auto rr = pr.arr;
	if (!n || !rr)
		return std::make_unique<HistogramRepa>(z);
	double f = 1.0 / z;
	auto ar = std::make_unique<HistogramRepa>();
	ar->dimension = n;
	ar->vectorVar = new std::size_t[n];
	auto vv1 = ar->vectorVar;
	ar->shape = new std::size_t[n];
	auto sh1 = ar->shape;
	std::size_t v = 1;
	std::size_t sz = 0;
	std::size_t* ii = new std::size_t[n];
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
	ar->arr = new double[v];
	auto rr1 = ar->arr;
	for (std::size_t j = 0; j < v; j++)
	{
		auto a = z;
		for (std::size_t i = 0; i < n; i++)
			a *= rr[xx[i] + ii[i]];
		rr1[j] = a;
		for (long long i = n - 1; i >= 0; i--)
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
	return ar;
}

// histogramRepaRedsListEntropy :: HistogramRepaRed -> [(Double,VariableRepa)]
std::unique_ptr<DoubleSizePairList> Alignment::histogramRepaRedsListEntropy(const HistogramRepaRed& pr)
{
	auto n = pr.dimension;
	auto vv = pr.vectorVar;
	auto sh = pr.shape;
	auto rr = pr.arr;
	auto ee = std::make_unique<DoubleSizePairList>();
	if (!n || !rr)
		return ee;
	ee->reserve(n);
	std::size_t j = 0;
	for (std::size_t i = 0; i < n; i++)
	{
		auto s = sh[i];
		double e = 0.0;
		for (std::size_t k = 0; k < s; k++)
		{
			auto a = rr[j];
			e -= a * log(a);
			j++;
		}
		if (e > 0.0)
			ee->push_back(DoubleSizePair(e, vv[i]));
	}
	return ee;
}



HistoryRepa::HistoryRepa() : _mapVarInt(0), dimension(0), vectorVar(0), shape(0), size(0), evient(false), arr(0)
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
		out << s;
	}
	out << "]," << z << "," << (hr.evient ? "true" : "false") << ",[";
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

void Alignment::HistoryRepa::transpose()
{
	if (!dimension || !size || !arr)
		return;
	auto n = dimension;
	auto z = size;
	auto arr1 = new unsigned char[z*n];
	if (evient)
	{
		std::size_t jn;
		for (std::size_t j = 0; j < z; j++)
		{
			jn = j*n;
			for (std::size_t i = 0; i < n; i++)
				arr1[i*z + j] = arr[jn + i];
		}
		evient = false;
	}
	else
	{
		std::size_t iz;
		for (std::size_t i = 0; i < n; i++)
		{
			iz = i*z;
			for (std::size_t j = 0; j < z; j++)
				arr1[j*n + i] = arr[iz + j];
		}
		evient = true;
	}
	auto arr2 = arr;
	arr = arr1;
	delete[] arr2;
	return;
}

void Alignment::HistoryRepa::reframe_u(const SizeSizeUMap& nn)
{
	auto n = dimension;
	auto vv = vectorVar;
	bool altered = false;
	for (std::size_t i = 0; i < n; i++)
	{
		auto it = nn.find(vv[i]);
		if (it != nn.end())
		{
			vv[i] = it->second;
			altered = true;
		}
	}
	if (!_mapVarInt && altered)
	{
		delete _mapVarInt;
		_mapVarInt = 0;
	}
}

// systemsHistoriesHistoryRepa_u :: System -> History -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::systemsHistoriesHistoryRepa_u(const System& uu, const SystemRepa& ur, const History& hh, unsigned char evient)
{
	auto& vvi = ur.mapVarSize();
	auto ww = historiesSetVar(hh);
	VarList ww1(ww->begin(), ww->end());
	auto hr = std::make_unique<HistoryRepa>();
	hr->dimension = ww1.size();
	auto n = hr->dimension;
	hr->vectorVar = new std::size_t[n];
	auto vv = hr->vectorVar;
	hr->shape = new std::size_t[n];
	auto sh = hr->shape;
	hr->size = hh.map_u().size();
	auto z = hr->size;
	std::vector<ValSizeUMap> mm(n);
	unsigned char ucmax = std::numeric_limits<unsigned char>::max();
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = ww1[i];
		vv[i] = vvi[v];
		auto xx = systemsVarsSetValue(uu, v);
		auto s = xx.size();
		if (s > ucmax + 1)
			throw std::out_of_range("systemsHistoriesHistoryRepa_u");
		sh[i] = s;
		auto& yy = mm[i];
		yy.reserve(s);
		std::size_t j = 0;
		for (auto& w : xx)
			yy.insert_or_assign(w, j++);
	}
	hr->evient = true;
	hr->arr = new unsigned char[z*n];
	auto rr = hr->arr;
	auto hm = sorted(hh.map_u());
	std::size_t j = 0;
	for (auto& is : hm)
	{
		auto& sm = is.second.map_u();
		for (std::size_t i = 0; i < n; i++)
		{
			rr[j] = (unsigned char)mm[i][sm.find(ww1[i])->second];
			j++;
		}
	}
	if (!evient)
		hr->transpose();
	return hr;
}


// systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History
std::unique_ptr<History> Alignment::systemsHistoryRepasHistory_u(const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
{
	auto& ivv = ur.listVarSizePair;
	auto n = hr.dimension;
	auto vv = hr.vectorVar;
	auto sh = hr.shape;
	auto z = hr.size;
	auto rr = hr.arr;
	std::vector<ValList> mm(n);
	VarList ww;
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = *ivv[vv[i]].first;
		ww.push_back(v);
		auto xx = systemsVarsSetValue(uu, v);
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
			ss.push_back(VarValPair(ww[i], mm[i][rr[hr.evient ? j*n + i : i*z + j]]));
		hm.insert_or_assign(Id((int)j + 1), State(ss));
	}
	return hh;
}

// eventsHistoryRepasHistoryRepaSelection_u :: [Int] -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::eventsHistoryRepasHistoryRepaSelection_u(std::size_t z1, const std::size_t* ll, const HistoryRepa& hr)
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
	hr1->shape = new std::size_t[n];
	auto sh1 = hr1->shape;
	for (std::size_t i = 0; i < n; i++)
	{
		vv1[i] = vv[i];
		sh1[i] = sh[i];
	}
	hr1->size = z1;
	hr1->evient = hr.evient;
	hr1->arr = new unsigned char[z1*n];
	auto rr1 = hr1->arr;
	std::size_t k = 0;
	if (hr.evient)
		for (std::size_t p = 0; p < z1; p++)
		{
			std::size_t jn = ll[p] * n;
			for (std::size_t i = 0; i < n; i++)
			{
				rr1[k] = rr[jn + i];
				k++;
			}
		}
	else
		for (std::size_t i = 0; i < n; i++)
		{
			std::size_t iz = i*z;
			for (std::size_t p = 0; p < z1; p++)
			{
				rr1[k] = rr[iz + ll[p]];
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
	if (ss.evient && hr.evient)
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
	else if (ss.evient && !hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			bool any = false;
			for (std::size_t k = 0; !any && k < y; k++)
			{
				std::size_t km = k*m;
				bool all = true;
				for (std::size_t i = 0; all && i < m; i++)
					all = rr1[km + i] == rr[pkk[i] * z + j];
				any = all;
			}
			if (any)
				ll.push_back(j);
		}
	else if (!ss.evient && hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			bool any = false;
			for (std::size_t k = 0; !any && k < y; k++)
			{
				bool all = true;
				for (std::size_t i = 0; all && i < m; i++)
					all = rr1[i*y + k] == rr[jn + pkk[i]];
				any = all;
			}
			if (any)
				ll.push_back(j);
		}
	else
		for (std::size_t j = 0; j < z; j++)
		{
			bool any = false;
			for (std::size_t k = 0; !any && k < y; k++)
			{
				bool all = true;
				for (std::size_t i = 0; all && i < m; i++)
					all = rr1[i*y + k] == rr[pkk[i] * z + j];
				any = all;
			}
			if (any)
				ll.push_back(j);
		}
	delete[] pkk;
	return hrsel(ll.size(), ll.data(), hr);
}

// setVarsHistoryRepasHistoryRepaReduced_u :: [VariableRepa] -> HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::setVarsHistoryRepasHistoryRepaReduced_u(std::size_t m, const std::size_t* kk, const HistoryRepa& hr)
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
	hr1->shape = new std::size_t[m];
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
	hr1->evient = hr.evient;
	hr1->arr = new unsigned char[z*m];
	auto rr1 = hr1->arr;
	std::size_t k = 0;
	if (hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			for (std::size_t i = 0; i < m; i++)
			{
				rr1[k] = rr[jn + pkk[i]];
				k++;
			}
		}
	else
		for (std::size_t i = 0; i < m; i++)
		{
			std::size_t iz = pkk[i] * z;
			for (std::size_t j = 0; j < z; j++)
			{
				rr1[k] = rr[iz + j];
				k++;
			}
		}
	delete[] pkk;
	return hr1;
}

// setVarsHistoryRepasReduce_u :: Double -> [VariableRepa] -> HistoryRepa -> HistogramRepa
std::unique_ptr<HistogramRepa> Alignment::setVarsHistoryRepasReduce_u(double f, std::size_t m, const std::size_t* kk, const HistoryRepa& hr)
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
	ar->shape = new std::size_t[m];
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
	if (m > 0 && hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			std::size_t k = rr[jn + pkk[0]];
			for (std::size_t i = 1; i < m; i++)
				k = skk[i] * k + rr[jn + pkk[i]];
			rr1[k] += f;
		}
	else if (m > 0 && !hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t k = rr[pkk[0] * z + j];
			for (std::size_t i = 1; i < m; i++)
				k = skk[i] * k + rr[pkk[i] * z + j];
			rr1[k] += f;
		}
	else
		rr1[0] = f*z;
	delete[] pkk;
	return ar;
}

// historyRepasRed :: HistoryRepa -> HistogramRepaRed
std::unique_ptr<HistogramRepaRed> Alignment::historyRepasRed(const HistoryRepa& hr)
{
	auto n = hr.dimension;
	auto vv = hr.vectorVar;
	auto sh = hr.shape;
	auto z = hr.size;
	auto rr = hr.arr;
	if (!n || !rr || !z)
		return std::make_unique<HistogramRepaRed>();
	double f = 1.0 / z;
	auto pr = std::make_unique<HistogramRepaRed>();
	pr->dimension = n;
	pr->vectorVar = new std::size_t[n];
	auto vv1 = pr->vectorVar;
	pr->shape = new std::size_t[n];
	auto sh1 = pr->shape;
	std::size_t sz = 0;
	std::size_t* xx = new std::size_t[n];
	for (std::size_t i = 0; i < n; i++)
	{
		vv1[i] = vv[i];
		auto s = sh[i];
		sh1[i] = s;
		xx[i] = sz;
		sz += s;
	}
	pr->arr = new double[sz];
	auto rr1 = pr->arr;
	for (std::size_t j = 0; j < sz; j++)
		rr1[j] = 0.0;
	if (hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			for (std::size_t i = 0; i < n; i++)
				rr1[xx[i] + rr[jn + i]] += f;
		}
	else
		for (std::size_t i = 0; i < n; i++)
		{
			std::size_t iz = i*z;
			for (std::size_t j = 0; j < z; j++)
				rr1[xx[i] + rr[iz + j]] += f;
		}
	delete[] xx;
	return pr;
}

// setVarsHistoryRepasRed_u :: Double -> [VariableRepa] -> HistoryRepa -> HistogramRepaRed
std::unique_ptr<HistogramRepaRed> Alignment::setVarsHistoryRepasRed_u(std::size_t m, const std::size_t* kk, const HistoryRepa& hr)
{
	auto n = hr.dimension;
	auto& mvv = hr.mapVarInt();
	auto sh = hr.shape;
	auto z = hr.size;
	auto rr = hr.arr;
	if (!n || !rr || !z)
		return std::make_unique<HistogramRepaRed>();
	double f = 1.0 / z;
	auto pr = std::make_unique<HistogramRepaRed>();
	pr->dimension = m;
	pr->vectorVar = new std::size_t[m];
	auto vv1 = pr->vectorVar;
	pr->shape = new std::size_t[m];
	auto sh1 = pr->shape;
	std::size_t sz = 0;
	std::size_t* pkk = new std::size_t[m];
	std::size_t* xx = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
	{
		auto x = kk[i];
		vv1[i] = x;
		auto p = mvv[x];
		pkk[i] = p;
		auto s = sh[p];
		sh1[i] = s;
		xx[i] = sz;
		sz += s;
	}
	pr->arr = new double[sz];
	auto rr1 = pr->arr;
	for (std::size_t j = 0; j < sz; j++)
		rr1[j] = 0.0;
	if (hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			for (std::size_t i = 0; i < m; i++)
				rr1[xx[i] + rr[jn + pkk[i]]] += f;
		}
	else
		for (std::size_t i = 0; i < m; i++)
		{
			std::size_t iz = pkk[i]*z;
			for (std::size_t j = 0; j < z; j++)
				rr1[xx[i] + rr[iz + j]] += f;
		}
	delete[] xx;
	delete[] pkk;
	return pr;
}


// vectorHistoryRepasConcat_u :: V.Vector HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::vectorHistoryRepasConcat_u(const HistoryRepaPtrList& ll)
{
	auto m = ll.size();
	auto hr1 = std::make_unique<HistoryRepa>();
	if (!m)
		return hr1;
	auto& hr0 = *ll[0];
	if (!hr0.arr)
		return hr1;
	auto n = hr0.dimension;
	auto vv = hr0.vectorVar;
	auto sh = hr0.shape;
	std::size_t z1 = 0;
	for (auto& hr : ll)
		z1 += hr->size;
	hr1->dimension = n;
	hr1->vectorVar = new std::size_t[n];
	auto vv1 = hr1->vectorVar;
	hr1->shape = new std::size_t[n];
	auto sh1 = hr1->shape;
	for (std::size_t i = 0; i < n; i++)
	{
		vv1[i] = vv[i];
		sh1[i] = sh[i];
	}
	hr1->size = z1;
	hr1->evient = true;
	hr1->arr = new unsigned char[z1*n];
	auto rr1 = hr1->arr;
	std::size_t y = 0;
	for (std::size_t k = 0; k < m; k++)
	{
		std::size_t yn = y*n;
		auto& hr = *ll[k];
		auto z = hr.size;
		auto rr = hr.arr;
		if (hr.evient)
			for (std::size_t j = 0; j < z; j++)
			{
				std::size_t jn = j*n;
				for (std::size_t i = 0; i < n; i++)
					rr1[yn + jn + i] = rr[jn + i];
			}
		else
			for (std::size_t i = 0; i < n; i++)
			{
				std::size_t iz = i*z;
				for (std::size_t j = 0; j < z; j++)
					rr1[yn + j*n + i] = rr[iz + j];
			}
		y += z;
	}
	if (!hr0.evient)
		hr1->transpose();
	return hr1;
}

// vectorHistoryRepasJoin_u :: V.Vector HistoryRepa -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::vectorHistoryRepasJoin_u(const HistoryRepaPtrList& ll)
{
	auto hr1 = std::make_unique<HistoryRepa>();
	if (!ll.size())
		return hr1;
	auto& hr0 = *ll[0];
	if (!hr0.arr)
		return hr1;
	auto z = hr0.size;
	std::size_t n1 = 0;
	for (auto& hr : ll)
		n1 += hr->dimension;
	hr1->dimension = n1;
	hr1->vectorVar = new std::size_t[n1];
	auto vv1 = hr1->vectorVar;
	hr1->shape = new std::size_t[n1];
	auto sh1 = hr1->shape;
	hr1->size = z;
	hr1->evient = false;
	hr1->arr = new unsigned char[z*n1];
	auto rr1 = hr1->arr;
	std::size_t k = 0;
	for (auto& hr : ll)
	{
		auto n = hr->dimension;
		auto vv = hr->vectorVar;
		auto sh = hr->shape;
		auto rr = hr->arr;
		for (std::size_t i = 0; i < n; i++)
		{
			vv1[k] = vv[i];
			sh1[k] = sh[i];
			std::size_t kz = k*z;
			std::size_t iz = i*z;
			for (std::size_t j = 0; j < z; j++)
				rr1[kz + j] = rr[hr->evient ? j*n + i : iz + j];
			k++;
		}
	}
	return hr1;
}


HistorySparse::HistorySparse() : size(0), vectorDimension(0), arr(0)
{
}

HistorySparse::~HistorySparse()
{
	delete[] arr;
	delete[] vectorDimension;
}

std::ostream& operator<<(std::ostream& out, const HistorySparse& hr)
{
	auto z = hr.size;
	out << "(" << z;
	if (z)
	{
		auto nn = hr.vectorDimension;
		auto rr = hr.arr;
		std::size_t k = 0;
		for (std::size_t j = 0; j < z; j++)
		{
			out << ",[";
			auto n = nn[j];
			for (std::size_t i = 0; i < n; i++)
			{
				if (i) out << ",";
				out << rr[k];
				k++;
			}
			out << "]";
		}
	}
	out << ")";
	return out;
}

// historyRepasHistorySparse :: HistoryRepa -> HistorySparse
std::unique_ptr<HistorySparse> Alignment::historyRepasHistorySparse(const HistoryRepa& hr)
{
	auto n = hr.dimension;
	auto vv = hr.vectorVar;
	auto z = hr.size;
	auto rr = hr.arr;
	std::map<std::size_t, std::size_t> ww;
	for (std::size_t i = 0; i < n; i++)
		ww.insert_or_assign(vv[i],i);
	auto hs = std::make_unique<HistorySparse>();
	hs->size = z;
	hs->vectorDimension = new std::size_t[z];
	auto dd = hs->vectorDimension;
	SizeList rr1;
	rr1.reserve(n*z);
	if (hr.evient)
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t jn = j*n;
			std::size_t k = 0;
			for (auto& p : ww)
			{
				std::size_t u = rr[jn + p.second];
				if (u)
				{
					rr1.push_back(p.first);
					k++;
				}
			}
			dd[j] = k;
		}
	else
		for (std::size_t j = 0; j < z; j++)
		{
			std::size_t k = 0;
			for (auto& p : ww)
			{
				std::size_t u = rr[p.second*z + j];
				if (u)
				{
					rr1.push_back(p.first);
					k++;
				}
			}
			dd[j] = k;
		}
	hs->arr = new std::size_t[rr1.size()];
	memcpy(hs->arr, rr1.data(), rr1.size() * sizeof(std::size_t));
	return hs;
}

// historySparsesHistoryRepa :: HistorySparse -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historySparsesHistoryRepa(const HistorySparse& hs)
{
	auto z = hs.size;
	auto dd = hs.vectorDimension;
	auto rr = hs.arr;
	std::set<std::size_t> ww;
	std::size_t m = 0;
	for (std::size_t j = 0; j < z; j++)
	{
		auto k = dd[j];
		for (std::size_t i = 0; i < k; i++)
		{
			ww.insert(rr[m]);
			m++;
		}
	}
	auto hr = std::make_unique<HistoryRepa>();
	hr->dimension = ww.size();
	auto n = hr->dimension;
	hr->vectorVar = new std::size_t[n];
	auto vv = hr->vectorVar;
	hr->shape = new std::size_t[n];
	auto sh = hr->shape;
	hr->size = z;
	{
		std::size_t i = 0;
		for (auto w : ww)
		{
			vv[i] = w;
			sh[i] = 2;
			i++;
		}
	}
	auto& mvv = hr->mapVarInt();
	hr->evient = true;
	hr->arr = new unsigned char[z*n];
	auto rr1 = hr->arr;
	memset(rr1, 0, z*n);
	m = 0;
	for (std::size_t j = 0; j < z; j++)
	{
		std::size_t jn = j*n;
		auto k = dd[j];
		for (std::size_t i = 0; i < k; i++)
		{
			rr1[jn + mvv[rr[m]]] = 1;
			m++;
		}
	}
	return hr;
}

// listSetIntsHistorySparse :: SizeSetList -> HistorySparse
std::unique_ptr<HistorySparse> Alignment::listSetIntsHistorySparse(const SizeSetList& ll)
{
	auto z = ll.size();
	auto hs = std::make_unique<HistorySparse>();
	hs->size = z;
	hs->vectorDimension = new std::size_t[z];
	auto dd = hs->vectorDimension;
	std::size_t m = 0;
	std::size_t j = 0;
	for (auto& qq : ll)
	{
		std::size_t k = qq.size();
		dd[j] = k;
		m += k;
		j++;
	}
	hs->arr = new std::size_t[m];	
	auto rr = hs->arr;
	m = 0;
	for (auto& qq : ll)
	{
		for (auto w : qq)
		{
			rr[m] = w;
			m++;
		}
	}
	return hs;
}

// historySparsesListSetInt :: HistorySparse -> SizeSetList
std::unique_ptr<SizeSetList> Alignment::historySparsesListSetInt(const HistorySparse& hs)
{
	auto z = hs.size;
	auto dd = hs.vectorDimension;
	auto rr = hs.arr;
	auto ll = std::make_unique<SizeSetList>(z);
	std::size_t m = 0;
	for (std::size_t j = 0; j < z; j++)
	{
		auto& qq = (*ll)[j];
		auto k = dd[j];
		for (std::size_t i = 0; i < k; i++)
		{
			qq.insert(rr[m]);
			m++;
		}
	}
	return ll;
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
		out << s;
	}
	out << "]," << u << ",[";
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

void Alignment::TransformRepa::reframe_u(const SizeSizeUMap& nn)
{
	auto n = dimension;
	auto vv = vectorVar;
	bool altered = false;
	for (std::size_t i = 0; i < n; i++)
	{
		auto it = nn.find(vv[i]);
		if (it != nn.end())
		{
			vv[i] = it->second;
			altered = true;
		}
	}
	{
		auto it = nn.find(derived);
		if (it != nn.end())
			derived = it->second;
	}
	if (!_mapVarInt && altered)
	{
		delete _mapVarInt;
		_mapVarInt = 0;
	}
}

// transformRepasTransformRepa :: TransformRepa -> TransformRepa
std::unique_ptr<TransformRepa> Alignment::transformRepasTransformRepa(const TransformRepa& tr)
{
	auto n = tr.dimension;
	auto vv = tr.vectorVar;
	auto sh = tr.shape;
	auto rr = tr.arr;
	auto tr1 = std::make_unique<TransformRepa>();
	if (!rr || !n)
		return tr1;
	tr1->derived = tr.derived;
	tr1->valency = tr.valency;
	tr1->dimension = n;
	tr1->vectorVar = new std::size_t[n];
	auto vv1 = tr1->vectorVar;
	tr1->shape = new std::size_t[n];
	auto sh1 = tr1->shape;
	std::size_t sz = 1;
	for (std::size_t i = 0; i < n; i++)
	{
		vv1[i] = vv[i];
		auto s = sh[i];
		sh1[i] = s;
		sz *= s;
	}
	tr1->arr = new unsigned char[sz];
	auto rr1 = tr1->arr;
	memcpy(rr1, rr, sz);
	return tr1;
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
	tr->shape = new std::size_t[n];
	auto sh = tr->shape;
	std::vector<ValSizeUMap> mm(n + 1);
	std::size_t sz = 1;
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = qq1[i];
		vv[i] = vvi[v];
		auto xx = systemsVarsSetValue(uu, v);
		auto s = xx.size();
		sh[i] = s;
		sz *= s;
		auto& yy = mm[i];
		yy.reserve(s);
		std::size_t j = 0;
		for (auto& x : xx)
			yy.insert_or_assign(x, j++);
	}
	unsigned char ucmax = std::numeric_limits<unsigned char>::max();
	{
		auto& v = *wit;
		tr->derived = vvi[v];
		auto xx = systemsVarsSetValue(uu, v);
		auto s = xx.size();
		if (s > ucmax + 1)
			throw std::out_of_range("systemsTransformsTransformRepa_u");
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
			j = sh[i] * j + mm[i][sm.find(qq1[i])->second];
		rr[j] = (unsigned char)mm[n][sm.find(*wit)->second];
	}
	return tr;
}

// systemsTransformRepasTransform_u :: System -> TransformRepa -> Transform
std::unique_ptr<Transform> Alignment::systemsTransformRepasTransform_u(const System& uu, const SystemRepa& ur, const TransformRepa& tr)
{
	auto& ivv = ur.listVarSizePair;
	auto n = tr.dimension;
	auto vv = tr.vectorVar;
	auto sh = tr.shape;
	auto rr = tr.arr;
	if (!tr.arr || !n)
		return std::make_unique<Transform>();
	std::size_t sz = 1;
	std::vector<ValList> mm(n + 1);
	VarList ww;
	for (std::size_t i = 0; i < n; i++)
	{
		auto& v = *ivv[vv[i]].first;
		ww.push_back(v);
		auto xx = systemsVarsSetValue(uu, v);
		auto s = xx.size();
		sz *= s;
		auto& yy = mm[i];
		yy.reserve(s);
		for (auto& x : xx)
			yy.push_back(x);
	}
	auto& w = *ivv[tr.derived].first;
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
		ss.reserve(n + 1);
		for (std::size_t i = 0; i < n; i++)
			ss.push_back(VarValPair(ww[i], mm[i][ii[i]]));
		ss.push_back(VarValPair(w, mm[n][rr[j]]));
		am.insert_or_assign(State(ss), Rational(1));
		for (long long k = n - 1; k >= 0; k--)
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

std::ostream& operator<<(std::ostream& out, const FudRepa& fr)
{
	out << "[";
	for (std::size_t i = 0; i < fr.layers.size(); i++)
	{
		out << (i ? "," : "") << "[";
		for (std::size_t j = 0; j < fr.layers[i].size(); j++)
			out << (j ? "," : "") << *fr.layers[i][j];
		out << "]";
	}
	out << "]";
	return out;
}

void Alignment::FudRepa::reframe_u(const SizeSizeUMap& nn)
{
	for (auto& ll : layers)
		for (auto& tr : ll)
			tr->reframe_u(nn);
}

// fudRepasFudRepa :: FudRepa -> FudRepa
std::unique_ptr<FudRepa> Alignment::fudRepasFudRepa(const FudRepa& fr)
{
	auto trcopy = transformRepasTransformRepa; 

	auto fr1 = std::make_unique<FudRepa>();
	fr1->layers.resize(fr.layers.size());
	for (std::size_t i = 0; i < fr.layers.size(); i++)
	{
		auto& ll = fr.layers[i];
		auto& ll1 = fr1->layers[i];
		ll1.reserve(ll.size());
		for (auto& tr : ll)
			ll1.push_back(std::move(trcopy(*tr)));
	}
	return fr1;
}

// setVariablesListTransformRepasFudRepa_u :: [VariableRepa] -> [TransformRepa] -> FudRepa
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

// fudRepasSize :: FudRepa -> Integer
std::size_t Alignment::fudRepasSize(const FudRepa& fr)
{
	std::size_t l = 0;
	for (auto& ll : fr.layers)
			l += ll.size();
	return l;
}

// fudRepasSetVar :: FudRepa -> [VariableRepa]
std::unique_ptr<SizeUSet> Alignment::fudRepasSetVar(const FudRepa& fr)
{
	std::size_t l = 0;
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
			l += tt->dimension + 1;
	auto vv = std::make_unique<SizeUSet>(l);
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
		{
			vv->insert(tt->derived);
			auto n = tt->dimension;
			auto xx = tt->vectorVar;
			for (std::size_t i = 0; i < n; i++)
				vv->insert(xx[i]);
		}
	return vv;
}

// fudRepasDerived :: FudRepa -> [VariableRepa]
std::unique_ptr<SizeUSet> Alignment::fudRepasDerived(const FudRepa& fr)
{
	auto vv = std::make_unique<SizeUSet>(fudRepasSize(fr));
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
			vv->insert(tt->derived);
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
		{
			auto n = tt->dimension;
			auto xx = tt->vectorVar;
			for (std::size_t i = 0; i < n; i++)
				vv->erase(xx[i]);
		}
	return vv;
}

// fudRepasUnderlying :: FudRepa -> [VariableRepa]
std::unique_ptr<SizeUSet> Alignment::fudRepasUnderlying(const FudRepa& fr)
{
	auto vv = std::make_unique<SizeUSet>(fudRepasSize(fr));
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
		{
			auto n = tt->dimension;
			auto xx = tt->vectorVar;
			for (std::size_t i = 0; i < n; i++)
				vv->insert(xx[i]);
		}
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
			vv->erase(tt->derived);
	return vv;
}

// fudRepasDefinitions :: FudRepa -> Map.Map VariableRepa TransformRepa
std::unique_ptr<SizeTransformRepaPtrMap> Alignment::fudRepasDefinitions(const FudRepa& fr)
{
	auto mm = std::make_unique<SizeTransformRepaPtrMap>();
	for (auto& ll : fr.layers)
		for (auto& tt : ll)
			mm->insert_or_assign(tt->derived, tt);
	return mm;
}

// fudRepasSetVarsDepends :: FudRepa -> [VariableRepa] -> FudRepa
std::unique_ptr<TransformRepaPtrList> Alignment::fudRepasSetVarsDepends(const FudRepa& fr, const SizeUSet& kk)
{
	auto mm = fudRepasDefinitions(fr);
	auto ll = std::make_unique<TransformRepaPtrList>();
	ll->reserve(mm->size());
	SizeUSet kk0(mm->size());
	kk0.insert(kk.begin(),kk.end());
	SizeUSet kk1(mm->size());
	SizeUSet kk2(mm->size());
	SizeUSet* kka = &kk0;
	SizeUSet* kkb = &kk1;
	while (kka->size())
	{
		for (auto& w : *kka)
		{
			if (kk2.find(w) != kk2.end())
				continue;
			auto it = mm->find(w);
			if (it == mm->end())
				continue;
			kk2.insert(w);
			auto& tr = it->second;
			ll->push_back(tr);
			auto n = tr->dimension;
			auto xx = tr->vectorVar;
			for (std::size_t i = 0; i < n; i++)
				kkb->insert(xx[i]);
		}
		SizeUSet* kkc = kka;
		kka = kkb;
		kkb = kkc;
		kkb->clear();
	}
	return ll;
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
		hr1->shape = new std::size_t[n1];
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
		hr1->evient = hr.evient;
		hr1->arr = new unsigned char[z*p];
		auto rr1 = hr1->arr;
		if (hr.evient)
			for (std::size_t j = 0; j < z; j++)
			{
				std::size_t jn = j*n;
				std::size_t jp = j*p;
				for (std::size_t i = 0; i < n; i++)
					rr1[jp + i] = rr[jn + i];
			}
		else
			memcpy(rr1, rr, z*n);
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
				if (hr.evient)
				{
					if (m > 0)
						for (std::size_t j = 0; j < z; j++)
						{
							std::size_t jp = j*p;
							std::size_t k = rr1[jp + pkk[0]];
							for (std::size_t i = 1; i < m; i++)
								k = sh[i] * k + rr1[jp + pkk[i]];
							rr1[jp + q] = ar[k];
						}
					else
						for (std::size_t j = 0; j < z; j++)
							rr1[j*p + q] = 0;
				}
				else
				{
					std::size_t qz = q*z;
					if (m > 0)
						for (std::size_t j = 0; j < z; j++)
						{
							std::size_t k = rr1[pkk[0] * z + j];
							for (std::size_t i = 1; i < m; i++)
								k = sh[i] * k + rr1[pkk[i] * z + j];
							rr1[qz + j] = ar[k];
						}
					else
						for (std::size_t j = 0; j < z; j++)
							rr1[qz + j] = 0;
				}
				q++;
			}
		delete[] pkk;
		return hr1;
	}

	void historyRepasFudRepasMultiply_up_layer(std::size_t tint, std::size_t mmax, SizeSizeUMap* mkk, std::size_t p, HistoryRepa* hr, std::size_t q, TransformRepaPtrList* ll, std::size_t t)
	{
		auto z = hr->size;
		auto rr1 = hr->arr;
		std::size_t* pkk = new std::size_t[mmax];
		for (std::size_t i = 0; i < ll->size(); i++)
			if (i % tint == t)
			{
				auto tr = (*ll)[i];
				auto m = tr->dimension;
				auto ww = tr->vectorVar;
				auto sh = tr->shape;
				auto ar = tr->arr;
				for (std::size_t i = 0; i < m; i++)
					pkk[i] = (*mkk)[ww[i]];
				if (hr->evient)
				{
					if (m > 0)
						for (std::size_t j = 0; j < z; j++)
						{
							std::size_t jp = j*p;
							std::size_t k = rr1[jp + pkk[0]];
							for (std::size_t i = 1; i < m; i++)
								k = sh[i] * k + rr1[jp + pkk[i]];
							rr1[jp + q+i] = ar[k];
						}
					else
						for (std::size_t j = 0; j < z; j++)
							rr1[j*p + q+i] = 0;
				}
				else
				{
					std::size_t qz = (q+i)*z;
					if (m > 0)
						for (std::size_t j = 0; j < z; j++)
						{
							std::size_t k = rr1[pkk[0] * z + j];
							for (std::size_t i = 1; i < m; i++)
								k = sh[i] * k + rr1[pkk[i] * z + j];
							rr1[qz + j] = ar[k];
						}
					else
						for (std::size_t j = 0; j < z; j++)
							rr1[qz + j] = 0;
				}
			}
		delete[] pkk;
	}

	// historyRepasFudRepasMultiply_up :: Int -> HistoryRepa -> FudRepa -> HistoryRepa
	std::unique_ptr<HistoryRepa> Alignment::historyRepasFudRepasMultiply_up(std::size_t tint, const HistoryRepa& hr, const FudRepa& fr)
	{
		auto layer = historyRepasFudRepasMultiply_up_layer;

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
		hr1->shape = new std::size_t[n1];
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
		hr1->evient = hr.evient;
		hr1->arr = new unsigned char[z*p];
		auto rr1 = hr1->arr;
		if (hr.evient)
			for (std::size_t j = 0; j < z; j++)
			{
				std::size_t jn = j*n;
				std::size_t jp = j*p;
				for (std::size_t i = 0; i < n; i++)
					rr1[jp + i] = rr[jn + i];
			}
		else
			memcpy(rr1, rr, z*n);
		auto q = n;
		for (auto& ll : fr.layers)
		{
			std::size_t tint1 = std::min(tint, ll.size());
			std::vector<std::thread> threads;
			threads.reserve(tint1);
			for (std::size_t t = 0; t < tint1; t++)
				threads.push_back(std::thread(layer, tint1, mmax, (SizeSizeUMap*)&mkk, p, hr1.get(), q, (TransformRepaPtrList*)&ll, t));
			for (auto& t : threads)
				t.join();
			q += ll.size();
		}
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
		hh.map_u().insert_or_assign(Id(1), ss);
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
		auto hh = systemsHistoryRepasHistory_u(uu, ur, hr);
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
		auto hh = systemsHistoryRepasHistory_u(uu, ur, hr);
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

std::ostream& operator<<(std::ostream& out, const ApplicationRepa& dr)
{
	out << "(" << dr.substrate << ",";
	if (dr.fud)
		out << *dr.fud;
	out << ",";
	if (dr.slices)
		out << *dr.slices;
	out << ")";
	return out;
}

void Alignment::ApplicationRepa::reframe_u(const SizeSizeUMap& nn)
{
	auto& vv = substrate;
	auto n = vv.size();
	for (std::size_t i = 0; i < n; i++)
	{
		auto it = nn.find(vv[i]);
		if (it != nn.end())
			vv[i] = it->second;
	}
	fud->reframe_u(nn);
	auto ll = treesPaths(*slices);
	bool altered = false;
	for (auto& pp : *ll)
		for (auto& v : pp)
		{
			auto it = nn.find(v);
			if (it != nn.end())
			{
				v = it->second;
				altered = true;
			}
		}
	if (altered)
		slices = std::move(pathsTree(*ll));
}

// applicationRepasApplicationRepa_u :: ApplicationRepa -> ApplicationRepa
std::unique_ptr<ApplicationRepa> Alignment::applicationRepasApplicationRepa_u(const ApplicationRepa& dr)
{
	auto dr1 = std::make_unique<ApplicationRepa>();
	dr1->substrate = dr.substrate;
	dr1->fud = std::move(fudRepasFudRepa(*dr.fud));
	dr1->slices = std::move(pathsTree(*treesPaths(*dr.slices)));
	return dr1;
}

// applicationRepaPairsJoin :: ApplicationRepa -> ApplicationRepa -> ApplicationRepa
std::unique_ptr<ApplicationRepa> Alignment::applicationRepaPairsJoin_u(const ApplicationRepa& dr1, const ApplicationRepa& dr2)
{
	auto llfr = setVariablesListTransformRepasFudRepa_u;
	auto frvars = fudRepasSetVar;
	auto frdep = fudRepasSetVarsDepends;

	auto dr3 = std::make_unique<ApplicationRepa>();
	dr3->substrate = dr1.substrate;
	dr3->slices = std::move(pathsTree(*treesPaths(*dr2.slices)));
	SizeUSet vv(dr1.substrate.begin(), dr1.substrate.end());
	auto sl = treesElements(*dr2.slices);
	SizeUSet ww(sl->begin(), sl->end());
	FudRepa fr;
	fr.layers.reserve(dr1.fud->layers.size() + dr2.fud->layers.size());
	fr.layers.insert(fr.layers.end(), dr1.fud->layers.begin(), dr1.fud->layers.end());
	fr.layers.insert(fr.layers.end(), dr2.fud->layers.begin(), dr2.fud->layers.end());
	auto vv1 = frvars(*dr1.fud);
	for (auto v : dr2.substrate)
		if (vv1->find(v) == vv1->end())
		{
			vv.insert(v);
			dr3->substrate.push_back(v);
		}
	dr3->fud = std::move(llfr(vv, *frdep(fr, ww)));
	return dr3;
}

std::size_t listVarsArrayHistoryEvientAlignedTop_u(
	std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t z1, std::size_t z2,
	std::size_t* ppww, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	long long x3;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t pj;
	std::size_t jn;
	std::size_t si;
	std::size_t sj;
	std::size_t u;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (ii = 0; ii<m - 1; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ij = ii + 1; ij<m; ij++)
		{
			pj = ppww[ij];
			sj = svv[pj];
			u = si*sj;
			if (u <= xmax)
			{
				s++;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				for (j = 0; j < z1; j++)
				{
					jn = n*j;
					aa[sj*phh1[jn + pi] + phh1[jn + pj]] += 1.0;
				}
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				for (j = 0; j < z2; j++)
				{
					jn = n*j;
					aa[sj*phh2[jn + pi] + phh2[jn + pj]] += f;
				}
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					x1 = zf*xx1[pi][i];
					x2 = zf*xx2[pi][i];
					for (j = 0; j<sj; j++)
					{
						a2 += alngam(x1*xx1[pj][j] + 1.0);
						b2 += alngam(x2*xx2[pj][j] + 1.0);
					}
				}
				if (t<omax)
				{
					tww1[t] = pi;
					tww2[t] = pj;
					ts1[t] = a1 - a2 - b1 + b2;
					ts2[t] = b2 - b1;
					ts3[t] = -(long long)u;
					t++;
					if (t == omax)
					{
						for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[i];
							if (t1>x1)
							{
								t1 = x1;
								t2 = ts2[i];
								t3 = ts3[i];
								tm = i;
							}
							else if (t1 == x1)
							{
								x2 = ts2[i];
								if (t2>x2)
								{
									t2 = x2;
									t3 = ts3[i];
									tm = i;
								}
								else if (t2 == x2)
								{
									x3 = ts3[i];
									if (t3>x3)
									{
										t3 = x3;
										tm = i;
									}
								}
							}
						}
					}
				}
				else
				{
					x1 = a1 - a2 - b1 + b2;
					x2 = b2 - b1;
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = pi;
						tww2[tm] = pj;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[i];
							if (t1>x1)
							{
								t1 = x1;
								t2 = ts2[i];
								t3 = ts3[i];
								tm = i;
							}
							else if (t1 == x1)
							{
								x2 = ts2[i];
								if (t2>x2)
								{
									t2 = x2;
									t3 = ts3[i];
									tm = i;
								}
								else if (t2 == x2)
								{
									x3 = ts3[i];
									if (t3>x3)
									{
										t3 = x3;
										tm = i;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	return t;
}

std::size_t listVarsArrayHistoryVarientAlignedTop_u(
	std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t z1, std::size_t z2,
	std::size_t* ppww, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	long long x3;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t pj;
	std::size_t qi;
	std::size_t qj;
	std::size_t si;
	std::size_t sj;
	std::size_t u;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (ii = 0; ii<m - 1; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ij = ii + 1; ij<m; ij++)
		{
			pj = ppww[ij];
			sj = svv[pj];
			u = si*sj;
			if (u <= xmax)
			{
				s++;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				qi = z1*pi;
				qj = z1*pj;
				for (j = 0; j < z1; j++)
					aa[sj*phh1[qi + j] + phh1[qj + j]] += 1.0;
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				qi = z2*pi;
				qj = z2*pj;
				for (j = 0; j < z2; j++)
					aa[sj*phh2[qi + j] + phh2[qj + j]] += f;
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					x1 = zf*xx1[pi][i];
					x2 = zf*xx2[pi][i];
					for (j = 0; j<sj; j++)
					{
						a2 += alngam(x1*xx1[pj][j] + 1.0);
						b2 += alngam(x2*xx2[pj][j] + 1.0);
					}
				}
				if (t<omax)
				{
					tww1[t] = pi;
					tww2[t] = pj;
					ts1[t] = a1 - a2 - b1 + b2;
					ts2[t] = b2 - b1;
					ts3[t] = -(long long)u;
					t++;
					if (t == omax)
					{
						for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[i];
							if (t1>x1)
							{
								t1 = x1;
								t2 = ts2[i];
								t3 = ts3[i];
								tm = i;
							}
							else if (t1 == x1)
							{
								x2 = ts2[i];
								if (t2>x2)
								{
									t2 = x2;
									t3 = ts3[i];
									tm = i;
								}
								else if (t2 == x2)
								{
									x3 = ts3[i];
									if (t3>x3)
									{
										t3 = x3;
										tm = i;
									}
								}
							}
						}
					}
				}
				else
				{
					x1 = a1 - a2 - b1 + b2;
					x2 = b2 - b1;
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = pi;
						tww2[tm] = pj;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[i];
							if (t1>x1)
							{
								t1 = x1;
								t2 = ts2[i];
								t3 = ts3[i];
								tm = i;
							}
							else if (t1 == x1)
							{
								x2 = ts2[i];
								if (t2>x2)
								{
									t2 = x2;
									t3 = ts3[i];
									tm = i;
								}
								else if (t2 == x2)
								{
									x3 = ts3[i];
									if (t3>x3)
									{
										t3 = x3;
										tm = i;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	return t;
}

void listVarsArrayHistoryVarientAlignedTop_up(
	std::size_t xmax, std::size_t omax, std::size_t tint, std::size_t n, std::size_t* svv, std::size_t m, std::size_t z1, std::size_t z2,
	std::size_t* ppww, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t h, std::size_t* tww1, std::size_t* tww2, double* ts1, std::size_t* t, std::size_t* s)
{
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	std::size_t tm;
	double x1;
	double x2;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t pj;
	std::size_t qi;
	std::size_t qj;
	std::size_t si;
	std::size_t sj;
	std::size_t u;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;
	std::size_t homax = h*omax;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (ii = 0; ii<m-1; ii++)
		if (ii % tint == h)
		{
			pi = ppww[ii];
			si = svv[pi];
			for (ij = ii + 1; ij<m; ij++)
			{
				pj = ppww[ij];
				sj = svv[pj];
				u = si*sj;
				if (u <= xmax)
				{
					s[h]++;
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					qi = z1*pi;
					qj = z1*pj;
					for (j = 0; j < z1; j++)
						aa[sj*phh1[qi + j] + phh1[qj + j]] += 1.0;
					for (a1 = 0.0, i = 0; i<u; i++)
						a1 += alngam(aa[i]);
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					qi = z2*pi;
					qj = z2*pj;
					for (j = 0; j < z2; j++)
						aa[sj*phh2[qi + j] + phh2[qj + j]] += f;
					for (b1 = 0.0, i = 0; i<u; i++)
						b1 += alngam(aa[i]);
					for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
					{
						x1 = zf*xx1[pi][i];
						x2 = zf*xx2[pi][i];
						for (j = 0; j<sj; j++)
						{
							a2 += alngam(x1*xx1[pj][j] + 1.0);
							b2 += alngam(x2*xx2[pj][j] + 1.0);
						}
					}
					if (t[h]<omax)
					{
						tww1[homax + t[h]] = pi;
						tww2[homax + t[h]] = pj;
						ts1[homax + t[h]] = a1 - a2 - b1 + b2;
						t[h]++;
						if (t[h] == omax)
						{
							for (t1 = ts1[homax], tm = 0, i = 1; i<omax; i++)
							{
								x1 = ts1[homax + i];
								if (t1>x1)
								{
									t1 = x1;
									tm = i;
								}
							}
						}
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						if (t1<x1)
						{
							tww1[homax + tm] = pi;
							tww2[homax + tm] = pj;
							ts1[homax + tm] = x1;
							for (t1 = ts1[homax], tm = 0, i = 1; i<omax; i++)
							{
								x1 = ts1[homax + i];
								if (t1>x1)
								{
									t1 = x1;
									tm = i;
								}
							}
						}
					}
				}
			}
		}
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
}

// parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->([(Double, [VariableRepa])],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t xmax, std::size_t omax, const SizeList& ww, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	std::size_t* tww1 = new std::size_t[omax];
	std::size_t* tww2 = new std::size_t[omax];
	double* ts1 = new double[omax];
	double* ts2 = new double[omax];
	long long* ts3 = new long long[omax];
	std::size_t s = 0;
	std::size_t t = 0;
	if (hh.evient)
		t = listVarsArrayHistoryEvientAlignedTop_u(xmax, omax, n, svv, m, z, zrr, pww, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	else
		t = listVarsArrayHistoryVarientAlignedTop_u(xmax, omax, n, svv, m, z, zrr, pww, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = 0; i < t; i++)
		qq->push_back(DoubleSizeListPair(ts1[i], SizeList{ vhh[tww1[i]],vhh[tww2[i]] }));
	delete[] ts3;
	delete[] ts2;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s);
}

// parametersSetVarsHistoryRepasSetSetVarsAlignedTop_up :: Integer -> Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->([(Double, [VariableRepa])],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsHistoryRepasSetSetVarsAlignedTop_up(std::size_t xmax, std::size_t omax, std::size_t tint, const SizeList& ww, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	if (hh.evient)
		throw std::out_of_range("parametersSetVarsHistoryRepasSetSetVarsAlignedTop_up");
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t tint1 = std::min(tint, m-1);
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	std::size_t* tww1 = new std::size_t[tint1*omax];
	std::size_t* tww2 = new std::size_t[tint1*omax];
	double* ts1 = new double[tint1*omax];
	std::size_t* s = new std::size_t[tint1];
	std::size_t* t = new std::size_t[tint1];
	std::vector<std::thread> threads;
	threads.reserve(tint1);
	for (std::size_t h = 0; h < tint1; h++)
	{
		s[h] = 0;
		t[h] = 0;
		threads.push_back(std::thread(listVarsArrayHistoryVarientAlignedTop_up, xmax, omax, tint1, n, svv, m, z, zrr, pww, phh1, pxx1, phh2, pxx2, h, tww1, tww2, ts1, t, s));
	}
	for (auto& h : threads)
		h.join();
	std::size_t s1 = 0;
	DoubleSizePairList qq1;
	qq1.reserve(tint1*omax);
	SizeListList qq2;
	qq2.reserve(tint1*omax);
	std::size_t j = 0;
	for (std::size_t h = 0; h < tint1; h++)
	{
		for (std::size_t i = 0; i < t[h]; i++)
		{
			qq1.push_back(DoubleSizePair(ts1[h*omax + i], j));
			qq2.push_back(SizeList{vhh[tww1[h*omax + i]],vhh[tww2[h*omax + i]]});
			j++;
		}
		s1 += s[h];
	}
	std::sort(qq1.begin(), qq1.end());
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = j - std::min(j, omax); i < j; i++)
		qq->push_back(DoubleSizeListPair(qq1[i].first,qq2[qq1[i].second]));
	delete[] s;
	delete[] t;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s1);
}

inline void incIndex(std::size_t n, std::size_t* svv, std::size_t* ivv)
{
	long long k;
	std::size_t y;

	for (k = n - 1; k >= 0; k--)
	{
		y = ivv[k] + 1;
		if (y == svv[k])
			ivv[k] = 0;
		else
		{
			ivv[k] = y;
			break;
		}
	}
}

inline std::size_t isdup(std::size_t e, std::size_t pi, std::size_t* pj, std::size_t qi, std::size_t* qj)
{
	std::size_t ok = 1;
	std::size_t i;
	std::size_t j;
	std::size_t pk;

	if (pi != qi)
	{
		for (i = 0, ok = 0; i<e; i++)
			if (pi == qj[i])
			{
				ok = 1;
				break;
			}
	}
	for (i = 0; ok && i<e; i++)
	{
		pk = pj[i];
		if (pk != qi)
		{
			for (j = 0, ok = 0; j<e; j++)
				if (pk == qj[j])
				{
					ok = 1;
					break;
				}
		}
	}
	return ok;
}

inline std::size_t hash(std::size_t n, std::size_t e, std::size_t pi, std::size_t* pj)
{
	std::size_t h;
	std::size_t i;
	std::size_t j;
	std::size_t mk;
	std::size_t pk;
	std::size_t pl;

	for (i = 0, h = 0, mk = n, pk = n; i <= e; i++, mk = pk, pk = n)
	{
		if (mk == n || mk<pi)
			pk = pi;
		for (j = 0; j<e; j++)
		{
			pl = pj[j];
			if (pl<pk && (mk == n || mk<pl))
				pk = pl;
		}
		h = n*h + pk;
	}
	return h;
}

std::size_t listVarsListTuplesArrayHistoryEvientsAlignedTop_u(
	unsigned char dense, std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e, std::size_t z1, std::size_t z2,
	std::size_t* ppww, std::size_t* ppdd, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	long long x3;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t jn;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;
	std::size_t ok;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ij = 0; ij<d; ij++)
		{
			pj = pdd[ij];
			for (k = 0, ok = 1, u1 = 1; k<e; k++)
			{
				pk = pj[k];
				if (pk == pi)
				{
					ok = 0;
					break;
				}
				sk = svv[pk];
				sj[k] = sk;
				u1 *= sk;
			}
			u = u1*si;
			if (ok && u <= xmax)
			{
				s++;
				x4 = hash(n, e, pi, pj);
				for (i = 0; i<t; i++)
					if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
					{
						ok = 0;
						break;
					}
				if (!ok)
					continue;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				pk = pj[0];
				sk = sj[0];
				for (j = 0; j<z1; j++)
				{
					jn = n*j;
					for (k = 1, a = sk*phh1[jn + pi] + phh1[jn + pk]; k<e; k++)
						a = sj[k] * a + phh1[jn + pj[k]];
					aa[a] += 1.0;
				}
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				for (j = 0; j<z2; j++)
				{
					jn = n*j;
					for (k = 1, a = sk*phh2[jn + pi] + phh2[jn + pk]; k<e; k++)
						a = sj[k] * a + phh2[jn + pj[k]];
					aa[a] += f;
				}
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (k = 0; k<e; k++)
					yj[k] = 0;
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					y1 = zf*xx1[pi][i];
					y2 = zf*xx2[pi][i];
					for (j = 0; j<u1; j++)
					{
						x1 = y1;
						x2 = y2;
						for (k = 0; k<e; k++)
						{
							pk = pj[k];
							yk = yj[k];
							x1 *= xx1[pk][yk];
							x2 *= xx2[pk][yk];
						}
						a2 += alngam(x1 + 1.0);
						b2 += alngam(x2 + 1.0);
						incIndex(e, sj, yj);
					}
				}
				if (t<omax)
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					tww1[t] = ii;
					tww2[t] = ij;
					ts1[t] = x1;
					ts2[t] = x2;
					ts3[t] = x3;
					ts4[t] = x4;
					t++;
					if (t == omax)
						findm = 1;
				}
				else
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = ii;
						tww2[tm] = ij;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						ts4[tm] = x4;
						findm = 1;
					}
				}
				if (findm)
				{
					for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
					{
						x1 = ts1[i];
						if (t1>x1)
						{
							t1 = x1;
							t2 = ts2[i];
							t3 = ts3[i];
							tm = i;
						}
						else if (t1 == x1)
						{
							x2 = ts2[i];
							if (t2>x2)
							{
								t2 = x2;
								t3 = ts3[i];
								tm = i;
							}
							else if (t2 == x2)
							{
								x3 = ts3[i];
								if (t3>x3)
								{
									t3 = x3;
									tm = i;
								}
							}
						}
					}
					findm = 0;
				}
			}
		}
	}
	delete[] yj;
	delete[] sj;
	delete[] ts4;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
	return t;
}

std::size_t listVarsListTuplesArrayHistoryVarientsAlignedTop_u(
	unsigned char dense, std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e, std::size_t z1, std::size_t z2,
	std::size_t* ppww, std::size_t* ppdd, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	long long x3;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t qi;
	std::size_t* qj = new std::size_t[e];
	std::size_t qk;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;
	std::size_t ok;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ij = 0; ij<d; ij++)
		{
			pj = pdd[ij];
			for (k = 0, ok = 1, u1 = 1; k<e; k++)
			{
				pk = pj[k];
				if (pk == pi)
				{
					ok = 0;
					break;
				}
				sk = svv[pk];
				sj[k] = sk;
				u1 *= sk;
			}
			u = u1*si;
			if (ok && u <= xmax)
			{
				s++;
				x4 = hash(n, e, pi, pj);
				for (i = 0; i<t; i++)
					if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
					{
						ok = 0;
						break;
					}
				if (!ok)
					continue;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				pk = pj[0];
				sk = sj[0];
				qi = z1*pi;
				qk = z1*pk;
				for (k = 1; k<e; k++)
					qj[k] = z1*pj[k];
				for (j = 0; j<z1; j++)
				{
					for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
						a = sj[k] * a + phh1[qj[k] + j];
					aa[a] += 1.0;
				}
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				qi = z2*pi;
				qk = z2*pk;
				for (k = 1; k<e; k++)
					qj[k] = z2*pj[k];
				for (j = 0; j<z2; j++)
				{
					for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
						a = sj[k] * a + phh2[qj[k] + j];
					aa[a] += f;
				}
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (k = 0; k<e; k++)
					yj[k] = 0;
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					y1 = zf*xx1[pi][i];
					y2 = zf*xx2[pi][i];
					for (j = 0; j<u1; j++)
					{
						x1 = y1;
						x2 = y2;
						for (k = 0; k<e; k++)
						{
							pk = pj[k];
							yk = yj[k];
							x1 *= xx1[pk][yk];
							x2 *= xx2[pk][yk];
						}
						a2 += alngam(x1 + 1.0);
						b2 += alngam(x2 + 1.0);
						incIndex(e, sj, yj);
					}
				}
				if (t<omax)
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					tww1[t] = ii;
					tww2[t] = ij;
					ts1[t] = x1;
					ts2[t] = x2;
					ts3[t] = x3;
					ts4[t] = x4;
					t++;
					if (t == omax)
						findm = 1;
				}
				else
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = ii;
						tww2[tm] = ij;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						ts4[tm] = x4;
						findm = 1;
					}
				}
				if (findm)
				{
					for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
					{
						x1 = ts1[i];
						if (t1>x1)
						{
							t1 = x1;
							t2 = ts2[i];
							t3 = ts3[i];
							tm = i;
						}
						else if (t1 == x1)
						{
							x2 = ts2[i];
							if (t2>x2)
							{
								t2 = x2;
								t3 = ts3[i];
								tm = i;
							}
							else if (t2 == x2)
							{
								x3 = ts3[i];
								if (t3>x3)
								{
									t3 = x3;
									tm = i;
								}
							}
						}
					}
					findm = 0;
				}
			}
		}
	}
	delete[] yj;
	delete[] sj;
	delete[] qj;
	delete[] ts4;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
	return t;
}

void listVarsListTuplesArrayHistoryVarientsAlignedTop_up(
	unsigned char dense, std::size_t xmax, std::size_t omax, std::size_t tint, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e, std::size_t z1, std::size_t z2,
	std::size_t* ppww, std::size_t* ppdd, unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t h, std::size_t* tww1, std::size_t* tww2, double* ts1, std::size_t* t, std::size_t* s)
{
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t qi;
	std::size_t* qj = new std::size_t[e];
	std::size_t qk;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;
	std::size_t ok;
	std::size_t homax = h*omax;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
		if (ii % tint == h)
		{
			pi = ppww[ii];
			si = svv[pi];
			for (ij = 0; ij<d; ij++)
			{
				pj = pdd[ij];
				for (k = 0, ok = 1, u1 = 1; k<e; k++)
				{
					pk = pj[k];
					if (pk == pi)
					{
						ok = 0;
						break;
					}
					sk = svv[pk];
					sj[k] = sk;
					u1 *= sk;
				}
				u = u1*si;
				if (ok && u <= xmax)
				{
					s[h]++;
					x4 = hash(n, e, pi, pj);
					for (i = 0; i<t[h]; i++)
						if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[homax + i]], pdd[tww2[homax + i]]))
						{
							ok = 0;
							break;
						}
					if (!ok)
						continue;
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					pk = pj[0];
					sk = sj[0];
					qi = z1*pi;
					qk = z1*pk;
					for (k = 1; k<e; k++)
						qj[k] = z1*pj[k];
					for (j = 0; j<z1; j++)
					{
						for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
							a = sj[k] * a + phh1[qj[k] + j];
						aa[a] += 1.0;
					}
					for (a1 = 0.0, i = 0; i<u; i++)
						a1 += alngam(aa[i]);
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					qi = z2*pi;
					qk = z2*pk;
					for (k = 1; k<e; k++)
						qj[k] = z2*pj[k];
					for (j = 0; j<z2; j++)
					{
						for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
							a = sj[k] * a + phh2[qj[k] + j];
						aa[a] += f;
					}
					for (b1 = 0.0, i = 0; i<u; i++)
						b1 += alngam(aa[i]);
					for (k = 0; k<e; k++)
						yj[k] = 0;
					for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
					{
						y1 = zf*xx1[pi][i];
						y2 = zf*xx2[pi][i];
						for (j = 0; j<u1; j++)
						{
							x1 = y1;
							x2 = y2;
							for (k = 0; k<e; k++)
							{
								pk = pj[k];
								yk = yj[k];
								x1 *= xx1[pk][yk];
								x2 *= xx2[pk][yk];
							}
							a2 += alngam(x1 + 1.0);
							b2 += alngam(x2 + 1.0);
							incIndex(e, sj, yj);
						}
					}
					if (t[h]<omax)
					{
						if (dense)
						{
							c = pow((double)u, 1.0 / ((double)(e + 1)));
							x1 = (a1 - a2 - b1 + b2) / c;
						}
						else
						{
							x1 = a1 - a2 - b1 + b2;
						}
						tww1[homax + t[h]] = ii;
						tww2[homax + t[h]] = ij;
						ts1[homax + t[h]] = x1;
						ts4[t[h]] = x4;
						t[h]++;
						if (t[h] == omax)
							findm = 1;
					}
					else
					{
						if (dense)
						{
							c = pow((double)u, 1.0 / ((double)(e + 1)));
							x1 = (a1 - a2 - b1 + b2) / c;
						}
						else
						{
							x1 = a1 - a2 - b1 + b2;
						}
						if (t1<x1)
						{
							tww1[homax + tm] = ii;
							tww2[homax + tm] = ij;
							ts1[homax + tm] = x1;
							ts4[tm] = x4;
							findm = 1;
						}
					}
					if (findm)
					{
						for (t1 = ts1[homax], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[homax + i];
							if (t1>x1)
							{
								t1 = x1;
								tm = i;
							}
						}
						findm = 0;
					}
				}
			}
		}
	delete[] yj;
	delete[] sj;
	delete[] qj;
	delete[] ts4;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
}


// parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa])],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u(std::size_t xmax, std::size_t omax, const SizeList& ww, const SizeListList& vdd, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	auto d = vdd.size();
	auto e = d ? vdd[0].size() : 0;
	std::size_t* pdd = new std::size_t[d*e];
	for (std::size_t i = 0; i < d; i++)
		for (std::size_t j = 0; j < e; j++)
			pdd[i*e + j] = mvv[vdd[i][j]];
	std::size_t* tww1 = new std::size_t[omax];
	std::size_t* tww2 = new std::size_t[omax];
	double* ts1 = new double[omax];
	double* ts2 = new double[omax];
	long long* ts3 = new long long[omax];
	std::size_t s = 0;
	std::size_t t = 0;
	if (hh.evient)
		t = listVarsListTuplesArrayHistoryEvientsAlignedTop_u(0, xmax, omax, n, svv, m, d, e, z, zrr, pww, pdd, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	else
		t = listVarsListTuplesArrayHistoryVarientsAlignedTop_u(0, xmax, omax, n, svv, m, d, e, z, zrr, pww, pdd, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = 0; i < t; i++)
	{
		SizeList ll;
		ll.push_back(ww[tww1[i]]);
		for (auto& jj : vdd[tww2[i]])
			ll.push_back(jj);
		qq->push_back(DoubleSizeListPair(ts1[i], ll));
	}
	delete[] ts3;
	delete[] ts2;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pdd;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s);
}

// parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_up :: Integer -> Integer -> Integer -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa])],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_up(std::size_t xmax, std::size_t omax, std::size_t tint, const SizeList& ww, const SizeListList& vdd, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	if (hh.evient)
		throw std::out_of_range("parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_up");
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t tint1 = std::min(tint, m);
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	auto d = vdd.size();
	auto e = d ? vdd[0].size() : 0;
	std::size_t* pdd = new std::size_t[d*e];
	for (std::size_t i = 0; i < d; i++)
		for (std::size_t j = 0; j < e; j++)
			pdd[i*e + j] = mvv[vdd[i][j]];
	std::size_t* tww1 = new std::size_t[tint1*omax];
	std::size_t* tww2 = new std::size_t[tint1*omax];
	double* ts1 = new double[tint1*omax];
	std::size_t* s = new std::size_t[tint1];
	std::size_t* t = new std::size_t[tint1];
	std::vector<std::thread> threads;
	threads.reserve(tint1);
	for (std::size_t h = 0; h < tint1; h++)
	{
		s[h] = 0;
		t[h] = 0;
		threads.push_back(std::thread(listVarsListTuplesArrayHistoryVarientsAlignedTop_up, 0, xmax, omax, tint1, n, svv, m, d, e, z, zrr, pww, pdd, phh1, pxx1, phh2, pxx2, h, tww1, tww2, ts1, t, s));
	}
	for (auto& h : threads)
		h.join();
	std::size_t s1 = 0;
	DoubleSizePairList qq1;
	qq1.reserve(tint1*omax);
	SizeListList qq2;
	qq2.reserve(tint1*omax);
	std::size_t j = 0;
	SizeSetSet qq3;
	for (std::size_t h = 0; h < tint1; h++)
	{
		for (std::size_t i = 0; i < t[h]; i++)
		{
			SizeList ll;
			ll.push_back(ww[tww1[h*omax + i]]);
			for (auto& jj : vdd[tww2[h*omax + i]])
				ll.push_back(jj);
			SizeSet ss(ll.begin(), ll.end());
			if (qq3.find(ss) == qq3.end())
			{
				qq1.push_back(DoubleSizePair(ts1[h*omax + i], j));
				qq2.push_back(ll);
				qq3.insert(ss);
				j++;
			}
		}
		s1 += s[h];
	}
	std::sort(qq1.begin(), qq1.end());
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = j - std::min(j, omax); i < j; i++)
		qq->push_back(DoubleSizeListPair(qq1[i].first, qq2[qq1[i].second]));
	delete[] s;
	delete[] t;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pdd;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s1);
}

inline std::size_t toIndexPerm(std::size_t n, std::size_t* ppp, std::size_t* svv, std::size_t* ivv)
{
	std::size_t k;
	std::size_t p;
	std::size_t a;

	for (k = 1, a = ivv[ppp[0]]; k<n; k++)
	{
		p = ppp[k];
		a = svv[p] * a + ivv[p];
	}
	return a;
}

void listListVarsArrayHistoryPairsPartition_u(
	std::size_t v, std::size_t n, std::size_t* svv, std::size_t* ppp, double* aa1, double* aa2,
	double* bb1, double* bb2)
{
	std::size_t i;
	std::size_t j;
	std::size_t* ivv = new std::size_t[n];

	for (i = 0; i<n; i++)
		ivv[i] = 0;

	for (j = 0; j<v; j++)
	{
		i = toIndexPerm(n, ppp, svv, ivv);
		bb1[i] = aa1[j];
		bb2[i] = aa2[j];
		incIndex(n, svv, ivv);
	}

	delete[] ivv;
}

void listListVarsArrayHistoryPairsPartitionIndependent_u(
	double z, std::size_t v, std::size_t n, std::size_t* svv, std::size_t m, std::size_t r,
	std::size_t* lyy, std::size_t* syy, std::size_t* pppp, double* aa1, double* aa2,
	double* bb1, double* bb2)
{
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t a;
	double f;
	double x1 = 0;
	double x2 = 0;
	std::size_t** ppp = new std::size_t*[m];
	double* pxx1 = new double[r];
	double* pxx2 = new double[r];
	double** xx1 = new double*[m];
	double** xx2 = new double*[m];
	std::size_t* ivv = new std::size_t[n];
	std::size_t* iyy = new std::size_t[m];

	for (k = 1, a = lyy[0], ppp[0] = pppp; k<m; k++)
	{
		ppp[k] = pppp + a;
		a += lyy[k];
	}

	for (k = 1, a = syy[0], xx1[0] = pxx1, xx2[0] = pxx2; k<m; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += syy[k];
	}

	for (i = 0; i<r; i++)
	{
		pxx1[i] = 0.0;
		pxx2[i] = 0.0;
	}

	for (i = 0; i<n; i++)
	{
		ivv[i] = 0;
	}

	for (j = 0; j<v; j++)
	{
		for (k = 0; k<m; k++)
		{
			i = toIndexPerm(lyy[k], ppp[k], svv, ivv);
			xx1[k][i] += aa1[j];
			xx2[k][i] += aa2[j];
		}
		incIndex(n, svv, ivv);
	}

	if (z != 1)
	{
		f = 1 / z;
		for (k = 0; k<m; k++)
		{
			a = syy[k];
			for (i = 0; i<a; i++)
			{
				xx1[k][i] *= f;
				xx2[k][i] *= f;
			}
		}
	}

	for (i = 0; i<m; i++)
	{
		iyy[i] = 0;
	}

	if (z != 1)
	{
		for (j = 0; j<v; j++)
		{
			for (k = 0, x1 = z, x2 = z; k<m; k++)
			{
				a = iyy[k];
				x1 *= xx1[k][a];
				x2 *= xx2[k][a];
			}
			bb1[j] = x1;
			bb2[j] = x2;
			incIndex(m, syy, iyy);
		}
	}
	else
	{
		for (j = 0; j<v; j++)
		{
			for (k = 1, a = iyy[0], x1 = xx1[0][a], x2 = xx2[0][a]; k<m; k++)
			{
				a = iyy[k];
				x1 *= xx1[k][a];
				x2 *= xx2[k][a];
			}
			bb1[j] = x1;
			bb2[j] = x2;
			incIndex(m, syy, iyy);
		}
	}
	delete[] iyy;
	delete[] ivv;
	delete[] xx2;
	delete[] xx1;
	delete[] pxx2;
	delete[] pxx1;
	delete[] ppp;
}



std::size_t listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
	std::size_t pmax, double z, std::size_t v, std::size_t n, std::size_t* svv, std::size_t q, double y1,
	std::size_t* qm, std::size_t* ql, std::size_t* qs, std::size_t* qp, double* aa1, double* aa2,
	std::size_t* tt)
{
	std::size_t t = 0;
	double* bb1 = new double[v];
	double* bb2 = new double[v];
	double* ts1 = new double[pmax];
	double* ts2 = new double[pmax];
	long long* ts3 = new long long[pmax];
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	long long x3;
	std::size_t p;
	std::size_t m;
	std::size_t r;
	std::size_t i;
	double a2;
	double b2;
	double c;

	for (p = 0; p < q; p++)
	{
		m = qm[p];
		c = pow((double)v, 1.0 / ((double)m));

		for (r = 0, i = 0; i<m; i++)
		{
			r += (qs + n*p)[i];
		}

		listListVarsArrayHistoryPairsPartitionIndependent_u(z, v, n, svv, m, r, ql + n*p, qs + n*p, qp + n*p, aa1, aa2, bb1, bb2);

		for (a2 = 0.0, b2 = 0.0, i = 0; i<v; i++)
		{
			a2 += alngam(bb1[i] + 1.0);
			b2 += alngam(bb2[i] + 1.0);
		}

		if (t < pmax)
		{
			tt[t] = p;
			ts1[t] = (y1 - a2 + b2) / c;
			ts2[t] = b2;
			ts3[t] = -(long long)m;
			t++;
			if (t == pmax)
			{
				for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i < pmax; i++)
				{
					x1 = ts1[i];
					if (t1 > x1)
					{
						t1 = x1;
						t2 = ts2[i];
						t3 = ts3[i];
						tm = i;
					}
					else if (t1 == x1)
					{
						x2 = ts2[i];
						if (t2 > x2)
						{
							t2 = x2;
							t3 = ts3[i];
							tm = i;
						}
						else if (t2 == x2)
						{
							x3 = ts3[i];
							if (t3 > x3)
							{
								t3 = x3;
								tm = i;
							}
						}
					}
				}
			}
		}
		else
		{
			x1 = (y1 - a2 + b2) / c;
			x2 = b2;
			x3 = -(long long)m;
			if (t1 < x1 || (t1 == x1 && t2 < x2) || (t1 == x1 && t2 == x2 && t3 < x3))
			{
				tt[tm] = p;
				ts1[tm] = x1;
				ts2[tm] = x2;
				ts3[tm] = x3;
				for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i < pmax; i++)
				{
					x1 = ts1[i];
					if (t1 > x1)
					{
						t1 = x1;
						t2 = ts2[i];
						t3 = ts3[i];
						tm = i;
					}
					else if (t1 == x1)
					{
						x2 = ts2[i];
						if (t2 > x2)
						{
							t2 = x2;
							t3 = ts3[i];
							tm = i;
						}
						else if (t2 == x2)
						{
							x3 = ts3[i];
							if (t3 > x3)
							{
								t3 = x3;
								tm = i;
							}
						}
					}
				}
			}
		}
	}
	delete[] ts3;
	delete[] ts2;
	delete[] ts1;
	delete[] bb2;
	delete[] bb1;
	return t;
}



// parametersHistogramRepaVecsSetTuplePartitionTopByM_u ::
//   Integer -> Integer -> Integer -> HistogramRepa -> HistogramRepa -> Double -> Double -> ([[[VariableRepa]]],Integer)
std::tuple<std::unique_ptr<SizeListListList>, std::size_t> Alignment::parametersHistogramRepaVecsSetTuplePartitionTopByM_u(std::size_t mmax, std::size_t umax, std::size_t pmax, const HistogramRepa& aa, const HistogramRepa& aarr, double z, double y1)
{
	auto n = aa.dimension;
	auto svv = aa.shape;
	std::size_t v = 1;
	for (std::size_t i = 0; i < n; i++)
		v *= svv[i];
	auto mm = std::unique_ptr<SizeListList>{ new SizeListList{ SizeList{ 0 } } };
	for (std::size_t i = 2; i <= n; i++)
	{
		auto mm1 = std::unique_ptr<SizeListList>{ new SizeListList{} };
		for (auto& xx : *mm)
			for (std::size_t j = 0; j < i; j++)
				if (j < mmax && j <= *max_element(xx.begin(), xx.end()) + 1)
				{
					SizeList yy(xx);
					yy.push_back(j);
					mm1->push_back(yy);
				}
		mm = std::move(mm1);
	}
	auto q0 = mm->size() - 1;
	SizeSizeSetMapList qq0(q0);
	for (std::size_t i = 0; i < q0; i++)
	{
		auto& nn = qq0[i];
		auto& ll = (*mm)[i + 1];
		for (std::size_t j = 0; j < n; j++)
			nn[ll[j]].insert(j);
	}
	SizeListList rr0(q0);
	for (std::size_t i = 0; i < q0; i++)
	{
		auto& nn = qq0[i];
		auto& ll = rr0[i];
		for (auto& cc : nn)
		{
			std::size_t sv = 1;
			for (auto& p : cc.second)
				sv *= svv[p];
			ll.push_back(sv);
		}
	}
	auto tt1 = std::make_unique<SizeListListList>();
	tt1->reserve(mmax*pmax);
	std::size_t q1 = 0;
	for (std::size_t m = 2; m <= mmax; m++)
	{
		SizeListListList qq;
		SizeListList rr;
		for (std::size_t i = 0; i < q0; i++)
		{
			auto& nn0 = qq0[i];
			if (nn0.size() == m)
			{
				auto& ll = rr0[i];
				bool all = true;
				for (std::size_t k = 0; all && k < m; k++)
					all = ll[k] <= umax;
				if (all)
				{
					rr.push_back(ll);
					SizeListList nn;
					nn.reserve(m);
					for (auto& cc0 : nn0)
						nn.push_back(SizeList(cc0.second.begin(), cc0.second.end()));
					qq.push_back(nn);
				}
			}
		}
		std::size_t q = qq.size();
		q1 += q;
		std::size_t* qm = new std::size_t[q];
		std::size_t* ql = new std::size_t[q*n];
		std::size_t* qs = new std::size_t[q*n];
		std::size_t* qp = new std::size_t[q*n];
		std::size_t j = 0;
		for (std::size_t i = 0; i < q; i++)
		{
			qm[i] = m;
			auto& nn = qq[i];
			auto& ll = rr[i];
			for (std::size_t k = 0; k < m; k++)
			{
				auto& cc = nn[k];
				std::size_t c = cc.size();
				ql[i*n + k] = c;
				qs[i*n + k] = ll[k];
				for (std::size_t d = 0; d < c; d++)
				{
					qp[j] = cc[d];
					j++;
				}
			}
			for (std::size_t k = m; k < n; k++)
			{
				ql[i*n + k] = 0;
				qs[i*n + k] = 0;
			}
		}
		std::size_t* tt = new std::size_t[pmax];
		std::size_t t = listListVarsArrayHistoryPairsSetTuplePartitionTop_u(pmax, z, v, n, svv, q, y1, qm, ql, qs, qp, aa.arr, aarr.arr, tt);
		for (std::size_t i = 0; i < t; i++)
			tt1->push_back(qq[tt[i]]);
		delete[] tt;
		delete[] qp;
		delete[] qs;
		delete[] ql;
		delete[] qm;
	}
	return std::tuple<std::unique_ptr<SizeListListList>, std::size_t>(std::move(tt1), q1);
}

inline std::size_t toIndex(std::size_t n, std::size_t* svv, std::size_t* ivv)
{
	std::size_t k;
	std::size_t a;

	for (k = 1, a = ivv[0]; k<n; k++)
		a = svv[k] * a + ivv[k];
	return a;
}

inline std::size_t toIndexInsert(std::size_t u, std::size_t r, std::size_t q, std::size_t n, std::size_t* svv, std::size_t* ivv)
{
	std::size_t k;
	std::size_t a;

	for (k = 0, a = 0; k<u; k++)
		a = svv[k] * a + ivv[k];
	a = r*a + q;
	for (; k<n; k++)
		a = svv[k] * a + ivv[k];
	return a;
}

std::size_t arrayHistoryPairsRollMax_u(
	std::size_t v, std::size_t n, std::size_t* svvy, std::size_t d, std::size_t nd,
	double* aay, double* aaxy, double* bby, double* bbxy,
	std::size_t* ppm)
{
	std::size_t srchd = 0;
	double fm;

	std::size_t* ivv = new std::size_t[n];
	std::size_t* syy = new std::size_t[n];
	std::size_t* svv = new std::size_t[n];
	std::size_t* szz = new std::size_t[n];
	double* aa = new double[v];
	double* aax = new double[v];
	double* bb = new double[v];
	double* bbx = new double[v];
	double* aaz = new double[v];
	double* aaxz = new double[v];
	double* bbz = new double[v];
	double* bbxz = new double[v];
	double* ff = new double[nd];
	std::size_t minv;
	std::size_t vc;
	std::size_t* ppc = new std::size_t[nd];
	double fc;
	std::size_t ww;
	std::size_t sw;
	std::size_t tw;
	double fw;
	std::size_t x;
	std::size_t y;
	std::size_t r;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t q;
	std::size_t w;
	std::size_t p;
	std::size_t s;
	std::size_t is;
	std::size_t t;
	std::size_t it;
	std::size_t m;
	std::size_t u;
	double f;
	double c;

	for (j = 0; j < v; j++)
	{
		aa[j] = aay[j];
		aax[j] = aaxy[j];
		bb[j] = bby[j];
		bbx[j] = bbxy[j];
	}

	for (i = 0; i < nd; i++)
		ff[i] = 0.0;

	for (minv = 1, i = 0; i < n; i++)
	{
		minv *= 2;
		svv[i] = svvy[i];
		for (j = 0; j < d; j++)
		{
			p = d*i + j;
			ppm[p] = j;
			ppc[p] = j;
		}
	}

	vc = v;

	for (i = 0; i < n; i++)
		ivv[i] = 0;
	for (fc = 0.0, i = 0; i < vc; i++)
	{
		f = alngam(aa[i] + 1.0) - alngam(aax[i] + 1.0)
			- alngam(bb[i] + 1.0) + alngam(bbx[i] + 1.0);
		fc += f;
		for (k = 0; k < n; k++)
			ff[d*k + ivv[k]] += f;
		incIndex(n, svv, ivv);
	}

	m = n - 1;

	for (q = 0; vc > minv; q++)
	{
		for (x = 0, w = 0; w < n; w++)
		{
			r = svv[w];
			if (r > 2)
			{
				for (i = 0; i < w; i++)
					syy[i] = svv[i];
				for (; i < m; i++)
					syy[i] = svv[i + 1];
				p = d * w;
				y = vc / r;
				c = 1.0 / pow((double)(y*(r - 1)), 1.0 / ((double)n));
				for (s = 1; s < r; s++)
					for (t = 0; t < s; t++)
					{
						for (i = 0; i < m; i++)
							ivv[i] = 0;
						for (f = fc - ff[p + s] - ff[p + t], i = 0; i < y; i++)
						{
							is = toIndexInsert(w, r, s, m, syy, ivv);
							it = toIndexInsert(w, r, t, m, syy, ivv);
							f += alngam(aa[is] + aa[it] + 1.0) - alngam(aax[is] + aax[it] + 1.0)
								- alngam(bb[is] + bb[it] + 1.0) + alngam(bbx[is] + bbx[it] + 1.0);
							incIndex(m, syy, ivv);
						}
						f *= c;
						srchd++;
						if ((x == 0) || (f > fw))
						{
							x++;
							ww = w;
							sw = s;
							tw = t;
							fw = f;
						}
					}
			}
		}
		for (i = 0; i < n; i++)
			szz[i] = svv[i];
		r = svv[ww] - 1;
		szz[ww] = r;
		vc /= r + 1;
		vc *= r;
		for (i = 0; i < nd; i++)
			ff[i] = 0.0;
		for (i = 0; i < n; i++)
			ivv[i] = 0;
		for (fc = 0.0, j = 0; j < vc; j++)
		{
			u = ivv[ww];
			if ((u < tw) || (u > tw && u < sw))
			{
				i = toIndex(n, svv, ivv);
				aaz[j] = aa[i];
				aaxz[j] = aax[i];
				bbz[j] = bb[i];
				bbxz[j] = bbx[i];
			}
			else if (u == tw)
			{
				i = toIndex(n, svv, ivv);
				aaz[j] = aa[i];
				aaxz[j] = aax[i];
				bbz[j] = bb[i];
				bbxz[j] = bbx[i];
				ivv[ww] = sw;
				i = toIndex(n, svv, ivv);
				aaz[j] += aa[i];
				aaxz[j] += aax[i];
				bbz[j] += bb[i];
				bbxz[j] += bbx[i];
				ivv[ww] = tw;
			}
			else
			{
				ivv[ww]++;
				i = toIndex(n, svv, ivv);
				aaz[j] = aa[i];
				aaxz[j] = aax[i];
				bbz[j] = bb[i];
				bbxz[j] = bbx[i];
				ivv[ww]--;
			}
			f = alngam(aaz[j] + 1.0) - alngam(aaxz[j] + 1.0)
				- alngam(bbz[j] + 1.0) + alngam(bbxz[j] + 1.0);
			fc += f;
			for (k = 0; k < n; k++)
				ff[d*k + ivv[k]] += f;
			incIndex(n, szz, ivv);
		}
		for (j = 0; j < vc; j++)
		{
			aa[j] = aaz[j];
			aax[j] = aaxz[j];
			bb[j] = bbz[j];
			bbx[j] = bbxz[j];
		}
		svv[ww] = r;
		for (j = 0; j < d; j++)
		{
			p = d*ww + j;
			u = ppc[p];
			if (u == sw)
				ppc[p] = tw;
			else if (u > sw)
				ppc[p] = u - 1;
		}
		if (q == 0 || fw > fm)
		{
			fm = fw;
			for (i = 0; i < nd; i++)
				ppm[i] = ppc[i];
		}
	}

	delete[] ppc;
	delete[] ff;
	delete[] bbxz;
	delete[] bbz;
	delete[] aaxz;
	delete[] aaz;
	delete[] bbx;
	delete[] bb;
	delete[] aax;
	delete[] aa;
	delete[] szz;
	delete[] svv;
	delete[] syy;
	delete[] ivv;
	return srchd;
}



// histogramRepaVecsRollMax :: [[VariableRepa]] -> HistogramRepa -> HistogramRepa -> Double -> Double -> ([[VariableRepa]],Integer)
std::tuple<std::unique_ptr<SizeListList>, std::size_t> Alignment::histogramRepaVecsRollMax(const SizeListList& pp, const HistogramRepa& aa, const HistogramRepa& aarr, double z)
{
	auto n = aa.dimension;
	auto svv = aa.shape;
	std::size_t v = 1;
	for (std::size_t i = 0; i < n; i++)
		v *= svv[i];
	auto m = pp.size();
	std::size_t* ppp = new std::size_t[n];
	std::size_t* lyy = new std::size_t[m];
	std::size_t* syy = new std::size_t[m];
	std::size_t r = 0;
	std::size_t d = 0;
	std::size_t j = 0;
	for (std::size_t k = 0; k < m; k++)
	{
		auto& cc = pp[k];
		auto e = cc.size();
		lyy[k] = e;
		std::size_t sz = 1;
		for (std::size_t c = 0; c < e; c++)
		{
			auto p = cc[c];
			ppp[j] = p;
			j++;
			sz *= svv[p];
		}
		syy[k] = sz;
		r += sz;
		if (d < sz)
			d = sz;
	}
	double* aay = new double[v];
	double* bby = new double[v];
	double* aaxy = new double[v];
	double* bbxy = new double[v];
	listListVarsArrayHistoryPairsPartition_u(v, n, svv, ppp, aa.arr, aarr.arr, aay, bby);
	listListVarsArrayHistoryPairsPartitionIndependent_u(z, v, n, svv, m, r, lyy, syy, ppp, aa.arr, aarr.arr, aaxy, bbxy);
	std::size_t* ppm = new std::size_t[m*d];
	auto s = arrayHistoryPairsRollMax_u(v, m, syy, d, m*d, aay, aaxy, bby, bbxy, ppm);
	auto tt1 = std::make_unique<SizeListList>();
	tt1->reserve(m);
	for (std::size_t k = 0; k < m; k++)
	{
		SizeList cc;
		auto e = syy[k];
		cc.reserve(e);
		for (std::size_t c = 0; c < e; c++)
			cc.push_back(ppm[d*k + c]);
		tt1->push_back(cc);
	}
	delete[] ppm;
	delete[] bbxy;
	delete[] aaxy;
	delete[] bby;
	delete[] aay;
	delete[] syy;
	delete[] lyy;
	delete[] ppp;
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(tt1), s);
}

std::size_t listVarsListTuplesArrayHistoryEvientsAlignedExcludeHiddenTop_u(
	unsigned char dense,
	std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e,
	std::size_t z1, std::size_t z2,
	std::size_t ccl, std::size_t* ppccd, std::size_t* ppccu,
	std::size_t* ppww, std::size_t* ppdd,
	unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ppccx = new std::size_t[ccl];
	std::size_t ccx;
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	long long x3;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t jn;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t h;
	std::size_t a;
	std::size_t ok;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ccx = 0, h = 0; h < ccl; h++)
			if (ppccu[h] == pi)
			{
				ppccx[ccx] = ppccd[h];
				ccx++;
			}
			else if (ppccd[h] == pi)
			{
				ppccx[ccx] = ppccu[h];
				ccx++;
			}
		for (ij = 0; ij<d; ij++)
		{
			pj = pdd[ij];
			for (k = 0, ok = 1, u1 = 1; k<e; k++)
			{
				pk = pj[k];
				if (pk == pi)
					ok = 0;
				for (h = 0; h<ccx; h++)
					if (pk == ppccx[h])
						ok = 0;
				if (!ok)
					break;
				sk = svv[pk];
				sj[k] = sk;
				u1 *= sk;
			}
			u = u1*si;
			if (ok && u <= xmax)
			{
				s++;
				x4 = hash(n, e, pi, pj);
				for (i = 0; i<t; i++)
					if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
					{
						ok = 0;
						break;
					}
				if (!ok)
					continue;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				pk = pj[0];
				sk = sj[0];
				for (j = 0; j<z1; j++)
				{
					jn = n*j;
					for (k = 1, a = sk*phh1[jn + pi] + phh1[jn + pk]; k<e; k++)
						a = sj[k] * a + phh1[jn + pj[k]];
					aa[a] += 1.0;
				}
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				for (j = 0; j<z2; j++)
				{
					jn = n*j;
					for (k = 1, a = sk*phh2[jn + pi] + phh2[jn + pk]; k<e; k++)
						a = sj[k] * a + phh2[jn + pj[k]];
					aa[a] += f;
				}
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (k = 0; k<e; k++)
					yj[k] = 0;
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					y1 = zf*xx1[pi][i];
					y2 = zf*xx2[pi][i];
					for (j = 0; j<u1; j++)
					{
						x1 = y1;
						x2 = y2;
						for (k = 0; k<e; k++)
						{
							pk = pj[k];
							yk = yj[k];
							x1 *= xx1[pk][yk];
							x2 *= xx2[pk][yk];
						}
						a2 += alngam(x1 + 1.0);
						b2 += alngam(x2 + 1.0);
						incIndex(e, sj, yj);
					}
				}
				if (t<omax)
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					tww1[t] = ii;
					tww2[t] = ij;
					ts1[t] = x1;
					ts2[t] = x2;
					ts3[t] = x3;
					ts4[t] = x4;
					t++;
					if (t == omax)
						findm = 1;
				}
				else
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = ii;
						tww2[tm] = ij;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						ts4[tm] = x4;
						findm = 1;
					}
				}
				if (findm)
				{
					for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
					{
						x1 = ts1[i];
						if (t1>x1)
						{
							t1 = x1;
							t2 = ts2[i];
							t3 = ts3[i];
							tm = i;
						}
						else if (t1 == x1)
						{
							x2 = ts2[i];
							if (t2>x2)
							{
								t2 = x2;
								t3 = ts3[i];
								tm = i;
							}
							else if (t2 == x2)
							{
								x3 = ts3[i];
								if (t3>x3)
								{
									t3 = x3;
									tm = i;
								}
							}
						}
					}
					findm = 0;
				}
			}
		}
	}
	delete[] yj;
	delete[] sj;
	delete[] ts4;
	delete[] ppccx;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
	return t;
}

std::size_t listVarsListTuplesArrayHistoryVarientsAlignedExcludeHiddenTop_u(
	unsigned char dense,
	std::size_t xmax, std::size_t omax, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e,
	std::size_t z1, std::size_t z2,
	std::size_t ccl, std::size_t* ppccd, std::size_t* ppccu,
	std::size_t* ppww, std::size_t* ppdd,
	unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t* tww1, std::size_t* tww2, double* ts1, double* ts2, long long* ts3, std::size_t& s)
{
	std::size_t t = 0;
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ppccx = new std::size_t[ccl];
	std::size_t ccx;
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	double t2;
	long long t3;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	long long x3;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t qi;
	std::size_t* qj = new std::size_t[e];
	std::size_t qk;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t h;
	std::size_t a;
	std::size_t ok;

	s = 0;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
	{
		pi = ppww[ii];
		si = svv[pi];
		for (ccx = 0, h = 0; h < ccl; h++)
			if (ppccu[h] == pi)
			{
				ppccx[ccx] = ppccd[h];
				ccx++;
			}
			else if (ppccd[h] == pi)
			{
				ppccx[ccx] = ppccu[h];
				ccx++;
			}
		for (ij = 0; ij<d; ij++)
		{
			pj = pdd[ij];
			for (k = 0, ok = 1, u1 = 1; k<e; k++)
			{
				pk = pj[k];
				if (pk == pi)
					ok = 0;
				for (h = 0; h<ccx; h++)
					if (pk == ppccx[h])
						ok = 0;
				if (!ok)
					break;
				sk = svv[pk];
				sj[k] = sk;
				u1 *= sk;
			}
			u = u1*si;
			if (ok && u <= xmax)
			{
				s++;
				x4 = hash(n, e, pi, pj);
				for (i = 0; i<t; i++)
					if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
					{
						ok = 0;
						break;
					}
				if (!ok)
					continue;
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				pk = pj[0];
				sk = sj[0];
				qi = z1*pi;
				qk = z1*pk;
				for (k = 1; k<e; k++)
					qj[k] = z1*pj[k];
				for (j = 0; j<z1; j++)
				{
					for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
						a = sj[k] * a + phh1[qj[k] + j];
					aa[a] += 1.0;
				}
				for (a1 = 0.0, i = 0; i<u; i++)
					a1 += alngam(aa[i]);
				for (i = 0; i<u; i++)
					aa[i] = 1.0;
				qi = z2*pi;
				qk = z2*pk;
				for (k = 1; k<e; k++)
					qj[k] = z2*pj[k];
				for (j = 0; j<z2; j++)
				{
					for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
						a = sj[k] * a + phh2[qj[k] + j];
					aa[a] += f;
				}
				for (b1 = 0.0, i = 0; i<u; i++)
					b1 += alngam(aa[i]);
				for (k = 0; k<e; k++)
					yj[k] = 0;
				for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
				{
					y1 = zf*xx1[pi][i];
					y2 = zf*xx2[pi][i];
					for (j = 0; j<u1; j++)
					{
						x1 = y1;
						x2 = y2;
						for (k = 0; k<e; k++)
						{
							pk = pj[k];
							yk = yj[k];
							x1 *= xx1[pk][yk];
							x2 *= xx2[pk][yk];
						}
						a2 += alngam(x1 + 1.0);
						b2 += alngam(x2 + 1.0);
						incIndex(e, sj, yj);
					}
				}
				if (t<omax)
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					tww1[t] = ii;
					tww2[t] = ij;
					ts1[t] = x1;
					ts2[t] = x2;
					ts3[t] = x3;
					ts4[t] = x4;
					t++;
					if (t == omax)
						findm = 1;
				}
				else
				{
					if (dense)
					{
						c = pow((double)u, 1.0 / ((double)(e + 1)));
						x1 = (a1 - a2 - b1 + b2) / c;
						x2 = (b2 - b1) / c;
					}
					else
					{
						x1 = a1 - a2 - b1 + b2;
						x2 = b2 - b1;
					}
					x3 = -(long long)u;
					if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
					{
						tww1[tm] = ii;
						tww2[tm] = ij;
						ts1[tm] = x1;
						ts2[tm] = x2;
						ts3[tm] = x3;
						ts4[tm] = x4;
						findm = 1;
					}
				}
				if (findm)
				{
					for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
					{
						x1 = ts1[i];
						if (t1>x1)
						{
							t1 = x1;
							t2 = ts2[i];
							t3 = ts3[i];
							tm = i;
						}
						else if (t1 == x1)
						{
							x2 = ts2[i];
							if (t2>x2)
							{
								t2 = x2;
								t3 = ts3[i];
								tm = i;
							}
							else if (t2 == x2)
							{
								x3 = ts3[i];
								if (t3>x3)
								{
									t3 = x3;
									tm = i;
								}
							}
						}
					}
					findm = 0;
				}
			}
		}
	}
	delete[] yj;
	delete[] sj;
	delete[] qj;
	delete[] ts4;
	delete[] ppccx;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
	return t;
}

void listVarsListTuplesArrayHistoryVarientsAlignedExcludeHiddenTop_up(
	unsigned char dense,
	std::size_t xmax, std::size_t omax, std::size_t tint, std::size_t n, std::size_t* svv, std::size_t m, std::size_t d, std::size_t e,
	std::size_t z1, std::size_t z2,
	std::size_t ccl, std::size_t* ppccd, std::size_t* ppccu,
	std::size_t* ppww, std::size_t* ppdd,
	unsigned char* phh1, double* pxx1, unsigned char* phh2, double* pxx2,
	std::size_t th, std::size_t* tww1, std::size_t* tww2, double* ts1, std::size_t* t, std::size_t* s)
{
	std::size_t findm = 0;
	std::size_t** pdd = new std::size_t*[d];
	double* aa = new double[xmax];
	double** xx1 = new double*[n];
	double** xx2 = new double*[n];
	std::size_t* ppccx = new std::size_t[ccl];
	std::size_t ccx;
	std::size_t* ts4 = new std::size_t[omax];
	double zf = (double)z1;
	double f = (double)z1 / (double)z2;
	double t1;
	std::size_t tm;
	double x1;
	double x2;
	double y1;
	double y2;
	std::size_t x4;
	double c;
	double a1;
	double a2;
	double b1;
	double b2;
	std::size_t ii;
	std::size_t ij;
	std::size_t pi;
	std::size_t* pj;
	std::size_t pk;
	std::size_t qi;
	std::size_t* qj = new std::size_t[e];
	std::size_t qk;
	std::size_t si;
	std::size_t* sj = new std::size_t[e];
	std::size_t* yj = new std::size_t[e];
	std::size_t sk;
	std::size_t yk;
	std::size_t u;
	std::size_t u1;
	std::size_t i;
	std::size_t j;
	std::size_t k;
	std::size_t h;
	std::size_t a;
	std::size_t ok;
	std::size_t thomax = th*omax;

	for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
	{
		xx1[k] = pxx1 + a;
		xx2[k] = pxx2 + a;
		a += svv[k];
	}

	for (k = 0; k<d; k++)
		pdd[k] = ppdd + e*k;

	for (ii = 0; ii<m; ii++)
		if (ii % tint == th)
		{
			pi = ppww[ii];
			si = svv[pi];
			for (ccx = 0, h = 0; h < ccl; h++)
				if (ppccu[h] == pi)
				{
					ppccx[ccx] = ppccd[h];
					ccx++;
				}
				else if (ppccd[h] == pi)
				{
					ppccx[ccx] = ppccu[h];
					ccx++;
				}
			for (ij = 0; ij<d; ij++)
			{
				pj = pdd[ij];
				for (k = 0, ok = 1, u1 = 1; k<e; k++)
				{
					pk = pj[k];
					if (pk == pi)
						ok = 0;
					for (h = 0; h<ccx; h++)
						if (pk == ppccx[h])
							ok = 0;
					if (!ok)
						break;
					sk = svv[pk];
					sj[k] = sk;
					u1 *= sk;
				}
				u = u1*si;
				if (ok && u <= xmax)
				{
					s[th]++;
					x4 = hash(n, e, pi, pj);
					for (i = 0; i<t[th]; i++)
						if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[thomax + i]], pdd[tww2[thomax + i]]))
						{
							ok = 0;
							break;
						}
					if (!ok)
						continue;
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					pk = pj[0];
					sk = sj[0];
					qi = z1*pi;
					qk = z1*pk;
					for (k = 1; k<e; k++)
						qj[k] = z1*pj[k];
					for (j = 0; j<z1; j++)
					{
						for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
							a = sj[k] * a + phh1[qj[k] + j];
						aa[a] += 1.0;
					}
					for (a1 = 0.0, i = 0; i<u; i++)
						a1 += alngam(aa[i]);
					for (i = 0; i<u; i++)
						aa[i] = 1.0;
					qi = z2*pi;
					qk = z2*pk;
					for (k = 1; k<e; k++)
						qj[k] = z2*pj[k];
					for (j = 0; j<z2; j++)
					{
						for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
							a = sj[k] * a + phh2[qj[k] + j];
						aa[a] += f;
					}
					for (b1 = 0.0, i = 0; i<u; i++)
						b1 += alngam(aa[i]);
					for (k = 0; k<e; k++)
						yj[k] = 0;
					for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
					{
						y1 = zf*xx1[pi][i];
						y2 = zf*xx2[pi][i];
						for (j = 0; j<u1; j++)
						{
							x1 = y1;
							x2 = y2;
							for (k = 0; k<e; k++)
							{
								pk = pj[k];
								yk = yj[k];
								x1 *= xx1[pk][yk];
								x2 *= xx2[pk][yk];
							}
							a2 += alngam(x1 + 1.0);
							b2 += alngam(x2 + 1.0);
							incIndex(e, sj, yj);
						}
					}
					if (t[th]<omax)
					{
						if (dense)
						{
							c = pow((double)u, 1.0 / ((double)(e + 1)));
							x1 = (a1 - a2 - b1 + b2) / c;
						}
						else
						{
							x1 = a1 - a2 - b1 + b2;
						}
						tww1[thomax + t[th]] = ii;
						tww2[thomax + t[th]] = ij;
						ts1[thomax + t[th]] = x1;
						ts4[t[th]] = x4;
						t[th]++;
						if (t[th] == omax)
							findm = 1;
					}
					else
					{
						if (dense)
						{
							c = pow((double)u, 1.0 / ((double)(e + 1)));
							x1 = (a1 - a2 - b1 + b2) / c;
						}
						else
						{
							x1 = a1 - a2 - b1 + b2;
						}
						if (t1<x1)
						{
							tww1[thomax + tm] = ii;
							tww2[thomax + tm] = ij;
							ts1[thomax + tm] = x1;
							ts4[tm] = x4;
							findm = 1;
						}
					}
					if (findm)
					{
						for (t1 = ts1[thomax], tm = 0, i = 1; i<omax; i++)
						{
							x1 = ts1[thomax + i];
							if (t1>x1)
							{
								t1 = x1;
								tm = i;
							}
						}
						findm = 0;
					}
				}
			}
		}
	delete[] yj;
	delete[] sj;
	delete[] qj;
	delete[] ts4;
	delete[] ppccx;
	delete[] xx2;
	delete[] xx1;
	delete[] aa;
	delete[] pdd;
}


// parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u :: Integer -> Integer -> [(Variable, Variable)] -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa]], Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u(std::size_t wmax, std::size_t omax, const SizeSizePairList& cc, const SizeList& ww, const SizeListList& vdd, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	auto d = vdd.size();
	auto e = d ? vdd[0].size() : 0;
	std::size_t* pdd = new std::size_t[d*e];
	for (std::size_t i = 0; i < d; i++)
		for (std::size_t j = 0; j < e; j++)
			pdd[i*e + j] = mvv[vdd[i][j]];
	auto ccl = cc.size();
	std::size_t* pccd = new std::size_t[ccl];
	std::size_t* pccu = new std::size_t[ccl];
	for (std::size_t i = 0; i < ccl; i++)
	{
		pccd[i] = mvv[cc[i].first];
		pccu[i] = mvv[cc[i].second];
	}
	std::size_t* tww1 = new std::size_t[omax];
	std::size_t* tww2 = new std::size_t[omax];
	double* ts1 = new double[omax];
	double* ts2 = new double[omax];
	long long* ts3 = new long long[omax];
	std::size_t s = 0;
	std::size_t t = 0;
	if (hh.evient)
		t = listVarsListTuplesArrayHistoryEvientsAlignedExcludeHiddenTop_u(1, wmax, omax, n, svv, m, d, e, z, zrr, ccl, pccd, pccu, pww, pdd, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	else
		t = listVarsListTuplesArrayHistoryVarientsAlignedExcludeHiddenTop_u(1, wmax, omax, n, svv, m, d, e, z, zrr, ccl, pccd, pccu, pww, pdd, phh1, pxx1, phh2, pxx2, tww1, tww2, ts1, ts2, ts3, s);
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = 0; i < t; i++)
	{
		SizeList ll;
		ll.push_back(ww[tww1[i]]);
		for (auto& jj : vdd[tww2[i]])
			ll.push_back(jj);
		qq->push_back(DoubleSizeListPair(ts1[i], ll));
	}
	delete[] ts3;
	delete[] ts2;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pccu;
	delete[] pccd;
	delete[] pdd;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s);
}


// parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_up :: Integer -> Integer -> Integer -> [(Variable, Variable)] -> [VariableRepa] -> [[VariableRepa]] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> ([(Double, [VariableRepa]], Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_up(std::size_t wmax, std::size_t omax, std::size_t tint, const SizeSizePairList& cc, const SizeList& ww, const SizeListList& vdd, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	if (hh.evient)
		throw std::out_of_range("parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_up");
	auto n = hh.dimension;
	auto vhh = hh.vectorVar;
	auto& mvv = hh.mapVarInt();
	auto svv = hh.shape;
	auto z = hh.size;
	auto phh1 = hh.arr;
	auto zrr = hhrr.size;
	auto pxx1 = hhx.arr;
	auto phh2 = hhrr.arr;
	auto pxx2 = hhrrx.arr;
	auto m = ww.size();
	std::size_t tint1 = std::min(tint, m);
	std::size_t* pww = new std::size_t[m];
	for (std::size_t i = 0; i < m; i++)
		pww[i] = mvv[ww[i]];
	auto d = vdd.size();
	auto e = d ? vdd[0].size() : 0;
	std::size_t* pdd = new std::size_t[d*e];
	for (std::size_t i = 0; i < d; i++)
		for (std::size_t j = 0; j < e; j++)
			pdd[i*e + j] = mvv[vdd[i][j]];
	auto ccl = cc.size();
	std::size_t* pccd = new std::size_t[ccl];
	std::size_t* pccu = new std::size_t[ccl];
	for (std::size_t i = 0; i < ccl; i++)
	{
		pccd[i] = mvv[cc[i].first];
		pccu[i] = mvv[cc[i].second];
	}
	std::size_t* tww1 = new std::size_t[tint1*omax];
	std::size_t* tww2 = new std::size_t[tint1*omax];
	double* ts1 = new double[tint1*omax];
	std::size_t* s = new std::size_t[tint1];
	std::size_t* t = new std::size_t[tint1];
	std::vector<std::thread> threads;
	threads.reserve(tint1);
	for (std::size_t h = 0; h < tint1; h++)
	{
		s[h] = 0;
		t[h] = 0;
		threads.push_back(std::thread(listVarsListTuplesArrayHistoryVarientsAlignedExcludeHiddenTop_up, 1, wmax, omax, tint1, n, svv, m, d, e, z, zrr, ccl, pccd, pccu, pww, pdd, phh1, pxx1, phh2, pxx2, h, tww1, tww2, ts1, t, s));
	}
	for (auto& h : threads)
		h.join();
	std::size_t s1 = 0;
	DoubleSizePairList qq1;
	qq1.reserve(tint1*omax);
	SizeListList qq2;
	qq2.reserve(tint1*omax);
	std::size_t j = 0;
	SizeSetSet qq3;
	for (std::size_t h = 0; h < tint1; h++)
	{
		for (std::size_t i = 0; i < t[h]; i++)
		{
			SizeList ll;
			ll.push_back(ww[tww1[h*omax + i]]);
			for (auto& jj : vdd[tww2[h*omax + i]])
				ll.push_back(jj);
			SizeSet ss(ll.begin(), ll.end());
			if (qq3.find(ss) == qq3.end())
			{
				qq1.push_back(DoubleSizePair(ts1[h*omax + i], j));
				qq2.push_back(ll);
				qq3.insert(ss);
				j++;
			}
		}
		s1 += s[h];
	}
	std::sort(qq1.begin(), qq1.end());
	auto qq = std::make_unique<DoubleSizeListPairList>();
	for (std::size_t i = j - std::min(j, omax); i < j; i++)
		qq->push_back(DoubleSizeListPair(qq1[i].first, qq2[qq1[i].second]));
	delete[] s;
	delete[] t;
	delete[] ts1;
	delete[] tww2;
	delete[] tww1;
	delete[] pccu;
	delete[] pccd;
	delete[] pdd;
	delete[] pww;
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(qq), s1);
}

DecompFudSlicedRepa::DecompFudSlicedRepa(std::size_t listFudRepaSizeExpectedA) : _mapVarInt(0), fudRepasSize(0), listFudRepaSizeExpected(listFudRepaSizeExpectedA)
{
	if (listFudRepaSizeExpected > 0)
		fuds.reserve(listFudRepaSizeExpected);
}

DecompFudSlicedRepa::~DecompFudSlicedRepa()
{
	delete _mapVarInt;
}

std::ostream& operator<<(std::ostream& out, const FudSlicedStruct& fr)
{
	out << "(" << fr.parent << ",[" << fr.children << "],[";	
	bool first = true;
	for (auto& tr : fr.fud)
	{
		if (first)
			first = false;
		else
			out << ",";		
		out << *tr;		
	}
	out << "])";
	return out;
}

std::ostream& operator<<(std::ostream& out, const DecompFudSlicedRepa& dr)
{
	out << "(" << dr.listFudRepaSizeExpected << ",[";	

	bool first = true;
	for (auto& fr : dr.fuds)
	{
		if (first)
			first = false;
		else
			out << ",";		
		out << fr;				
	}
	out << "],"  << dr.fudRepasSize << ")";
	return out;
}

SizeSizeUMap& Alignment::DecompFudSlicedRepa::mapVarInt() const
{
	if (!_mapVarInt)
	{
		auto fudsSize = fuds.size();
		const_cast<DecompFudSlicedRepa*>(this)->_mapVarInt = new SizeSizeUMap(listFudRepaSizeExpected > fudsSize ? listFudRepaSizeExpected : fudsSize);
		for (std::size_t i = 0; i < fudsSize; i++)
		{	
			const_cast<DecompFudSlicedRepa*>(this)->_mapVarInt->insert_or_assign(fuds[i].parent, i);			
		}
	}
	return *_mapVarInt;
}

std::unique_ptr<DecompFudSlicedRepa> Alignment::applicationRepasDecompFudSlicedRepa_u(const ApplicationRepa& er)
{	
	auto dr = std::make_unique<DecompFudSlicedRepa>(treesSize(*er.slices));
	dr->fuds.resize(1);
	auto& vi = dr->mapVarInt();
	SizeSizeTreePairList nodesA {SizeSizeTreePair(0,er.slices)};
	SizeSizeTreePairList nodesB;
	while (nodesA.size())
	{
		for (auto sl : nodesA)
		{
			auto& fs = dr->fuds[vi[sl.first]];
			fs.children.reserve(sl.second->_list.size());
			for (auto& pp : sl.second->_list)
				fs.children.push_back(pp.first);
			TransformRepaPtrList kk;
			SizeSet vv(fs.children.begin(),fs.children.end());
			for (auto ll = er.fud->layers.rbegin(); ll != er.fud->layers.rend(); ll++) 
				for (auto& tr : *ll)
				{
					auto w = tr->derived;
					if (vv.find(w) != vv.end())
					{
						kk.push_back(tr);
						auto n = tr->dimension;
						auto xx = tr->vectorVar;
						for (std::size_t i = 0; i < n; i++)
						{
							auto x = xx[i];
							if (x != fs.parent)
								vv.insert(x);
						}					
					}
				}
			fs.fud.insert(fs.fud.end(),kk.rbegin(),kk.rend());
			for (auto& pp : sl.second->_list)
			{
				if (pp.second && pp.second->_list.size())
				{
					dr->fuds.resize(dr->fuds.size()+1);
					vi[pp.first] = dr->fuds.size()-1;
					dr->fuds.back().parent = pp.first;
					nodesB.push_back(pp);
				}
			}			
		}		
		nodesA = nodesB;
		nodesB.clear();
	}
	return dr;
}

std::unique_ptr<ApplicationRepa> Alignment::decompFudSlicedRepasApplicationRepa_u(const DecompFudSlicedRepa& dr)
{
	auto er = std::make_unique<ApplicationRepa>();
	if (!dr.fuds.size())
		return er;		
	er->fud = std::make_shared<FudRepa>();
	er->slices = std::make_shared<SizeTree>();
	std::vector<std::shared_ptr<SizeTree>> nodes;
	nodes.push_back(er->slices);
	std::map<std::size_t, std::size_t> mm;
	auto& vi = dr.mapVarInt();
	for (auto& fs : dr.fuds)
	{
		auto ss = nodes[mm[fs.parent]];	
		for (auto tr : fs.fud)
			er->fud->layers.push_back(TransformRepaPtrList {tr});
		for (auto sl : fs.children)
		{
			auto tt = std::make_shared<SizeTree>();
			ss->_list.push_back(SizeSizeTreePair(sl, tt));
			nodes.push_back(tt);
			mm[sl] = nodes.size()-1;
		}
	}
	return er;
}