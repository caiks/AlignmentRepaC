#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentApprox.h"
#include "AlignmentAeson.h"
#include "AlignmentRepa.h"
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/filereadstream.h"
#include <iomanip>
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <thread>
#include <chrono>
#include <ctime>

using namespace Alignment;
namespace js = rapidjson;
using namespace std;

void main()
{
    if (false)
    {
	auto suit = Variable("suit");
	auto rank = Variable("rank");

	auto vv = VarList{ suit,rank };

	HistogramRepa hr;
	hr.vectorVar = vv;

	cout << "hr.vectorVar" << endl
	    << hr.vectorVar << endl << endl;

	cout << "hr.mapVarInt()" << endl
	    << hr.mapVarInt() << endl << endl;
    }

    if (false)
    {
	auto regcart = histogramRegularCartesian_u;
	auto regsing = histogramRegularUnitSingleton_u;
	auto regdiag = histogramRegularUnitDiagonal_u;
	auto sys = histogramsSystemImplied;
	auto aarr = systemsHistogramsHistogramRepa_u;
	auto rraa = systemsHistogramRepasHistogram_u;
	auto arred = [](const HistogramRepa& aa, const VarList& vv)
	{
	    return setVarsHistogramRepasReduce_u(vv, aa);
	};

	auto aa = regdiag(2, 2);
	cout << "aa = regdiag(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	auto ar = aarr(*sys(*aa),*aa);
	cout << "ar = aarr(sys(aa),aa)" << endl;
	cout << "ar.vectorVar" << endl
	    << ar->vectorVar << endl << endl;
	cout << "ar.shape" << endl
	    << ar->shape << endl << endl;
	cout << "ar.arr" << endl
	    << ar->arr << endl << endl;
	cout << "rraa(sys(aa),ar)" << endl
	    << *rraa(*sys(*aa), *ar) << endl << endl;

	auto br = arred(*ar, VarList{ Variable(2),Variable(1) });
	cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{ Variable(2) });
	cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{});
	cout << "br = arred(ar, VarList{})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;

	aa = regcart(2, 2);
	cout << "aa = regcart(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	ar = aarr(*sys(*aa), *aa);
	cout << "ar = aarr(sys(aa),aa)" << endl;
	cout << "ar.vectorVar" << endl
	    << ar->vectorVar << endl << endl;
	cout << "ar.shape" << endl
	    << ar->shape << endl << endl;
	cout << "ar.arr" << endl
	    << ar->arr << endl << endl;
	cout << "rraa(sys(aa),ar)" << endl
	    << *rraa(*sys(*aa), *ar) << endl << endl;

	br = arred(*ar, VarList{ Variable(2),Variable(1) });
	cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{ Variable(2) });
	cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{});
	cout << "br = arred(ar, VarList{})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;

	aa = regsing(2, 2);
	cout << "aa = regsing(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	ar = aarr(*sys(*aa), *aa);
	cout << "ar = aarr(sys(aa),aa)" << endl;
	cout << "ar.vectorVar" << endl
	    << ar->vectorVar << endl << endl;
	cout << "ar.shape" << endl
	    << ar->shape << endl << endl;
	cout << "ar.arr" << endl
	    << ar->arr << endl << endl;
	cout << "rraa(sys(aa),ar)" << endl
	    << *rraa(*sys(*aa), *ar) << endl << endl;

	br = arred(*ar, VarList{ Variable(2),Variable(1) });
	cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{ Variable(2) });
	cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	br = arred(*ar, VarList{});
	cout << "br = arred(ar, VarList{})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
    }

    if (false)
    {
	auto uvars = systemsSetVar;
	auto cart = systemsSetVarsSetStateCartesian_u;
	typedef std::pair<int, ValList> IntValListPair;
	typedef std::vector<IntValListPair> IntValListPairList;
	auto llhh = [](const VarList& vv, const IntValListPairList& ee)
	{
	    std::vector<IdStatePair> ii;
	    for (auto& pp : ee)
	    {
		auto i = pp.first;
		auto& ll = pp.second;
		auto jj = std::vector<VarValPair>();
		for (int j = 0; j < ll.size(); j++)
		    jj.push_back(VarValPair(vv[j], ll[j]));
		ii.push_back(IdStatePair(Id(i), *listsState(jj)));
	    }
	    return listsHistory_u(ii);
	};
	auto hhll = historiesList;
	auto hvars = historiesSetVar;
	auto hsize = historiesSize;
	auto hred = [](const History& hh, const VarUSet& vv)
	{
	    return setVarsHistoriesReduce(vv, hh);
	};
	auto hhaa = historiesHistogram;
	auto aahh = histogramsHistory_u;
	auto aall = histogramsList;
	auto vars = histogramsSetVar;
	auto size = histogramsSize;
	auto trim = histogramsTrim;
	auto unit = setStatesHistogramUnit_u;
	auto norm = [](const Histogram& aa)
	{
	    return histogramsResize(1, aa);
	};
	auto ared = [](const Histogram& aa, const VarUSet& vv)
	{
	    return setVarsHistogramsReduce(vv, aa);
	};
	auto ind = histogramsIndependent;
	auto aarr = systemsHistogramsHistogramRepa_u;
	auto rraa = systemsHistogramRepasHistogram_u;
	auto arred = [](const HistogramRepa& aa, const VarList& vv)
	{
	    return setVarsHistogramRepasReduce_u(vv, aa);
	};

	auto pressure = Variable("pressure");
	auto cloud = Variable("cloud");
	auto wind = Variable("wind");
	auto rain = Variable("rain");
	auto low = Value("low");
	auto medium = Value("medium");
	auto high = Value("high");
	auto none = Value("none");
	auto light = Value("light");
	auto heavy = Value("heavy");
	auto strong = Value("strong");
	auto uu = listsSystem_u(std::vector<VarValSetPair>{
	    VarValSetPair(pressure, ValSet{ low,medium,high }),
		VarValSetPair(cloud, ValSet{ none,light,heavy }),
		VarValSetPair(wind, ValSet{ none,light,strong }),
		VarValSetPair(rain, ValSet{ none,light,heavy })});
	auto hh = llhh(VarList{ pressure, cloud, wind, rain }, IntValListPairList{
	    IntValListPair(1, ValList{ high, none, none, none }),
	    IntValListPair(2, ValList{ medium, light, none, light }),
	    IntValListPair(3, ValList{ high, none, light, none }),
	    IntValListPair(4, ValList{ low, heavy, strong, heavy }),
	    IntValListPair(5, ValList{ low, none, light, light }),
	    IntValListPair(6, ValList{ medium, none, light, light }),
	    IntValListPair(7, ValList{ low, heavy, light, heavy }),
	    IntValListPair(8, ValList{ high, none, light, none }),
	    IntValListPair(9, ValList{ medium, light, strong, heavy }),
	    IntValListPair(10, ValList{ medium, light, light, light }),
	    IntValListPair(11, ValList{ high, light, light, heavy }),
	    IntValListPair(12, ValList{ medium, none, none, none }),
	    IntValListPair(13, ValList{ medium, light, none, none }),
	    IntValListPair(14, ValList{ high, light, strong, light }),
	    IntValListPair(15, ValList{ medium, none, light, light }),
	    IntValListPair(16, ValList{ low, heavy, strong, heavy }),
	    IntValListPair(17, ValList{ low, heavy, light, heavy }),
	    IntValListPair(18, ValList{ high, none, none, none }),
	    IntValListPair(19, ValList{ low, light, none, light }),
	    IntValListPair(20, ValList{ high, none, none, none }) });
	auto aa = hhaa(*hh);

	cout << "rpln(aall(aa))" << endl;
	rpln(cout, sorted(*aall(*aa))); cout << endl;

	auto ar = aarr(*uu, *aa);
	cout << "ar = aarr(uu,aa)" << endl;
	cout << "ar.vectorVar" << endl
	    << ar->vectorVar << endl << endl;
	cout << "ar.shape" << endl
	    << ar->shape << endl << endl;
	cout << "ar.arr" << endl
	    << ar->arr << endl << endl;

	cout << "rpln(aall(trim(rraa(uu,ar))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ar))))); cout << endl;

	auto br = arred(*ar, VarList{ pressure, rain, cloud, wind });
	cout << "br = arred(ar, VarList{ pressure, rain, cloud, wind })" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;
	br = arred(*ar, VarList{ wind, cloud, rain });
	cout << "br = arred(ar, VarList{ wind, cloud, rain})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;
	br = arred(*ar, VarList{ rain, wind  });
	cout << "br = arred(ar, VarList{ rain, wind})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;
	br = arred(*ar, VarList{ rain });
	cout << "br = arred(ar, VarList{ rain})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;
	br = arred(*ar, VarList{});
	cout << "br = arred(ar, VarList{})" << endl;
	cout << "br.vectorVar" << endl
	    << br->vectorVar << endl << endl;
	cout << "br.shape" << endl
	    << br->shape << endl << endl;
	cout << "br.arr" << endl
	    << br->arr << endl << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;

    }

    if (false)
    {
	auto suit = Variable("suit");
	auto rank = Variable("rank");

	auto vv = VarList{ suit,rank };

	HistogramRepaVec hr;
	hr.vectorVar = vv;

	cout << "hr.vectorVar" << endl
	    << hr.vectorVar << endl << endl;

	cout << "hr.mapVarInt()" << endl
	    << hr.mapVarInt() << endl << endl;
    }

    if (true)
    {
	auto regcart = histogramRegularCartesian_u;
	auto regsing = histogramRegularUnitSingleton_u;
	auto regdiag = histogramRegularUnitDiagonal_u;
	auto sys = histogramsSystemImplied;
	auto araa = systemsHistogramRepasHistogram_u;
	auto aahr = [](const System& uu, const Histogram& aa)
	{
	    return systemsHistoriesHistoryRepa_u(uu, *histogramsHistory_u(aa));
	};
	auto hraa = [](const System& uu, const HistoryRepa& hr)
	{
	    return historiesHistogram(*systemsHistoryRepasHistory_u(uu,hr));
	};
	auto hrhh = [](const System& uu, const HistoryRepa& hr)
	{
	    return systemsHistoryRepasHistory_u(uu, hr);
	};
	auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
	{
	    return eventsHistoryRepasHistoryRepaSelection_u(ll, hr);
	};
	auto hrhrred = [](const HistoryRepa& hr, const VarList& kk)
	{
	    return setVarsHistoryRepasHistoryRepaReduced_u(kk, hr);
	};
	auto hrred = [](const HistoryRepa& hr, const VarList& kk)
	{
	    return setVarsHistoryRepasReduce_u(1.0,kk, hr);
	};

	auto aa = regdiag(2, 2);
	cout << "aa = regdiag(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	auto hr = aahr(*sys(*aa), *aa);
	cout << "hr = aahr(sys(aa),aa)" << endl;
	cout << "hr.vectorVar" << endl
	    << hr->vectorVar << endl << endl;
	cout << "hr.size" << endl
	    << hr->size << endl << endl;
	cout << "hr.shape" << endl
	    << hr->shape << endl << endl;
//	cout << "hr.arr" << endl
//	    << hr->arr << endl << endl;
	cout << "hraa(sys(aa),hr)" << endl
	    << *hraa(*sys(*aa), *hr) << endl << endl;

	aa = regcart(2, 2);
	cout << "aa = regcart(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	hr = aahr(*sys(*aa), *aa);
	cout << "hr = aahr(sys(aa),aa)" << endl;
	cout << "hr.vectorVar" << endl
	    << hr->vectorVar << endl << endl;
	cout << "hr.size" << endl
	    << hr->size << endl << endl;
	cout << "hr.shape" << endl
	    << hr->shape << endl << endl;
//	cout << "hr.arr" << endl
//	    << hr->arr << endl << endl;
	cout << "hraa(sys(aa),hr)" << endl
	    << *hraa(*sys(*aa), *hr) << endl << endl;

	aa = regsing(2, 2);
	cout << "aa = regsing(2, 2)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	hr = aahr(*sys(*aa), *aa);
	cout << "hr = aahr(sys(aa),aa)" << endl;
	cout << "hr.vectorVar" << endl
	    << hr->vectorVar << endl << endl;
	cout << "hr.size" << endl
	    << hr->size << endl << endl;
	cout << "hr.shape" << endl
	    << hr->shape << endl << endl;
//	cout << "hr.arr" << endl
//	    << hr->arr << endl << endl;
	cout << "hraa(sys(aa),hr)" << endl
	    << *hraa(*sys(*aa), *hr) << endl << endl;

	aa = regcart(257,1);
	cout << "aa = regcart(257,1)" << endl;
	cout << "aa" << endl
	    << *aa << endl << endl;

	try
	{
	    hr = aahr(*sys(*aa), *aa);
	}
	catch (const std::out_of_range& e)
	{
	    cout << "caught out_of_range: " << e.what() << endl << endl;
	}

	aa = regcart(3, 2);
	cout << "aa = regcart(3, 2)" << endl;

	hr = aahr(*sys(*aa), *aa);
	cout << "hr = aahr(sys(aa),aa)" << endl;
	cout << "hrhh(sys(aa),hr)" << endl
	    << *hrhh(*sys(*aa), *hr) << endl << endl;
	cout << "hrhh(sys(aa),hrsel(hr,SizeList{}))" << endl
	    << *hrhh(*sys(*aa), *hrsel(*hr, SizeList{})) << endl << endl;
	cout << "hrhh(sys(aa),hrsel(hr,SizeList{0}))" << endl
	    << *hrhh(*sys(*aa), *hrsel(*hr, SizeList{0})) << endl << endl;
	cout << "hrhh(sys(aa),hrsel(hr,SizeList{0,4,8}))" << endl
	    << *hrhh(*sys(*aa), *hrsel(*hr, SizeList{ 0,4,8 })) << endl << endl;

	cout << "hraa(sys(aa),hrhrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
	    << *hraa(*sys(*aa), *hrhrred(*hr, VarList{ Variable(2),Variable(1) })) << endl << endl;
	cout << "hraa(sys(aa),hrhrred(hr, VarList{ Variable(1) }))" << endl
	    << *hraa(*sys(*aa), *hrhrred(*hr, VarList{ Variable(1) })) << endl << endl;
	cout << "hraa(sys(aa),hrhrred(hr, VarList{ }))" << endl
	    << *hraa(*sys(*aa), *hrhrred(*hr, VarList{ })) << endl << endl;

	cout << "araa(sys(aa),hrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
	    << *araa(*sys(*aa), *hrred(*hr, VarList{ Variable(2),Variable(1) })) << endl << endl;
	cout << "araa(sys(aa),hrred(hr, VarList{ Variable(1) }))" << endl
	    << *araa(*sys(*aa), *hrred(*hr, VarList{ Variable(1) })) << endl << endl;
	cout << "araa(sys(aa),hrred(hr, VarList{  }))" << endl
	    << *araa(*sys(*aa), *hrred(*hr, VarList{ })) << endl << endl;

    }

    if (true)
    {
	auto uvars = systemsSetVar;
	auto cart = systemsSetVarsSetStateCartesian_u;
	typedef std::pair<int, ValList> IntValListPair;
	typedef std::vector<IntValListPair> IntValListPairList;
	auto llhh = [](const VarList& vv, const IntValListPairList& ee)
	{
	    std::vector<IdStatePair> ii;
	    for (auto& pp : ee)
	    {
		auto i = pp.first;
		auto& ll = pp.second;
		auto jj = std::vector<VarValPair>();
		for (int j = 0; j < ll.size(); j++)
		    jj.push_back(VarValPair(vv[j], ll[j]));
		ii.push_back(IdStatePair(Id(i), *listsState(jj)));
	    }
	    return listsHistory_u(ii);
	};
	auto hhll = historiesList;
	auto hvars = historiesSetVar;
	auto hsize = historiesSize;
	auto hred = [](const History& hh, const VarUSet& vv)
	{
	    return setVarsHistoriesReduce(vv, hh);
	};
	auto hhaa = historiesHistogram;
	auto aahh = histogramsHistory_u;
	auto aall = histogramsList;
	auto vars = histogramsSetVar;
	auto size = histogramsSize;
	auto trim = histogramsTrim;
	auto unit = setStatesHistogramUnit_u;
	auto norm = [](const Histogram& aa)
	{
	    return histogramsResize(1, aa);
	};
	auto ared = [](const Histogram& aa, const VarUSet& vv)
	{
	    return setVarsHistogramsReduce(vv, aa);
	};
	auto ind = histogramsIndependent;
	auto aarr = systemsHistogramsHistogramRepa_u;
	auto rraa = systemsHistogramRepasHistogram_u;
	auto arred = [](const HistogramRepa& aa, const VarList& vv)
	{
	    return setVarsHistogramRepasReduce_u(vv, aa);
	};
	auto aahr = [](const System& uu, const Histogram& aa)
	{
	    return systemsHistoriesHistoryRepa_u(uu, *histogramsHistory_u(aa));
	};
	auto hhhr = [](const System& uu, const History& hh)
	{
	    return systemsHistoriesHistoryRepa_u(uu, hh);
	};
	auto hraa = [](const System& uu, const HistoryRepa& hr)
	{
	    return historiesHistogram(*systemsHistoryRepasHistory_u(uu, hr));
	};
	auto hrhh = [](const System& uu, const HistoryRepa& hr)
	{
	    return systemsHistoryRepasHistory_u(uu, hr);
	};
	auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
	{
	    return eventsHistoryRepasHistoryRepaSelection_u(ll, hr);
	};
	auto hrhrred = [](const HistoryRepa& hr, const VarList& kk)
	{
	    return setVarsHistoryRepasHistoryRepaReduced_u(kk, hr);
	};
	auto hrred = [](const HistoryRepa& hr, const VarList& kk)
	{
	    return setVarsHistoryRepasReduce_u(1.0, kk, hr);
	};

	auto pressure = Variable("pressure");
	auto cloud = Variable("cloud");
	auto wind = Variable("wind");
	auto rain = Variable("rain");
	auto low = Value("low");
	auto medium = Value("medium");
	auto high = Value("high");
	auto none = Value("none");
	auto light = Value("light");
	auto heavy = Value("heavy");
	auto strong = Value("strong");
	auto uu = listsSystem_u(std::vector<VarValSetPair>{
	    VarValSetPair(pressure, ValSet{ low,medium,high }),
		VarValSetPair(cloud, ValSet{ none,light,heavy }),
		VarValSetPair(wind, ValSet{ none,light,strong }),
		VarValSetPair(rain, ValSet{ none,light,heavy })});
	auto hh = llhh(VarList{ pressure, cloud, wind, rain }, IntValListPairList{
	    IntValListPair(1, ValList{ high, none, none, none }),
	    IntValListPair(2, ValList{ medium, light, none, light }),
	    IntValListPair(3, ValList{ high, none, light, none }),
	    IntValListPair(4, ValList{ low, heavy, strong, heavy }),
	    IntValListPair(5, ValList{ low, none, light, light }),
	    IntValListPair(6, ValList{ medium, none, light, light }),
	    IntValListPair(7, ValList{ low, heavy, light, heavy }),
	    IntValListPair(8, ValList{ high, none, light, none }),
	    IntValListPair(9, ValList{ medium, light, strong, heavy }),
	    IntValListPair(10, ValList{ medium, light, light, light }),
	    IntValListPair(11, ValList{ high, light, light, heavy }),
	    IntValListPair(12, ValList{ medium, none, none, none }),
	    IntValListPair(13, ValList{ medium, light, none, none }),
	    IntValListPair(14, ValList{ high, light, strong, light }),
	    IntValListPair(15, ValList{ medium, none, light, light }),
	    IntValListPair(16, ValList{ low, heavy, strong, heavy }),
	    IntValListPair(17, ValList{ low, heavy, light, heavy }),
	    IntValListPair(18, ValList{ high, none, none, none }),
	    IntValListPair(19, ValList{ low, light, none, light }),
	    IntValListPair(20, ValList{ high, none, none, none }) });
	cout << "rpln(hhll(hh))" << endl;
	rpln(cout, sorted(*hhll(*hh))); cout << endl;

	auto hr = hhhr(*uu, *hh);
	cout << "hr = hhhr(uu,hh)" << endl;
	cout << "hraa(uu,hr)" << endl
	    << *hraa(*uu,*hr) << endl << endl;

	cout << "hrhh(uu,hrsel(hr,SizeList{}))" << endl
	    << *hrhh(*uu, *hrsel(*hr, SizeList{})) << endl << endl;
	cout << "hrhh(uu,hrsel(hr,SizeList{0}))" << endl
	    << *hrhh(*uu, *hrsel(*hr, SizeList{ 0 })) << endl << endl;
	cout << "hrhh(uu,hrsel(hr,SizeList{0,4,8}))" << endl
	    << *hrhh(*uu, *hrsel(*hr, SizeList{ 0,4,8 })) << endl << endl;

	auto hr1 = hrhrred(*hr, VarList{ pressure, rain, cloud, wind });
	cout << "hr1 = hrhrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
	cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
	rpln(cout, sorted(*aall(*trim(*hraa(*uu,*hr1))))); cout << endl;
	hr1 = hrhrred(*hr, VarList{ wind, cloud, rain });
	cout << "hr1 = hrhrred(hr, VarList{ wind, cloud, rain })" << endl;
	cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
	rpln(cout, sorted(*aall(*trim(*hraa(*uu, *hr1))))); cout << endl;
	hr1 = hrhrred(*hr, VarList{ rain, wind });
	cout << "hr1 = hrhrred(hr, VarList{ rain, wind })" << endl;
	cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
	rpln(cout, sorted(*aall(*trim(*hraa(*uu, *hr1))))); cout << endl;
	hr1 = hrhrred(*hr, VarList{ rain });
	cout << "hr1 = hrhrred(hr, VarList{ rain })" << endl;
	cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
	rpln(cout, sorted(*aall(*trim(*hraa(*uu, *hr1))))); cout << endl;
	hr1 = hrhrred(*hr, VarList{ });
	cout << "hr1 = hrhrred(hr, VarList{ })" << endl;
	cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
	rpln(cout, sorted(*aall(*trim(*hraa(*uu, *hr1))))); cout << endl;

	auto br = hrred(*hr, VarList{ pressure, rain, cloud, wind });
	cout << "br = hrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu,*br))))); cout << endl;
	br = hrred(*hr, VarList{ wind, cloud, rain });
	cout << "br = hrred(hr, VarList{ wind, cloud, rain})" << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu,*br))))); cout << endl;
	br = hrred(*hr, VarList{ rain, wind });
	cout << "br = hrred(hr, VarList{ rain, wind})" << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu,*br))))); cout << endl;
	br = hrred(*hr, VarList{ rain });
	cout << "br = hrred(hr, VarList{ rain})" << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;
	br = hrred(*hr, VarList{});
	cout << "br = hrred(hr, VarList{})" << endl;
	cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
	rpln(cout, sorted(*aall(*trim(*rraa(*uu, *br))))); cout << endl;

    }

    if (true)
    {
	auto lluu = listsSystem_u;
	auto cart = systemsSetVarsSetStateCartesian_u;
	auto llss = listsState;
	auto sys = histogramsSystemImplied;
	auto unit = setStatesHistogramUnit_u;
	auto aall = histogramsList;
	auto size = histogramsSize;
	auto resize = histogramsResize;
	auto norm = [](const Histogram& aa)
	{
	    return histogramsResize(1, aa);
	};
	auto add = pairHistogramsAdd_u;
	auto scalar = histogramScalar_u;
	auto regsing = histogramRegularUnitSingleton_u;
	auto regdiag = histogramRegularUnitDiagonal_u;
	auto regcart = histogramRegularCartesian_u;
	auto ared = [](const Histogram& aa, const VarUSet& vv)
	{
	    return setVarsHistogramsReduce(vv, aa);
	};
	auto llhh = [llss](const VarList& vv, const IntValListPairList& ee)
	{
	    std::vector<IdStatePair> ii;
	    for (auto& pp : ee)
	    {
		auto i = pp.first;
		auto& ll = pp.second;
		auto jj = std::vector<VarValPair>();
		for (int j = 0; j < ll.size(); j++)
		    jj.push_back(VarValPair(vv[j], ll[j]));
		ii.push_back(IdStatePair(Id(i), *llss(jj)));
	    }
	    return listsHistory_u(ii);
	};
	auto hhaa = historiesHistogram;
	auto trans = histogramsSetVarsTransform;
	auto ttaa = transformsHistogram;
	auto und = transformsUnderlying;
	auto der = transformsDerived;
	auto tmul = [](const Histogram& aa, const Transform& tt)
	{
	    return transformsHistogramsApply(tt, aa);
	};
	auto tttr = systemsTransformsTransformRepa_u;

	auto suit = Variable("suit");
	auto rank = Variable("rank");
	auto vv = VarUSet{ suit,rank };
	auto hearts = Value("hearts");
	auto clubs = Value("clubs");
	auto diamonds = Value("diamonds");
	auto spades = Value("spades");
	auto wws = ValSet{ hearts,clubs,diamonds,spades };
	auto jack = Value("J");
	auto queen = Value("Q");
	auto king = Value("K");
	auto ace = Value("A");
	auto wwr = ValSet{ jack,queen,king,ace };
	for (int i = 2; i <= 10; i++)
	    wwr.insert(Value(i));
	auto uu = lluu(std::vector<VarValSetPair>{VarValSetPair(suit, wws), VarValSetPair(rank, wwr)});

	cout << "uu" << endl
	    << *uu << endl << endl;

	cout << "vv" << endl
	    << sorted(vv) << endl << endl;

	auto aa = unit(*cart(*uu, vv));
	cout << "rpln(aall(aa))" << endl;
	rpln(cout, sorted(*aall(*aa))); cout << endl;

	auto colour = Variable("colour");
	auto red = Value("red");
	auto black = Value("black");

	auto xx = hhaa(*llhh(VarList{ suit, colour }, IntValListPairList{
	    IntValListPair(1, ValList{ hearts, red }),
	    IntValListPair(2, ValList{ clubs, black }),
	    IntValListPair(3, ValList{ diamonds, red }),
	    IntValListPair(4, ValList{ spades, black }) }));

	cout << "rpln(aall(xx))" << endl;
	rpln(cout, sorted(*aall(*xx))); cout << endl;

	auto ww = VarUSet{ colour };

	auto tt = trans(*xx, ww);

	cout << "trans(xx,ww)" << endl
	    << *tt << endl << endl;

	cout << "rpln(aall(ttaa(tt)))" << endl;
	rpln(cout, sorted(*aall(ttaa(*tt)))); cout << endl;

	cout << "und(tt)" << endl
	    << sorted(*und(*tt)) << endl << endl;

	cout << "der(tt)" << endl
	    << sorted(der(*tt)) << endl << endl;

	cout << "rpln(aall(tmul(aa, tt)))" << endl;
	rpln(cout, sorted(*aall(*tmul(*aa, *tt)))); cout << endl;

	auto uu1 = sys(tt->histogram_u());
	auto tr = tttr(*uu1, *tt);
	cout << "tr = tttr(uu1,tt)" << endl;
	cout << "tr.vectorVar" << endl
	    << tr->vectorVar << endl << endl;
	cout << "tr.derived" << endl
	    << *tr->derived << endl << endl;
	cout << "tr.valency" << endl
	    << (std::size_t)(tr->valency) << endl << endl;
	cout << "tr.shape" << endl
	    << tr->shape << endl << endl;
//	cout << "tr.arr" << endl
//	    << tr->arr << endl << endl;

    }

}
