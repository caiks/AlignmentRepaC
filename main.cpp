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

    if (false)
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
	auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
	{
	    return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
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

	auto bb = regcart(1,1);
	cout << "bb = regcart(1,1)" << endl;
	auto ss = aahr(*sys(*bb), *bb);
	cout << "ss = aahr(sys(bb),bb)" << endl;
	cout << "hrhh(sys(aa),hrhrsel(hr,ss))" << endl
	    << *hrhh(*sys(*aa),*hrhrsel(*hr, *ss)) << endl << endl;
	bb = regcart(2, 1);
	cout << "bb = regcart(2,1)" << endl;
	ss = aahr(*sys(*bb), *bb);
	cout << "ss = aahr(sys(bb),bb)" << endl;
	cout << "hrhh(sys(aa),hrhrsel(hr,ss))" << endl
	    << *hrhh(*sys(*aa), *hrhrsel(*hr, *ss)) << endl << endl;
	bb = regcart(3, 1);
	cout << "bb = regcart(3,1)" << endl;
	ss = aahr(*sys(*bb), *bb);
	cout << "ss = aahr(sys(bb),bb)" << endl;
	cout << "hrhh(sys(aa),hrhrsel(hr,ss))" << endl
	    << *hrhh(*sys(*aa), *hrhrsel(*hr, *ss)) << endl << endl;

	hr = aahr(*sys(*aa), *aa);
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

    if (false)
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
	auto trtt = systemsTransformRepasTransform_u;

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
	cout << "trtt(uu1,tr)" << endl
	    << *trtt(*uu1,*tr) << endl << endl;

    }

    if (false)
    {
	auto uvars = systemsSetVar;
	auto uunion = pairSystemsUnion;
	auto cart = systemsSetVarsSetStateCartesian_u;
	auto llss = listsState;
	auto vol = systemsSetVarsVolume_u;
	auto ssplit = [](const VarUSet& vv, const Histogram& aa)
	{
	    return setVarsSetStatesSplit(vv, *histogramsStates(aa));
	};
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
	auto hhll = historiesList;
	auto hvars = historiesSetVar;
	auto hsize = historiesSize;
	auto hred = [](const History& hh, const VarUSet& vv)
	{
	    return setVarsHistoriesReduce(vv, hh);
	};
	auto hhaa = historiesHistogram;
	auto aahh = histogramsHistory_u;
	auto ind = histogramsIndependent;
	auto ent = histogramsEntropy;
	auto lent = [size, ent, ared](const VarUSet& vv, const Histogram& aa)
	{
	    return size(aa).getDouble() * (ent(aa) - ent(*ared(aa, vv)));
	};
	auto algn = histogramsAlignment;
	auto trans = [](std::unique_ptr<Histogram>& xx, const VarUSet& ww)
	{
	    return std::make_shared<Transform>(std::move(xx), ww);
	};
	auto ttaa = transformsHistogram;
	auto und = transformsUnderlying;
	auto der = transformsDerived;
	auto tmul = [](const Histogram& aa, const Transform& tt)
	{
	    return transformsHistogramsApply(tt, aa);
	};
	auto lltt = [llss, trans](const VarList& kk, const VarList& ww, const ValListList& qq)
	{
	    VarList vv(kk.begin(), kk.end());
	    vv.insert(vv.end(), ww.begin(), ww.end());
	    std::vector<StateRationalPair> ii;
	    for (auto& ll : qq)
	    {
		auto jj = std::vector<VarValPair>();
		for (int j = 0; j < ll.size(); j++)
		    jj.push_back(VarValPair(vv[j], ll[j]));
		ii.push_back(StateRationalPair(*llss(jj), 1));
	    }
	    return trans(std::make_unique<Histogram>(ii), VarUSet(ww.begin(), ww.end()));
	};
	auto llff = setTransformsFud_u;
	auto fhis = fudsSetHistogram;
	auto fvars = fudsSetVar;
	auto fder = fudsDerived;
	auto fund = fudsUnderlying;
	auto fftt = fudsTransform;
	auto fsys = fudsSystemImplied;
	auto dep = fudsSetVarsDepends_u;
	auto hraa = [](const System& uu, const HistoryRepa& hr)
	{
	    return historiesHistogram(*systemsHistoryRepasHistory_u(uu, hr));
	};
	auto hhhr = [](const System& uu, const History& hh)
	{
	    return systemsHistoriesHistoryRepa_u(uu, hh);
	};
	auto hrhh = [](const System& uu, const HistoryRepa& hr)
	{
	    return systemsHistoryRepasHistory_u(uu, hr);
	};
	auto fffr = systemsFudsFudRepa_u;
	auto frff = systemsFudRepasFud_u;
	auto frmul = historyRepasFudRepasMultiply_u;

	auto pressure = Variable("pressure");
	auto cloud = Variable("cloud");
	auto wind = Variable("wind");
	auto rain = Variable("rain");
	auto cloud_and_wind = Variable("cloud_and_wind");
	auto cloud_and_pressure = Variable("cloud_and_pressure");
	auto cloud_wind_pressure = Variable("cloud_wind_pressure");
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
	auto vv = uvars(*uu);
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

	cout << "size(aa)" << endl
	    << size(*aa) << endl << endl;

	auto ttcw = lltt(VarList{ cloud, wind }, VarList{ cloud_and_wind }, ValListList{
	    ValList{ none, none, none },
	    ValList{ none, light, light },
	    ValList{ none, strong, light },
	    ValList{ light, none, light },
	    ValList{ light, light, light },
	    ValList{ light, strong, light },
	    ValList{ heavy, none, strong },
	    ValList{ heavy, light, strong },
	    ValList{ heavy, strong, strong } });

	auto ttcp = lltt(VarList{ cloud, pressure }, VarList{ cloud_and_pressure }, ValListList{
	    ValList{ none, high, none },
	    ValList{ none, medium, light },
	    ValList{ none, low, light },
	    ValList{ light, high, light },
	    ValList{ light, medium, light },
	    ValList{ light, low, light },
	    ValList{ heavy, high, strong },
	    ValList{ heavy, medium, strong },
	    ValList{ heavy, low, strong } });

	auto ff = llff(TransformPtrList{ ttcw, ttcp });

	cout << "fder(ff)" << endl
	    << sorted(*fder(*ff)) << endl << endl;

	cout << "fund(ff)" << endl
	    << sorted(*fund(*ff)) << endl << endl;

	cout << "fsys(ff)" << endl
	    << *fsys(*ff) << endl << endl;

	cout << "der(fftt(ff))" << endl
	    << sorted(der(*fftt(*ff))) << endl << endl;

	cout << "und(fftt(ff))" << endl
	    << sorted(*und(*fftt(*ff))) << endl << endl;

	auto ttcwp = lltt(VarList{ cloud_and_wind,cloud_and_pressure }, VarList{ cloud_wind_pressure }, ValListList{
	    ValList{ none, none, none },
	    ValList{ none, light, none },
	    ValList{ none, strong, none },
	    ValList{ light, none, none },
	    ValList{ light, light, light },
	    ValList{ light, strong, light },
	    ValList{ strong, none, none },
	    ValList{ strong, light, light },
	    ValList{ strong, strong, strong } });

	auto gg = llff(TransformPtrList{ ttcw, ttcp, ttcwp });

	cout << "fder(gg)" << endl
	    << sorted(*fder(*gg)) << endl << endl;

	cout << "fund(gg)" << endl
	    << sorted(*fund(*gg)) << endl << endl;

	cout << "fsys(gg)" << endl
	    << *fsys(*gg) << endl << endl;

	cout << "der(fftt(gg))" << endl
	    << sorted(der(*fftt(*gg))) << endl << endl;

	cout << "und(fftt(gg))" << endl
	    << sorted(*und(*fftt(*gg))) << endl << endl;

	auto uu1 = fsys(*gg);
	auto fr = fffr(*uu1, *ff);
	cout << "fr = fffr(uu1,ff)" << endl;
	for (std::size_t i = 0; i < fr->layers.size(); i++)
	{
	    cout << "layer " << i << " :";
	    for (std::size_t j = 0; j < fr->layers[i].size(); j++)
		cout << " " << *fr->layers[i][j]->derived;
	    cout << endl;
	}
	cout << "frff(uu1, fr)" << endl
	    << *frff(*uu1, *fr) << endl << endl;

	fr = fffr(*uu1, *gg);
	cout << "fr = fffr(uu1,gg)" << endl;
	for (std::size_t i = 0; i < fr->layers.size(); i++)
	{
	    cout << "layer " << i << " :";
	    for (std::size_t j = 0; j < fr->layers[i].size(); j++)
		cout << " " << *fr->layers[i][j]->derived;
	    cout << endl;
	}
	cout << "frff(uu1, fr)" << endl
	    << *frff(*uu1, *fr) << endl << endl;

	auto hr = hhhr(*uu, *hh);
	cout << "hr = hhhr(uu,hh)" << endl;
	cout << "hrhh(uu,hr)" << endl
	    << *hrhh(*uu, *hr) << endl << endl;

	auto uu2 = uunion(*uu, *uu1);
	auto hr1 = frmul(*hr, *fr);
	cout << "hr1 = frmul(hr,fr)" << endl;
	cout << "hrhh(uu2,hr1)" << endl
	    << *hrhh(*uu2,*hr1) << endl << endl;
	cout << "rpln(aall(hraa(uu2,hr1)))" << endl;
	rpln(cout, sorted(*aall(*hraa(*uu2,*hr1)))); cout << endl;
    }

    if (false)
    {
	auto lluu = listsSystem_u;
	auto cart = systemsSetVarsSetStateCartesian_u;
	auto llss = listsState;
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
	auto trans = [](std::unique_ptr<Histogram>& xx, const VarUSet& ww)
	{
	    return std::make_shared<Transform>(std::move(xx), ww);
	};
	auto ttaa = transformsHistogram;
	auto und = transformsUnderlying;
	auto der = transformsDerived;
	auto tmul = [](const Histogram& aa, const Transform& tt)
	{
	    return transformsHistogramsApply(tt, aa);
	};
	auto lltt = [llss, trans](const VarList& kk, const VarList& ww, const ValListList& qq)
	{
	    VarList vv(kk.begin(), kk.end());
	    vv.insert(vv.end(), ww.begin(), ww.end());
	    std::vector<StateRationalPair> ii;
	    for (auto& ll : qq)
	    {
		auto jj = std::vector<VarValPair>();
		for (int j = 0; j < ll.size(); j++)
		    jj.push_back(VarValPair(vv[j], ll[j]));
		ii.push_back(StateRationalPair(*llss(jj), 1));
	    }
	    return trans(std::make_unique<Histogram>(ii), VarUSet(ww.begin(), ww.end()));
	};
	auto llff = setTransformsFud_u;
	auto fhis = fudsSetHistogram;
	auto fvars = fudsSetVar;
	auto fder = fudsDerived;
	auto fund = fudsUnderlying;
	auto fftt = fudsTransform;
	auto fsys = fudsSystemImplied;
	auto dep = fudsSetVarsDepends_u;
	typedef std::pair<VarValPairList, TransformPtrList> VarValPairListTransformPtrListPair;
	typedef std::vector<VarValPairListTransformPtrListPair> VarValPairListTransformPtrListPairList;
	typedef std::vector<VarValPairListTransformPtrListPairList> VarValPairListTransformPtrListPairListList;
	auto lldf = [llss, llff](const VarValPairListTransformPtrListPairListList& zz)
	{
	    auto jj = std::vector<std::vector<StatePtrFudPtrPair>>();
	    jj.reserve(zz.size());
	    for (auto& ll : zz)
	    {
		auto kk = std::vector<StatePtrFudPtrPair>();
		kk.reserve(ll.size());
		for (auto& pp : ll)
		{
		    StatePtrFudPtrPair qq(std::move(llss(pp.first)), std::move(llff(pp.second)));
		    kk.push_back(qq);
		}
		jj.push_back(kk);
	    }
	    auto tt = std::make_unique<DecompFud>();
	    tt->tree_u() = *pathsTree(jj);
	    return tt;
	};
	auto dfll = [](const DecompFud& df)
	{
	    return treesPaths(df.tree_u());
	};
	auto dfund = decompFudsUnderlying;
	auto dfff = decompFudsFud;
	auto dfdr = systemsDecompFudsDecompFudRepa_u;
	auto drdf = systemsDecompFudRepasDecompFud_u;

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

	auto colour = Variable("colour");
	auto red = Value("red");
	auto black = Value("black");

	auto ttc = lltt(VarList{ suit }, VarList{ colour }, ValListList{
	    ValList{ hearts, red },
	    ValList{ clubs, black },
	    ValList{ diamonds, red },
	    ValList{ spades, black } });

	auto pip_or_face = Variable("pip_or_face");
	auto pip = Value("pip");
	auto face = Value("face");

	auto xxt = ValListList{
	    ValList{ ace, pip },
	    ValList{ king, face },
	    ValList{ queen, face },
	    ValList{ jack, face } };
	for (int i = 2; i <= 10; i++)
	    xxt.push_back(ValList{ Value(i),pip });

	auto ttt = lltt(VarList{ rank }, VarList{ pip_or_face }, xxt);

	auto odd_pip = Variable("odd_pip");
	auto yes = Value("yes");
	auto no = Value("no");

	auto xxop = ValListList{
	    ValList{ ace, yes },
	    ValList{ king, no },
	    ValList{ queen, no },
	    ValList{ jack, no } };
	for (int i = 2; i <= 10; i += 2)
	    xxop.push_back(ValList{ Value(i),no });
	for (int i = 3; i <= 9; i += 2)
	    xxop.push_back(ValList{ Value(i),yes });

	auto ttop = lltt(VarList{ rank }, VarList{ odd_pip }, xxop);

	auto df = lldf(VarValPairListTransformPtrListPairListList{
	    VarValPairListTransformPtrListPairList{
	    VarValPairListTransformPtrListPair(VarValPairList{},TransformPtrList{ ttop }),
	    VarValPairListTransformPtrListPair(VarValPairList{ VarValPair(odd_pip, no) },TransformPtrList{ ttc, ttt }) },
	    VarValPairListTransformPtrListPairList{
	    VarValPairListTransformPtrListPair(VarValPairList{},TransformPtrList{ ttop }),
	    VarValPairListTransformPtrListPair(VarValPairList{ VarValPair(odd_pip, yes) },TransformPtrList{ ttc }) } });

	cout << "df" << endl
	    << *df << endl << endl;

	cout << "treesSize(dfzz(df))" << endl
	    << treesSize(df->tree_u()) << endl << endl;

	cout << "dfund(df)" << endl
	    << sorted(*dfund(*df)) << endl << endl;

	cout << "rpln(dfll(df))" << endl;
	rpln(cout, *dfll(*df)); cout << endl;

	auto uu1 = fsys(*dfff(*df));
	auto dr = dfdr(*uu1, *df);
	cout << "dr = dfdr(uu1,df)" << endl;
	cout << "drdf(uu1,dr)" << endl
	    << *drdf(*uu1,*dr) << endl << endl;
    }

    if (true)
    {
	auto fsys = fudsSystemImplied;
	auto dfund = decompFudsUnderlying;
	auto dfff = decompFudsFud;
	auto fvars = fudsSetVar;
	auto dfdr = systemsDecompFudsDecompFudRepa_u;
	auto drdf = systemsDecompFudRepasDecompFud_u;

	ifstream istrm("C:/zzz/caiks/NISTPy-master/NIST_model2.json");

	auto start = chrono::system_clock::now();
	auto df = persistentsDecompFud(istrm);
	auto end = chrono::system_clock::now();
	cout << "df = persistentsDecompFud(istrm) " << ((chrono::duration<double>)(end-start)).count() << "s" << endl;

	cout << "len(fvars(dfff(df)))" << endl
	    << fvars(*dfff(*df))->size() << endl << endl;

	auto uu = fsys(*dfff(*df));
	cout << "len(uu)" << endl
	    << uu->map_u().size() << endl << endl;

	start = chrono::system_clock::now();
	auto dr = dfdr(*uu,*df);
	end = chrono::system_clock::now();
	cout << "dr = dfdr(*uu,*df) " << ((chrono::duration<double>)(end-start)).count() << "s" << endl;

	start = chrono::system_clock::now();
	auto df1 = drdf(*uu,*dr);
	end = chrono::system_clock::now();
	cout << "df1 = drdf(*uu,*dr) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

	cout << "len(fvars(dfff(df1)))" << endl
	    << fvars(*dfff(*df1))->size() << endl << endl;


    }


}
