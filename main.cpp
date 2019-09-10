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

    if (true)
    {
	auto regcart = histogramRegularCartesian_u;
	auto regsing = histogramRegularUnitSingleton_u;
	auto regdiag = histogramRegularUnitDiagonal_u;
	auto sys = histogramsSystemImplied;
	auto aarr = systemsHistogramsHistogramRepa_u;
	auto rraa = systemsHistogramRepasHistogram_u;

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

    }

}
