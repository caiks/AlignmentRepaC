#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentApprox.h"
#include "AlignmentAeson.h"
#include "AlignmentRepa.h"
#include "AlignmentAesonRepa.h"
#include "AlignmentRandomRepa.h"
#include "AlignmentPracticableRepa.h"
#include "AlignmentPracticableIORepa.h"
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

#define ECHO(x) cout << #x << endl; x
#define EVAL(x) cout << #x << ": " << (x) << endl
#define EVALL(x) cout << #x << ": " << endl << (x) << endl
#define TRUTH(x) cout << #x << ": " << ((x) ? "true" : "false") << endl


int main(int argc, char **argv)
{
	if (false)
	{
		auto regcart = histogramRegularCartesian_u;
		auto regsing = histogramRegularUnitSingleton_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
		auto sys = histogramsSystemImplied;
		auto uuur = systemsSystemRepa;
		auto uruu = systemsRepasSystem;

		auto aa = regdiag(2, 2);
		cout << "aa = regdiag(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		auto uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;
		cout << "ur.mapVarSize()" << endl
			<< ur->mapVarSize() << endl << endl;

		ur->listVarSizePair.push_back(VarSizePair(std::make_shared<Variable>(Variable(98)), 3));
		ur->listVarSizePair.push_back(VarSizePair(std::make_shared<Variable>(Variable(99)), 5));

		cout << "ur" << endl
			<< *ur << endl << endl;
		cout << "ur.mapVarSize()" << endl
			<< ur->mapVarSize() << endl << endl;

		uruu(*ur,*uu);
		cout << "uruu(ur,uu)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		aa = regcart(257, 1);
		cout << "aa = regcart(257,1)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		try
		{
			ur = uuur(*uu);
		}
		catch (const std::out_of_range& e)
		{
			cout << "caught out_of_range: " << e.what() << endl << endl;
		}
	}

	if (false)
	{
		auto regcart = histogramRegularCartesian_u;
		auto regsing = histogramRegularUnitSingleton_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
		auto sys = histogramsSystemImplied;
		auto uuur = systemsSystemRepa;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto arred = [](const HistogramRepa& ar, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistogramRepasReduce_u(m, kk1.data(), ar);
		};

		auto aa = regdiag(2, 2);
		cout << "aa = regdiag(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		auto uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		auto ar = aarr(*uu,*ur,*aa);
		cout << "ar = aarr(uu,ur,aa)" << endl;
		cout << "ar" << endl
			<< *ar << endl << endl;

		auto aa1 = rraa(*uu,*ur,*ar);
		cout << "aa1 = rraa(uu,ur,ar)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		auto br = arred(*ar, *ur, VarList{ Variable(2),Variable(1) });
		cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		auto bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{ Variable(2) });
		cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{});
		cout << "br = arred(ar, VarList{})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		aa = regcart(2, 2);
		cout << "aa = regcart(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		ar = aarr(*uu, *ur, *aa);
		cout << "ar = aarr(uu,ur,aa)" << endl;
		cout << "ar" << endl
			<< *ar << endl << endl;

		aa1 = rraa(*uu, *ur, *ar);
		cout << "aa1 = rraa(uu,ur,ar)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		br = arred(*ar, *ur, VarList{ Variable(2),Variable(1) });
		cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{ Variable(2) });
		cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{});
		cout << "br = arred(ar, VarList{})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		aa = regsing(2, 2);
		cout << "aa = regsing(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		ar = aarr(*uu, *ur, *aa);
		cout << "ar = aarr(uu,ur,aa)" << endl;
		cout << "ar" << endl
			<< *ar << endl << endl;

		aa1 = rraa(*uu, *ur, *ar);
		cout << "aa1 = rraa(uu,ur,ar)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		br = arred(*ar, *ur, VarList{ Variable(2),Variable(1) });
		cout << "br = arred(ar, VarList{ Variable(2),Variable(1) })" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{ Variable(2) });
		cout << "br = arred(ar, VarList{ Variable(2)})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;

		br = arred(*ar, *ur, VarList{});
		cout << "br = arred(ar, VarList{})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "bb1" << endl
			<< *bb1 << endl << endl;
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
		auto uuur = systemsSystemRepa;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto arred = [](const HistogramRepa& ar, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistogramRepasReduce_u(m, kk1.data(), ar);
		};
		auto arpr = histogramRepasRed;
		auto prar = histogramRepaRedsIndependent;
		auto arent = histogramRepasEntropy;

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

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		auto ar = aarr(*uu, *ur, *aa);
		cout << "ar = aarr(uu,ur,aa)" << endl;
		cout << "ar" << endl
			<< *ar << endl << endl;
		cout << "20 * arent(ar)" << endl
			<< 20 * arent(*ar) << endl << endl;

		auto aa1 = rraa(*uu, *ur, *ar);
		cout << "aa1 = rraa(uu,ur,ar)" << endl;
		cout << "rpln(aall(trim(aa1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*aa1)))); cout << endl;

		auto pr = arpr(20, *ar);
		cout << "pr = arpr(20,ar)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		auto xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;

		auto br = arred(*ar, *ur, VarList{ pressure, rain, cloud, wind });
		cout << "br = arred(ar, VarList{ pressure, rain, cloud, wind })" << endl;
		cout << "br" << endl
			<< *br << endl << endl;
		cout << "20 * arent(br)" << endl
			<< 20 * arent(*br) << endl << endl;

		auto bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "rpln(aall(trim(bb1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*bb1)))); cout << endl;

		pr = arpr(20, *br);
		cout << "pr = arpr(20,br)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;

		br = arred(*ar, *ur, VarList{ wind, cloud, rain });
		cout << "br = arred(ar, VarList{ wind, cloud, rain})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;
		cout << "20 * arent(br)" << endl
			<< 20 * arent(*br) << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "rpln(aall(trim(bb1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*bb1)))); cout << endl;

		pr = arpr(20, *br);
		cout << "pr = arpr(20,br)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;

		br = arred(*ar, *ur, VarList{ wind, rain });
		cout << "br = arred(ar, VarList{ wind, rain})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;
		cout << "20 * arent(br)" << endl
			<< 20 * arent(*br) << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "rpln(aall(trim(bb1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*bb1)))); cout << endl;

		pr = arpr(20, *br);
		cout << "pr = arpr(20,br)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;

		br = arred(*ar, *ur, VarList{ rain });
		cout << "br = arred(ar, VarList{ rain})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;
		cout << "20 * arent(br)" << endl
			<< 20 * arent(*br) << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "rpln(aall(trim(bb1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*bb1)))); cout << endl;

		pr = arpr(20, *br);
		cout << "pr = arpr(20,br)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;

		br = arred(*ar, *ur, VarList{ });
		cout << "br = arred(ar, VarList{})" << endl;
		cout << "br" << endl
			<< *br << endl << endl;
		cout << "20 * arent(br)" << endl
			<< 20 * arent(*br) << endl << endl;

		bb1 = rraa(*uu, *ur, *br);
		cout << "bb1 = rraa(uu,ur,br)" << endl;
		cout << "rpln(aall(trim(bb1)))" << endl;
		rpln(cout, sorted(*aall(*trim(*bb1)))); cout << endl;

		pr = arpr(20, *br);
		cout << "pr = arpr(20,br)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		xr = prar(20, *pr);
		cout << "xr = prar(20,pr)" << endl;
		cout << "xr" << endl
			<< *xr << endl << endl;
		cout << "20 * arent(xr)" << endl
			<< 20 * arent(*xr) << endl << endl;
	}

//	if (false)
//	{
//		auto suit = Variable("suit");
//		auto rank = Variable("rank");
//
//		auto vv = VarList{ suit,rank };
//
//		HistogramRepaVec hr;
//		hr.vectorVar = vv;
//
//		cout << "hr.vectorVar" << endl
//			<< hr.vectorVar << endl << endl;
//
//		cout << "hr.mapVarInt()" << endl
//			<< hr.mapVarInt() << endl << endl;
//	}

	if (false)
	{
		auto regcart = histogramRegularCartesian_u;
		auto regsing = histogramRegularUnitSingleton_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
		auto sys = histogramsSystemImplied;
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu,ur,hr));
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
		{
			return eventsHistoryRepasHistoryRepaSelection_u(ll.size(), (std::size_t*)ll.data(), hr);
		};		
		auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
		{
			return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
		};
		auto hrhrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasHistoryRepaReduced_u(m, kk1.data(), hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};

		auto aa = regdiag(2, 2);
		cout << "aa = regdiag(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		auto uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		auto hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		auto aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regcart(2, 2);
		cout << "aa = regcart(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regsing(2, 2);
		cout << "aa = regsing(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regcart(3, 2);
		cout << "aa = regcart(3, 2)" << endl;

		uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		auto hh = hrhh(*uu, *ur, *hr);
		cout << "hh = hraa(uu,hr)" << endl;
		cout << "hh" << endl
			<< *hh << endl << endl;

		auto hr1 = hrsel(*hr, SizeList{});
		cout << "hr1 = hrsel(hr, SizeList{})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		auto hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		hr1 = hrsel(*hr, SizeList{ 0 });
		cout << "hr1 = hrsel(hr, SizeList{0})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		hr1 = hrsel(*hr, SizeList{ 0,4,8 });
		cout << "hr1 = hrsel(hr, SizeList{0,4,8})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		auto bb = regcart(1,1);
		cout << "bb = regcart(1,1)" << endl;
		auto ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu,*ur,*hrhrsel(*hr,*ss)) << endl << endl;
		bb = regcart(2, 1);
		cout << "bb = regcart(2,1)" << endl;
		ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu, *ur, *hrhrsel(*hr, *ss)) << endl << endl;
		bb = regcart(3, 1);
		cout << "bb = regcart(3,1)" << endl;
		ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu, *ur, *hrhrsel(*hr, *ss)) << endl << endl;

		cout << "hraa(uu,hrhrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{ Variable(2),Variable(1) })) << endl << endl;
		cout << "hraa(uu,hrhrred(hr, VarList{ Variable(1) }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{ Variable(1) })) << endl << endl;
		cout << "hraa(uu,hrhrred(hr, VarList{ }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{ })) << endl << endl;

		cout << "araa(uu,hrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{ Variable(2),Variable(1) })) << endl << endl;
		cout << "araa(uu,hrred(hr, VarList{ Variable(1) }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{ Variable(1) })) << endl << endl;
		cout << "araa(uu,hrred(hr, VarList{ }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{ })) << endl << endl;
	}

	if (false)
	{
		auto regcart = histogramRegularCartesian_u;
		auto regsing = histogramRegularUnitSingleton_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
		auto sys = histogramsSystemImplied;
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa));
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
		{
			return eventsHistoryRepasHistoryRepaSelection_u(ll.size(), (std::size_t*)ll.data(), hr);
		};
		auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
		{
			return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
		};
		auto hrhrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasHistoryRepaReduced_u(m, kk1.data(), hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};

		auto aa = regdiag(2, 2);
		cout << "aa = regdiag(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		auto uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		auto hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		auto aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regcart(2, 2);
		cout << "aa = regcart(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regsing(2, 2);
		cout << "aa = regsing(2, 2)" << endl;
		cout << "aa" << endl
			<< *aa << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		aa1 = hraa(*uu, *ur, *hr);
		cout << "aa1 = hraa(uu,ur,hr)" << endl;
		cout << "aa1" << endl
			<< *aa1 << endl << endl;

		aa = regcart(3, 2);
		cout << "aa = regcart(3, 2)" << endl;

		uu = sys(*aa);
		cout << "uu = sys(aa)" << endl;
		cout << "uu" << endl
			<< *uu << endl << endl;

		ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;
		cout << "ur" << endl
			<< *ur << endl << endl;

		hr = aahr(*uu, *ur, *aa);
		cout << "hr = aarr(uu,ur,aa)" << endl;
		cout << "hr" << endl
			<< *hr << endl << endl;

		auto hh = hrhh(*uu, *ur, *hr);
		cout << "hh = hraa(uu,hr)" << endl;
		cout << "hh" << endl
			<< *hh << endl << endl;

		auto hr1 = hrsel(*hr, SizeList{});
		cout << "hr1 = hrsel(hr, SizeList{})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		auto hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		hr1 = hrsel(*hr, SizeList{ 0 });
		cout << "hr1 = hrsel(hr, SizeList{0})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		hr1 = hrsel(*hr, SizeList{ 0,4,8 });
		cout << "hr1 = hrsel(hr, SizeList{0,4,8})" << endl;
		cout << "hr1" << endl
			<< *hr1 << endl << endl;

		hh1 = hrhh(*uu, *ur, *hr1);
		cout << "hh1 = hraa(uu,hr1)" << endl;
		cout << "hh1" << endl
			<< *hh1 << endl << endl;

		auto bb = regcart(1, 1);
		cout << "bb = regcart(1,1)" << endl;
		auto ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu, *ur, *hrhrsel(*hr, *ss)) << endl << endl;
		bb = regcart(2, 1);
		cout << "bb = regcart(2,1)" << endl;
		ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu, *ur, *hrhrsel(*hr, *ss)) << endl << endl;
		bb = regcart(3, 1);
		cout << "bb = regcart(3,1)" << endl;
		ss = aahr(*uu, *ur, *bb);
		cout << "ss = aahr(uu,bb)" << endl;
		cout << "hrhh(uu,hrhrsel(hr,ss))" << endl
			<< *hrhh(*uu, *ur, *hrhrsel(*hr, *ss)) << endl << endl;

		cout << "hraa(uu,hrhrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{ Variable(2),Variable(1) })) << endl << endl;
		cout << "hraa(uu,hrhrred(hr, VarList{ Variable(1) }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{ Variable(1) })) << endl << endl;
		cout << "hraa(uu,hrhrred(hr, VarList{ }))" << endl
			<< *hraa(*uu, *ur, *hrhrred(*hr, *ur, VarList{})) << endl << endl;

		cout << "araa(uu,hrred(hr, VarList{ Variable(2),Variable(1) }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{ Variable(2),Variable(1) })) << endl << endl;
		cout << "araa(uu,hrred(hr, VarList{ Variable(1) }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{ Variable(1) })) << endl << endl;
		cout << "araa(uu,hrred(hr, VarList{ }))" << endl
			<< *araa(*uu, *ur, *hrred(*hr, *ur, VarList{})) << endl << endl;
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
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
		{
			return eventsHistoryRepasHistoryRepaSelection_u(ll.size(), (std::size_t*)ll.data(), hr);
		};
		auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
		{
			return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
		};
		auto hrhrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasHistoryRepaReduced_u(m, kk1.data(), hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto prents = histogramRepaRedsListEntropy;
		auto hrpr1 = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasRed_u(m, kk1.data(), hr);
		};
		auto hrshuffle = historyRepasShuffle_u;
		auto cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u;

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

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;

		auto hr = hhhr(*uu, *ur, *hh);
		cout << "hr = hhhr(uu,hh)" << endl;
		cout << "hraa(uu,hr)" << endl
			<< *hraa(*uu, *ur, *hr) << endl << endl;

		cout << "hrhh(uu,hrsel(hr,SizeList{}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{})) << endl << endl;
		cout << "hrhh(uu,hrsel(hr,SizeList{0}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{ 0 })) << endl << endl;
		cout << "hrhh(uu,hrsel(hr,SizeList{0,4,8}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{ 0,4,8 })) << endl << endl;

		auto pr = hrpr(*hr);
		cout << "pr = hrpr(hr)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;
		auto pr1 = hrpr1(*hr, *ur, VarList{ rain, pressure, cloud, wind });
		cout << "pr1 = hrpr1(*hr, *ur, VarList{ rain, pressure, cloud, wind })" << endl
			<< *pr1 << endl << endl;
		cout << "prents(pr1)" << endl
			<< *prents(*pr1) << endl << endl;
		pr1 = hrpr1(*hr, *ur, VarList{ pressure, cloud, wind });
		cout << "pr1 = hrpr1(*hr, *ur, VarList{ pressure, cloud, wind })" << endl
			<< *pr1 << endl << endl;
		cout << "prents(pr1)" << endl
			<< *prents(*pr1) << endl << endl;
		pr1 = hrpr1(*hr, *ur, VarList{ cloud, wind });
		cout << "pr1 = hrpr1(*hr, *ur, VarList{ cloud, wind })" << endl
			<< *pr1 << endl << endl;
		cout << "prents(pr1)" << endl
			<< *prents(*pr1) << endl << endl;
		pr1 = hrpr1(*hr, *ur, VarList{ wind });
		cout << "pr1 = hrpr1(*hr, *ur, VarList{ wind })" << endl
			<< *pr1 << endl << endl;
		cout << "prents(pr1)" << endl
			<< *prents(*pr1) << endl << endl;
		pr1 = hrpr1(*hr, *ur, VarList{ });
		cout << "pr1 = hrpr1(*hr, *ur, VarList{ })" << endl
			<< *pr1 << endl << endl;
		cout << "prents(pr1)" << endl
			<< *prents(*pr1) << endl << endl;

		auto hr1 = hrhrred(*hr, *ur, VarList{ pressure, rain, cloud, wind });
		cout << "hr1 = hrhrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ wind, cloud, rain });
		cout << "hr1 = hrhrred(hr, VarList{ wind, cloud, rain })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ rain, wind });
		cout << "hr1 = hrhrred(hr, VarList{ rain, wind })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ rain });
		cout << "hr1 = hrhrred(hr, VarList{ rain })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ });
		cout << "hr1 = hrhrred(hr, VarList{ })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		auto br = hrred(*hr, *ur, VarList{ pressure, rain, cloud, wind });
		cout << "br = hrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ wind, cloud, rain });
		cout << "br = hrred(hr, VarList{ wind, cloud, rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ rain, wind });
		cout << "br = hrred(hr, VarList{ rain, wind})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ rain });
		cout << "br = hrred(hr, VarList{ rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{});
		cout << "br = hrred(hr, VarList{})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		auto hrs = hrshuffle(*hr,7);
		cout << "hrs = hrshuffle(7,hr)" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hrs))))); cout << endl;

		br = hrred(*hrs, *ur, VarList{ rain });
		cout << "br = hrred(hrs, VarList{ rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		br = hrred(*hrs, *ur, VarList{});
		cout << "br = hrred(hrs, VarList{})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		{
			VarList ww{ pressure, rain, cloud, wind };
			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);
			pr = hrpr(*hr);
			auto prs = hrpr(*hrs);
			std::size_t xmax = 81;
			std::size_t omax = 5;
			auto xx = cross(xmax, omax, ww1, *hr, *pr, *hrs, *prs);
			auto s = std::get<1>(xx);
			cout << "s" << endl
				<< s << endl << endl;
			auto tt = sorted(*std::get<0>(xx));
			for (auto& mm : tt)
				cout << mm.first << "," << mm.second << endl;
		}
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
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa));
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
		{
			return eventsHistoryRepasHistoryRepaSelection_u(ll.size(), (std::size_t*)ll.data(), hr);
		};
		auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
		{
			return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
		};
		auto hrhrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasHistoryRepaReduced_u(m, kk1.data(), hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto hrconcat = vectorHistoryRepasConcat_u;
		auto cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u;

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

		auto ur = uuur(*uu);
		cout << "ur = uuur(*uu)" << endl;

		auto hr = hhhr(*uu, *ur, *hh);
		cout << "hr = hhhr(uu,hh)" << endl;
		cout << "hraa(uu,hr)" << endl
			<< *hraa(*uu, *ur, *hr) << endl << endl;

		cout << "hrhh(uu,hrsel(hr,SizeList{}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{})) << endl << endl;
		cout << "hrhh(uu,hrsel(hr,SizeList{0}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{ 0 })) << endl << endl;
		cout << "hrhh(uu,hrsel(hr,SizeList{0,4,8}))" << endl
			<< *hrhh(*uu, *ur, *hrsel(*hr, SizeList{ 0,4,8 })) << endl << endl;

		auto pr = hrpr(*hr);
		cout << "pr = hrpr(hr)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		auto hr1 = hrhrred(*hr, *ur, VarList{ pressure, rain, cloud, wind });
		cout << "hr1 = hrhrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ wind, cloud, rain });
		cout << "hr1 = hrhrred(hr, VarList{ wind, cloud, rain })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ rain, wind });
		cout << "hr1 = hrhrred(hr, VarList{ rain, wind })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{ rain });
		cout << "hr1 = hrhrred(hr, VarList{ rain })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		hr1 = hrhrred(*hr, *ur, VarList{});
		cout << "hr1 = hrhrred(hr, VarList{ })" << endl;
		cout << "rpln(aall(trim(hraa(uu,hr1))))" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hr1))))); cout << endl;
		pr = hrpr(*hr1);
		cout << "pr = hrpr(hr1)" << endl;
		cout << "pr" << endl
			<< *pr << endl << endl;

		auto br = hrred(*hr, *ur, VarList{ pressure, rain, cloud, wind });
		cout << "br = hrred(hr, VarList{ pressure, rain, cloud, wind })" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ wind, cloud, rain });
		cout << "br = hrred(hr, VarList{ wind, cloud, rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ rain, wind });
		cout << "br = hrred(hr, VarList{ rain, wind})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{ rain });
		cout << "br = hrred(hr, VarList{ rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;
		br = hrred(*hr, *ur, VarList{});
		cout << "br = hrred(hr, VarList{})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		auto hrs = hrshuffle(*hr, 7);
		cout << "hrs = hrshuffle(7,hr)" << endl;
		rpln(cout, sorted(*aall(*trim(*hraa(*uu, *ur, *hrs))))); cout << endl;

		br = hrred(*hrs, *ur, VarList{ rain });
		cout << "br = hrred(hrs, VarList{ rain})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		br = hrred(*hrs, *ur, VarList{});
		cout << "br = hrred(hrs, VarList{})" << endl;
		cout << "rpln(aall(trim(rraa(uu,br))))" << endl;
		rpln(cout, sorted(*aall(*trim(*rraa(*uu, *ur, *br))))); cout << endl;

		{
			VarList ww{ pressure, rain, cloud, wind };
			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);
			pr = hrpr(*hr);
			auto prs = hrpr(*hrs);
			std::size_t xmax = 81;
			std::size_t omax = 5;
			auto xx = cross(xmax, omax, ww1, *hr, *pr, *hrs, *prs);
			auto s = std::get<1>(xx);
			cout << "s" << endl
				<< s << endl << endl;
			auto tt = sorted(*std::get<0>(xx));
			for (auto& mm : tt)
				cout << mm.first << "," << mm.second << endl;
		}

		HistoryRepaPtr hr0 = std::move(hr);
		auto hr2 = hrconcat(HistoryRepaPtrList{ hr0,hr0,hr0 });
		cout << "hr2 = hrconcat(hr,hr,hr)" << endl;
		cout << "hraa(uu,hr2)" << endl
			<< *hraa(*uu, *ur, *hr2) << endl << endl;
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
		auto uuur = systemsSystemRepa;
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
		cout << "uu1 = sys(tt->histogram_u())" << endl;
		cout << "uu1" << endl
			<< *uu1 << endl << endl;

		auto ur1 = uuur(*uu1);
		cout << "ur1 = uuur(*uu1)" << endl;
		cout << "ur1" << endl
			<< *ur1 << endl << endl;

		auto tr = tttr(*uu1, *ur1, *tt);
		cout << "tr = tttr(uu1,tt)" << endl;
		cout << "tr" << endl
			<< *tr << endl << endl;

		cout << "trtt(uu1,tr)" << endl
			<< *trtt(*uu1, *ur1, *tr) << endl << endl;

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
			return std::make_shared<Transform>(xx, ww);
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
			auto aa = std::make_unique<Histogram>(ii);
			return trans(aa, VarUSet(ww.begin(), ww.end()));
		};
		auto llff = setTransformsFud_u;
		auto fhis = fudsSetHistogram;
		auto fvars = fudsSetVar;
		auto fder = fudsDerived;
		auto fund = fudsUnderlying;
		auto fftt = fudsTransform;
		auto fsys = fudsSystemImplied;
		auto dep = fudsSetVarsDepends_u;
		auto uuur = systemsSystemRepa;
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh);
		};
		auto hrhh = systemsHistoryRepasHistory_u;
		auto llfr = setVariablesListTransformRepasFudRepa_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;
		auto frmul = historyRepasFudRepasMultiply_u;
		auto frvars = fudRepasSetVar;
		auto frder = fudRepasDerived;
		auto frund = fudRepasUnderlying;
		auto frdep = fudRepasSetVarsDepends;

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
		auto uu2 = uunion(*uu, *uu1);
		auto ur2 = uuur(*uu2);
		auto fr = fffr(*uu1, *ur2, *ff);
		cout << "fr = fffr(uu1,ff)" << endl;
		for (std::size_t i = 0; i < fr->layers.size(); i++)
		{
			cout << "layer " << i << " :";
			for (std::size_t j = 0; j < fr->layers[i].size(); j++)
				cout << " " << fr->layers[i][j]->derived;
			cout << endl;
		}
		cout << "fr" << endl
			<< *fr << endl << endl;
		cout << "frder(fr)" << endl
			<< sorted(*frder(*fr)) << endl << endl;
		cout << "frund(fr)" << endl
			<< sorted(*frund(*fr)) << endl << endl;
		cout << "frvars(fr)" << endl
			<< sorted(*frvars(*fr)) << endl << endl;
		cout << "frff(uu1, fr)" << endl
			<< *frff(*uu1, *ur2, *fr) << endl << endl;
		SizeUSet vv1;
		for (auto& v : *vv)
			vv1.insert(ur2->mapVarSize()[v]);
		SizeUSet kk;
		cout << "frdep(fr,{})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.insert(1);
		cout << "frdep(fr,{1})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.insert(2);
		cout << "frdep(fr,{1,2})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;

		fr = fffr(*uu1, *ur2, *gg);
		cout << "fr = fffr(uu1,gg)" << endl;
		for (std::size_t i = 0; i < fr->layers.size(); i++)
		{
			cout << "layer " << i << " :";
			for (std::size_t j = 0; j < fr->layers[i].size(); j++)
				cout << " " << fr->layers[i][j]->derived;
			cout << endl;
		}
		cout << "fr" << endl
			<< *fr << endl << endl;
		cout << "frder(fr)" << endl
			<< sorted(*frder(*fr)) << endl << endl;
		cout << "frund(fr)" << endl
			<< sorted(*frund(*fr)) << endl << endl;
		cout << "frvars(fr)" << endl
			<< sorted(*frvars(*fr)) << endl << endl;
		cout << "frff(uu1, fr)" << endl
			<< *frff(*uu1, *ur2, *fr) << endl << endl;
		kk.clear();
		cout << "frdep(fr,{})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.insert(1);
		cout << "frdep(fr,{1})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.insert(2);
		cout << "frdep(fr,{1,2})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.insert(3);
		cout << "frdep(fr,{1,2,3})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;
		kk.clear();
		kk.insert(3);
		cout << "frdep(fr,{3})" << endl
			<< *llfr(vv1, *frdep(*fr, kk)) << endl << endl;

		auto hr = hhhr(*uu, *ur2, *hh);
		cout << "hr = hhhr(uu,hh)" << endl;
		cout << "hrhh(uu,hr)" << endl
			<< *hrhh(*uu, *ur2, *hr) << endl << endl;

		auto hr1 = frmul(*hr, *fr);
		cout << "hr1 = frmul(hr,fr)" << endl;
		cout << "hrhh(uu2,hr1)" << endl
			<< *hrhh(*uu2, *ur2, *hr1) << endl << endl;
		cout << "rpln(aall(hraa(uu2,hr1)))" << endl;
		rpln(cout, sorted(*aall(*hraa(*uu2, *ur2, *hr1)))); cout << endl;

		hr->transpose();
		cout << "hr.transpose()" << endl;
		hr1 = frmul(*hr, *fr);
		cout << "hr1 = frmul(hr,fr)" << endl;
		cout << "hrhh(uu2,hr1)" << endl
			<< *hrhh(*uu2, *ur2, *hr1) << endl << endl;
		cout << "rpln(aall(hraa(uu2,hr1)))" << endl;
		rpln(cout, sorted(*aall(*hraa(*uu2, *ur2, *hr1)))); cout << endl;
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
			return std::make_shared<Transform>(xx, ww);
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
			auto aa = std::make_unique<Histogram>(ii);
			return trans(aa, VarUSet(ww.begin(), ww.end()));
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
		auto uuur = systemsSystemRepa;
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
		auto ur1 = uuur(*uu1);
		auto dr = dfdr(*uu1, *ur1, *df);
		cout << "dr = dfdr(uu1,df)" << endl;
		cout << "drdf(uu1,dr)" << endl
			<< *drdf(*uu1, *ur1, *dr) << endl << endl;
	}

	if (false)
	{
		auto fsys = fudsSystemImplied;
		auto dfund = decompFudsUnderlying;
		auto dfff = decompFudsFud;
		auto fvars = fudsSetVar;
		auto uuur = systemsSystemRepa;
		auto dfdr = systemsDecompFudsDecompFudRepa_u;
		auto drdf = systemsDecompFudRepasDecompFud_u;

		ifstream istrm("C:/zzz/caiks/NISTPy-master/NIST_model2.json");

		auto start = chrono::system_clock::now();
		auto df = persistentsDecompFud(istrm);
		auto end = chrono::system_clock::now();
		cout << "df = persistentsDecompFud(istrm) " << ((chrono::duration<double>)(end-start)).count() << "s" << endl;

		auto ff = dfff(*df);
		cout << "len(fvars(dfff(df)))" << endl
			<< fvars(*ff)->size() << endl << endl;

		auto uu = fsys(*ff);
		cout << "len(uu)" << endl
			<< uu->map_u().size() << endl << endl;

		start = chrono::system_clock::now();
		std::unique_ptr<DecompFudRepa> dr;
		try
		{
			auto ur = uuur(*uu);
			dr = dfdr(*uu, *ur, *df);
			end = chrono::system_clock::now();
			cout << "dr = dfdr(*uu,*df) " << ((chrono::duration<double>)(end-start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto df1 = drdf(*uu, *ur, *dr);
			end = chrono::system_clock::now();
			cout << "df1 = drdf(*uu,*dr) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			auto ff1 = dfff(*df1);
			cout << "len(fvars(dfff(df1)))" << endl
				<< fvars(*ff1)->size() << endl << endl;
		}
		catch (const std::exception& e)
		{
			cout << "caught exception: " << e.what() << endl << endl;
		}
	}

	if (false)
	{
		auto fsys = fudsSystemImplied;
		auto dfund = decompFudsUnderlying;
		auto dfff = decompFudsFud;
		auto fvars = fudsSetVar;
		auto uuur = systemsSystemRepa;
		auto dfdr = systemsDecompFudsDecompFudRepa_u;
		auto drdf = systemsDecompFudRepasDecompFud_u;

		ifstream istrm("C:/zzz/caiks/NISTPy-master/NIST_model24.json");

		auto start = chrono::system_clock::now();
		auto df = persistentsDecompFud(istrm);
		auto end = chrono::system_clock::now();
		cout << "df = persistentsDecompFud(istrm) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		auto ff = dfff(*df);
		cout << "len(fvars(dfff(df)))" << endl
			<< fvars(*ff)->size() << endl << endl;

		auto uu = fsys(*ff);
		cout << "len(uu)" << endl
			<< uu->map_u().size() << endl << endl;

		start = chrono::system_clock::now();
		std::unique_ptr<DecompFudRepa> dr;
		try
		{
			auto ur = uuur(*uu);
			dr = dfdr(*uu, *ur, *df);
			end = chrono::system_clock::now();
			cout << "dr = dfdr(*uu,*df) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto df1 = drdf(*uu, *ur, *dr);
			end = chrono::system_clock::now();
			cout << "df1 = drdf(*uu,*dr) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			auto ff1 = dfff(*df1);
			cout << "len(fvars(dfff(df1)))" << endl
				<< fvars(*ff1)->size() << endl << endl;
		}
		catch (const std::exception& e)
		{
			cout << "caught exception: " << e.what() << endl << endl;
		}
		/*
		df = persistentsDecompFud(istrm) 11.435s
		len(fvars(dfff(df)))
		5966

		len(uu)
		5966

		dr = dfdr(*uu,*df) 4.29008s
		df1 = drdf(*uu,*dr) 4.21208s
		len(fvars(dfff(df1)))
		5966
		*/
	}

	if (false)
	{
		auto uvars = systemsSetVar;
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
		auto trim = histogramsTrim;
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
			return std::make_shared<Transform>(xx, ww);
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
			auto aa = std::make_unique<Histogram>(ii);
			return trans(aa, VarUSet(ww.begin(), ww.end()));
		};
		auto llff = setTransformsFud_u;
		auto fhis = fudsSetHistogram;
		auto fvars = fudsSetVar;
		auto fder = fudsDerived;
		auto fund = fudsUnderlying;
		auto fftt = fudsTransform;
		auto fsys = fudsSystemImplied;
		auto dep = fudsSetVarsDepends_u;
		auto uuur = systemsSystemRepa;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa),1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = systemsHistoryRepasHistory_u;
		auto tttr = systemsTransformsTransformRepa_u;
		auto trtt = systemsTransformRepasTransform_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;

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

		StrVarPtrMap m;

		std::stringstream str;
		systemsPersistent(*uu, str);
		auto uu1 = persistentsSystem(str, m);
		auto ur1 = uuur(*uu1);
		{
			auto hr = hhhr(*uu1, *ur1, *hh);
			cout << "hr = hhhr(uu1,hh)" << endl;
			rpln(cout, sorted(*aall(*trim(*hraa(*uu1, *ur1, *hr))))); cout << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "historyRepasPersistent(hr,out)" << endl;
			historyRepasPersistent(*hr,out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto hr2 = persistentsHistoryRepa(in);
			cout << "hr2 = persistentsHistoryRepa(in)" << endl;
			rpln(cout, sorted(*aall(*trim(*hraa(*uu1, *ur1, *hr2))))); cout << endl;
			in.close();
		}

		auto uu2 = fsys(*gg);
		auto ur2 = uuur(*uu2);
		{
			auto tr = tttr(*uu2, *ur2, *ttcw);
			cout << "tr = tttr(uu2,ttcw)" << endl;
			cout << "trtt(uu2,tr)" << endl
				<< *trtt(*uu2, *ur2, *tr) << endl << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "transformRepasPersistent(tr,out)" << endl;
			transformRepasPersistent(*tr,out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto tr2 = persistentsTransformRepa(in);
			cout << "tr2 = persistentsTransformRepa(in)" << endl;
			cout << "trtt(uu2,tr2)" << endl
				<< *trtt(*uu2, *ur2, *tr2) << endl << endl;
			in.close();
		}

		{
			auto fr = fffr(*uu2, *ur2, *gg);
			cout << "fr = fffr(uu2,gg)" << endl;
			cout << "frff(uu2, fr)" << endl
				<< *frff(*uu2, *ur2, *fr) << endl << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "fudRepasPersistent(fr,out)" << endl;
			fudRepasPersistent(*fr, out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto fr2 = persistentsFudRepa(in);
			cout << "fr2 = persistentsFudRepa(in)" << endl;
			cout << "frff(uu2, fr2)" << endl
				<< *frff(*uu2, *ur2, *fr2) << endl << endl;
			in.close();
		}
	}

	if (false)
	{
		auto uvars = systemsSetVar;
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
		auto trim = histogramsTrim;
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
			return std::make_shared<Transform>(xx, ww);
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
			auto aa = std::make_unique<Histogram>(ii);
			return trans(aa, VarUSet(ww.begin(), ww.end()));
		};
		auto llff = setTransformsFud_u;
		auto fhis = fudsSetHistogram;
		auto fvars = fudsSetVar;
		auto fder = fudsDerived;
		auto fund = fudsUnderlying;
		auto fftt = fudsTransform;
		auto fsys = fudsSystemImplied;
		auto dep = fudsSetVarsDepends_u;
		auto uuur = systemsSystemRepa;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa));
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh);
		};
		auto hrhh = systemsHistoryRepasHistory_u;
		auto tttr = systemsTransformsTransformRepa_u;
		auto trtt = systemsTransformRepasTransform_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;

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

		StrVarPtrMap m;

		std::stringstream str;
		systemsPersistent(*uu, str);
		auto uu1 = persistentsSystem(str, m);
		auto ur1 = uuur(*uu1);
		{
			auto hr = hhhr(*uu1, *ur1, *hh);
			cout << "hr = hhhr(uu1,hh)" << endl;
			rpln(cout, sorted(*aall(*trim(*hraa(*uu1, *ur1, *hr))))); cout << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "historyRepasPersistent(hr,out)" << endl;
			historyRepasPersistent(*hr, out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto hr2 = persistentsHistoryRepa(in);
			cout << "hr2 = persistentsHistoryRepa(in)" << endl;
			rpln(cout, sorted(*aall(*trim(*hraa(*uu1, *ur1, *hr2))))); cout << endl;
			in.close();
		}

		auto uu2 = fsys(*gg);
		auto ur2 = uuur(*uu2);
		{
			auto tr = tttr(*uu2, *ur2, *ttcw);
			cout << "tr = tttr(uu2,ttcw)" << endl;
			cout << "trtt(uu2,tr)" << endl
				<< *trtt(*uu2, *ur2, *tr) << endl << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "transformRepasPersistent(tr,out)" << endl;
			transformRepasPersistent(*tr, out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto tr2 = persistentsTransformRepa(in);
			cout << "tr2 = persistentsTransformRepa(in)" << endl;
			cout << "trtt(uu2,tr2)" << endl
				<< *trtt(*uu2, *ur2, *tr2) << endl << endl;
			in.close();
		}

		{
			auto fr = fffr(*uu2, *ur2, *gg);
			cout << "fr = fffr(uu2,gg)" << endl;
			cout << "frff(uu2, fr)" << endl
				<< *frff(*uu2, *ur2, *fr) << endl << endl;

			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			cout << "fudRepasPersistent(fr,out)" << endl;
			fudRepasPersistent(*fr, out); cout << endl;
			out.close();

			std::ifstream in(filename, std::ios::binary);
			auto fr2 = persistentsFudRepa(in);
			cout << "fr2 = persistentsFudRepa(in)" << endl;
			cout << "frff(uu2, fr2)" << endl
				<< *frff(*uu2, *ur2, *fr2) << endl << endl;
			in.close();
		}
	}

	if (false)
	{
		auto fsys = fudsSystemImplied;
		auto dfund = decompFudsUnderlying;
		auto dfff = decompFudsFud;
		auto fvars = fudsSetVar;
		auto uuur = systemsSystemRepa;
		auto dfdr = systemsDecompFudsDecompFudRepa_u;
		auto drdf = systemsDecompFudRepasDecompFud_u;

		ifstream istrm("C:/zzz/caiks/NISTPy-master/NIST_model2.json");

		auto start = chrono::system_clock::now();
		auto df = persistentsDecompFud(istrm);
		auto end = chrono::system_clock::now();
		cout << "df = persistentsDecompFud(istrm) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		auto ff = dfff(*df);
		cout << "len(fvars(dfff(df)))" << endl
			<< fvars(*ff)->size() << endl << endl;

		auto uu = fsys(*ff);
		cout << "len(uu)" << endl
			<< uu->map_u().size() << endl << endl;

		auto ur = uuur(*uu);

		std::unique_ptr<DecompFudRepa> dr;
		try
		{
			StrVarPtrMap m;

			start = chrono::system_clock::now();
			std::stringstream str;
			systemsPersistent(*uu, str);
			end = chrono::system_clock::now();
			cout << "systemsPersistent(*uu, str) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto uu1 = persistentsSystem(str, m);
			end = chrono::system_clock::now();
			cout << "persistentsSystem(str, m) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			cout << "len(uu1)" << endl
				<< uu1->map_u().size() << endl << endl;

			start = chrono::system_clock::now();
			dr = dfdr(*uu, *ur, *df);
			end = chrono::system_clock::now();
			cout << "dr = dfdr(*uu,*df) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			decompFudRepasPersistent(*dr, out); cout << endl;
			out.close();
			end = chrono::system_clock::now();
			cout << "decompFudRepasPersistent(dr,out) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			std::ifstream in(filename, std::ios::binary);
			auto dr2 = persistentsDecompFudRepa(in);
			in.close();
			end = chrono::system_clock::now();
			cout << "dr2 = persistentsDecompFudRepa(in) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto df2 = drdf(*uu, *ur, *dr2);
			end = chrono::system_clock::now();
			cout << "df2 = drdf(*uu,*dr2) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			auto ff2 = dfff(*df2);
			cout << "len(fvars(dfff(df2)))" << endl
				<< fvars(*ff2)->size() << endl << endl;
		}
		catch (const std::exception& e)
		{
			cout << "caught exception: " << e.what() << endl << endl;
		}
		/*
		df = persistentsDecompFud(istrm) 0.171603s
		len(fvars(dfff(df)))
		515

		len(uu)
		515

		systemsPersistent(*uu, str) 0s
		persistentsSystem(str, m) 0.0156003s
		len(uu1)
		515

		dr = dfdr(*uu,*df) 0.0312006s

		decompFudRepasPersistent(dr,out) 0.0156003s
		dr2 = persistentsDecompFudRepa(in, m) 0s
		df2 = drdf(*uu,*dr2) 0.0624012s
		len(fvars(dfff(df2)))
		515
		*/
	}

	if (false)
	{
		auto fsys = fudsSystemImplied;
		auto dfund = decompFudsUnderlying;
		auto dfff = decompFudsFud;
		auto fvars = fudsSetVar;
		auto uuur = systemsSystemRepa;
		auto dfdr = systemsDecompFudsDecompFudRepa_u;
		auto drdf = systemsDecompFudRepasDecompFud_u;

		ifstream istrm("C:/zzz/caiks/NISTPy-master/NIST_model24.json");

		auto start = chrono::system_clock::now();
		auto df = persistentsDecompFud(istrm);
		auto end = chrono::system_clock::now();
		cout << "df = persistentsDecompFud(istrm) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		auto ff = dfff(*df);
		cout << "len(fvars(dfff(df)))" << endl
			<< fvars(*ff)->size() << endl << endl;

		auto uu = fsys(*ff);
		cout << "len(uu)" << endl
			<< uu->map_u().size() << endl << endl;

		auto ur = uuur(*uu);

		std::unique_ptr<DecompFudRepa> dr;
		try
		{
			StrVarPtrMap m;

			start = chrono::system_clock::now();
			std::stringstream str;
			systemsPersistent(*uu, str);
			end = chrono::system_clock::now();
			cout << "systemsPersistent(*uu, str) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto uu1 = persistentsSystem(str, m);
			end = chrono::system_clock::now();
			cout << "persistentsSystem(str, m) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			cout << "len(uu1)" << endl
				<< uu1->map_u().size() << endl << endl;

			start = chrono::system_clock::now();
			dr = dfdr(*uu, *ur, *df);
			end = chrono::system_clock::now();
			cout << "dr = dfdr(*uu,*df) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			std::string filename = "test.bin";
			std::ofstream out(filename, std::ios::binary);
			decompFudRepasPersistent(*dr, out); cout << endl;
			out.close();
			end = chrono::system_clock::now();
			cout << "decompFudRepasPersistent(dr,out) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			std::ifstream in(filename, std::ios::binary);
			auto dr2 = persistentsDecompFudRepa(in);
			in.close();
			end = chrono::system_clock::now();
			cout << "dr2 = persistentsDecompFudRepa(in) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			start = chrono::system_clock::now();
			auto df2 = drdf(*uu, *ur, *dr2);
			end = chrono::system_clock::now();
			cout << "df2 = drdf(*uu,*dr2) " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

			auto ff2 = dfff(*df2);
			cout << "len(fvars(dfff(df2)))" << endl
				<< fvars(*ff2)->size() << endl << endl;
		}
		catch (const std::exception& e)
		{
			cout << "caught exception: " << e.what() << endl << endl;
		}
		/*
		df = persistentsDecompFud(istrm) 11.0398s
		len(fvars(dfff(df)))
		5966

		len(uu)
		5966

		systemsPersistent(*uu, str) 0.0468009s
		persistentsSystem(str, m) 0.124802s
		len(uu1)
		5966

		dr = dfdr(*uu,*df) 4.35248s

		decompFudRepasPersistent(dr,out) 0.0468009s
		dr2 = persistentsDecompFudRepa(in, m) 0.0468009s
		df2 = drdf(*uu,*dr2) 4.24328s
		len(fvars(dfff(df2)))
		5966
		*/
	}

	if (false)
	{
		cout << "sizeof(std::size_t): " << sizeof(std::size_t) << endl;
		cout << "sizeof(Variable): " << sizeof(Variable) << endl;
	}

	if (false)
	{
		std::size_t mmax = 2;
		std::size_t n = 8;
		SizeListList mm{ SizeList{0} };
		for (std::size_t i = 2; i <= n; i++)
		{
			SizeListList mm1;
			for (auto& xx : mm)
				for (std::size_t j = 0; j < i; j++)
					if (j < mmax && j <= *max_element(xx.begin(), xx.end()) + 1)
					{
						SizeList yy(xx);
						yy.push_back(j);
						mm1.push_back(yy);
					}
			mm = mm1;
		}
		cout << "mm" << endl
			<< mm << endl << endl;
		cout << "mm.size()" << endl
			<< mm.size() << endl << endl;
	}

	if (false)
	{
		std::size_t mmax = 4;
		std::size_t n = 8;
		SizeListList mm{ SizeList{ 0 } };
		for (std::size_t i = 2; i <= n; i++)
		{
			SizeListList mm1;
			for (auto& xx : mm)
				for (std::size_t j = 0; j < i; j++)
					if (j < mmax && j <= *max_element(xx.begin(), xx.end()) + 1)
					{
						SizeList yy(xx);
						yy.push_back(j);
						mm1.push_back(yy);
					}
			mm = mm1;
		}
		cout << "mm" << endl
			<< mm << endl << endl;
		cout << "mm.size()" << endl
			<< mm.size() << endl << endl;
	}


	if (false)
	{
		std::size_t mmax = 2;
		std::size_t n = 8;
		auto start = chrono::system_clock::now();
		for (std::size_t k = 0; k < 1000; k++)
		{
			SizeListList mm{ SizeList{ 0 } };
			for (std::size_t i = 2; i <= n; i++)
			{
				SizeListList mm1;
				for (auto& xx : mm)
					for (std::size_t j = 0; j < i; j++)
						if (j < mmax && j <= *max_element(xx.begin(), xx.end()) + 1)
						{
							SizeList yy(xx);
							yy.push_back(j);
							mm1.push_back(yy);
						}
				mm = mm1;
			}
		}
		auto end = chrono::system_clock::now();
		cout << "partition " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		// partition 0.202803s
	}

	if (false)
	{
		std::size_t mmax = 4;
		std::size_t n = 8;
		auto start = chrono::system_clock::now();
		for (std::size_t k = 0; k < 1000; k++)
		{
			SizeListList mm{ SizeList{ 0 } };
			for (std::size_t i = 2; i <= n; i++)
			{
				SizeListList mm1;
				for (auto& xx : mm)
					for (std::size_t j = 0; j < i; j++)
						if (j < mmax && j <= *max_element(xx.begin(), xx.end()) + 1)
						{
							SizeList yy(xx);
							yy.push_back(j);
							mm1.push_back(yy);
						}
				mm = mm1;
			}
		}
		auto end = chrono::system_clock::now();
		cout << "partition " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		// partition 3.02644s
	}

	if (false)
	{
		std::size_t mmax = 4;
		std::size_t n = 8;
		auto start = chrono::system_clock::now();
		for (std::size_t k = 0; k < 1000; k++)
		{
			auto mm = std::unique_ptr<SizeListList>{new SizeListList{ SizeList{ 0 } } };
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
		}
		auto end = chrono::system_clock::now();
		cout << "partition " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		// partition 2.26203s
	}

	if (false)
	{
		std::size_t mmax = 4;
		std::size_t n = 8;

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
		std::vector<std::map<std::size_t, std::set<std::size_t>>> qq(mm->size() - 1);
		for (std::size_t i = 0; i < mm->size() - 1; i++)
		{
			auto& pp = qq[i];
			auto& ll = (*mm)[i + 1];
			for (std::size_t j = 0; j < n; j++)
				pp[ll[j]].insert(j);
		}
		cout << "qq" << endl
			<< qq << endl << endl;
		cout << "qq.size()" << endl
			<< qq.size() << endl << endl;
	}

	if (false)
	{
		std::size_t mmax = 2;
		std::size_t n = 8;
		auto start = chrono::system_clock::now();
		for (std::size_t k = 0; k < 1000; k++)
		{
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
			std::vector<std::map<std::size_t, std::set<std::size_t>>> qq(mm->size() - 1);
			for (std::size_t i = 0; i < mm->size() - 1; i++)
			{
				auto& pp = qq[i];
				auto& ll = (*mm)[i + 1];
				for (std::size_t j = 0; j < n; j++)
					pp[ll[j]].insert(j);
			}
		}
		auto end = chrono::system_clock::now();
		cout << "partition " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		// partition 0.499206s
	}

	if (false)
	{
		std::size_t mmax = 4;
		std::size_t n = 8;
		auto start = chrono::system_clock::now();
		for (std::size_t k = 0; k < 1000; k++)
		{
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
			std::vector<std::map<std::size_t, std::set<std::size_t>>> qq(mm->size() - 1);
			for (std::size_t i = 0; i < mm->size() - 1; i++)
			{
				auto& pp = qq[i];
				auto& ll = (*mm)[i + 1];
				for (std::size_t j = 0; j < n; j++)
					pp[ll[j]].insert(j);
			}
		}
		auto end = chrono::system_clock::now();
		cout << "partition " << ((chrono::duration<double>)(end - start)).count() << "s" << endl;

		// partition 12.4958s
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
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrsel = [](const HistoryRepa& hr, const SizeList& ll)
		{
			return eventsHistoryRepasHistoryRepaSelection_u(ll.size(), (std::size_t*)ll.data(), hr);
		};
		auto hrhrsel = [](const HistoryRepa& hr, const HistoryRepa& ss)
		{
			return historyRepasHistoryRepasHistoryRepaSelection_u(ss, hr);
		};
		auto hrhrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasHistoryRepaReduced_u(m, kk1.data(), hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto parter = parametersSystemsPartitionerMaxRollByMRepa_ui;

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

		auto ur = uuur(*uu);

		auto hr = hhhr(*uu, *ur, *hh);
		VarList ww{ pressure, rain, cloud, wind };
		{
			auto hrs = hrshuffle(*hr, 7);

			std::size_t mmax = 2;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 1;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i?",":"") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

		{
			auto hrs = hrshuffle(*hr, 7);

			std::size_t mmax = 4;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 1;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i ? "," : "") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

		{
			auto hrs = hrshuffle(*hr, 7);

			std::size_t mmax = 4;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 2;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i ? "," : "") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

		{
			auto hrs = hrshuffle(*hr, 17);

			std::size_t mmax = 2;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 1;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i ? "," : "") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

		{
			auto hrs = hrshuffle(*hr, 17);

			std::size_t mmax = 4;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 1;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i ? "," : "") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

		{
			auto hrs = hrshuffle(*hr, 17);

			std::size_t mmax = 4;
			std::size_t umax = 3 * 3 * 3;
			std::size_t pmax = 2;

			auto& vvi = ur->mapVarSize();
			std::size_t m = ww.size();
			SizeList ww1;
			for (std::size_t i = 0; i < m; i++)
				ww1.push_back(vvi[ww[i]]);

			auto t = parter(mmax, umax, pmax, ww1, *hr, *hrs);
			auto xx = std::move(std::get<0>(t));
			auto s = std::get<1>(t);
			cout << "steps: " << s << endl;
			for (auto& nn : *xx)
			{
				for (std::size_t i = 0; i < nn.size(); i++)
				{
					VarSet cc;
					for (std::size_t k = 0; k < nn[i].size(); k++)
						cc.insert(ww[nn[i][k]]);
					cout << (i ? "," : "") << cc;
				}
				cout << endl;
			}
			cout << endl;
		}

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
		auto uuur = systemsSystemRepa;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const VarList& kk)
		{
			auto& vvi = ur.mapVarSize();
			std::size_t m = kk.size();
			SizeList kk1;
			for (std::size_t i = 0; i < m; i++)
				kk1.push_back(vvi[kk[i]]);
			return setVarsHistoryRepasReduce_u(1.0, m, kk1.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
		auto parter = parametersSystemsPartitionerMaxRollByMRepa_ui;
		auto roller = parametersRollerMaximumRollExcludedSelfRepa_i;

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

		auto ur = uuur(*uu);

		auto hr = hhhr(*uu, *ur, *hh);
		VarList ww{ pressure, rain, cloud, wind };

		auto hrs = hrshuffle(*hr, 7);

		std::size_t mmax = 4;
		std::size_t umax = 3 * 3 * 3;
		std::size_t pmax = 2;

		auto z = (double)hr->size;
		auto ar = hrred(*hr, *ur, ww);
		auto ars = hrred(*hrs, *ur, ww);
		double y1 = ar->facLn() - ars->facLn();
		auto t = rrvqqy(mmax, umax, pmax, *ar, *ars, z, y1);
		auto xx = std::move(std::get<0>(t));
		auto s = std::get<1>(t);
		cout << "steps: " << s << endl;
		for (auto& nn : *xx)
		{
			cout << "partition " << nn << endl;
			auto t2 = roller(nn, *ar, *ars, z);
			auto xx2 = std::move(std::get<0>(t2));
			auto s2 = std::get<1>(t2);
			cout << "roll steps: " << s2 << endl;
			cout << "roll " << *xx2 << endl;
		}
		/*
		steps: 14
		partition [[0,2,3],[1]]
		roll steps: 3353
		roll [[[0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,1,1]]]
		partition [[0,3],[1,2]]
		roll steps: 440
		roll [[[0,0,1,1,1,1,1,0,0],[0,1,0,1,0,0,1,1,1]]]
		partition [[0,3],[1],[2]]
		roll steps: 171
		roll [[[0,0,1,2,1,2,1,0,0],[0,1,0],[0,1,1]]]
		partition [[0],[1],[2,3]]
		roll steps: 185
		roll [[[0,1,1],[0,1,2],[0,0,0,1,1,0,1,1,1]]]
		partition [[0],[1],[2],[3]]
		roll steps: 30
		roll [[[0,1,0],[0,1,1],[0,1,1],[0,0,1]]]
		*/
	}

	if (false)
	{
		auto sysreg = systemRegular_u;
		auto sys = histogramsSystemImplied;
		auto uunion = pairSystemsUnion;
		auto uvars = systemsSetVar;
		auto cart = systemsSetVarsSetStateCartesian_u;
		auto reframe = histogramsMapVarsFrame_u;
		auto mul = pairHistogramsMultiply;
		auto add = pairHistogramsAdd_u;
		auto resize = histogramsResize;
		auto scalar = histogramScalar_u;
		auto regpivot = histogramRegularUnitPivot_u;
		auto regcart = histogramRegularCartesian_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
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
		auto cdaa = [](const SizeListList& ll)
		{
			auto llss = listsState;
			auto llaa = listsHistogram_u;
			std::vector<StateRationalPair> qq;
			qq.reserve(ll.size());
			for (auto& mm : ll)
			{
				std::vector<VarValPair> ss;
				ss.reserve(mm.size());
				for (int i = 0; i < mm.size(); i++)
					ss.push_back(VarValPair(Variable(i+1), Value((int)mm[i])));
				qq.push_back(StateRationalPair(*llss(ss), 1));
			}
			return llaa(qq);
		};
		auto cdtt = [cdaa](const SizeList& pp, const SizeListList& ll)
		{
			auto aa = cdaa(ll);
			auto vv = histogramsSetVar(*aa);
			VarVarUMap nn;
			std::size_t i = 0;
			for (auto& v : *vv)
			{
				if (i < pp.size())
					nn.insert_or_assign(v,Variable((int)pp[i]));
				i++;
			}
			auto bb = histogramsMapVarsFrame_u(*aa, nn);
			VarUSet ww;
			if (pp.size())
				ww.insert(Variable((int)pp.back()));
			return histogramsSetVarsTransform(*bb, ww);
		};
		auto llff = setTransformsFud_u;
		auto fsys = fudsSystemImplied;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto uuur = systemsSystemRepa;
		auto uruu = systemsRepasSystem;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const SizeList& kk)
		{
			return setVarsHistoryRepasReduce_u(1.0, kk.size(), kk.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto frmul = historyRepasFudRepasMultiply_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;
		auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui;
		auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
		auto parter = parametersSystemsPartitionerMaxRollByMRepa_ui;
		auto roller = parametersRollerMaximumRollExcludedSelfRepa_i;
		auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;
		auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u;

		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(900,*mul(*regpivot(3, 3), *reframe(*regcart(1, 1), VarVarUMap{ { Variable(1), Variable(4) } })));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = add(*hraa(*uu, *ur, *hrshuffle(*hr, 1)), *hraa(*uu, *ur, *hrshuffle(*hr, 2)));
			auto hrs = aahr(*uu, *ur, *rr);
			auto pr = hrpr(*hr);
			auto prs = hrpr(*hrs);

			auto t = tupler(3*3*3, 100, 1000, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			std::vector<VarSet> pp;
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3 * 3, 100, 2, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3 * 3, 100, 1, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3 * 3*3, 100, 1000, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3 * 3 * 3, 100, 2, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3 * 3 * 3, 100, 1, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;

			t = tupler(3 * 3, 100, 1000, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3, 100, 2, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
			t = tupler(3 * 3, 100, 1, 1, vv, FudRepa(), *hr, *pr, *hrs, *prs);
			pp.clear();
			for (auto& kk : *std::get<0>(t))
			{
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				pp.push_back(qq);
			}
			cout << pp << endl << endl;
		}

		cout << endl;
		{
				auto uu = sysreg(3, 4);
				auto ur = uuur(*uu);
				auto vv = SizeList{ 0,1,2,3 };
				auto aa = resize(150, *regdiag(3, 4));
				auto hr = aahr(*uu, *ur, *aa);
				auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
				aa = add(*resize(150, *regdiag(3, 4)), *rr);
				hr = aahr(*uu, *ur, *aa);
				auto hrs = hrshuffle(*hr, 1);

				auto t = parter(4, 3*3*3, 20, vv, *hr, *hrs);
				for (auto& nn : *std::get<0>(t))
				{
					std::vector<VarSet> pp;
					for (auto& kk : nn)
					{
						VarSet qq;
						for (std::size_t i = 0; i < kk.size(); i++)
							qq.insert(*(ur->listVarSizePair[kk[i]]).first);
						pp.push_back(qq);
					}
					cout << pp << endl;
				}
				cout << endl;

				t = parter(4, 3 * 3 * 3, 1, vv, *hr, *hrs);
				for (auto& nn : *std::get<0>(t))
				{
					std::vector<VarSet> pp;
					for (auto& kk : nn)
					{
						VarSet qq;
						for (std::size_t i = 0; i < kk.size(); i++)
							qq.insert(*(ur->listVarSizePair[kk[i]]).first);
						pp.push_back(qq);
					}
					cout << pp << endl;
				}
		}

		cout << endl;
		{
			auto aa = resize(900, *mul(*regpivot(3, 3), *reframe(*regcart(1, 1), VarVarUMap{ { Variable(1), Variable(4) } })));
			auto uu = sys(*aa);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = add(*hraa(*uu, *ur, *hrshuffle(*hr, 1)), *hraa(*uu, *ur, *hrshuffle(*hr, 2)));
			auto hrs = aahr(*uu, *ur, *rr);

			auto z = (double)hr->size;
			auto ar = hrred(*hr, *ur, vv);
			auto ars = hrred(*hrs, *ur, vv);
			double y1 = ar->facLn() - ars->facLn();
			auto t = rrvqqy(4, 3 * 3 * 3, 1, *ar, *ars, z, y1);
			for (auto& nn : *std::get<0>(t))
			{
				cout << "partition " << nn << endl;
				auto t2 = roller(nn, *ar, *ars, z);
				cout << "roll " << *std::get<0>(t2) << endl;
			}
			cout << endl;
			t = rrvqqy(4, 3 * 3 * 3, 20, *ar, *ars, z, y1);
			for (auto& nn : *std::get<0>(t))
			{
				cout << "partition " << nn << endl;
				auto t2 = roller(nn, *ar, *ars, z);
				cout << "roll " << *std::get<0>(t2) << endl;
			}
		}

		cout << endl;
		{
			cout << *cdaa(SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) << endl;
			cout << *cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) << endl;
			cout << *cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) << endl;
		}

		cout << endl;
		{
			auto uu = sysreg(3, 3);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2};
			auto aa = mul(*scalar(81), *regpivot(3, 3));
			auto hr = aahr(*uu, *ur, *aa);
			auto hrs = hrshuffle(*hr, 17);
			auto pr = hrpr(*hr);
			auto prs = hrpr(*hrs);

			{
				auto t = deriveder(3 * 3 * 3, 100, FudRepa(), *hr, *pr, *hrs, *prs);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,2 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 3,17 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,2 } }), cdtt(SizeList{ 3,17 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 3,17 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 3,18 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 100, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 3,17 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 3,18 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 1, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}

			{
				auto ff = llff(TransformPtrList{ cdtt(SizeList{ 1,15 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 2,16 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 15,17 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }), cdtt(SizeList{ 16,18 }, SizeListList{ SizeList{ 1,1 },SizeList{ 2,2 }, SizeList{ 3,3 } }) });
				auto uu1 = uunion(*uu, *fsys(*ff));
				auto ur1 = uuur(*uu1);
				auto fr = fffr(*uu1, *ur1, *ff);
				auto hr1 = frmul(*hr, *fr);
				auto hrs1 = frmul(*hrs, *fr);
				auto pr1 = hrpr(*hr1);
				auto prs1 = hrpr(*hrs1);
				auto t = deriveder(3 * 3 * 3, 1, *fr, *hr1, *pr1, *hrs1, *prs1);
				std::vector<VarSet> pp;
				for (auto& kk : *std::get<0>(t))
				{
					VarSet qq;
					for (std::size_t i = 0; i < kk.second.size(); i++)
						qq.insert(*(ur1->listVarSizePair[kk.second[i]]).first);
					pp.push_back(qq);
				}
				cout << pp << endl << endl;
			}



		}

		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(900, *mul(*regpivot(3, 3), *reframe(*regcart(1, 1), VarVarUMap{ { Variable(1), Variable(4) } })));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = add(*hraa(*uu, *ur, *hrshuffle(*hr, 1)), *hraa(*uu, *ur, *hrshuffle(*hr, 2)));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, vv, *hr, *hrs, 1, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			>>> layerer
			>>> layer	   fud: 1  layer: 1
			substrate cardinality: 4
			fud cardinality: 0
			tupler  searched: 5	 rate: 26615.7
			tupler 0.000187859s
			tuple cardinality: 1
			max tuple algn: 586.053
			layer cardinality: 2
			parter  searched: 3	 rate: 44532.2
			parter 6.7367e-05s
			roller  searched: 158   rate: 815977
			roller 0.000193633s
			der vars algn density: 154.389
			dervarser	   searched: 2	 rate: 45977
			dervarser 4.35e-05s
			application 0.11146s
			<<< layer 0.119519s
			>>> layer	   fud: 1  layer: 2
			substrate cardinality: 4
			fud cardinality: 2
			tupler  searched: 25	rate: 54254.2
			tupler 0.000460794s
			tuple cardinality: 1
			max tuple algn: 471.452
			layer cardinality: 2
			parter  searched: 3	 rate: 57302.2
			parter 5.2354e-05s
			roller  searched: 34	rate: 664075
			roller 5.1199e-05s
			der vars algn density: 154.389
			dervarser	   searched: 4	 rate: 41563.2
			dervarser 9.6239e-05s
			application 0.107381s
			<<< layer 0.11536s
			<<< layerer 0.251153s
			a = 154.389
			kk = {<<1,1>,1>,<<1,1>,2>}
			ff = {
			({({(1,1),(2,1),(<<1,1>,1>,0)},1 % 1),({(1,1),(2,2),(<<1,1>,1>,0)},1 % 1),({(1,1),(2,3),(<<1,1>,1>,1)},1 % 1),({(1,2),(2,1),(<<1,1>,1>,1)},1 % 1),({(1,2),(2,2),(<<1,1>,1>,1)},1 % 1),({(1,2),(2,3),(<<1,1>,1>,1)},1 % 1),({(1,3),(2,1),(<<1,1>,1>,1)},1 % 1),({(1,3),(2,2),(<<1,1>,1>,1)},1 % 1),({(1,3),(2,3),(<<1,1>,1>,1)},1 % 1)},{<<1,1>,1>}),
			({({(3,1),(<<1,1>,2>,0)},1 % 1),({(3,2),(<<1,1>,2>,1)},1 % 1),({(3,3),(<<1,1>,2>,1)},1 % 1)},{<<1,1>,2>})}
*/
		}
		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(150, *regdiag(3, 4));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			aa = add(*resize(1500, *regdiag(3, 4)), *rr);
			hr = aahr(*uu, *ur, *aa);
			rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, vv, *hr, *hrs, 1, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			>>> layerer
			>>> layer	   fud: 1  layer: 1
			substrate cardinality: 4
			fud cardinality: 0
			tupler  searched: 10	rate: 39121.9
			tupler 0.000255611s
			tuple cardinality: 1
			max tuple algn: 2853.7
			layer cardinality: 2
			parter  searched: 3	 rate: 39359.2
			parter 7.6221e-05s
			roller  searched: 143   rate: 530673
			roller 0.000269469s
			der vars algn density: 475.058
			dervarser	   searched: 2	 rate: 41232.9
			dervarser 4.8505e-05s
			application 0.10771s
			<<< layer 0.115672s
			>>> layer	   fud: 1  layer: 2
			substrate cardinality: 4
			fud cardinality: 2
			tupler  searched: 29	rate: 37367.6
			tupler 0.000776073s
			tuple cardinality: 1
			max tuple algn: 2200
			layer cardinality: 2
			parter  searched: 3	 rate: 42124.8
			parter 7.1217e-05s
			roller  searched: 143   rate: 555262
			roller 0.000257536s
			der vars algn density: 883.399
			dervarser	   searched: 10	rate: 39720.1
			dervarser 0.000251762s
			application 0.106361s
			<<< layer 0.118203s
			>>> layer	   fud: 1  layer: 3
			substrate cardinality: 4
			fud cardinality: 4
			tupler  searched: 53	rate: 41607.1
			tupler 0.00127382s
			tuple cardinality: 1
			max tuple algn: 1382.39
			layer cardinality: 2
			parter  searched: 3	 rate: 25220.3
			parter 0.000118952s
			roller  searched: 143   rate: 600114
			roller 0.000238288s
			der vars algn density: 883.399
			dervarser	   searched: 14	rate: 40319
			dervarser 0.000347231s
			application 0.106129s
			<<< layer 0.11497s
			<<< layerer 0.373039s
			a = 883.399
			kk = {<<1,1>,1>,<<1,1>,2>,<<1,2>,1>}
			ff = {
			({({(2,1),(3,1),(<<1,1>,1>,0)},1 % 1),({(2,1),(3,2),(<<1,1>,1>,1)},1 % 1),({(2,1),(3,3),(<<1,1>,1>,1)},1 % 1),({(2,2),(3,1),(<<1,1>,1>,1)},1 % 1),({(2,2),(3,2),(<<1,1>,1>,1)},1 % 1),({(2,2),(3,3),(<<1,1>,1>,1)},1 % 1),({(2,3),(3,1),(<<1,1>,1>,1)},1 % 1),({(2,3),(3,2),(<<1,1>,1>,1)},1 % 1),({(2,3),(3,3),(<<1,1>,1>,2)},1 % 1)},{<<1,1>,1>}),
			({({(4,1),(<<1,1>,2>,0)},1 % 1),({(4,2),(<<1,1>,2>,1)},1 % 1),({(4,3),(<<1,1>,2>,2)},1 % 1)},{<<1,1>,2>}),
			({({(1,1),(4,1),(<<1,2>,1>,0)},1 % 1),({(1,1),(4,2),(<<1,2>,1>,1)},1 % 1),({(1,1),(4,3),(<<1,2>,1>,1)},1 % 1),({(1,2),(4,1),(<<1,2>,1>,1)},1 % 1),({(1,2),(4,2),(<<1,2>,1>,1)},1 % 1),({(1,2),(4,3),(<<1,2>,1>,1)},1 % 1),({(1,3),(4,1),(<<1,2>,1>,1)},1 % 1),({(1,3),(4,2),(<<1,2>,1>,1)},1 % 1),({(1,3),(4,3),(<<1,2>,1>,2)},1 % 1)},{<<1,2>,1>}),({({(<<1,1>,1>,0),(<<1,2>,2>,0)},1 % 1),({(<<1,1>,1>,1),(<<1,2>,2>,1)},1 % 1),({(<<1,1>,1>,2),(<<1,2>,2>,2)},1 % 1)},{<<1,2>,2>})}
			*/
		}
	}

	if (false)
	{
		auto sysreg = systemRegular_u;
		auto sys = histogramsSystemImplied;
		auto uunion = pairSystemsUnion;
		auto uvars = systemsSetVar;
		auto cart = systemsSetVarsSetStateCartesian_u;
		auto reframe = histogramsMapVarsFrame_u;
		auto mul = pairHistogramsMultiply;
		auto add = pairHistogramsAdd_u;
		auto resize = histogramsResize;
		auto scalar = histogramScalar_u;
		auto regpivot = histogramRegularUnitPivot_u;
		auto regcart = histogramRegularCartesian_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
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
		auto cdaa = [](const SizeListList& ll)
		{
			auto llss = listsState;
			auto llaa = listsHistogram_u;
			std::vector<StateRationalPair> qq;
			qq.reserve(ll.size());
			for (auto& mm : ll)
			{
				std::vector<VarValPair> ss;
				ss.reserve(mm.size());
				for (int i = 0; i < mm.size(); i++)
					ss.push_back(VarValPair(Variable(i + 1), Value((int)mm[i])));
				qq.push_back(StateRationalPair(*llss(ss), 1));
			}
			return llaa(qq);
		};
		auto cdtt = [cdaa](const SizeList& pp, const SizeListList& ll)
		{
			auto aa = cdaa(ll);
			auto vv = histogramsSetVar(*aa);
			VarVarUMap nn;
			std::size_t i = 0;
			for (auto& v : *vv)
			{
				if (i < pp.size())
					nn.insert_or_assign(v, Variable((int)pp[i]));
				i++;
			}
			auto bb = histogramsMapVarsFrame_u(*aa, nn);
			VarUSet ww;
			if (pp.size())
				ww.insert(Variable((int)pp.back()));
			return histogramsSetVarsTransform(*bb, ww);
		};
		auto llff = setTransformsFud_u;
		auto fsys = fudsSystemImplied;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto uuur = systemsSystemRepa;
		auto uruu = systemsRepasSystem;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const SizeList& kk)
		{
			return setVarsHistoryRepasReduce_u(1.0, kk.size(), kk.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto frmul = historyRepasFudRepasMultiply_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;
		auto applicationer = parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa;

		{
			auto uu = sysreg(3, 5);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3,4 };
			auto aa = resize(150, *regdiag(3, 5));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			aa = add(*resize(1500, *regdiag(3, 5)), *rr);
			hr = aahr(*uu, *ur, *aa);
			hr->transpose();
			auto dr = applicationer(27, 3, 81, 2, 2, 2, 9, 1, 1, 2, 17, vv, *hr, *ur);
			cout << "dr = " << *dr << endl;
		}

		cout << endl;
		{
			auto uu = sysreg(3, 5);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3,4 };
			auto aa = resize(150, *regdiag(3, 5));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			aa = add(*resize(1500, *regdiag(3, 5)), *rr);
			hr = aahr(*uu, *ur, *aa);
			hr->transpose();
			auto dr = applicationer(27, 3, 81, 2, 2, 2, 9, 1, 15, 2, 17, vv, *hr, *ur);
			cout << "dr = " << *dr << endl;
			cout << "treesPaths(dr.slices)" << endl;
			rpln(cout, *treesPaths(*dr->slices)); cout << endl;

			{
				std::string filename = "test.bin";
				std::ofstream out(filename, std::ios::binary);
				systemRepasPersistent(*ur, out); cout << endl;
				out.close();
				cout << "systemRepasPersistent(ur,out) " << endl;
				cout << "ur = " << *ur << endl;

				StrVarPtrMap m;
				std::ifstream in(filename, std::ios::binary);
				auto ur2 = persistentsSystemRepa(in,m);
				in.close();
				cout << "ur2 = persistentsSystemRepa(in) " << endl;
				cout << "ur2 = " << *ur2 << endl;
			}

			{
				std::string filename = "test.bin";
				std::ofstream out(filename, std::ios::binary);
				applicationRepasPersistent(*dr, out); cout << endl;
				out.close();
				cout << "applicationRepasPersistent(dr,out) " << endl;

				std::ifstream in(filename, std::ios::binary);
				auto dr2 = persistentsApplicationRepa(in);
				in.close();
				cout << "dr2 = persistentsApplicationRepa(in) " << endl;
				cout << "dr2 = " << *dr2 << endl;
			}
		}

	}

	if (true)
	{
		auto regcart = histogramRegularCartesian_u;
		auto regsing = histogramRegularUnitSingleton_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
		auto sys = histogramsSystemImplied;
		auto uuur = systemsSystemRepa;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hrhs = historyRepasHistorySparse;
		auto hshr = historySparsesHistoryRepa;
		auto llhs = listSetIntsHistorySparse;
		auto hsll = historySparsesListSetInt;
		
		auto aa = regdiag(2, 3);
		EVALL(*aa);

		auto uu = sys(*aa);
		auto ur = uuur(*uu);
		auto hr = aahr(*uu, *ur, *aa);
		EVALL(*hr);

		auto hs = hrhs(*hr);
		EVALL(*hs);

		EVALL(*hshr(*hs));
		
		EVALL(*hsll(*hs));
		
		EVALL(*llhs(*hsll(*hs)));

		std::string filename = "test.bin";
		std::ofstream out(filename, std::ios::binary);
		ECHO(historySparsesPersistent(*hs, out));
		out.close();

		std::ifstream in(filename, std::ios::binary);
		ECHO(auto hs2 = persistentsHistorySparse(in));
		in.close();
		EVALL(*hs2);
		EVALL(*hshr(*hs2));
	}

	if (false)
	{
		auto sysreg = systemRegular_u;
		auto sys = histogramsSystemImplied;
		auto uunion = pairSystemsUnion;
		auto uvars = systemsSetVar;
		auto cart = systemsSetVarsSetStateCartesian_u;
		auto reframe = histogramsMapVarsFrame_u;
		auto mul = pairHistogramsMultiply;
		auto add = pairHistogramsAdd_u;
		auto resize = histogramsResize;
		auto scalar = histogramScalar_u;
		auto regpivot = histogramRegularUnitPivot_u;
		auto regcart = histogramRegularCartesian_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
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
		auto cdaa = [](const SizeListList& ll)
		{
			auto llss = listsState;
			auto llaa = listsHistogram_u;
			std::vector<StateRationalPair> qq;
			qq.reserve(ll.size());
			for (auto& mm : ll)
			{
				std::vector<VarValPair> ss;
				ss.reserve(mm.size());
				for (int i = 0; i < mm.size(); i++)
					ss.push_back(VarValPair(Variable(i+1), Value((int)mm[i])));
				qq.push_back(StateRationalPair(*llss(ss), 1));
			}
			return llaa(qq);
		};
		auto cdtt = [cdaa](const SizeList& pp, const SizeListList& ll)
		{
			auto aa = cdaa(ll);
			auto vv = histogramsSetVar(*aa);
			VarVarUMap nn;
			std::size_t i = 0;
			for (auto& v : *vv)
			{
				if (i < pp.size())
					nn.insert_or_assign(v,Variable((int)pp[i]));
				i++;
			}
			auto bb = histogramsMapVarsFrame_u(*aa, nn);
			VarUSet ww;
			if (pp.size())
				ww.insert(Variable((int)pp.back()));
			return histogramsSetVarsTransform(*bb, ww);
		};
		auto llff = setTransformsFud_u;
		auto fsys = fudsSystemImplied;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto uuur = systemsSystemRepa;
		auto uruu = systemsRepasSystem;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const SizeList& kk)
		{
			return setVarsHistoryRepasReduce_u(1.0, kk.size(), kk.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto frmul = historyRepasFudRepasMultiply_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;
		auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui;
		auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
		auto parter = parametersSystemsPartitionerMaxRollByMRepa_ui;
		auto roller = parametersRollerMaximumRollExcludedSelfRepa_i;
		auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;
		auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u;

		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(900, *mul(*regpivot(3, 3), *reframe(*regcart(1, 1), VarVarUMap{ { Variable(1), Variable(4) } })));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = add(*hraa(*uu, *ur, *hrshuffle(*hr, 1)), *hraa(*uu, *ur, *hrshuffle(*hr, 2)));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, vv, *hr, *hrs, 1, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			>>> layerer
			>>> layer       fud: 1  layer: 1
			substrate cardinality: 4
			fud cardinality: 0
			tupler  searched: 5     rate: 17743.1
			tupler 0.0002818s
			tuple cardinality: 1
			max tuple algn: 588.469
			layer cardinality: 2
			parter  searched: 3     rate: 9216.59
			parter 0.0003255s
			roller  searched: 158   rate: 728111
			roller 0.000217s
			der vars algn density: 154.34
			dervarser       searched: 2     rate: 124224
			dervarser 1.61e-05s
			application 0.0004421s
			<<< layer 0.0012853s
			>>> layer       fud: 1  layer: 2
			substrate cardinality: 4
			fud cardinality: 2
			tupler  searched: 25    rate: 173370
			tupler 0.0001442s
			tuple cardinality: 1
			max tuple algn: 607.304
			layer cardinality: 1
			parter  searched: 3     rate: 303030
			parter 9.9e-06s
			roller  searched: 34    rate: 2.85714e+06
			roller 1.19e-05s
			der vars algn density: 154.34
			dervarser       searched: 3     rate: 160428
			dervarser 1.87e-05s
			application 0.0001674s
			<<< layer 0.0003704s
			<<< layerer 0.0017068s
			a = 154.34
			kk = {<<<0,1>,1>,1>,<<<0,1>,1>,2>}
			ff = {({({(1,1),(3,1),(<<<0,1>,1>,1>,0)},1 % 1),({(1,1),(3,2),(<<<0,1>,1>,1>,0)},1 % 1),({(1,1),(3,3),(<<<0,1>,1>,1>,0)},1 % 1),({(1,2),(3,1),(<<<0,1>,1>,1>,1)},1 % 1),({(1,2),(3,2),(<<<0,1>,1>,1>,1)},1 % 1),({(1,2),(3,3),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,1),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,2),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,3),(<<<0,1>,1>,1>,1)},1 % 1)},{<<<0,1>,1>,1>}),({({(2,1),(<<<0,1>,1>,2>,0)},1 % 1),({(2,2),(<<<0,1>,1>,2>,1)},1 % 1),({(2,3),(<<<0,1>,1>,2>,1)},1 % 1)},{<<<0,1>,1>,2>})}
*/
		}
		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(150, *regdiag(3, 4));
			auto hr = aahr(*uu, *ur, *aa);
			// auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 17));
			aa = add(*resize(1500, *regdiag(3, 4)), *rr);
			hr = aahr(*uu, *ur, *aa);
			// rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			rr = hraa(*uu, *ur, *hrshuffle(*hr, 17));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, vv, *hr, *hrs, 1, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			>>> layerer
			>>> layer       fud: 1  layer: 1
			substrate cardinality: 4
			fud cardinality: 0
			tupler  searched: 10    rate: 130548
			tupler 7.66e-05s
			tuple cardinality: 1
			max tuple algn: 2882.52
			layer cardinality: 1
			parter  searched: 3     rate: 209790
			parter 1.43e-05s
			roller  searched: 143   rate: 1.63242e+06
			roller 8.76e-05s
			no der vars sets
			dervarser       searched: 0     rate: 0
			dervarser 3.5e-06s
			application 0.0001355s
			<<< layer 0.0003209s
			<<< layerer 0.0003933s
			ff = {}
			*/
		}
	}
	
	if (false)
	{
		auto sysreg = systemRegular_u;
		auto sys = histogramsSystemImplied;
		auto uunion = pairSystemsUnion;
		auto uvars = systemsSetVar;
		auto cart = systemsSetVarsSetStateCartesian_u;
		auto reframe = histogramsMapVarsFrame_u;
		auto mul = pairHistogramsMultiply;
		auto add = pairHistogramsAdd_u;
		auto resize = histogramsResize;
		auto scalar = histogramScalar_u;
		auto regpivot = histogramRegularUnitPivot_u;
		auto regcart = histogramRegularCartesian_u;
		auto regdiag = histogramRegularUnitDiagonal_u;
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
		auto cdaa = [](const SizeListList& ll)
		{
			auto llss = listsState;
			auto llaa = listsHistogram_u;
			std::vector<StateRationalPair> qq;
			qq.reserve(ll.size());
			for (auto& mm : ll)
			{
				std::vector<VarValPair> ss;
				ss.reserve(mm.size());
				for (int i = 0; i < mm.size(); i++)
					ss.push_back(VarValPair(Variable(i+1), Value((int)mm[i])));
				qq.push_back(StateRationalPair(*llss(ss), 1));
			}
			return llaa(qq);
		};
		auto cdtt = [cdaa](const SizeList& pp, const SizeListList& ll)
		{
			auto aa = cdaa(ll);
			auto vv = histogramsSetVar(*aa);
			VarVarUMap nn;
			std::size_t i = 0;
			for (auto& v : *vv)
			{
				if (i < pp.size())
					nn.insert_or_assign(v,Variable((int)pp[i]));
				i++;
			}
			auto bb = histogramsMapVarsFrame_u(*aa, nn);
			VarUSet ww;
			if (pp.size())
				ww.insert(Variable((int)pp.back()));
			return histogramsSetVarsTransform(*bb, ww);
		};
		auto llff = setTransformsFud_u;
		auto fsys = fudsSystemImplied;
		auto aarr = systemsHistogramsHistogramRepa_u;
		auto rraa = systemsHistogramRepasHistogram_u;
		auto uuur = systemsSystemRepa;
		auto uruu = systemsRepasSystem;
		auto araa = systemsHistogramRepasHistogram_u;
		auto aahr = [](const System& uu, const SystemRepa& ur, const Histogram& aa)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, *histogramsHistory_u(aa), 1);
		};
		auto hraa = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return historiesHistogram(*systemsHistoryRepasHistory_u(uu, ur, hr));
		};
		auto hhhr = [](const System& uu, const SystemRepa& ur, const History& hh)
		{
			return systemsHistoriesHistoryRepa_u(uu, ur, hh, 1);
		};
		auto hrhh = [](const System& uu, const SystemRepa& ur, const HistoryRepa& hr)
		{
			return systemsHistoryRepasHistory_u(uu, ur, hr);
		};
		auto hrred = [](const HistoryRepa& hr, const SystemRepa& ur, const SizeList& kk)
		{
			return setVarsHistoryRepasReduce_u(1.0, kk.size(), kk.data(), hr);
		};
		auto hrpr = historyRepasRed;
		auto hrshuffle = historyRepasShuffle_u;
		auto frmul = historyRepasFudRepasMultiply_u;
		auto fffr = systemsFudsFudRepa_u;
		auto frff = systemsFudRepasFud_u;
		auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui;
		auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
		auto parter = parametersSystemsPartitionerMaxRollByMRepa_ui;
		auto roller = parametersRollerMaximumRollExcludedSelfRepa_i;
		auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;
		auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestLogIORepa_up;

		auto log = [](const std::string& str)
		{
			std::cout << str << std::endl;
			return;
		};
		
		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(900, *mul(*regpivot(3, 3), *reframe(*regcart(1, 1), VarVarUMap{ { Variable(1), Variable(4) } })));
			auto hr = aahr(*uu, *ur, *aa);
			auto rr = add(*hraa(*uu, *ur, *hrshuffle(*hr, 1)), *hraa(*uu, *ur, *hrshuffle(*hr, 2)));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, 1, vv, *hr, *hrs, 1, log, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			layerer status: begin
			layer   status: begin   fud: 1  layer: 1        substrate cardinality: 4        fud cardinality: 0
			tupler  searched: 5     rate: 11739.8   time 0.0004259s tuple cardinality: 1
			layer   max tuple algn: 588.469 layer cardinality: 2
			parter  searched: 3     rate: 126582    time 2.37e-05s
			roller  searched: 158   rate: 3.62385e+06       time 4.36e-05s
			layer   der vars algn density: 154.34
			dervarser       searched: 2     rate: 9350.16   time 0.0002139s
			application     time 0.0001708s
			layer   status: end     time 0.0008824s
			layer   status: begin   fud: 1  layer: 2        substrate cardinality: 4        fud cardinality: 2
			tupler  searched: 25    rate: 54229.9   time 0.000461s  tuple cardinality: 1
			layer   max tuple algn: 607.304 layer cardinality: 1
			parter  searched: 3     rate: 252101    time 1.19e-05s
			roller  searched: 34    rate: 2.5e+06   time 1.36e-05s
			layer   der vars algn density: 154.34
			dervarser       searched: 3     rate: 14910.5   time 0.0002012s
			application     time 0.0001902s
			layer   status: end     time 0.0008948s
			layerer status: end     time 0.0018519s
			a = 154.34
			kk = {<<<0,1>,1>,1>,<<<0,1>,1>,2>}
			ff = {({({(1,1),(3,1),(<<<0,1>,1>,1>,0)},1 % 1),({(1,1),(3,2),(<<<0,1>,1>,1>,0)},1 % 1),({(1,1),(3,3),(<<<0,1>,1>,1>,0)},1 % 1),({(1,2),(3,1),(<<<0,1>,1>,1>,1)},1 % 1),({(1,2),(3,2),(<<<0,1>,1>,1>,1)},1 % 1),({(1,2),(3,3),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,1),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,2),(<<<0,1>,1>,1>,1)},1 % 1),({(1,3),(3,3),(<<<0,1>,1>,1>,1)},1 % 1)},{<<<0,1>,1>,1>}),({({(2,1),(<<<0,1>,1>,2>,0)},1 % 1),({(2,2),(<<<0,1>,1>,2>,1)},1 % 1),({(2,3),(<<<0,1>,1>,2>,1)},1 % 1)},{<<<0,1>,1>,2>})}
			*/
		}
		cout << endl;
		{
			auto uu = sysreg(3, 4);
			auto ur = uuur(*uu);
			auto vv = SizeList{ 0,1,2,3 };
			auto aa = resize(150, *regdiag(3, 4));
			auto hr = aahr(*uu, *ur, *aa);
			// auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			auto rr = hraa(*uu, *ur, *hrshuffle(*hr, 17));
			aa = add(*resize(1500, *regdiag(3, 4)), *rr);
			hr = aahr(*uu, *ur, *aa);
			// rr = hraa(*uu, *ur, *hrshuffle(*hr, 1));
			rr = hraa(*uu, *ur, *hrshuffle(*hr, 17));
			auto hrs = aahr(*uu, *ur, *rr);
			hr->transpose();
			hrs->transpose();
			auto t = layerer(27, 3, 27, 2, 2, 2, 9, 1, 1, vv, *hr, *hrs, 1, log, *ur);
			auto& fr = std::get<0>(t);
			auto& mm = std::get<1>(t);
			if (mm->size())
			{
				auto& a = mm->back().first;
				auto& kk = mm->back().second;
				VarSet qq;
				for (std::size_t i = 0; i < kk.size(); i++)
					qq.insert(*(ur->listVarSizePair[kk[i]]).first);
				cout << "a = " << a << endl;
				cout << "kk = " << qq << endl;
			}
			uruu(*ur, *uu);
			auto ff = frff(*uu, *ur, *fr);
			cout << "ff = " << *ff << endl;
			/*
			layerer status: begin
			layer   status: begin   fud: 1  layer: 1        substrate cardinality: 4        fud cardinality: 0
			tupler  searched: 10    rate: 27367.3   time 0.0003654s tuple cardinality: 1
			layer   max tuple algn: 2882.52 layer cardinality: 1
			parter  searched: 3     rate: 223881    time 1.34e-05s
			roller  searched: 143   rate: 2.31392e+06       time 6.18e-05s
			layer   no der vars sets
			dervarser       searched: 0     rate: 0 time 0.0001175s
			application     time 0.0001846s
			layer   status: end     time 0.0007611s
			layerer status: end     time 0.0008309s
			ff = {}
			*/
		}
	}
	
	return 0;
}
