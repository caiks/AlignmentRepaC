#include "AlignmentPracticableIORepa.h"
#include <chrono>
#include <ctime>

using namespace Alignment;

typedef std::chrono::duration<double> sec;
typedef std::chrono::system_clock clk;

// parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> HistogramRepaRed -> HistoryRepa-> HistogramRepaRed -> Integer ->
//   IO(SystemRepa, FudRepa, [(Double, [VariableRepa]])
std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> Alignment::parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, const SizeList& vv, const HistoryRepa& hr, const HistoryRepa& hrs, std::size_t f, SystemRepa& ur)
{
    auto hrred = [](double f, const HistoryRepa& hr, const SizeList& kk)
    {
	return setVarsHistoryRepasReduce_u(f, kk.size(), kk.data(), hr);
    };
    auto hrpr = historyRepasRed;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui;
    auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
    auto roller = parametersRollerMaximumRollExcludedSelfRepa_i;
    auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;

    auto t0 = clk::now();
    std::cout << ">>> layerer" << std::endl;
    auto& llu = ur.listVarSizePair;
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    auto z = (double)hr.size;
    auto zr = (double)hrs.size;
    auto vf = std::make_shared<Variable>(f);
    auto fr = std::make_unique<FudRepa>();
    auto mm = std::make_unique<DoubleSizeListPairList>();
    auto hr1 = frmul(hr, *fr);
    auto hrs1 = frmul(hrs, *fr);
    std::size_t l = 1;
    while (l <= lmax)
    {
	auto t4 = clk::now();
	std::cout << ">>> layer\tfud: " << f << "\tlayer: " << l << std::endl;
	auto vl = std::make_shared<Variable>(l);
	auto vfl = std::make_shared<Variable>(vf,vl);
	std::size_t b = 1;
	auto pr1 = hrpr(*hr1);
	auto prs1 = hrpr(*hrs1);
	auto t1 = clk::now();
	std::cout << ">>> tupler" << std::endl;
	std::cout << "substrate cardinality: " << vv.size() << std::endl;
	std::cout << "fud cardinality: " << fudRepasSize(*fr) << std::endl;
	auto tt2 = tupler(xmax, omax, bmax, mmax, vv, *fr, *hr1, *pr1, *hrs1, *prs1);
	auto& x2 = std::get<0>(tt2);
	auto s2 = std::get<1>(tt2);
	std::cout << "tuple cardinality: " << x2->size() << std::endl;
	auto t2 = clk::now();
	std::cout << "tupler\tsearched: " << s2 << "\trate: " << s2 / ((sec)(t2 - t1)).count() << std::endl;
	std::cout << "<<< tupler " << ((sec)(t2 - t1)).count() << "s" << std::endl;
	if (!x2->size())
	    break;
	std::vector<SizeSet> qq;
	qq.reserve(x2->size() * mmax * pmax);
	TransformRepaPtrList ll;
	ll.reserve(x2->size() * mmax * pmax);
	std::map<std::string, double> time;
	std::map<std::string, std::size_t> steps;
	auto mark = clk::now();
	for (auto& kk : *x2)
	{
	    auto ar = hrred(1.0, *hr1, kk);
	    auto ars = hrred(z/zr, *hrs1, kk);
	    double y1 = ar->facLn() - ars->facLn();
	    mark = clk::now();
	    auto tt3 = rrvqqy(mmax, umax, pmax, *ar, *ars, z, y1);
	    time["parter"] += ((sec)(clk::now() - mark)).count();
	    steps["parter"] += std::get<1>(tt3);
	    auto& x3 = std::get<0>(tt3);
	    for (auto& nn : *x3)
	    {
		mark = clk::now();
		auto tt4 = roller(nn, *ar, *ars, z);
		time["roller"] += ((sec)(clk::now() - mark)).count();
		steps["roller"] += std::get<1>(tt4);
		auto& x4 = std::get<0>(tt4);
		for (auto& pp : *x4)
		{
		    auto m = nn.size();
		    for (std::size_t i = 0; i < m; i++)
		    {
			auto& cc = nn[i];
			auto e = cc.size();
			SizeList jj;
			jj.reserve(e);
			for (std::size_t j = 0; j < e; j++)
			    jj.push_back(kk[cc[j]]);
			SizeSet jj1(jj.begin(), jj.end());
			bool dup = false;
			for (auto& jj2 : qq)
			    if (jj2 == jj1)
			    {
				dup = true;
				break;
			    }
			if (!dup)
			{
			    qq.push_back(jj1);
			    auto tr = std::make_shared<TransformRepa>();
			    tr->dimension = e;
			    tr->vectorVar = new std::size_t[e];
			    auto ww = tr->vectorVar;
			    tr->shape = new std::size_t[e];
			    auto sh = tr->shape;
			    for (std::size_t j = 0; j < e; j++)
			    {
				auto& v = jj[j];
				ww[j] = v;
				sh[j] = llu[v].second;
			    }
			    auto& rr1 = pp[i];
			    auto sz = rr1.size();
			    std::size_t s = 1;
			    tr->arr = new unsigned char[sz];
			    auto rr = tr->arr;
			    for (std::size_t j = 0; j < sz; j++)
			    {
				auto u = rr1[j];
				rr[j] = u;
				if (u > s)
				    s = u;
			    }
			    s++;
			    if (s > ucmax + 1)
				throw std::out_of_range("parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u");
			    tr->valency = s;
			    auto vb = std::make_shared<Variable>(b);
			    auto vflb = std::make_shared<Variable>(vfl, vb);
			    llu.push_back(VarSizePair(*vflb,s));
			    auto w = llu.size()-1;
			    tr->derived = w;
			    ll.push_back(tr);
			    b++;
			}
		    }
		}
	    }
	}
	l++;
    }
    std::cout << "<<< layerer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;

    return std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>>(std::move(fr), std::move(mm));
}



