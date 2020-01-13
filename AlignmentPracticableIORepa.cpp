#include "AlignmentRandomRepa.h"
#include "AlignmentApprox.h"
#include "AlignmentPracticableIORepa.h"
#include <chrono>
#include <ctime>
#include <cmath>
#include <thread>

using namespace Alignment;

typedef std::chrono::duration<double> sec;
typedef std::chrono::high_resolution_clock clk;

const double repaRounding = 1e-6;

// parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> HistoryRepa-> Integer ->
//   IO (SystemRepa, FudRepa, [(Double, [VariableRepa]])
std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> Alignment::parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, const SizeList& vv, const HistoryRepa& hr, const HistoryRepa& hrs, std::size_t f, SystemRepa& ur)
{
    auto hrred = [](double f, const HistoryRepa& hr, const SizeList& kk)
    {
	return setVarsHistoryRepasReduce_u(f, kk.size(), kk.data(), hr);
    };
    auto hrpr = historyRepasRed;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui;
    auto parter = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
    auto roller = histogramRepaVecsRollMax;
    auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;

    auto t0 = clk::now();
    std::cout << ">>> layerer" << std::endl;
    auto& llu = ur.listVarSizePair;
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    auto z = (double)hr.size;
    auto zr = (double)hrs.size;
    auto vd = std::make_shared<Variable>(0);
    auto vf = std::make_shared<Variable>((int)f);
    auto vdf = std::make_shared<Variable>(vd, vf);
    auto fr = std::make_unique<FudRepa>();
    auto mm = std::make_unique<DoubleSizeListPairList>();
    auto hr1 = frmul(hr, *fr);
    auto hrs1 = frmul(hrs, *fr);
    auto pr1 = hrpr(*hr1);
    auto prs1 = hrpr(*hrs1);
    std::size_t l = 1;
    bool layering = true;
    while (layering && l <= lmax)
    {
	std::map<std::string, double> time;
	std::map<std::string, std::size_t> steps;
	auto t1 = clk::now();
	layering = false;
	std::cout << ">>> layer\tfud: " << f << "\tlayer: " << l << std::endl;
	auto vl = std::make_shared<Variable>((int)l);
	auto vfl = std::make_shared<Variable>(vdf,vl);
	std::size_t b = 1;
	std::cout << "substrate cardinality: " << vv.size() << std::endl;
	std::cout << "fud cardinality: " << fudRepasSize(*fr) << std::endl;
	auto mark = clk::now();
	auto tt2 = tupler(xmax, omax, bmax, mmax, vv, *fr, *hr1, *pr1, *hrs1, *prs1);
	time["tupler"] += ((sec)(clk::now() - mark)).count();
	steps["tupler"] += std::get<1>(tt2);
	std::cout << "tupler\tsearched: " << steps["tupler"] << "\trate: " << ((double)steps["tupler"]) / time["tupler"] << std::endl;
	std::cout << "tupler " << time["tupler"] << "s" << std::endl;
	auto& x2 = std::get<0>(tt2);
	std::cout << "tuple cardinality: " << x2->size() << std::endl;
	if (!x2->size())
	    break;
	std::vector<SizeSet> qq;
	qq.reserve(x2->size() * mmax * pmax);
	FudRepa gr;
	gr.layers.push_back(TransformRepaPtrList());
	auto& ll = gr.layers.back();
	ll.reserve(x2->size() * mmax * pmax);
	double ymax = 0.0;
	for (auto& kk : *x2)
	{
	    auto ar = hrred(1.0, *hr1, kk);
	    auto ars = hrred(z/zr, *hrs1, kk);
	    double y1 = ar->facLn() - ars->facLn();
	    if (!ll.size() || y1 > ymax)
		ymax = y1;
	    mark = clk::now();
	    auto tt3 = parter(mmax, umax, pmax, *ar, *ars, z, y1);
	    time["parter"] += ((sec)(clk::now() - mark)).count();
	    steps["parter"] += std::get<1>(tt3);
	    auto& x3 = std::get<0>(tt3);
	    for (auto& nn : *x3)
	    {
		mark = clk::now();
		auto m = nn.size();
		auto tt4 = roller(nn, *ar, *ars, z);
		time["roller"] += ((sec)(clk::now() - mark)).count();
		steps["roller"] += std::get<1>(tt4);
		auto& x4 = std::get<0>(tt4);
		if (x4->size() != m)
		    continue;
		bool any = false;
		for (auto& rr1 : *x4)
		{
		    auto sz = rr1.size();
		    std::size_t s = 1;
		    for (std::size_t j = 0; j < sz; j++)
		    {
			auto u = rr1[j];
			if (u > s)
			    s = u;
		    }
		    if (s < sz - 1)
		    {
			any = true;
			break;
		    }
		}
		if (!any)
		    continue;
		for (std::size_t i = 0; i < m; i++)
		{
		    auto& cc = nn[i];
		    auto e = cc.size();
		    SizeList jj;
		    jj.reserve(e);
		    for (std::size_t j = 0; j < e; j++)
			jj.push_back(kk[cc[j]]);
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
		    auto& rr1 = (*x4)[i];
		    auto sz = rr1.size();
		    std::size_t s = 0;
		    tr->arr = new unsigned char[sz];
		    auto rr = tr->arr;
		    for (std::size_t j = 0; j < sz; j++)
		    {
			auto u = rr1[j];
			rr[j] = (unsigned char)u;
			if (u > s)
			    s = u;
		    }
		    s++;
		    if (s == 1)
			continue;
		    if (e == 1 && s == sz)
			continue;
		    if (s > ucmax + 1)
			throw std::out_of_range("parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u");
		    SizeSet jj1(jj.begin(), jj.end());
		    bool dup = false;
		    for (auto& jj2 : qq)
			if (jj2 == jj1)
			{
			    dup = true;
			    break;
			}
		    if (dup)
			continue;
		    qq.push_back(jj1);
		    tr->valency = s;
		    auto vb = std::make_shared<Variable>((int)b);
		    auto vflb = std::make_shared<Variable>(vfl, vb);
		    llu.push_back(VarSizePair(vflb,s));
		    auto w = llu.size() - 1;
		    tr->derived = w;
		    ll.push_back(tr);
		    b++;
		}
	    }
	}
	std::cout << "max tuple algn: " << ymax << std::endl;
	std::cout << "layer cardinality: " << ll.size() << std::endl;
	std::cout << "parter\tsearched: " << steps["parter"] << "\trate: " << ((double)steps["parter"]) / time["parter"] << std::endl;
	std::cout << "parter " << time["parter"] << "s" << std::endl;
	std::cout << "roller\tsearched: " << steps["roller"] << "\trate: " << ((double)steps["roller"]) / time["roller"] << std::endl;
	std::cout << "roller " << time["roller"] << "s" << std::endl;
	if (ll.size())
	{
	    hr1 = frmul(*hr1, gr);
	    hrs1 = frmul(*hrs1, gr);
	    pr1 = hrpr(*hr1);
	    prs1 = hrpr(*hrs1);
	    fr->layers.push_back(ll);
	    mark = clk::now();
	    auto tt5 = deriveder(wmax, omax, *fr, *hr1, *pr1, *hrs1, *prs1);
	    time["dervarser"] += ((sec)(clk::now() - mark)).count();
	    steps["dervarser"] += std::get<1>(tt5);
	    auto& mm1 = std::get<0>(tt5);
	    if (mm1->size())
		std::cout << "der vars algn density: " << mm1->back().first << std::endl;
	    else
		std::cout << "no der vars sets" << std::endl;
	    std::cout << "dervarser\tsearched: " << steps["dervarser"] << "\trate: " << ((double)steps["dervarser"]) / time["dervarser"] << std::endl;
	    std::cout << "dervarser " << time["dervarser"] << "s" << std::endl;
	    if (mm1->size() && (!mm->size() || mm1->back().first > mm->back().first + repaRounding))
	    {
		mm = std::move(mm1);
		layering = true;
	    }
	    else
		fr->layers.pop_back();
	}
	time["application"] = ((sec)(clk::now() - t1)).count() - time["tupler"] - time["parter"] - time["roller"] - time["dervarser"];
	std::cout << "application " << time["application"] << "s" << std::endl;
	std::cout << "<<< layer " << ((sec)(clk::now() - t1)).count() << "s" << std::endl;
	l++;
    }
    std::cout << "<<< layerer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>>(std::move(fr), std::move(mm));
}

// parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_up ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> HistoryRepa-> Integer ->
//   IO (SystemRepa, FudRepa, [(Double, [VariableRepa]])
std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>> Alignment::parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_up(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t tint, const SizeList& vv, const HistoryRepa& hr, const HistoryRepa& hrs, std::size_t f, SystemRepa& ur)
{
    auto hrred = [](double f, const HistoryRepa& hr, const SizeList& kk)
    {
	return setVarsHistoryRepasReduce_u(f, kk.size(), kk.data(), hr);
    };
    auto hrpr = historyRepasRed;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto tupler = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_uip;
    auto parter = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;
    auto roller = histogramRepaVecsRollMax;
    auto deriveder = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui;

    auto t0 = clk::now();
    std::cout << ">>> layerer" << std::endl;
    auto& llu = ur.listVarSizePair;
    unsigned char ucmax = std::numeric_limits<unsigned char>::max();
    auto z = (double)hr.size;
    auto zr = (double)hrs.size;
    auto vd = std::make_shared<Variable>(0);
    auto vf = std::make_shared<Variable>((int)f);
    auto vdf = std::make_shared<Variable>(vd, vf);
    auto fr = std::make_unique<FudRepa>();
    auto mm = std::make_unique<DoubleSizeListPairList>();
    auto hr1 = frmul(hr, *fr);
    auto hrs1 = frmul(hrs, *fr);
    auto pr1 = hrpr(*hr1);
    auto prs1 = hrpr(*hrs1);
    std::size_t l = 1;
    bool layering = true;
    while (layering && l <= lmax)
    {
	std::map<std::string, double> time;
	std::map<std::string, std::size_t> steps;
	auto t1 = clk::now();
	layering = false;
	std::cout << ">>> layer\tfud: " << f << "\tlayer: " << l << std::endl;
	auto vl = std::make_shared<Variable>((int)l);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	std::size_t b = 1;
	std::cout << "substrate cardinality: " << vv.size() << std::endl;
	std::cout << "fud cardinality: " << fudRepasSize(*fr) << std::endl;
	auto mark = clk::now();
	auto tt2 = tupler(xmax, omax, bmax, mmax, tint, vv, *fr, *hr1, *pr1, *hrs1, *prs1);
	time["tupler"] += ((sec)(clk::now() - mark)).count();
	steps["tupler"] += std::get<1>(tt2);
	std::cout << "tupler\tsearched: " << steps["tupler"] << "\trate: " << ((double)steps["tupler"]) / time["tupler"] << std::endl;
	std::cout << "tupler " << time["tupler"] << "s" << std::endl;
	auto& x2 = std::get<0>(tt2);
	std::cout << "tuple cardinality: " << x2->size() << std::endl;
	if (!x2->size())
	    break;
	std::vector<SizeSet> qq;
	qq.reserve(x2->size() * mmax * pmax);
	FudRepa gr;
	gr.layers.push_back(TransformRepaPtrList());
	auto& ll = gr.layers.back();
	ll.reserve(x2->size() * mmax * pmax);
	double ymax = 0.0;
	for (auto& kk : *x2)
	{
	    auto ar = hrred(1.0, *hr1, kk);
	    auto ars = hrred(z / zr, *hrs1, kk);
	    double y1 = ar->facLn() - ars->facLn();
	    if (!ll.size() || y1 > ymax)
		ymax = y1;
	    mark = clk::now();
	    auto tt3 = parter(mmax, umax, pmax, *ar, *ars, z, y1);
	    time["parter"] += ((sec)(clk::now() - mark)).count();
	    steps["parter"] += std::get<1>(tt3);
	    auto& x3 = std::get<0>(tt3);
	    for (auto& nn : *x3)
	    {
		mark = clk::now();
		auto m = nn.size();
		auto tt4 = roller(nn, *ar, *ars, z);
		time["roller"] += ((sec)(clk::now() - mark)).count();
		steps["roller"] += std::get<1>(tt4);
		auto& x4 = std::get<0>(tt4);
		if (x4->size() != m)
		    continue;
		bool any = false;
		for (auto& rr1 : *x4)
		{
		    auto sz = rr1.size();
		    std::size_t s = 1;
		    for (std::size_t j = 0; j < sz; j++)
		    {
			auto u = rr1[j];
			if (u > s)
			    s = u;
		    }
		    if (s < sz - 1)
		    {
			any = true;
			break;
		    }
		}
		if (!any)
		    continue;
		for (std::size_t i = 0; i < m; i++)
		{
		    auto& cc = nn[i];
		    auto e = cc.size();
		    SizeList jj;
		    jj.reserve(e);
		    for (std::size_t j = 0; j < e; j++)
			jj.push_back(kk[cc[j]]);
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
		    auto& rr1 = (*x4)[i];
		    auto sz = rr1.size();
		    std::size_t s = 0;
		    tr->arr = new unsigned char[sz];
		    auto rr = tr->arr;
		    for (std::size_t j = 0; j < sz; j++)
		    {
			auto u = rr1[j];
			rr[j] = (unsigned char)u;
			if (u > s)
			    s = u;
		    }
		    s++;
		    if (s == 1)
			continue;
		    if (e == 1 && s == sz)
			continue;
		    if (s > ucmax + 1)
			throw std::out_of_range("parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u");
		    SizeSet jj1(jj.begin(), jj.end());
		    bool dup = false;
		    for (auto& jj2 : qq)
			if (jj2 == jj1)
			{
			    dup = true;
			    break;
			}
		    if (dup)
			continue;
		    qq.push_back(jj1);
		    tr->valency = s;
		    auto vb = std::make_shared<Variable>((int)b);
		    auto vflb = std::make_shared<Variable>(vfl, vb);
		    llu.push_back(VarSizePair(vflb, s));
		    auto w = llu.size() - 1;
		    tr->derived = w;
		    ll.push_back(tr);
		    b++;
		}
	    }
	}
	std::cout << "max tuple algn: " << ymax << std::endl;
	std::cout << "layer cardinality: " << ll.size() << std::endl;
	std::cout << "parter\tsearched: " << steps["parter"] << "\trate: " << ((double)steps["parter"]) / time["parter"] << std::endl;
	std::cout << "parter " << time["parter"] << "s" << std::endl;
	std::cout << "roller\tsearched: " << steps["roller"] << "\trate: " << ((double)steps["roller"]) / time["roller"] << std::endl;
	std::cout << "roller " << time["roller"] << "s" << std::endl;
	if (ll.size())
	{
	    hr1 = frmul(*hr1, gr);
	    hrs1 = frmul(*hrs1, gr);
	    pr1 = hrpr(*hr1);
	    prs1 = hrpr(*hrs1);
	    fr->layers.push_back(ll);
	    mark = clk::now();
	    auto tt5 = deriveder(wmax, omax, *fr, *hr1, *pr1, *hrs1, *prs1);
	    time["dervarser"] += ((sec)(clk::now() - mark)).count();
	    steps["dervarser"] += std::get<1>(tt5);
	    auto& mm1 = std::get<0>(tt5);
	    if (mm1->size())
		std::cout << "der vars algn density: " << mm1->back().first << std::endl;
	    else
		std::cout << "no der vars sets" << std::endl;
	    std::cout << "dervarser\tsearched: " << steps["dervarser"] << "\trate: " << ((double)steps["dervarser"]) / time["dervarser"] << std::endl;
	    std::cout << "dervarser " << time["dervarser"] << "s" << std::endl;
	    if (mm1->size() && (!mm->size() || mm1->back().first > mm->back().first + repaRounding))
	    {
		mm = std::move(mm1);
		layering = true;
	    }
	    else
		fr->layers.pop_back();
	}
	time["application"] = ((sec)(clk::now() - t1)).count() - time["tupler"] - time["parter"] - time["roller"] - time["dervarser"];
	std::cout << "application " << time["application"] << "s" << std::endl;
	std::cout << "<<< layer " << ((sec)(clk::now() - t1)).count() << "s" << std::endl;
	l++;
    }
    std::cout << "<<< layerer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return std::tuple<std::unique_ptr<FudRepa>, std::unique_ptr<DoubleSizeListPairList>>(std::move(fr), std::move(mm));
}

// parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t seed, const SizeList& vv, const HistoryRepa& hr, SystemRepa& ur)
{
    return parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa(wmax, lmax, xmax, 0, omax, bmax, mmax, umax, pmax, fmax, mult, 0, seed, vv, FudRepa(), hr, 0, ur);
}

// parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsFudRepasHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t seed, const SizeList& vv, const FudRepa& er, const HistoryRepa& hr, SystemRepa& ur)
{
    return parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa(wmax, lmax, xmax, 0, omax, bmax, mmax, umax, pmax, fmax, mult, 0, seed, vv, er, hr, 0, ur);
}

// parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa(std::size_t wmax, std::size_t lmax, std::size_t xmax, double znnmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t smin, std::size_t seed, const SizeList& vv, const FudRepa& er, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto hrred = setVarsHistoryRepasReduce_u;
    auto hrhrred = setVarsHistoryRepasHistoryRepaReduced_u;
    auto hrpr = setVarsHistoryRepasRed_u;
    auto prents = histogramRepaRedsListEntropy;
    auto hrconcat = vectorHistoryRepasConcat_u;
    auto hrshuffle = historyRepasShuffle_u;
    auto llfr = setVariablesListTransformRepasFudRepa_u;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto frdep = fudRepasSetVarsDepends;
    auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u;

    auto t0 = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    SizeUSet vv00(hr0.dimension);
    for (std::size_t i = 0; i < hr0.dimension; i++)
	vv00.insert(hr0.vectorVar[i]);
    SizeUSet vv1(vv.begin(), vv.end());
    auto& llu = ur.listVarSizePair;
    auto z = hr0.size;
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax*(lmax + 1));
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    auto mark = clk::now();
    std::cout << ">>> applier " << std::endl;
    std::shared_ptr<HistoryRepa> hr2 = hrhrred(vv.size(), vv.data(), *frmul(hr0, er));
    std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
    mark = clk::now();
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	std::size_t mult1 = smin/z + 1 > mult ? smin/z + 1 : mult;
	HistoryRepaPtrList qq;
	qq.reserve(mult1);
	for (std::size_t i = 1; i <= mult1; i++)
	    qq.push_back(hrshuffle(hr0, (unsigned int)(seed + i*z)));
	auto hrs0 = hrconcat(qq);
	qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> substrater " << std::endl;
	auto hr = hr2;
	SizeList vv0;
	std::unique_ptr<FudRepa> er0;
	auto nmax = (std::size_t)std::sqrt(znnmax / (double)(z + mult1*z));
	if (nmax > bmax)
	{
	    auto ee = prents(*hrpr(vv.size(), vv.data(), *hr));
	    auto m = ee->size();
	    if (m > nmax)
	    {
		std::sort(ee->begin(), ee->end());
		vv0.reserve(nmax);
		for (std::size_t i = m - nmax; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    else
	    {
		vv0.reserve(m);
		for (std::size_t i = 0; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    if (vv0.size() < vv.size())
	    {
		SizeUSet vv01(vv0.begin(), vv0.end());
		er0 = llfr(vv00, *frdep(er, vv01));
	    }
	}
	time["substrater"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< substrater " << time["substrater"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	std::unique_ptr<HistoryRepa> hrs;
	if (er0)
	{
	    hr = frmul(hr0, *er0);
	    hrs = frmul(*hrs0, *er0);
	}
	else
	{
	    hrs = hrhrred(vv.size(), vv.data(), *frmul(*hrs0, er));
	}
	hrs0.reset();
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, nmax > bmax ? vv0 : vv, *hr, *hrs, f, ur);
	    fr = std::move(std::get<0>(t));
	    mm = std::move(std::get<1>(t));
	}
	catch (const std::out_of_range& e)
	{
	    std::cout << "out of range exception: " << e.what() << std::endl;
	    fr.reset();
	    mm.reset();
	}
	if (!mm || !mm->size())
	{
	    std::cout << "no fud" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto& a = mm->back().first;
	auto& kk = mm->back().second;
	auto m = kk.size();
	if (m < 2 || a <= repaRounding)
	{
	    std::cout << "no algn" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "derived cardinality: " << m << std::endl;
	std::cout << "derived algn density: " << a << std::endl;
	std::cout << "derived impl bi-valency percent: " << 100.0 * (std::exp(a / (double)z / (double)(m - 1)) - 1.0) << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(*hr, *gr));
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = 1;
	auto skk = ar->shape;
	auto rr0 = ar->arr;
	for (std::size_t i = 0; i < m; i++)
	    sz *= skk[i];
	sl.reserve(sz);
	ll.reserve(sz);
	bool remainder = false;
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    if (rr0[i] <= 0.0)
	    {
		remainder = true;
		continue;
	    }
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m;
	    tr->vectorVar = new std::size_t[m];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m];
	    auto sh = tr->shape;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j] = kk[j];
		sh[j] = skk[j];
	    }
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	if (remainder)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m;
	    tr->vectorVar = new std::size_t[m];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m];
	    auto sh = tr->shape;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j] = kk[j];
		sh[j] = skk[j];
	    }
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = rr0[j] <= 0.0 ? 1 : 0;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.insert(dr->fud->layers.end(), gr->layers.begin(), gr->layers.end());
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	auto hr = frmul(*hr2, *dr->fud);
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto n = hr->dimension;
	auto& mvv = hr->mapVarInt();
	auto rr = hr->arr;
	auto nn = treesLeafNodes(*dr->slices);
	SizeSizePairList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t a = 0;
	    auto pk = mvv[v];
	    if (hr->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[j*n + pk];
		    if (u)
			a++;
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[pk*z + j];
		    if (u)
			a++;
		}
	    if (a > 1)
		zs.push_back(SizeSizePair(a, v));
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = zs.back().first;
	auto v = zs.back().second;
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v].first << std::endl;
	SizeList ev;
	ev.reserve(z2);
	{
	    auto pk = mvv[v];
	    if (hr->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[j*n + pk];
		    if (u)
			ev.push_back(j);
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[pk*z + j];
		    if (u)
			ev.push_back(j);
		}
	}
	hr = hrsel(ev.size(), ev.data(), *hr2);
	auto hr1 = hrsel(ev.size(), ev.data(), hr0);
	ev.clear();
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	std::size_t mult1 = smin / z2 + 1 > mult ? smin / z2 + 1 : mult;
	HistoryRepaPtrList qq;
	qq.reserve(mult1);
	for (std::size_t i = 1; i <= mult1; i++)
	    qq.push_back(hrshuffle(*hr1, (unsigned int)(seed + i*z2)));
	auto hrs1 = hrconcat(qq);
	qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> substrater " << std::endl;
	SizeList vv0;
	std::unique_ptr<FudRepa> er0;
	auto nmax = (std::size_t)std::sqrt(znnmax / (double)(z2 + mult1*z2));
	if (nmax > bmax)
	{
	    auto ee = prents(*hrpr(vv.size(), vv.data(), *hr));
	    auto m = ee->size();
	    if (m > nmax)
	    {
		std::sort(ee->begin(), ee->end());
		vv0.reserve(nmax);
		for (std::size_t i = m - nmax; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    else
	    {
		vv0.reserve(m);
		for (std::size_t i = 0; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    if (vv0.size() < vv.size())
	    {
		SizeUSet vv01(vv0.begin(), vv0.end());
		er0 = llfr(vv00, *frdep(er, vv01));
	    }
	}
	time["substrater"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< substrater " << time["substrater"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	std::unique_ptr<HistoryRepa> hrs;
	if (er0)
	{
	    hr = frmul(*hr1, *er0);
	    hrs = frmul(*hrs1, *er0);
	}
	else
	{
	    hrs = hrhrred(vv.size(), vv.data(), *frmul(*hrs1, er));
	}
	hr1.reset();
	hrs1.reset();
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, nmax > bmax ? vv0 : vv, *hr, *hrs, f + 1, ur);
	    fr = std::move(std::get<0>(t));
	    mm = std::move(std::get<1>(t));
	}
	catch (const std::out_of_range& e)
	{
	    std::cout << "out of range exception: " << e.what() << std::endl;
	    fr.reset();
	    mm.reset();
	}
	if (!mm || !mm->size())
	{
	    ig.insert(v);
	    std::cout << "no fud" << std::endl;
	    continue;
	}
	auto& a = mm->back().first;
	auto& kk = mm->back().second;
	auto m = kk.size();
	if (m < 2 || a <= repaRounding)
	{
	    ig.insert(v);
	    std::cout << "no algn" << std::endl;
	    continue;
	}
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "derived cardinality: " << m << std::endl;
	std::cout << "derived algn density: " << a << std::endl;
	std::cout << "derived impl bi-valency percent: " << 100.0 * (std::exp(a / (double)z2 / (double)(m - 1)) - 1.0) << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(*hr, *gr));
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = 1;
	auto skk = ar->shape;
	auto rr0 = ar->arr;
	for (std::size_t i = 0; i < m; i++)
	    sz *= skk[i];
	sl.reserve(sz);
	ll.reserve(sz);
	bool remainder = false;
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    if (rr0[i] <= 0.0)
	    {
		remainder = true;
		continue;
	    }
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m + 1;
	    tr->vectorVar = new std::size_t[m + 1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m + 1];
	    auto sh = tr->shape;
	    ww[0] = v;
	    sh[0] = 2;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j + 1] = kk[j];
		sh[j + 1] = skk[j];
	    }
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	if (remainder)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m + 1;
	    tr->vectorVar = new std::size_t[m + 1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m + 1];
	    auto sh = tr->shape;
	    ww[0] = v;
	    sh[0] = 2;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j + 1] = kk[j];
		sh[j + 1] = skk[j];
	    }
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = j >= sz && rr0[j - sz] <= 0.0 ? 1 : 0;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.insert(dr->fud->layers.end(), gr->layers.begin(), gr->layers.end());
	dr->fud->layers.push_back(ll);
	for (auto& p : *nn)
	    if (p.first == v)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}

// parametersSystemsHistoryRepasApplicationerMaxRollByMExcludedSelfHighestFmaxIORepa_p ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer ->
//   Integer -> Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsFudRepasHistoryRepasApplicationerSubstrateEntropyMaxRollByMExcludedSelfHighestFmaxIORepa_p(std::size_t wmax, std::size_t lmax, std::size_t xmax, double znnmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t umax, std::size_t pmax, std::size_t fmax, std::size_t mult, std::size_t smin, std::size_t seed, std::size_t tint, const SizeList& vv, const FudRepa& er, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto hrred = setVarsHistoryRepasReduce_u;
    auto hrhrred = setVarsHistoryRepasHistoryRepaReduced_u;
    auto hrpr = setVarsHistoryRepasRed_u;
    auto prents = histogramRepaRedsListEntropy;
    auto hrconcat = vectorHistoryRepasConcat_u;
    auto hrshuffle = historyRepasShuffle_u;
    auto llfr = setVariablesListTransformRepasFudRepa_u;
    auto frmul = historyRepasFudRepasMultiply_up;
    auto frdep = fudRepasSetVarsDepends;
    auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_up;

    auto t0 = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    SizeUSet vv00(hr0.dimension);
    for (std::size_t i = 0; i < hr0.dimension; i++)
	vv00.insert(hr0.vectorVar[i]);
    SizeUSet vv1(vv.begin(), vv.end());
    auto& llu = ur.listVarSizePair;
    auto z = hr0.size;
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax*(lmax + 1));
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    auto mark = clk::now();
    std::cout << ">>> applier " << std::endl;
    std::shared_ptr<HistoryRepa> hr2 = hrhrred(vv.size(), vv.data(), *frmul(tint, hr0, er));
    std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
    mark = clk::now();
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	std::size_t mult1 = smin / z + 1 > mult ? smin / z + 1 : mult;
	HistoryRepaPtrList qq;
	qq.reserve(mult1);
	for (std::size_t i = 1; i <= mult1; i++)
	    qq.push_back(hrshuffle(hr0, (unsigned int)(seed + i*z)));
	auto hrs0 = hrconcat(qq);
	qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> substrater " << std::endl;
	auto hr = hr2;
	SizeList vv0;
	std::unique_ptr<FudRepa> er0;
	auto nmax = (std::size_t)std::sqrt(znnmax / (double)(z + mult1*z));
	if (nmax > bmax)
	{
	    auto ee = prents(*hrpr(vv.size(), vv.data(), *hr));
	    auto m = ee->size();
	    if (m > nmax)
	    {
		std::sort(ee->begin(), ee->end());
		vv0.reserve(nmax);
		for (std::size_t i = m - nmax; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    else
	    {
		vv0.reserve(m);
		for (std::size_t i = 0; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    if (vv0.size() < vv.size())
	    {
		SizeUSet vv01(vv0.begin(), vv0.end());
		er0 = llfr(vv00, *frdep(er, vv01));
	    }
	}
	time["substrater"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< substrater " << time["substrater"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	std::unique_ptr<HistoryRepa> hrs;
	if (er0)
	{
	    hr = frmul(tint, hr0, *er0);
	    hrs = frmul(tint, *hrs0, *er0);
	}
	else
	{
	    hrs = hrhrred(vv.size(), vv.data(), *frmul(tint, *hrs0, er));
	}
	hrs0.reset();
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, tint, nmax > bmax ? vv0 : vv, *hr, *hrs, f, ur);
	    fr = std::move(std::get<0>(t));
	    mm = std::move(std::get<1>(t));
	}
	catch (const std::out_of_range& e)
	{
	    std::cout << "out of range exception: " << e.what() << std::endl;
	    fr.reset();
	    mm.reset();
	}
	if (!mm || !mm->size())
	{
	    std::cout << "no fud" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto& a = mm->back().first;
	auto& kk = mm->back().second;
	auto m = kk.size();
	if (m < 2 || a <= repaRounding)
	{
	    std::cout << "no algn" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "derived cardinality: " << m << std::endl;
	std::cout << "derived algn density: " << a << std::endl;
	std::cout << "derived impl bi-valency percent: " << 100.0 * (std::exp(a / (double)z / (double)(m - 1)) - 1.0) << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(tint, *hr, *gr));
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = 1;
	auto skk = ar->shape;
	auto rr0 = ar->arr;
	for (std::size_t i = 0; i < m; i++)
	    sz *= skk[i];
	sl.reserve(sz);
	ll.reserve(sz);
	bool remainder = false;
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    if (rr0[i] <= 0.0)
	    {
		remainder = true;
		continue;
	    }
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m;
	    tr->vectorVar = new std::size_t[m];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m];
	    auto sh = tr->shape;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j] = kk[j];
		sh[j] = skk[j];
	    }
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	if (remainder)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m;
	    tr->vectorVar = new std::size_t[m];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m];
	    auto sh = tr->shape;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j] = kk[j];
		sh[j] = skk[j];
	    }
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = rr0[j] <= 0.0 ? 1 : 0;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.insert(dr->fud->layers.end(), gr->layers.begin(), gr->layers.end());
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	auto hr = frmul(tint, *hr2, *dr->fud);
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto n = hr->dimension;
	auto& mvv = hr->mapVarInt();
	auto rr = hr->arr;
	auto nn = treesLeafNodes(*dr->slices);
	SizeSizePairList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t a = 0;
	    auto pk = mvv[v];
	    if (hr->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[j*n + pk];
		    if (u)
			a++;
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[pk*z + j];
		    if (u)
			a++;
		}
	    if (a > 1)
		zs.push_back(SizeSizePair(a, v));
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = zs.back().first;
	auto v = zs.back().second;
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v].first << std::endl;
	SizeList ev;
	ev.reserve(z2);
	{
	    auto pk = mvv[v];
	    if (hr->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[j*n + pk];
		    if (u)
			ev.push_back(j);
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[pk*z + j];
		    if (u)
			ev.push_back(j);
		}
	}
	hr = hrsel(ev.size(), ev.data(), *hr2);
	auto hr1 = hrsel(ev.size(), ev.data(), hr0);
	ev.clear();
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	std::size_t mult1 = smin / z2 + 1 > mult ? smin / z2 + 1 : mult;
	HistoryRepaPtrList qq;
	qq.reserve(mult1);
	for (std::size_t i = 1; i <= mult1; i++)
	    qq.push_back(hrshuffle(*hr1, (unsigned int)(seed + i*z2)));
	auto hrs1 = hrconcat(qq);
	qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> substrater " << std::endl;
	SizeList vv0;
	std::unique_ptr<FudRepa> er0;
	auto nmax = (std::size_t)std::sqrt(znnmax / (double)(z2 + mult1*z2));
	if (nmax > bmax)
	{
	    auto ee = prents(*hrpr(vv.size(), vv.data(), *hr));
	    auto m = ee->size();
	    if (m > nmax)
	    {
		std::sort(ee->begin(), ee->end());
		vv0.reserve(nmax);
		for (std::size_t i = m - nmax; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    else
	    {
		vv0.reserve(m);
		for (std::size_t i = 0; i < m; i++)
		    vv0.push_back((*ee)[i].second);
	    }
	    if (vv0.size() < vv.size())
	    {
		SizeUSet vv01(vv0.begin(), vv0.end());
		er0 = llfr(vv00, *frdep(er, vv01));
	    }
	}
	time["substrater"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< substrater " << time["substrater"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> applier " << std::endl;
	std::unique_ptr<HistoryRepa> hrs;
	if (er0)
	{
	    hr = frmul(tint, *hr1, *er0);
	    hrs = frmul(tint, *hrs1, *er0);
	}
	else
	{
	    hrs = hrhrred(vv.size(), vv.data(), *frmul(tint, *hrs1, er));
	}
	hr1.reset();
	hrs1.reset();
	time["applier"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, tint, nmax > bmax ? vv0 : vv, *hr, *hrs, f + 1, ur);
	    fr = std::move(std::get<0>(t));
	    mm = std::move(std::get<1>(t));
	}
	catch (const std::out_of_range& e)
	{
	    std::cout << "out of range exception: " << e.what() << std::endl;
	    fr.reset();
	    mm.reset();
	}
	if (!mm || !mm->size())
	{
	    ig.insert(v);
	    std::cout << "no fud" << std::endl;
	    continue;
	}
	auto& a = mm->back().first;
	auto& kk = mm->back().second;
	auto m = kk.size();
	if (m < 2 || a <= repaRounding)
	{
	    ig.insert(v);
	    std::cout << "no algn" << std::endl;
	    continue;
	}
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "derived cardinality: " << m << std::endl;
	std::cout << "derived algn density: " << a << std::endl;
	std::cout << "derived impl bi-valency percent: " << 100.0 * (std::exp(a / (double)z2 / (double)(m - 1)) - 1.0) << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(tint, *hr, *gr));
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = 1;
	auto skk = ar->shape;
	auto rr0 = ar->arr;
	for (std::size_t i = 0; i < m; i++)
	    sz *= skk[i];
	sl.reserve(sz);
	ll.reserve(sz);
	bool remainder = false;
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    if (rr0[i] <= 0.0)
	    {
		remainder = true;
		continue;
	    }
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m + 1;
	    tr->vectorVar = new std::size_t[m + 1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m + 1];
	    auto sh = tr->shape;
	    ww[0] = v;
	    sh[0] = 2;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j + 1] = kk[j];
		sh[j + 1] = skk[j];
	    }
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	if (remainder)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = m + 1;
	    tr->vectorVar = new std::size_t[m + 1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m + 1];
	    auto sh = tr->shape;
	    ww[0] = v;
	    sh[0] = 2;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j + 1] = kk[j];
		sh[j + 1] = skk[j];
	    }
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = j >= sz && rr0[j - sz] <= 0.0 ? 1 : 0;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	for (std::size_t i = 0; i < gr->layers.size(); i++)
	{
	    if (i < dr->fud->layers.size())
		dr->fud->layers[i].insert(dr->fud->layers[i].end(), gr->layers[i].begin(), gr->layers[i].end());
	    else
		dr->fud->layers.push_back(gr->layers[i]);
	}
	{
	    bool found = false;
	    std::size_t i = 0;
	    for (i = 0; !found && i < dr->fud->layers.size(); i++)
		for (auto tr : dr->fud->layers[i])
		    if (tr->derived == v)
		    {
			found = true;
			break;
		    }
	    i = std::max(i+1, gr->layers.size());
	    if (i < dr->fud->layers.size())
		dr->fud->layers[i].insert(dr->fud->layers[i].end(), ll.begin(), ll.end());
	    else
		dr->fud->layers.push_back(ll);
	}
	for (auto& p : *nn)
	    if (p.first == v)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}


typedef std::tuple<double, std::size_t, std::size_t> DoubleSizeSizeTuple;
typedef std::vector<DoubleSizeSizeTuple> DoubleSizeSizeTupleList;

// parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u ::
//  Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u(std::size_t fmax, const SizeList& vv, std::size_t l, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto frmul = historyRepasFudRepasMultiply_u;

    auto t0 = clk::now();
    auto mark = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    auto n = vv.size();
    auto& llu = ur.listVarSizePair;
    auto n0 = hr0.dimension;
    auto z = hr0.size;
    auto& mvv = hr0.mapVarInt();
    auto sh0 = hr0.shape;
    auto rr0 = hr0.arr;
    std::size_t ms = 2;
    std::size_t ms1 = 2;
    for (auto v : vv)
    {
	auto s = sh0[mvv[v]];
	if (s > ms)
	    ms = s;
	if (v != l && s > ms1)
	    ms1 = s;
    }
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax);
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    HistoryRepa hr1;
    hr1.evient = false;
    std::size_t m1 = n + 1;
    std::size_t n1 = m1 + fmax*ms1;
    hr1.dimension = n1;
    hr1.vectorVar = new std::size_t[n1];
    auto vv1 = hr1.vectorVar;
    hr1.shape = new std::size_t[n1];
    auto sh1 = hr1.shape;
    auto pp1 = new std::size_t[m1];
    vv1[0] = l;
    pp1[0] = mvv[l];
    auto ls = sh0[pp1[0]];
    sh1[0] = ls;
    for (std::size_t i = 0; i < n; i++)
    {
	auto v = vv[i];
	vv1[i+1] = v;
	auto p = mvv[v];
	pp1[i+1] = p;
	sh1[i+1] = sh0[p];
    }
    hr1.arr = new unsigned char[z*n1];
    auto rr1 = hr1.arr;
    if (hr0.evient)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn0 = j*n0;
	    for (std::size_t i = 0; i < m1; i++)
		rr1[i*z + j] = rr0[jn0 + pp1[i]];
	}
    else
	for (std::size_t i = 0; i < m1; i++)
	{
	    std::size_t iz0 = pp1[i]*z;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
		rr1[iz + j] = rr0[iz0 + j];
	}
    delete pp1;
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	double g = 1.0 / (double)z;
	auto ee0 = new std::size_t[ls];
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 0;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t w = rr1[j];
	    ee0[w]++;
	}
	double e0 = 0.0;
	for (std::size_t k = 0; k < ls; k++)
	{
	    double a = (double)ee0[k] * g;
	    if (a > 0.0)
		e0 -= a * log(a);
	}
	delete ee0;
	if (e0 <= repaRounding)
	{
	    std::cout << "no label entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	double e = e0;
	std::size_t i1 = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (std::size_t i = 1; i <= n; i++)
	{
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 0;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 0;
	    auto iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		std::size_t w = rr1[j];
		ee1[w*s + u]++;
		ee2[u]++;
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
	    {
		double a = (double)ee1[k] * g;
		if (a > 0.0)
		    e1 += a * log(a);
	    }
	    for (std::size_t k = 0; k < s; k++)
	    {
		double a = (double)ee2[k] * g;
		if (a > 0.0)
		    e2 += a * log(a);
	    }
	    if (vv1[i] != l && e > e2 - e1)
	    {
		e = e2 - e1;
		i1 = i;
	    }
	}
	delete ee2;
	delete ee1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    std::cout << "no conditional entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "sized entropy label : " << (double)z * e0 << std::endl;
	std::cout << "sized conditional entropy: " << (double)z * e << std::endl;
	std::cout << "sized entropy decrease: " << (double)z * (e0 - e) << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 1;
	    tr->vectorVar = new std::size_t[1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[1];
	    auto sh = tr->shape;
	    ww[0] = v1;
	    sh[0] = sz;
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto ee0 = new std::size_t[ls];
	auto nn = treesLeafNodes(*dr->slices);
	DoubleSizeSizeTupleList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t i = 1+n;
	    while (i < m1)
	    {
		if (vv1[i] == v)
		    break;
		i++;
	    }
	    if (i==m1)
		continue;
	    std::size_t a = 0;
	    for (std::size_t k = 0; k < ls; k++)
		ee0[k] = 0;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		if (u)
		{
		    a++;
		    std::size_t w = rr1[j];
		    ee0[w]++;
		}
	    }
	    if (a > 1)
	    {
		double g = 1.0 / (double)a;
		double e = 0.0;
		for (std::size_t k = 0; k < ls; k++)
		{
		    double b = (double)ee0[k] * g;
		    if (b > 0.0)
			e -= b * log(b);
		}
		if ((double)a * e > repaRounding)
		    zs.push_back(DoubleSizeSizeTuple((double)a * e, a, i));
	    }
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = std::get<1>(zs.back());
	auto i2 = std::get<2>(zs.back());
	auto v2 = vv1[i2];
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v2].first << std::endl;
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	double g = 1.0 / (double)z2;
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 0;
	auto i2z = i2*z;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t t = rr1[i2z + j];
	    if (t)
	    {
		std::size_t w = rr1[j];
		ee0[w]++;
	    }
	}
	double e0 = 0.0;
	for (std::size_t k = 0; k < ls; k++)
	{
	    double a = (double)ee0[k] * g;
	    if (a > 0.0)
		e0 -= a * log(a);
	}
	delete ee0;
	double e = e0;
	std::size_t i1 = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (std::size_t i = 1; i <= n; i++)
	{
	    auto iz = i*z;
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 0;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 0;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t t = rr1[i2z + j];
		if (t)
		{
		    std::size_t u = rr1[iz + j];
		    std::size_t w = rr1[j];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
	    {
		double a = (double)ee1[k] * g;
		if (a > 0.0)
		    e1 += a * log(a);
	    }
	    for (std::size_t k = 0; k < s; k++)
	    {
		double a = (double)ee2[k] * g;
		if (a > 0.0)
		    e2 += a * log(a);
	    }
	    if (vv1[i] != l && e > e2 - e1)
	    {
		e = e2 - e1;
		i1 = i;
	    }
	}
	delete ee2;
	delete ee1;
	if (e >= e0 - repaRounding)
	{
	    ig.insert(v2);
	    std::cout << "no conditional entropy" << std::endl;
	    continue;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "sized entropy label : " << (double)z2 * e0 << std::endl;
	std::cout << "sized conditional entropy: " << (double)z2 * e << std::endl;
	std::cout << "sized entropy decrease: " << (double)z2 * (e0 - e) << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 2;
	    tr->vectorVar = new std::size_t[2];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[2];
	    auto sh = tr->shape;
	    ww[0] = v2;
	    sh[0] = 2;
	    ww[1] = v1;
	    sh[1] = sz;
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i2z + j];
		k = sz*k + rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	for (auto& p : *nn)
	    if (p.first == v2)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}

// parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u_1 ::
//  Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerConditionalFmaxIORepa_u_1(std::size_t fmax, const SizeList& vv, std::size_t l, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto frmul = historyRepasFudRepasMultiply_u;

    auto t0 = clk::now();
    auto mark = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    auto& llu = ur.listVarSizePair;
    auto z = hr0.size;
    auto& mvv = hr0.mapVarInt();
    auto sh = hr0.shape;
    std::size_t ms = 2;
    for (auto v : vv)
    {
	auto s = sh[mvv[v]];
	if (s > ms)
	    ms = s;
    }
    auto pl = mvv[l];
    auto ls = sh[pl];
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax);
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	auto n = hr0.dimension;
	auto plz = pl*z;
	auto rr = hr0.arr;
	double g = 1.0 / (double)z;
	auto ee0 = new std::size_t[ls];
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 0;
	if (hr0.evient)
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t w = rr[j*n + pl];
		ee0[w]++;
	    }
	else
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t w = rr[plz + j];
		ee0[w]++;
	    }
	double e0 = 0.0;
	for (std::size_t k = 0; k < ls; k++)
	{
	    double a = (double)ee0[k] * g;
	    if (a > 0.0)
		e0 -= a * log(a);
	}
	delete ee0;
	if (e0 <= repaRounding)
	{
	    std::cout << "no label entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	double e = e0;
	std::size_t ew = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (auto v : vv)
	{
	    auto p = mvv[v];
	    auto pz = p*z;
	    auto s = sh[p];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 0;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 0;
	    if (hr0.evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    auto jn = j*n;
		    std::size_t u = rr[jn + p];
		    std::size_t w = rr[jn + pl];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr[pz + j];
		    std::size_t w = rr[plz + j];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
	    {
		double a = (double)ee1[k] * g;
		if (a > 0.0)
		    e1 += a * log(a);
	    }
	    for (std::size_t k = 0; k < s; k++)
	    {
		double a = (double)ee2[k] * g;
		if (a > 0.0)
		    e2 += a * log(a);
	    }
	    if (v != l && e > e2 - e1)
	    {
		e = e2 - e1;
		ew = v;
	    }
	}
	delete ee2;
	delete ee1;
	if (e >= e0 - repaRounding)
	{
	    std::cout << "no conditional entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "sized entropy label : " << (double)z * e0 << std::endl;
	std::cout << "sized conditional entropy: " << (double)z * e << std::endl;
	std::cout << "sized entropy decrease: " << (double)z * (e0 - e) << std::endl;
	std::cout << "entropy variable: " << *llu[ew].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh[mvv[ew]];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 1;
	    tr->vectorVar = new std::size_t[1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[1];
	    auto sh = tr->shape;
	    ww[0] = ew;
	    sh[0] = sz;
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	auto nn = treesLeafNodes(*dr->slices);
	auto ee0 = new std::size_t[ls];
	std::unique_ptr<HistoryRepa> hr;
	std::size_t z2 = z;
	std::size_t v2 = 0;
	{
	    mark = clk::now();
	    std::cout << ">>> applier " << std::endl;
	    hr = frmul(hr0, *dr->fud);
	    time["applier"] = ((sec)(clk::now() - mark)).count();
	    std::cout << "<<< applier " << time["applier"] << "s" << std::endl;
	    mark = clk::now();
	    std::cout << ">>> slicer " << std::endl;
	    auto n = hr->dimension;
	    auto& mvv1 = hr->mapVarInt();
	    auto rr = hr->arr;
	    DoubleSizeSizeTupleList zs;
	    zs.reserve(nn->size());
	    for (auto& p : *nn)
	    {
		auto v = p.first;
		if (ig.find(v) != ig.end())
		    continue;
		std::size_t a = 0;
		auto pk = mvv1[v];
		for (std::size_t k = 0; k < ls; k++)
		    ee0[k] = 0;
		if (hr->evient)
		    for (std::size_t j = 0; j < z; j++)
		    {
			std::size_t u = rr[j*n + pk];
			if (u)
			{
			    a++;
			    std::size_t w = rr[j*n + pl];
			    ee0[w]++;
			}
		    }
		else
		    for (std::size_t j = 0; j < z; j++)
		    {
			std::size_t u = rr[pk*z + j];
			if (u)
			{
			    a++;
			    std::size_t w = rr[pl*z + j];
			    ee0[w]++;
			}
		    }
		if (a > 1)
		{
		    double g = 1.0 / (double)a;
		    double e = 0.0;
		    for (std::size_t k = 0; k < ls; k++)
		    {
			double b = (double)ee0[k] * g;
			if (b > 0.0)
			    e -= b * log(b);
		    }
		    if ((double)a * e > repaRounding)
			zs.push_back(DoubleSizeSizeTuple((double)a * e, a, v));
		}
	    }
	    if (!zs.size())
	    {
		std::cout << "no slices" << std::endl;
		break;
	    }
	    std::sort(zs.begin(), zs.end());
	    z2 = std::get<1>(zs.back());
	    v2 = std::get<2>(zs.back());
	    std::cout << "slice size: " << z2 << std::endl;
	    std::cout << "slice variable: " << *llu[v2].first << std::endl;
	    SizeList ev;
	    ev.reserve(z2);
	    {
		auto pk = mvv1[v2];
		if (hr->evient)
		    for (std::size_t j = 0; j < z; j++)
		    {
			std::size_t u = rr[j*n + pk];
			if (u)
			    ev.push_back(j);
		    }
		else
		    for (std::size_t j = 0; j < z; j++)
		    {
			std::size_t u = rr[pk*z + j];
			if (u)
			    ev.push_back(j);
		    }
	    }
	    hr = hrsel(ev.size(), ev.data(), hr0);
	    time["slicer"] = ((sec)(clk::now() - mark)).count();
	    std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	}
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	auto n = hr->dimension;
	auto plz = pl*z2;
	auto rr = hr->arr;
	double g = 1.0 / (double)z2;
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 0;
	if (hr->evient)
	    for (std::size_t j = 0; j < z2; j++)
	    {
		std::size_t w = rr[j*n + pl];
		ee0[w]++;
	    }
	else
	    for (std::size_t j = 0; j < z2; j++)
	    {
		std::size_t w = rr[plz + j];
		ee0[w]++;
	    }
	double e0 = 0.0;
	for (std::size_t k = 0; k < ls; k++)
	{
	    double a = (double)ee0[k] * g;
	    if (a > 0.0)
		e0 -= a * log(a);
	}
	delete ee0;
	double e = e0;
	std::size_t ew = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (auto v : vv)
	{
	    auto p = mvv[v];
	    auto pz = p*z2;
	    auto s = sh[p];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 0;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 0;
	    if (hr->evient)
		for (std::size_t j = 0; j < z2; j++)
		{
		    auto jn = j*n;
		    std::size_t u = rr[jn + p];
		    std::size_t w = rr[jn + pl];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    else
		for (std::size_t j = 0; j < z2; j++)
		{
		    std::size_t u = rr[pz + j];
		    std::size_t w = rr[plz + j];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
	    {
		double a = (double)ee1[k] * g;
		if (a > 0.0)
		    e1 += a * log(a);
	    }
	    for (std::size_t k = 0; k < s; k++)
	    {
		double a = (double)ee2[k] * g;
		if (a > 0.0)
		    e2 += a * log(a);
	    }
	    if (v != l && e > e2 - e1)
	    {
		e = e2 - e1;
		ew = v;
	    }
	}
	delete ee2;
	delete ee1;
	if (e >= e0 - repaRounding)
	{
	    ig.insert(v2);
	    std::cout << "no conditional entropy" << std::endl;
	    continue;
	}
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "sized entropy label : " << (double)z2 * e0 << std::endl;
	std::cout << "sized conditional entropy: " << (double)z2 * e << std::endl;
	std::cout << "sized entropy decrease: " << (double)z2 * (e0 - e) << std::endl;
	std::cout << "entropy variable: " << *llu[ew].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh[mvv[ew]];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 2;
	    tr->vectorVar = new std::size_t[2];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[2];
	    auto sh = tr->shape;
	    ww[0] = v2;
	    sh[0] = 2;
	    ww[1] = ew;
	    sh[1] = sz;
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.push_back(ll);
	for (auto& p : *nn)
	    if (p.first == v2)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}

// parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_u ::
//  Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_u(std::size_t fmax, const SizeList& vv, std::size_t l, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto frmul = historyRepasFudRepasMultiply_u;

    auto t0 = clk::now();
    auto mark = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    auto n = vv.size();
    auto& llu = ur.listVarSizePair;
    auto n0 = hr0.dimension;
    auto z = hr0.size;
    auto& mvv = hr0.mapVarInt();
    auto sh0 = hr0.shape;
    auto rr0 = hr0.arr;
    std::size_t ms = 2;
    std::size_t ms1 = 2;
    for (auto v : vv)
    {
	auto s = sh0[mvv[v]];
	if (s > ms)
	    ms = s;
	if (v != l && s > ms1)
	    ms1 = s;
    }
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax);
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    HistoryRepa hr1;
    hr1.evient = false;
    std::size_t m1 = n + 1;
    std::size_t n1 = m1 + fmax*ms1;
    hr1.dimension = n1;
    hr1.vectorVar = new std::size_t[n1];
    auto vv1 = hr1.vectorVar;
    hr1.shape = new std::size_t[n1];
    auto sh1 = hr1.shape;
    auto pp1 = new std::size_t[m1];
    vv1[0] = l;
    pp1[0] = mvv[l];
    auto ls = sh0[pp1[0]];
    sh1[0] = ls;
    for (std::size_t i = 0; i < n; i++)
    {
	auto v = vv[i];
	vv1[i + 1] = v;
	auto p = mvv[v];
	pp1[i + 1] = p;
	sh1[i + 1] = sh0[p];
    }
    hr1.arr = new unsigned char[z*n1];
    auto rr1 = hr1.arr;
    if (hr0.evient)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn0 = j*n0;
	    for (std::size_t i = 0; i < m1; i++)
		rr1[i*z + j] = rr0[jn0 + pp1[i]];
	}
    else
	for (std::size_t i = 0; i < m1; i++)
	{
	    std::size_t iz0 = pp1[i] * z;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
		rr1[iz + j] = rr0[iz0 + j];
	}
    delete pp1;
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	auto ee0 = new std::size_t[ls];
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t w = rr1[j];
	    ee0[w]++;
	}
	double e0 = alngam((double)z + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	if (e0 <= repaRounding)
	{
	    std::cout << "no label entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	double e = e0;
	std::size_t i1 = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (std::size_t i = 1; i <= n; i++)
	{
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 1;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 1;
	    auto iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		std::size_t w = rr1[j];
		ee1[w*s + u]++;
		ee2[u]++;
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
		e1 += alngam((double)ee1[k]);
	    for (std::size_t k = 0; k < s; k++)
		e2 += alngam((double)ee2[k]);
	    if (vv1[i] != l && (e > e2 - e1))
	    {
		e = e2 - e1;
		i1 = i;
	    }
	}
	delete ee2;
	delete ee1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    std::cout << "no conditional entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 1;
	    tr->vectorVar = new std::size_t[1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[1];
	    auto sh = tr->shape;
	    ww[0] = v1;
	    sh[0] = sz;
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto ee0 = new std::size_t[ls];
	auto nn = treesLeafNodes(*dr->slices);
	DoubleSizeSizeTupleList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t i = 1 + n;
	    while (i < m1)
	    {
		if (vv1[i] == v)
		    break;
		i++;
	    }
	    if (i == m1)
		continue;
	    std::size_t a = 0;
	    for (std::size_t k = 0; k < ls; k++)
		ee0[k] = 1;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		if (u)
		{
		    a++;
		    std::size_t w = rr1[j];
		    ee0[w]++;
		}
	    }
	    if (a > 1)
	    {
		double e = alngam((double)a + 1.0);
		for (std::size_t k = 0; k < ls; k++)
		    e -= alngam((double)ee0[k]);
		if (e > repaRounding)
		    zs.push_back(DoubleSizeSizeTuple(e, a, i));
	    }
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = std::get<1>(zs.back());
	auto i2 = std::get<2>(zs.back());
	auto v2 = vv1[i2];
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v2].first << std::endl;
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	auto i2z = i2*z;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t t = rr1[i2z + j];
	    if (t)
	    {
		std::size_t w = rr1[j];
		ee0[w]++;
	    }
	}
	double e0 = alngam((double)z2 + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	double e = e0;
	std::size_t i1 = 0;
	auto ee1 = new std::size_t[ms*ls];
	auto ee2 = new std::size_t[ms];
	for (std::size_t i = 1; i <= n; i++)
	{
	    auto iz = i*z;
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 1;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 1;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t t = rr1[i2z + j];
		if (t)
		{
		    std::size_t u = rr1[iz + j];
		    std::size_t w = rr1[j];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
		e1 += alngam((double)ee1[k]);
	    for (std::size_t k = 0; k < s; k++)
		e2 += alngam((double)ee2[k]);
	    if (vv1[i] != l && e > e2 - e1)
	    {
		e = e2 - e1;
		i1 = i;
	    }
	}
	delete ee2;
	delete ee1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    ig.insert(v2);
	    std::cout << "no conditional entropy" << std::endl;
	    continue;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 2;
	    tr->vectorVar = new std::size_t[2];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[2];
	    auto sh = tr->shape;
	    ww[0] = v2;
	    sh[0] = 2;
	    ww[1] = v1;
	    sh[1] = sz;
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i2z + j];
		k = sz*k + rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	for (auto& p : *nn)
	    if (p.first == v2)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}


void parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_root(std::size_t tint, std::size_t n, std::size_t z, std::size_t* vv1, std::size_t* sh1, unsigned char* rr1, std::size_t ms, std::size_t l, std::size_t ls, std::size_t t, double* e, std::size_t* i1)
{
    auto ee1 = new std::size_t[ms*ls];
    auto ee2 = new std::size_t[ms];
    for (std::size_t i = 1; i <= n; i++)
	if (i % tint == t)
	{
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 1;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 1;
	    auto iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		std::size_t w = rr1[j];
		ee1[w*s + u]++;
		ee2[u]++;
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
		e1 += alngam((double)ee1[k]);
	    for (std::size_t k = 0; k < s; k++)
		e2 += alngam((double)ee2[k]);
	    if (vv1[i] != l && (e[t] > e2 - e1))
	    {
		e[t] = e2 - e1;
		i1[t] = i;
	    }
	}
    delete ee2;
    delete ee1;
}

void parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_child(std::size_t tint, std::size_t n, std::size_t z, std::size_t i2z, std::size_t* vv1, std::size_t* sh1, unsigned char* rr1, std::size_t ms, std::size_t l, std::size_t ls, std::size_t t, double* e, std::size_t* i1)
{
    auto ee1 = new std::size_t[ms*ls];
    auto ee2 = new std::size_t[ms];
    for (std::size_t i = 1; i <= n; i++)
	if (i % tint == t)
	{
	    auto iz = i*z;
	    auto s = sh1[i];
	    for (std::size_t k = 0; k < s*ls; k++)
		ee1[k] = 1;
	    for (std::size_t k = 0; k < s; k++)
		ee2[k] = 1;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t q = rr1[i2z + j];
		if (q)
		{
		    std::size_t u = rr1[iz + j];
		    std::size_t w = rr1[j];
		    ee1[w*s + u]++;
		    ee2[u]++;
		}
	    }
	    double e1 = 0.0;
	    double e2 = 0.0;
	    for (std::size_t k = 0; k < s*ls; k++)
		e1 += alngam((double)ee1[k]);
	    for (std::size_t k = 0; k < s; k++)
		e2 += alngam((double)ee2[k]);
	    if (vv1[i] != l && (e[t] > e2 - e1))
	    {
		e[t] = e2 - e1;
		i1[t] = i;
	    }
	}
    delete ee2;
    delete ee1;
}

// parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up ::
//  Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up(std::size_t fmax, std::size_t tint, const SizeList& vv, std::size_t l, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto root = parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_root;
    auto child = parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_child;

    auto t0 = clk::now();
    auto mark = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    auto n = vv.size();
    std::size_t tint1 = std::min(tint, n);
    auto& llu = ur.listVarSizePair;
    auto n0 = hr0.dimension;
    auto z = hr0.size;
    auto& mvv = hr0.mapVarInt();
    auto sh0 = hr0.shape;
    auto rr0 = hr0.arr;
    std::size_t ms = 2;
    std::size_t ms1 = 2;
    for (auto v : vv)
    {
	auto s = sh0[mvv[v]];
	if (s > ms)
	    ms = s;
	if (v != l && s > ms1)
	    ms1 = s;
    }
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax);
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    HistoryRepa hr1;
    hr1.evient = false;
    std::size_t m1 = n + 1;
    std::size_t n1 = m1 + fmax*ms1;
    hr1.dimension = n1;
    hr1.vectorVar = new std::size_t[n1];
    auto vv1 = hr1.vectorVar;
    hr1.shape = new std::size_t[n1];
    auto sh1 = hr1.shape;
    auto pp1 = new std::size_t[m1];
    vv1[0] = l;
    pp1[0] = mvv[l];
    auto ls = sh0[pp1[0]];
    sh1[0] = ls;
    for (std::size_t i = 0; i < n; i++)
    {
	auto v = vv[i];
	vv1[i + 1] = v;
	auto p = mvv[v];
	pp1[i + 1] = p;
	sh1[i + 1] = sh0[p];
    }
    hr1.arr = new unsigned char[z*n1];
    auto rr1 = hr1.arr;
    if (hr0.evient)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn0 = j*n0;
	    for (std::size_t i = 0; i < m1; i++)
		rr1[i*z + j] = rr0[jn0 + pp1[i]];
	}
    else
	for (std::size_t i = 0; i < m1; i++)
	{
	    std::size_t iz0 = pp1[i] * z;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
		rr1[iz + j] = rr0[iz0 + j];
	}
    delete pp1;
    std::size_t f = 1;
    SizeSizeTreePairList nn;
    SizeSizeTreePairList nn1;
    nn1.reserve(fmax*ms1);
    {
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	auto ee0 = new std::size_t[ls];
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t w = rr1[j];
	    ee0[w]++;
	}
	double e0 = alngam((double)z + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	if (e0 <= repaRounding)
	{
	    std::cout << "no label entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto ee = new double[tint1];
	auto ii1 = new std::size_t[tint1];
	std::vector<std::thread> threads;
	threads.reserve(tint1);
	for (std::size_t t = 0; t < tint1; t++)
	{
	    ee[t] = e0;
	    ii1[t] = 0;
	    threads.push_back(std::thread(root, tint1, n, z, vv1, sh1, rr1, ms, l, ls, t, ee, ii1));
	}
	for (auto& t : threads)
	    t.join();
	double e = e0;
	std::size_t i1 = 0;
	for (std::size_t t = 0; t < tint1; t++)
	    if (e > ee[t])
	    {
		e = ee[t];
		i1 = ii1[t];
	    }
	delete ee;
	delete ii1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    std::cout << "no conditional entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 1;
	    tr->vectorVar = new std::size_t[1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[1];
	    auto sh = tr->shape;
	    ww[0] = v1;
	    sh[0] = sz;
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	{
	    SizeSizeTreePair q(s, std::make_shared<SizeTree>());
	    nn.push_back(q);
	    nn1.push_back(q);
	    dr->slices->_list.push_back(q);
	}
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    DoubleSizeSizeTupleList zs;
    zs.reserve(fmax*ms1);
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	if (zs.size())
	    zs.pop_back();
	auto ee0 = new std::size_t[ls];
	for (auto& p : nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t i = 1 + n;
	    while (i < m1)
	    {
		if (vv1[i] == v)
		    break;
		i++;
	    }
	    if (i == m1)
		continue;
	    std::size_t a = 0;
	    for (std::size_t k = 0; k < ls; k++)
		ee0[k] = 1;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		if (u)
		{
		    a++;
		    std::size_t w = rr1[j];
		    ee0[w]++;
		}
	    }
	    if (a > 1)
	    {
		double e = alngam((double)a + 1.0);
		for (std::size_t k = 0; k < ls; k++)
		    e -= alngam((double)ee0[k]);
		if (e > repaRounding)
		    zs.push_back(DoubleSizeSizeTuple(e, a, i));
	    }
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = std::get<1>(zs.back());
	auto i2 = std::get<2>(zs.back());
	auto v2 = vv1[i2];
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v2].first << std::endl;
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	auto i2z = i2*z;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t t = rr1[i2z + j];
	    if (t)
	    {
		std::size_t w = rr1[j];
		ee0[w]++;
	    }
	}
	double e0 = alngam((double)z2 + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	auto ee = new double[tint1];
	auto ii1 = new std::size_t[tint1];
	std::vector<std::thread> threads;
	threads.reserve(tint1);
	for (std::size_t t = 0; t < tint1; t++)
	{
	    ee[t] = e0;
	    ii1[t] = 0;
	    threads.push_back(std::thread(child, tint1, n, z, i2z, vv1, sh1, rr1, ms, l, ls, t, ee, ii1));
	}
	for (auto& t : threads)
	    t.join();
	double e = e0;
	std::size_t i1 = 0;
	for (std::size_t t = 0; t < tint1; t++)
	    if (e > ee[t])
	    {
		e = ee[t];
		i1 = ii1[t];
	    }
	delete ee;
	delete ii1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    ig.insert(v2);
	    std::cout << "no conditional entropy" << std::endl;
	    continue;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 2;
	    tr->vectorVar = new std::size_t[2];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[2];
	    auto sh = tr->shape;
	    ww[0] = v2;
	    sh[0] = 2;
	    ww[1] = v1;
	    sh[1] = sz;
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i2z + j];
		k = sz*k + rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	nn.clear();
	for (auto& p : nn1)
	    if (p.first == v2)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		{
		    SizeSizeTreePair q(s, std::make_shared<SizeTree>());
		    nn.push_back(q);
		    nn1.push_back(q);
		    p.second->_list.push_back(q);
		}
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}

// parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_1 ::
//  Integer -> Integer ->
//   [VariableRepa] -> HistoryRepa -> Integer ->
//   IO (SystemRepa, ApplicationRepa)
std::unique_ptr<ApplicationRepa> Alignment::parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_1(std::size_t fmax, std::size_t tint, const SizeList& vv, std::size_t l, const HistoryRepa& hr0, int d, SystemRepa& ur)
{
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto root = parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_root;
    auto child = parametersSystemsHistoryRepasApplicationerCondMultinomialFmaxIORepa_up_child;

    auto t0 = clk::now();
    auto mark = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    auto n = vv.size();
    auto& llu = ur.listVarSizePair;
    auto n0 = hr0.dimension;
    auto z = hr0.size;
    auto& mvv = hr0.mapVarInt();
    auto sh0 = hr0.shape;
    auto rr0 = hr0.arr;
    std::size_t ms = 2;
    std::size_t ms1 = 2;
    for (auto v : vv)
    {
	auto s = sh0[mvv[v]];
	if (s > ms)
	    ms = s;
	if (v != l && s > ms1)
	    ms1 = s;
    }
    auto vd = std::make_shared<Variable>(d);
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax);
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    HistoryRepa hr1;
    hr1.evient = false;
    std::size_t m1 = n + 1;
    std::size_t n1 = m1 + fmax*ms1;
    hr1.dimension = n1;
    hr1.vectorVar = new std::size_t[n1];
    auto vv1 = hr1.vectorVar;
    hr1.shape = new std::size_t[n1];
    auto sh1 = hr1.shape;
    auto pp1 = new std::size_t[m1];
    vv1[0] = l;
    pp1[0] = mvv[l];
    auto ls = sh0[pp1[0]];
    sh1[0] = ls;
    for (std::size_t i = 0; i < n; i++)
    {
	auto v = vv[i];
	vv1[i + 1] = v;
	auto p = mvv[v];
	pp1[i + 1] = p;
	sh1[i + 1] = sh0[p];
    }
    hr1.arr = new unsigned char[z*n1];
    auto rr1 = hr1.arr;
    if (hr0.evient)
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t jn0 = j*n0;
	    for (std::size_t i = 0; i < m1; i++)
		rr1[i*z + j] = rr0[jn0 + pp1[i]];
	}
    else
	for (std::size_t i = 0; i < m1; i++)
	{
	    std::size_t iz0 = pp1[i] * z;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
		rr1[iz + j] = rr0[iz0 + j];
	}
    delete pp1;
    std::size_t f = 1;
    {
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	auto ee0 = new std::size_t[ls];
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t w = rr1[j];
	    ee0[w]++;
	}
	double e0 = alngam((double)z + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	if (e0 <= repaRounding)
	{
	    std::cout << "no label entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto ee = new double[tint];
	auto ii1 = new std::size_t[tint];
	std::vector<std::thread> threads;
	threads.reserve(tint);
	for (std::size_t t = 0; t < tint; t++)
	{
	    ee[t] = e0;
	    ii1[t] = 0;
	    threads.push_back(std::thread(root, tint, n, z, vv1, sh1, rr1, ms, l, ls, t, ee, ii1));
	}
	for (auto& t : threads)
	    t.join();
	double e = e0;
	std::size_t i1 = 0;
	for (std::size_t t = 0; t < tint; t++)
	    if (e > ee[t])
	    {
		e = ee[t];
		i1 = ii1[t];
	    }
	delete ee;
	delete ii1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    std::cout << "no conditional entropy" << std::endl;
	    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	    return dr;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 1;
	    tr->vectorVar = new std::size_t[1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[1];
	    auto sh = tr->shape;
	    ww[0] = v1;
	    sh[0] = sz;
	    tr->arr = new unsigned char[sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < sz; j++)
		rr[j] = 0;
	    rr[i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	dr->slices->_list.reserve(sz);
	for (auto& s : sl)
	    dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    SizeSet ig;
    while (f < fmax)
    {
	mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto ee0 = new std::size_t[ls];
	auto nn = treesLeafNodes(*dr->slices);
	DoubleSizeSizeTupleList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t i = 1 + n;
	    while (i < m1)
	    {
		if (vv1[i] == v)
		    break;
		i++;
	    }
	    if (i == m1)
		continue;
	    std::size_t a = 0;
	    for (std::size_t k = 0; k < ls; k++)
		ee0[k] = 1;
	    std::size_t iz = i*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t u = rr1[iz + j];
		if (u)
		{
		    a++;
		    std::size_t w = rr1[j];
		    ee0[w]++;
		}
	    }
	    if (a > 1)
	    {
		double e = alngam((double)a + 1.0);
		for (std::size_t k = 0; k < ls; k++)
		    e -= alngam((double)ee0[k]);
		if (e > repaRounding)
		    zs.push_back(DoubleSizeSizeTuple(e, a, i));
	    }
	}
	if (!zs.size())
	{
	    std::cout << "no slices" << std::endl;
	    break;
	}
	std::sort(zs.begin(), zs.end());
	auto z2 = std::get<1>(zs.back());
	auto i2 = std::get<2>(zs.back());
	auto v2 = vv1[i2];
	std::cout << "slice size: " << z2 << std::endl;
	std::cout << "slice variable: " << *llu[v2].first << std::endl;
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> searcher " << std::endl;
	for (std::size_t k = 0; k < ls; k++)
	    ee0[k] = 1;
	auto i2z = i2*z;
	for (std::size_t j = 0; j < z; j++)
	{
	    std::size_t t = rr1[i2z + j];
	    if (t)
	    {
		std::size_t w = rr1[j];
		ee0[w]++;
	    }
	}
	double e0 = alngam((double)z2 + 1.0);
	for (std::size_t k = 0; k < ls; k++)
	    e0 -= alngam((double)ee0[k]);
	delete ee0;
	auto ee = new double[tint];
	auto ii1 = new std::size_t[tint];
	std::vector<std::thread> threads;
	threads.reserve(tint);
	for (std::size_t t = 0; t < tint; t++)
	{
	    ee[t] = e0;
	    ii1[t] = 0;
	    threads.push_back(std::thread(child, tint, n, z, i2z, vv1, sh1, rr1, ms, l, ls, t, ee, ii1));
	}
	for (auto& t : threads)
	    t.join();
	double e = e0;
	std::size_t i1 = 0;
	for (std::size_t t = 0; t < tint; t++)
	    if (e > ee[t])
	    {
		e = ee[t];
		i1 = ii1[t];
	    }
	delete ee;
	delete ii1;
	if (i1 == 0 || e >= e0 - repaRounding)
	{
	    ig.insert(v2);
	    std::cout << "no conditional entropy" << std::endl;
	    continue;
	}
	auto v1 = vv1[i1];
	time["searcher"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< searcher " << time["searcher"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> transer " << std::endl;
	f++;
	std::cout << "fud: " << f << std::endl;
	std::cout << "fud slice size: " << z2 << std::endl;
	std::cout << "sized entropy label : " << e0 << std::endl;
	std::cout << "sized conditional entropy: " << e << std::endl;
	std::cout << "sized entropy decrease: " << e0 - e << std::endl;
	std::cout << "entropy variable: " << *llu[v1].first << std::endl;
	auto vf = std::make_shared<Variable>((int)f);
	auto vdf = std::make_shared<Variable>(vd, vf);
	auto vfl = std::make_shared<Variable>(vdf, vl);
	SizeList sl;
	TransformRepaPtrList ll;
	std::size_t sz = sh1[i1];
	sl.reserve(sz);
	ll.reserve(sz);
	std::size_t b = 1;
	for (std::size_t i = 0; i < sz; i++)
	{
	    auto tr = std::make_shared<TransformRepa>();
	    tr->dimension = 2;
	    tr->vectorVar = new std::size_t[2];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[2];
	    auto sh = tr->shape;
	    ww[0] = v2;
	    sh[0] = 2;
	    ww[1] = v1;
	    sh[1] = sz;
	    tr->arr = new unsigned char[2 * sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2 * sz; j++)
		rr[j] = 0;
	    rr[sz + i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>((int)b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	    vv1[m1] = w;
	    sh1[m1] = 2;
	    std::size_t i1z = i1*z;
	    std::size_t m1z = m1*z;
	    for (std::size_t j = 0; j < z; j++)
	    {
		std::size_t k = rr1[i2z + j];
		k = sz*k + rr1[i1z + j];
		rr1[m1z + j] = rr[k];
	    }
	    m1++;
	}
	dr->fud->layers.push_back(ll);
	for (auto& p : *nn)
	    if (p.first == v2)
	    {
		p.second->_list.reserve(sz);
		for (auto& s : sl)
		    p.second->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));
		break;
	    }
	std::cout << "fud slice cardinality: " << ll.size() << std::endl;
	time["transer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< transer " << time["transer"] << "s" << std::endl;
    }
    std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
    return dr;
}


