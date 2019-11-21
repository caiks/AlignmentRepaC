#include "AlignmentRandomRepa.h"
#include "AlignmentPracticableIORepa.h"
#include <chrono>
#include <ctime>
#include <cmath>

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
    auto vf = std::make_shared<Variable>(f);
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
	auto vl = std::make_shared<Variable>(l);
	auto vfl = std::make_shared<Variable>(vf,vl);
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
			rr[j] = u;
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
		    auto vb = std::make_shared<Variable>(b);
		    auto vflb = std::make_shared<Variable>(vfl, vb);
		    llu.push_back(VarSizePair(*vflb,s));
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
    auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
    auto hrred = setVarsHistoryRepasReduce_u;
    auto hrconcat = vectorHistoryRepasConcat_u;
    auto hrshuffle = historyRepasShuffle_u;
    auto llfr = setVariablesListTransformRepasFudRepa_u;
    auto frmul = historyRepasFudRepasMultiply_u;
    auto frdep = fudsSetVarsDepends;
    auto layerer = parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u;

    auto t0 = clk::now();
    std::map<std::string, double> time;
    std::cout << ">>> applicationer" << std::endl;
    SizeUSet vv1(vv.begin(), vv.end());
    auto& llu = ur.listVarSizePair;
    auto z = hr.size;
    auto vl = std::make_shared<Variable>("s");
    auto dr = std::make_unique<ApplicationRepa>();
    dr->substrate = vv;
    dr->fud = std::make_shared<FudRepa>();
    dr->slices = std::make_shared<SizeTree>();
    dr->fud->layers.reserve(fmax*(lmax+1));
    if (!z)
    {
	std::cout << "empty history" << std::endl;
	std::cout << "<<< applicationer " << ((sec)(clk::now() - t0)).count() << "s" << std::endl;
	return dr;
    }
    std::size_t f = 1;
    {
	auto mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	HistoryRepaPtrList qq;
	qq.reserve(mult);
	for (std::size_t i = 1; i <= mult; i++)
	    qq.push_back(std::move(hrshuffle(hr, seed + i*z)));
	auto hrs = hrconcat(qq);
        qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, vv, hr, *hrs, f, ur);
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
	std::cout << "derived algn density per size: " << a / (double)z << std::endl;
	std::cout << "derived impl bi-valency percent: " << 50.0 * std::exp (a / (double)z / (double)(m -1)) << std::endl;
	auto vf = std::make_shared<Variable>(f);
	auto vfl = std::make_shared<Variable>(vf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(hr, *gr));
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
	    auto vb = std::make_shared<Variable>(b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(*vflb, 2));
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
	    auto vb = std::make_shared<Variable>(b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(*vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.insert(dr->fud->layers.end(), fr->layers.begin(), fr->layers.end());
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
	auto mark = clk::now();
	std::cout << ">>> slicer " << std::endl;
	auto hr1 = frmul(hr, *dr->fud);
	auto n1 = hr1->dimension;
	auto& mvv1 = hr1->mapVarInt();
	auto rr1 = hr1->arr;
	auto nn = treesLeafNodes(*dr->slices);
	SizeSizePairList zs;
	zs.reserve(nn->size());
	for (auto& p : *nn)
	{
	    auto v = p.first;
	    if (ig.find(v) != ig.end())
		continue;
	    std::size_t a = 0;
	    auto pk = mvv1[v];
	    if (hr1->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr1[j*n1 + pk];
		    if (u)
			a++;
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr1[pk*z + j];
		    if (u)
			a++;
		}
	    if (a > 1)
		zs.push_back(SizeSizePair(a,v));
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
	std::cout << "slice variable: " << llu[v] << std::endl;
	SizeList ev;
	ev.reserve(z2);
	{
	    auto pk = mvv1[v];
	    if (hr1->evient)
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr1[j*n1 + pk];
		    if (u)
			ev.push_back(j);
		}
	    else
		for (std::size_t j = 0; j < z; j++)
		{
		    std::size_t u = rr1[pk*z + j];
		    if (u)
			ev.push_back(j);
		}
	}
	auto hr2 = hrsel(ev.size(), ev.data(), hr);
	time["slicer"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< slicer " << time["slicer"] << "s" << std::endl;
	mark = clk::now();
	std::cout << ">>> shuffler " << std::endl;
	HistoryRepaPtrList qq;
	qq.reserve(mult);
	for (std::size_t i = 1; i <= mult; i++)
	    qq.push_back(std::move(hrshuffle(*hr2, seed + i*z2)));
	auto hr2s = hrconcat(qq);
        qq.clear();
	time["shuffler"] = ((sec)(clk::now() - mark)).count();
	std::cout << "<<< shuffler " << time["shuffler"] << "s" << std::endl;
	std::unique_ptr<FudRepa> fr;
	std::unique_ptr<DoubleSizeListPairList> mm;
	try
	{
	    auto t = layerer(wmax, lmax, xmax, omax, bmax, mmax, umax, pmax, vv, *hr2, *hr2s, f+1, ur);
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
	std::cout << "derived algn density per size: " << a / (double)z2 << std::endl;
	std::cout << "derived impl bi-valency percent: " << 50.0 * std::exp (a / (double)z2 / (double)(m -1)) << std::endl;
	auto vf = std::make_shared<Variable>(f);
	auto vfl = std::make_shared<Variable>(vf, vl);
	SizeUSet kk1(kk.begin(), kk.end());
	auto gr = llfr(vv1, *frdep(*fr, kk1));
	auto ar = hrred(1.0, m, kk.data(), *frmul(*hr2, *gr));
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
	    tr->dimension = m+1;
	    tr->vectorVar = new std::size_t[m+1];
	    auto ww = tr->vectorVar;
	    tr->shape = new std::size_t[m+1];
	    auto sh = tr->shape;
	    ww[0] = v;
	    sh[0] = 2;
	    for (std::size_t j = 0; j < m; j++)
	    {
		ww[j+1] = kk[j];
		sh[j+1] = skk[j];
	    }
	    tr->arr = new unsigned char[2*sz];
	    auto rr = tr->arr;
	    for (std::size_t j = 0; j < 2*sz; j++)
		rr[j] = 0;
	    rr[sz+i] = 1;
	    tr->valency = 2;
	    auto vb = std::make_shared<Variable>(b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(*vflb, 2));
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
	    auto vb = std::make_shared<Variable>(b++);
	    auto vflb = std::make_shared<Variable>(vfl, vb);
	    llu.push_back(VarSizePair(*vflb, 2));
	    auto w = llu.size() - 1;
	    tr->derived = w;
	    sl.push_back(w);
	    ll.push_back(tr);
	}
	dr->fud->layers.insert(dr->fud->layers.end(), fr->layers.begin(), fr->layers.end());
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

