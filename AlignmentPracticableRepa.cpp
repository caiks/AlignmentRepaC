﻿#include "AlignmentPracticableRepa.h"

using namespace Alignment;

// parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui ::
//   Integer -> Integer -> Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([[VariableRepa]],Integer)
std::tuple<std::unique_ptr<SizeListList>, std::size_t> Alignment::parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, const SizeList& vv, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto frvars = fudRepasSetVar;
	auto frder = fudRepasDerived;
	auto cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u;
	auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u;

	SizeList vv1;
	vv1.reserve(vv.size());
	{
		auto n = hh.dimension;
		auto& mvv = hh.mapVarInt();
		auto sh = hh.shape;
		auto pxx1 = hhx.arr;
		std::size_t sz = 0;
		SizeList xx;
		xx.reserve(n);
		for (std::size_t i = 0; i < n; i++)
		{
			xx.push_back(sz);
			sz += sh[i];
		}
		for (auto& w : vv)
		{
			auto p = mvv[w];
			auto e = sh[p];
			bool first = false;
			bool second = false;
			for (std::size_t i = 0; !second && i < e; i++)
			{
				bool nz = pxx1[xx[p] + i] > 0.0;
				second = second || (first && nz);
				first = first || nz;
			}
			if (second)
				vv1.push_back(w);
		}
	}
	std::size_t s = 0;
	auto xx1 = std::make_unique<SizeListList>();
	if (vv1.size() < 2)
		return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
	DoubleSizeListPairList xx;
	xx.reserve(omax*10);
	if (!fr.layers.size())
	{
		auto t = cross(xmax, omax, vv1, hh, hhx, hhrr, hhrrx);
		auto& xc = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xc->begin(), xc->end());
		xx1->reserve(xc->size());
		for (auto& pp : *xc)
			xx1->push_back(pp.second);
	}
	else
	{
		auto vvf = frvars(fr);
		vv1.reserve(vv1.size() + vvf->size());
		vv1.insert(vv1.end(), vvf->begin(), vvf->end());
		auto wwf = frder(fr);
		xx1->reserve(wwf->size());
		for (auto& w : *wwf)
			xx1->push_back(SizeList{ w });
	}
	while (xx1->size())
	{
		auto t = append(xmax, omax, vv1, *xx1, hh, hhx, hhrr, hhrrx);
		auto& xa = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xa->begin(), xa->end());
		xx1->clear();
		xx1->reserve(xa->size());
		for (auto& pp : *xa)
			xx1->push_back(pp.second);
	}
	if (!xx.size())
		return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
	std::sort(xx.begin(), xx.end());
	xx1->reserve(bmax / mmax);
	std::size_t start = xx.size() - 1;
	std::size_t end = xx.size() > bmax/mmax ? xx.size() - bmax/mmax : 0;
	for (long long i = start; i >= (long long)end; i--)
		xx1->push_back(xx[i].second);
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
}

// parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_uip ::
//   Integer -> Integer -> Integer -> Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([[VariableRepa]],Integer)
std::tuple<std::unique_ptr<SizeListList>, std::size_t> Alignment::parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_uip(std::size_t xmax, std::size_t omax, std::size_t bmax, std::size_t mmax, std::size_t tint, const SizeList& vv, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto frvars = fudRepasSetVar;
	auto frder = fudRepasDerived;
	auto cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_up;
	auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_up;

	SizeList vv1;
	vv1.reserve(vv.size());
	{
		auto n = hh.dimension;
		auto& mvv = hh.mapVarInt();
		auto sh = hh.shape;
		auto pxx1 = hhx.arr;
		std::size_t sz = 0;
		SizeList xx;
		xx.reserve(n);
		for (std::size_t i = 0; i < n; i++)
		{
			xx.push_back(sz);
			sz += sh[i];
		}
		for (auto& w : vv)
		{
			auto p = mvv[w];
			auto e = sh[p];
			bool first = false;
			bool second = false;
			for (std::size_t i = 0; !second && i < e; i++)
			{
				bool nz = pxx1[xx[p] + i] > 0.0;
				second = second || (first && nz);
				first = first || nz;
			}
			if (second)
				vv1.push_back(w);
		}
	}
	std::size_t s = 0;
	auto xx1 = std::make_unique<SizeListList>();
	if (vv1.size() < 2)
		return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
	DoubleSizeListPairList xx;
	xx.reserve(omax * 10);
	if (!fr.layers.size())
	{
		auto t = cross(xmax, omax, tint, vv1, hh, hhx, hhrr, hhrrx);
		auto& xc = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xc->begin(), xc->end());
		xx1->reserve(xc->size());
		for (auto& pp : *xc)
			xx1->push_back(pp.second);
	}
	else
	{
		auto vvf = frvars(fr);
		vv1.reserve(vv1.size() + vvf->size());
		vv1.insert(vv1.end(), vvf->begin(), vvf->end());
		auto wwf = frder(fr);
		xx1->reserve(wwf->size());
		for (auto& w : *wwf)
			xx1->push_back(SizeList{ w });
	}
	while (xx1->size())
	{
		auto t = append(xmax, omax, tint, vv1, *xx1, hh, hhx, hhrr, hhrrx);
		auto& xa = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xa->begin(), xa->end());
		xx1->clear();
		xx1->reserve(xa->size());
		for (auto& pp : *xa)
			xx1->push_back(pp.second);
	}
	if (!xx.size())
		return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
	std::sort(xx.begin(), xx.end());
	xx1->reserve(bmax / mmax);
	std::size_t start = xx.size() - 1;
	std::size_t end = xx.size() > bmax / mmax ? xx.size() - bmax / mmax : 0;
	for (long long i = start; i >= (long long)end; i--)
		xx1->push_back(xx[i].second);
	return std::tuple<std::unique_ptr<SizeListList>, std::size_t>(std::move(xx1), s);
}

// parametersSystemsPartitionerMaxRollByMRepa_ui ::
//   Integer -> Integer -> Integer -> [VariableRepa] -> HistoryRepa -> HistoryRepa -> 
//   ([[[VariableRepa]]],Integer)
std::tuple<std::unique_ptr<SizeListListList>, std::size_t> Alignment::parametersSystemsPartitionerMaxRollByMRepa_ui(std::size_t mmax, std::size_t umax, std::size_t pmax, const SizeList& kk, const HistoryRepa& hh, const HistoryRepa& hhrr)
{
	auto hrred = [](double f, const HistoryRepa& hr, const SizeList& kk)
	{
		return setVarsHistoryRepasReduce_u(f, kk.size(), kk.data(), hr);
	};
	auto rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u;

	auto z = (double)hh.size;
	auto zr = (double)hhrr.size;
	auto aa = hrred(1.0, hh, kk);
	auto aarr = hrred(z/zr, hhrr, kk);
	double y1 = aa->facLn() - aarr->facLn();
	return rrvqqy(mmax, umax, pmax, *aa, *aarr, z, y1);
}


// parametersRollerMaximumRollExcludedSelfRepa_i ::
//   [[VariableRepa]] -> -> HistogramRepa -> -> HistogramRepa ->
//   ([[[Int]]],Integer)
std::tuple<std::unique_ptr<SizeListListList>, std::size_t> Alignment::parametersRollerMaximumRollExcludedSelfRepa_i(const SizeListList& pp, const HistogramRepa& aa, const HistogramRepa& aarr, double z)
{
	auto t = histogramRepaVecsRollMax(pp, aa, aarr, z);
	auto tt = std::move(std::get<0>(t));
	auto s = std::get<1>(t);
	auto tt1 = std::make_unique<SizeListListList>();
	bool any = false;
	for (auto& vv : *tt)
		if (*max_element(vv.begin(), vv.end()) < vv.size() - 1)
		{
			tt1->push_back(*tt);
			break;
		}
	return std::tuple<std::unique_ptr<SizeListListList>, std::size_t>(std::move(tt1), s);
}

// parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui ::
//   Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([(Double, [VariableRepa]],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(std::size_t wmax, std::size_t omax, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto fsize = fudRepasSize;
	auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u;

	auto l = fsize(fr);
	SizeUSet wwf;
	wwf.reserve(l);
	SizeList yy;
	yy.reserve(l);
	SizeSizeSetMap mm;
	bool first = true;
	for (auto& ll : fr.layers)
	{
		for (auto& tt : ll)
		{
			auto& w = tt->derived;
			yy.push_back(w);
			wwf.insert(w);
			auto& nn = mm[w];
			if (!first)
			{
				auto m = tt->dimension;
				auto& zz = tt->vectorVar;
				for (std::size_t i = 0; i < m; i++)
				{
					auto& v = zz[i];
					wwf.erase(v);
					nn.insert(v);
					auto it = mm.find(v);
					if (it != mm.end())
						nn.insert(it->second.begin(), it->second.end());
				}
			}
		}
		first = false;
	}
	SizeSizePairList cc;
	for (auto& pp : mm)
	{
		auto& w = pp.first;
		for (auto& v : pp.second)
			cc.push_back(SizeSizePair(w, v));
	}
	std::size_t s = 0;
	DoubleSizeListPairList xx;
	xx.reserve(omax * 10);
	SizeListList xx1;
	xx1.reserve(wwf.size());
	for (auto& w : wwf)
		xx1.push_back(SizeList{ w });
	while (xx1.size())
	{
		auto t = append(wmax, omax, cc, yy, xx1, hh, hhx, hhrr, hhrrx);
		auto& xa = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xa->begin(), xa->end());
		xx1.clear();
		xx1.reserve(xa->size());
		for (auto& pp : *xa)
			xx1.push_back(pp.second);
	}
	std::sort(xx.begin(), xx.end());
	auto xx2 = std::make_unique<DoubleSizeListPairList>();
	xx2->reserve(1);
	if (xx.size())
		xx2->push_back(xx.back());
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(xx2), s);
}

// parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_uip ::
//   Integer -> Integer -> Integer -> [VariableRepa] -> Fud ->
//   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->
//   ([(Double, [VariableRepa]],Integer)
std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t> Alignment::parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_uip(std::size_t wmax, std::size_t omax, std::size_t tint, const FudRepa& fr, const HistoryRepa& hh, const HistogramRepaRed& hhx, const HistoryRepa& hhrr, const HistogramRepaRed& hhrrx)
{
	auto fsize = fudRepasSize;
	auto append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_up;

	auto l = fsize(fr);
	SizeUSet wwf;
	wwf.reserve(l);
	SizeList yy;
	yy.reserve(l);
	SizeSizeSetMap mm;
	bool first = true;
	for (auto& ll : fr.layers)
	{
		for (auto& tt : ll)
		{
			auto& w = tt->derived;
			yy.push_back(w);
			wwf.insert(w);
			auto& nn = mm[w];
			if (!first)
			{
				auto m = tt->dimension;
				auto& zz = tt->vectorVar;
				for (std::size_t i = 0; i < m; i++)
				{
					auto& v = zz[i];
					wwf.erase(v);
					nn.insert(v);
					auto it = mm.find(v);
					if (it != mm.end())
						nn.insert(it->second.begin(), it->second.end());
				}
			}
		}
		first = false;
	}
	SizeSizePairList cc;
	for (auto& pp : mm)
	{
		auto& w = pp.first;
		for (auto& v : pp.second)
			cc.push_back(SizeSizePair(w, v));
	}
	std::size_t s = 0;
	DoubleSizeListPairList xx;
	xx.reserve(omax * 10);
	SizeListList xx1;
	xx1.reserve(wwf.size());
	for (auto& w : wwf)
		xx1.push_back(SizeList{ w });
	while (xx1.size())
	{
		auto t = append(wmax, omax, tint, cc, yy, xx1, hh, hhx, hhrr, hhrrx);
		auto& xa = std::get<0>(t);
		s += std::get<1>(t);
		xx.insert(xx.end(), xa->begin(), xa->end());
		xx1.clear();
		xx1.reserve(xa->size());
		for (auto& pp : *xa)
			xx1.push_back(pp.second);
	}
	std::sort(xx.begin(), xx.end());
	auto xx2 = std::make_unique<DoubleSizeListPairList>();
	xx2->reserve(1);
	if (xx.size())
		xx2->push_back(xx.back());
	return std::tuple<std::unique_ptr<DoubleSizeListPairList>, std::size_t>(std::move(xx2), s);
}

// systemsHistoryRepasApplicationsHistoryHistoryPartitionedRepa_u ::
//   HistoryRepa -> ApplicationRepa ->
//   (SystemRepa, HistoryRepa)
std::unique_ptr<HistoryRepa> Alignment::systemsHistoryRepasApplicationsHistoryHistoryPartitionedRepa_u(const HistoryRepa& hr, const ApplicationRepa& dr, SystemRepa& ur)
{
	auto frmul = historyRepasFudRepasMultiply_u;

	auto& llu1 = ur.listVarSizePair;	
	auto hr2 = frmul(hr, *dr.fud);
	
	auto sl = treesLeafElements(*dr.slices);
	auto m = sl->size();
	auto z = hr2->size;
	auto& mvv2 = hr2->mapVarInt();
	SizeList pvv2;
	for (auto v : *sl)
		pvv2.push_back(mvv2[v]);
	auto hr3 = std::make_unique<HistoryRepa>();
	std::size_t n3 = 2;
	std::size_t s0 = (m-1)/256 + 1;
	std::size_t s1 = m > 256 ? 256 : m;
	llu1.push_back(VarSizePair(std::make_shared<Variable>("partition0"), s0));
	llu1.push_back(VarSizePair(std::make_shared<Variable>("partition1"), s1));
	std::size_t v0 = llu1.size() - 2;
	std::size_t v1 = llu1.size() - 1;
	hr3->dimension = n3;
	hr3->vectorVar = new std::size_t[n3];
	auto vv3 = hr3->vectorVar;
	hr3->shape = new std::size_t[n3];
	auto sh3 = hr3->shape;
	vv3[0] = v0;
	sh3[0] = s0;
	vv3[1] = v1;
	sh3[1] = s1;
	hr3->size = z;
	hr3->evient = false;
	hr3->arr = new unsigned char[z*n3];
	auto rr2 = hr2->arr;
	auto rr3 = hr3->arr;
	for (std::size_t j = 0; j < z; j++)
		for (std::size_t i = 0; i < m; i++)
		{
			std::size_t u = rr2[pvv2[i]*z + j];
			if (u)
			{
				rr3[j] = (unsigned char)(i/256);
				rr3[z+j] = (unsigned char)(i%256);
				break;
			}
		}
	return hr3;
}

