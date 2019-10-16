﻿#include "AlignmentAesonRepa.h"

#include <iostream>
#include <sstream>

using namespace Alignment;

// historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
void Alignment::historyRepasPersistent(const HistoryRepa& hr, std::ostream& out)
{
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto sh = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
    for (std::size_t i = 0; i < n; i++)
    {
	out.write(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&sh[i]), 1);
    }
    out.write(reinterpret_cast<char*>(&z), sizeof(std::size_t));
    out.write(reinterpret_cast<char*>(hr.arr), z*n);
}

// persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::persistentsHistoryRepa(std::istream& in, StrVarPtrMap& vm)
{
    auto hr = std::make_unique<HistoryRepa>();
    std::size_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
    hr->dimension = n;
    hr->vectorVar = new std::size_t[n];
    auto vv = hr->vectorVar;
    hr->shape = new unsigned char[n];
    auto sh = hr->shape;
    for (std::size_t i = 0; i < n; i++)
    {
	in.read(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
	in.read(reinterpret_cast<char*>(&sh[i]), 1);
    }
    std::size_t z;
    in.read(reinterpret_cast<char*>(&z), sizeof(std::size_t));
    hr->size = z;
    hr->arr = new unsigned char[z*n];
    in.read(reinterpret_cast<char*>(hr->arr), z*n);
    return hr;
}

// transformRepasPersistent :: TransformRepa -> TransformRepaPersistent
void Alignment::transformRepasPersistent(const TransformRepa& tr, std::ostream& out)
{
    auto n = tr.dimension;
    auto vv = tr.vectorVar;
    auto w = tr.derived;
    auto u = tr.valency;
    auto sh = tr.shape;
    auto rr = tr.arr;
    out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	out.write(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
	auto s = sh[i];
	sz *= s;
	out.write(reinterpret_cast<char*>(&s), 1);
    }
    out.write(reinterpret_cast<char*>(&w), sizeof(std::size_t));
    out.write(reinterpret_cast<char*>(&u), 1);
    out.write(reinterpret_cast<char*>(tr.arr), sz);
}

// persistentsTransformRepa :: TransformRepaPersistent -> Maybe TransformRepa
std::unique_ptr<TransformRepa> Alignment::persistentsTransformRepa(std::istream& in, StrVarPtrMap& vm)
{
    auto tr = std::make_unique<TransformRepa>();
    std::size_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
    tr->dimension = n;
    tr->vectorVar = new std::size_t[n];
    auto vv = tr->vectorVar;
    tr->shape = new unsigned char[n];
    auto sh = tr->shape;
    std::size_t sz = 1;
    for (std::size_t i = 0; i < n; i++)
    {
	in.read(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
	std::size_t s;
	in.read(reinterpret_cast<char*>(&s), 1);
	sz *= s;
	sh[i] = s;
    }
    std::size_t w;
    in.read(reinterpret_cast<char*>(&w), sizeof(std::size_t));
    tr->derived = w;
    unsigned char u;
    in.read(reinterpret_cast<char*>(&u), 1);
    tr->valency = u;
    tr->arr = new unsigned char[sz];
    in.read(reinterpret_cast<char*>(tr->arr), sz);
    return tr;
}

// fudRepasPersistent :: FudRepa -> FudRepaPersistent
void Alignment::fudRepasPersistent(const FudRepa& fr, std::ostream& out)
{
    auto& ll = fr.layers;
    std::size_t l = ll.size();
    out.write(reinterpret_cast<char*>(&l), sizeof l);
    for (auto& mm : ll)
    {
	std::size_t m = mm.size();
	out.write(reinterpret_cast<char*>(&m), sizeof m);
	for (auto& tr : mm)
	    transformRepasPersistent(*tr,out);
    }
}

// persistentsFudRepa :: FudRepaPersistent -> Maybe FudRepa
std::unique_ptr<FudRepa> Alignment::persistentsFudRepa(std::istream& in, StrVarPtrMap& vm)
{
    auto fr = std::make_unique<FudRepa>();
    auto& ll = fr->layers;
    std::size_t l;
    in.read(reinterpret_cast<char*>(&l), sizeof l);
    ll.resize(l);
    for (auto& mm : ll)
    {
	std::size_t m;
	in.read(reinterpret_cast<char*>(&m), sizeof m);
	for (std::size_t i = 0; i < m; i++)
	{
	    auto tr = persistentsTransformRepa(in, vm);
	    mm.push_back(std::move(tr));
	}
    }
    return fr;
}

// historyRepaPtrFudRepaPtrPairTreesPersistent :: Tree HistoryRepaPtrFudRepaPtrPair -> DecompFudRepaPersistent
void historyRepaPtrFudRepaPtrPairTreesPersistent(const Tree<HistoryRepaPtrFudRepaPtrPair>& zz, std::ostream& out)
{
    auto& ll = zz._list;
    std::size_t l = ll.size();
    out.write(reinterpret_cast<char*>(&l), sizeof l);
    for (auto& pp : ll)
    {
	historyRepasPersistent(*pp.first._state, out);
	fudRepasPersistent(*pp.first._fud, out);
	historyRepaPtrFudRepaPtrPairTreesPersistent(*pp.second, out);
    }
}

// decompFudRepasPersistent :: DecompFudRepa -> DecompFudRepaPersistent
void Alignment::decompFudRepasPersistent(const DecompFudRepa& dr, std::ostream& out)
{
    auto& zz = dr.tree;
    auto& ll = zz._list;
    std::size_t l = ll.size();
    out.write(reinterpret_cast<char*>(&l), sizeof l);
    for (auto& pp : ll)
    {
	historyRepasPersistent(*pp.first._state, out);
	fudRepasPersistent(*pp.first._fud, out);
	historyRepaPtrFudRepaPtrPairTreesPersistent(*pp.second, out);
    }
}

typedef std::shared_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> HistoryRepaPtrFudRepaPtrPairTreePtr;
typedef std::pair<HistoryRepaPtrFudRepaPtrPair, HistoryRepaPtrFudRepaPtrPairTreePtr> HistoryRepaPtrFudRepaPtrPairTreePtrPair;

// persistentsHistoryRepaPtrFudRepaPtrPairTree :: Tree HistoryRepaPtrFudRepaPtrPair -> Maybe DecompFudRepa
std::unique_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> persistentsHistoryRepaPtrFudRepaPtrPairTree(std::istream& in, StrVarPtrMap& vm)
{
    auto zz = std::make_unique<Tree<HistoryRepaPtrFudRepaPtrPair>>();
    auto& ll = zz->_list;
    std::size_t l;
    in.read(reinterpret_cast<char*>(&l), sizeof l);
    ll.reserve(l);
    for (std::size_t i = 0; i < l; i++)
    {
	auto hr = persistentsHistoryRepa(in, vm);
	auto fr = persistentsFudRepa(in, vm);
	HistoryRepaPtrFudRepaPtrPair mr(std::move(hr), std::move(fr));
	auto zr = persistentsHistoryRepaPtrFudRepaPtrPairTree(in, vm);
	ll.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mr, std::move(zr)));
    }
    return zz;
}

// persistentsDecompFudRepa :: DecompFudRepaPersistent -> Maybe DecompFudRepa
std::unique_ptr<DecompFudRepa> Alignment::persistentsDecompFudRepa(std::istream& in, StrVarPtrMap& vm)
{
    auto dr = std::make_unique<DecompFudRepa>();
    auto& zz = dr->tree;
    auto& ll = zz._list;
    std::size_t l;
    in.read(reinterpret_cast<char*>(&l), sizeof l);
    ll.reserve(l);
    for (std::size_t i = 0; i < l; i++)
    {
	auto hr = persistentsHistoryRepa(in, vm);
	auto fr = persistentsFudRepa(in, vm);
	HistoryRepaPtrFudRepaPtrPair mr(std::move(hr), std::move(fr));
	auto zr = persistentsHistoryRepaPtrFudRepaPtrPairTree(in, vm);
	ll.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mr, std::move(zr)));
    }
    return dr;
}


