#include "AlignmentAesonRepa.h"

#include <iostream>
#include <sstream>
#include <cstring>

using namespace Alignment;

// systemRepasPersistent :: SystemRepa -> SystemRepaPersistent
void Alignment::systemRepasPersistent(const SystemRepa& ur, std::ostream& out)
{
	auto& llu = ur.listVarSizePair;
	auto n = llu.size();
	out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	for (std::size_t i = 0; i < n; i++)
	{
		auto& p = llu[i];
		std::ostringstream sout;
		sout << *p.first;
		std::string s(sout.str());
		std::size_t l = s.size();
		out.write(reinterpret_cast<char*>(&l), sizeof(std::size_t));
		out.write(reinterpret_cast<char*>((char*)s.data()), l);
		std::size_t u = p.second;
		out.write(reinterpret_cast<char*>(&u), sizeof(std::size_t));
	}
}

// persistentsSystemRepa :: SystemRepaPersistent -> Maybe SystemRepa
std::unique_ptr<SystemRepa> Alignment::persistentsSystemRepa(std::istream& in, StrVarPtrMap& m)
{
	auto ur = std::make_unique<SystemRepa>();
	auto& llu = ur->listVarSizePair;
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	llu.reserve(n);
	for (std::size_t i = 0; i < n; i++)
	{
		std::size_t l;
		in.read(reinterpret_cast<char*>(&l), sizeof(std::size_t));
		std::string s(l,' ');
		in.read(reinterpret_cast<char*>((char*)s.data()), l);
		std::size_t u;
		in.read(reinterpret_cast<char*>(&u), sizeof(std::size_t));
		llu.push_back(VarSizePair(stringsVariable(s, m), u));
	}
	return ur;
}

// historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
void Alignment::historyRepasPersistent(const HistoryRepa& hr, std::ostream& out)
{
	auto n = hr.dimension;
	auto vv = hr.vectorVar;
	auto sh = hr.shape;
	auto z = hr.size;
	auto evient = hr.evient;
	auto rr = hr.arr;
	out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	for (std::size_t i = 0; i < n; i++)
	{
		out.write(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
		out.write(reinterpret_cast<char*>(&sh[i]), sizeof(std::size_t));
	}
	out.write(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&evient), 1);
	out.write(reinterpret_cast<char*>(hr.arr), z*n);
}

// persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::persistentsHistoryRepa(std::istream& in)
{
	auto hr = std::make_unique<HistoryRepa>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	hr->dimension = n;
	if (n)
	{
		hr->vectorVar = new std::size_t[n];
		auto vv = hr->vectorVar;
		hr->shape = new std::size_t[n];
		auto sh = hr->shape;
		for (std::size_t i = 0; i < n; i++)
		{
			in.read(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
			in.read(reinterpret_cast<char*>(&sh[i]), sizeof(std::size_t));
		}		
	}
	std::size_t z;
	in.read(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	in.read(reinterpret_cast<char*>(&hr->evient), 1);
	hr->size = z;
	if (z && n)
	{
		hr->arr = new unsigned char[z*n];
		in.read(reinterpret_cast<char*>(hr->arr), z*n);		
	}
	return hr;
}

// historyRepasPersistentInitial :: HistoryRepa -> Int -> HistoryRepaPersistentInitial
// assumes evient
void Alignment::historyRepasPersistentInitial(const HistoryRepa& hr, std::size_t y, std::ostream& out)
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
		out.write(reinterpret_cast<char*>(&sh[i]), sizeof(std::size_t));
	}
	out.write(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&y), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(hr.arr), y*n);
}

// persistentInitialsHistoryRepa :: HistoryRepaPersistentInitial -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::persistentInitialsHistoryRepa(std::istream& in)
{
	auto hr = std::make_unique<HistoryRepa>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	hr->dimension = n;
	if (n)
	{
		hr->vectorVar = new std::size_t[n];
		auto vv = hr->vectorVar;
		hr->shape = new std::size_t[n];
		auto sh = hr->shape;
		for (std::size_t i = 0; i < n; i++)
		{
			in.read(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
			in.read(reinterpret_cast<char*>(&sh[i]), sizeof(std::size_t));
		}		
	}
	std::size_t z;
	in.read(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	std::size_t y;
	in.read(reinterpret_cast<char*>(&y), sizeof(std::size_t));
	hr->size = z;
	hr->evient = true;
	if (z && n)
	{
		hr->arr = new unsigned char[z*n];
		in.read(reinterpret_cast<char*>(hr->arr), y*n);
	}
	return hr;
}

// historySparsesPersistent :: HistorySparse -> HistorySparsePersistent
void Alignment::historySparsesPersistent(const HistorySparse& hs, std::ostream& out)
{
	auto z = hs.size;
	auto dd = hs.vectorDimension;
	auto rr = hs.arr;
	std::size_t m = 0;
	for (std::size_t j = 0; j < z; j++)
	{
		auto k = dd[j];
		out.write(reinterpret_cast<char*>(&k), sizeof(std::size_t));
		for (std::size_t i = 0; i < k; i++)
		{
			out.write(reinterpret_cast<char*>(&rr[m]), sizeof(std::size_t));
			m++;
		}
	}
}

// persistentsHistorySparse :: HistorySparsePersistent -> Maybe HistorySparse
std::unique_ptr<HistorySparse> Alignment::persistentsHistorySparse(std::istream& in)
{
	std::size_t z = 0;
	SizeList dd1;
	SizeList rr1;
	while (true)
	{
		std::size_t k = 0;
		in.read(reinterpret_cast<char*>(&k), sizeof(std::size_t));
		if (in.eof())
			break;
		z++;
		dd1.push_back(k);
		for (std::size_t i = 0; i < k; i++)
		{
			std::size_t w = 0;
			in.read(reinterpret_cast<char*>(&w), sizeof(std::size_t));
			rr1.push_back(w);
		}
	}
	auto hs = std::make_unique<HistorySparse>();
	hs->size = z;
	if (z)
	{
		hs->vectorDimension = new std::size_t[z];
		auto dd = hs->vectorDimension;
		std::memcpy(dd, dd1.data(), dd1.size() * sizeof(std::size_t));
		hs->arr = new std::size_t[rr1.size()];
		std::memcpy(hs->arr, rr1.data(), rr1.size() * sizeof(std::size_t));		
	}
	return hs;
}

// historySparseArraysPersistent :: HistorySparseArray -> HistorySparseArrayPersistent
void Alignment::historySparseArraysPersistent(const HistorySparseArray& hr, std::ostream& out)
{
	auto n = hr.capacity;
	auto z = hr.size;
	auto rr = hr.arr;
	out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(hr.arr), z*n*sizeof(std::size_t));
}

// persistentsHistorySparseArray :: HistorySparseArrayPersistent -> Maybe HistorySparseArray
std::unique_ptr<HistorySparseArray> Alignment::persistentsHistorySparseArray(std::istream& in)
{
	auto hr = std::make_unique<HistorySparseArray>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	hr->capacity = n;
	std::size_t z;
	in.read(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	hr->size = z;
	if (z && n)
	{
		hr->arr = new std::size_t[z*n];
		in.read(reinterpret_cast<char*>(hr->arr), z*n*sizeof(std::size_t));		
	}
	return hr;
}

// historySparseArraysPersistentInitial :: HistorySparseArray -> HistorySparseArrayPersistentInitial
void Alignment::historySparseArraysPersistentInitial(const HistorySparseArray& hr, std::size_t y, std::ostream& out)
{
	auto n = hr.capacity;
	auto z = hr.size;
	auto rr = hr.arr;
	out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&y), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(hr.arr), y*n*sizeof(std::size_t));
}

// persistentInitialsHistorySparseArray :: HistorySparseArrayPersistentInitial -> Maybe HistorySparseArray
std::unique_ptr<HistorySparseArray> Alignment::persistentInitialsHistorySparseArray(std::istream& in)
{
	auto hr = std::make_unique<HistorySparseArray>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	hr->capacity = n;
	std::size_t z;
	in.read(reinterpret_cast<char*>(&z), sizeof(std::size_t));
	std::size_t y;
	in.read(reinterpret_cast<char*>(&y), sizeof(std::size_t));
	hr->size = z;
	if (z && n)
	{
		hr->arr = new std::size_t[z*n];
		in.read(reinterpret_cast<char*>(hr->arr), y*n*sizeof(std::size_t));
	}
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
		out.write(reinterpret_cast<char*>(&s), sizeof(std::size_t));
	}
	out.write(reinterpret_cast<char*>(&w), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(&u), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>(tr.arr), sz);
}

// persistentsTransformRepa :: TransformRepaPersistent -> Maybe TransformRepa
std::unique_ptr<TransformRepa> Alignment::persistentsTransformRepa(std::istream& in)
{
	auto tr = std::make_unique<TransformRepa>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	tr->dimension = n;
	if (n)
	{
		tr->vectorVar = new std::size_t[n];
		auto vv = tr->vectorVar;
		tr->shape = new std::size_t[n];
		auto sh = tr->shape;
		std::size_t sz = 1;
		for (std::size_t i = 0; i < n; i++)
		{
			in.read(reinterpret_cast<char*>(&vv[i]), sizeof(std::size_t));
			unsigned char s;
			in.read(reinterpret_cast<char*>(&s), sizeof(std::size_t));
			sz *= s;
			sh[i] = s;
		}
		std::size_t w;
		in.read(reinterpret_cast<char*>(&w), sizeof(std::size_t));
		tr->derived = w;
		std::size_t u;
		in.read(reinterpret_cast<char*>(&u), sizeof(std::size_t));
		tr->valency = u;
		tr->arr = new unsigned char[sz];
		in.read(reinterpret_cast<char*>(tr->arr), sz);		
	}
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
std::unique_ptr<FudRepa> Alignment::persistentsFudRepa(std::istream& in)
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
			auto tr = persistentsTransformRepa(in);
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
std::unique_ptr<Tree<HistoryRepaPtrFudRepaPtrPair>> persistentsHistoryRepaPtrFudRepaPtrPairTree(std::istream& in)
{
	auto zz = std::make_unique<Tree<HistoryRepaPtrFudRepaPtrPair>>();
	auto& ll = zz->_list;
	std::size_t l;
	in.read(reinterpret_cast<char*>(&l), sizeof l);
	ll.reserve(l);
	for (std::size_t i = 0; i < l; i++)
	{
		auto hr = persistentsHistoryRepa(in);
		auto fr = persistentsFudRepa(in);
		HistoryRepaPtrFudRepaPtrPair mr(std::move(hr), std::move(fr));
		auto zr = persistentsHistoryRepaPtrFudRepaPtrPairTree(in);
		ll.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mr, std::move(zr)));
	}
	return zz;
}

// persistentsDecompFudRepa :: DecompFudRepaPersistent -> Maybe DecompFudRepa
std::unique_ptr<DecompFudRepa> Alignment::persistentsDecompFudRepa(std::istream& in)
{
	auto dr = std::make_unique<DecompFudRepa>();
	auto& zz = dr->tree;
	auto& ll = zz._list;
	std::size_t l;
	in.read(reinterpret_cast<char*>(&l), sizeof l);
	ll.reserve(l);
	for (std::size_t i = 0; i < l; i++)
	{
		auto hr = persistentsHistoryRepa(in);
		auto fr = persistentsFudRepa(in);
		HistoryRepaPtrFudRepaPtrPair mr(std::move(hr), std::move(fr));
		auto zr = persistentsHistoryRepaPtrFudRepaPtrPairTree(in);
		ll.push_back(HistoryRepaPtrFudRepaPtrPairTreePtrPair(mr, std::move(zr)));
	}
	return dr;
}

void sizeTreesPersistent(const SizeTree& zz, std::ostream& out)
{
	auto& ll = zz._list;
	std::size_t l = ll.size();
	out.write(reinterpret_cast<char*>(&l), sizeof l);
	for (auto& pp : ll)
	{
		std::size_t x = pp.first;
		out.write(reinterpret_cast<char*>(&x), sizeof(std::size_t));
		if (pp.second)
			sizeTreesPersistent(*pp.second, out);
		else
			sizeTreesPersistent(SizeTree(), out);
	}
}

// applicationRepasPersistent :: ApplicationRepa -> ApplicationRepaPersistent
void Alignment::applicationRepasPersistent(const ApplicationRepa& dr, std::ostream& out)
{
	auto n = dr.substrate.size();
	out.write(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	out.write(reinterpret_cast<char*>((std::size_t*)dr.substrate.data()), n * sizeof(std::size_t));
	if (dr.fud)
		fudRepasPersistent(*dr.fud, out);
	else 
		fudRepasPersistent(FudRepa(), out);
	if (dr.slices)
		sizeTreesPersistent(*dr.slices, out);
	else
		sizeTreesPersistent(SizeTree(), out);
}

typedef std::pair<std::size_t, std::shared_ptr<SizeTree>> SizeSizeTreePtrPair;

std::unique_ptr<SizeTree> persistentsSizeTree(std::istream& in)
{
	auto zz = std::make_unique<SizeTree>();
	auto& ll = zz->_list;
	std::size_t l;
	in.read(reinterpret_cast<char*>(&l), sizeof l);
	ll.reserve(l);
	for (std::size_t i = 0; i < l; i++)
	{
		std::size_t x;
		in.read(reinterpret_cast<char*>(&x), sizeof(std::size_t));
		auto zr = persistentsSizeTree(in);
		ll.push_back(SizeSizeTreePtrPair(x, std::move(zr)));
	}
	return zz;
}

// persistentsApplicationRepa :: ApplicationRepaPersistent -> Maybe ApplicationRepa
std::unique_ptr<ApplicationRepa> Alignment::persistentsApplicationRepa(std::istream& in)
{
	auto dr = std::make_unique<ApplicationRepa>();
	std::size_t n;
	in.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));
	dr->substrate.resize(n);
	in.read(reinterpret_cast<char*>((std::size_t*)dr->substrate.data()), n * sizeof(std::size_t));
	dr->fud = std::move(persistentsFudRepa(in));
	dr->slices = std::move(persistentsSizeTree(in));
	return dr;
}


// decompFudSlicedRepasPersistent :: DecompFudSlicedRepa -> DecompFudSlicedRepaPersistent
void Alignment::decompFudSlicedRepasPersistent(const DecompFudSlicedRepa& dr, std::ostream& out)
{
	out.write(reinterpret_cast<char*>((std::size_t*)&dr.listFudRepaSizeExpected), sizeof(std::size_t));
	auto& ll = dr.fuds;
	std::size_t l = ll.size();
	out.write(reinterpret_cast<char*>(&l), sizeof(std::size_t));
	for (auto& fs : ll)
	{
		out.write(reinterpret_cast<char*>((std::size_t*)&fs.parent), sizeof(std::size_t));
		std::size_t c = fs.children.size();
		out.write(reinterpret_cast<char*>(&c), sizeof(std::size_t));
		for (auto sl : fs.children)
			out.write(reinterpret_cast<char*>(&sl), sizeof(std::size_t));
		std::size_t m = fs.fud.size();
		out.write(reinterpret_cast<char*>(&m), sizeof m);
		for (auto& tr : fs.fud)
			transformRepasPersistent(*tr,out);
	}
}

// persistentsDecompFudSlicedRepa :: DecompFudSlicedRepaPersistent -> Maybe DecompFudSlicedRepa
std::unique_ptr<DecompFudSlicedRepa> Alignment::persistentsDecompFudSlicedRepa(std::istream& in)
{
	auto dr = std::make_unique<DecompFudSlicedRepa>();
	in.read(reinterpret_cast<char*>(&dr->listFudRepaSizeExpected), sizeof(std::size_t));
	auto& ll = dr->fuds;
	std::size_t l;
	in.read(reinterpret_cast<char*>(&l), sizeof(std::size_t));
	ll.resize(l);
	for (auto& fs : ll)
	{
		in.read(reinterpret_cast<char*>(&fs.parent), sizeof(std::size_t));
		std::size_t c;
		in.read(reinterpret_cast<char*>(&c), sizeof(std::size_t));
		for (std::size_t i = 0; i < c; i++)
		{
			std::size_t sl;
			in.read(reinterpret_cast<char*>(&sl), sizeof(std::size_t));
			fs.children.push_back(sl);
		}		
		std::size_t m;
		in.read(reinterpret_cast<char*>(&m), sizeof(std::size_t));
		for (std::size_t i = 0; i < m; i++)
		{
			auto tr = persistentsTransformRepa(in);
			fs.fud.push_back(std::move(tr));
		}
		dr->fudRepasSize += m;
	}
	return dr;
}
