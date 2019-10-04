#include "AlignmentAesonRepa.h"

#include <iostream>
#include <sstream>

using namespace Alignment;

// historyRepasPersistent :: HistoryRepa -> HistoryRepaPersistent
void Alignment::historyRepasPersistent(const HistoryRepa& hr, std::ostream& out)
{
    auto& vv = hr.vectorVar;
    auto n = vv.size();
    out.write(reinterpret_cast<char*>(&n), sizeof n);
    for (std::size_t i = 0; i < n; i++)
    {
	std::stringstream str1;
	str1 << vv[i];
	auto s = str1.str();
	auto l = s.size();
	out.write(reinterpret_cast<char*>(&l), sizeof l);
	out << s;
    }
    auto z = hr.size;
    out.write(reinterpret_cast<char*>(&z), sizeof z);
    for (std::size_t i = 0; i < n; i++)
    {
	auto s = hr.shape[i];
	out.write(reinterpret_cast<char*>(&s), sizeof s);
    }
    out.write(reinterpret_cast<char*>(hr.arr), z*n);
}

// persistentsHistoryRepa :: HistoryRepaPersistent -> Maybe HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::persistentsHistoryRepa(std::istream& in, StrVarPtrMap& vm)
{
    auto hr = std::make_unique<HistoryRepa>();
    auto& vv = hr->vectorVar;
    std::size_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof n);
    for (std::size_t i = 0; i < n; i++)
    {
	std::size_t l;
	in.read(reinterpret_cast<char*>(&l), sizeof l);
	std::string s(l,'\0');
	in.read(&s[0],l);
	vv.push_back(*stringsVariable(s, vm));
    }
    std::size_t z;
    in.read(reinterpret_cast<char*>(&z), sizeof z);
    hr->size = z;
    auto& sh = hr->shape;
    for (std::size_t i = 0; i < n; i++)
    {
	std::size_t s;
	in.read(reinterpret_cast<char*>(&s), sizeof s);
	sh.push_back(s);
    }
    hr->arr = new unsigned char[z*n];
    in.read(reinterpret_cast<char*>(hr->arr), z*n);
    return hr;
}
