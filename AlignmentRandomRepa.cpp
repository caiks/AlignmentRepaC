#include "AlignmentRandomRepa.h"

#include <stdlib.h>

using namespace Alignment;

// historyRepasShuffle_u :: HistoryRepa -> Int -> HistoryRepa
std::unique_ptr<HistoryRepa> Alignment::historyRepasShuffle_u(const HistoryRepa& hr, unsigned s)
{
    auto n = hr.dimension;
    auto vv = hr.vectorVar;
    auto sh = hr.shape;
    auto z = hr.size;
    auto rr = hr.arr;
    auto hr1 = std::make_unique<HistoryRepa>();
    if (!rr)
	return hr1;
    hr1->dimension = n;
    hr1->vectorVar = new std::size_t[n];
    auto vv1 = hr1->vectorVar;
    hr1->shape = new unsigned char[n];
    auto sh1 = hr1->shape;
    for (std::size_t i = 0; i < n; i++)
    {
	vv1[i] = vv[i];
	sh1[i] = sh[i];
    }
    hr1->size = z;
    hr1->arr = new unsigned char[z*n];
    auto rr1 = hr1->arr;
    memcpy(rr1, rr, z*n);
    srand(s);
    for (std::size_t i = 0; i < n; i++)
	for (std::size_t j = z-1; j != 0; j--)
	{
	    auto k = rand() % (j+1);
	    auto x = rr1[k*n + i];
	    rr1[k*n + i] = rr1[j*n + i];
	    rr1[j*n + i] = x;
	}
    return hr1;
}


