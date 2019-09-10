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



    }


}
