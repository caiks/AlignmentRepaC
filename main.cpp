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
    if (true)
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

}
