#include "ReadAscii.h"
#include "Streams.h"

using namespace std;

string GetOmicronFilePattern(const string aChannelName, const int aTimeStart, const int aTimeEnd);
string GetOmicronFilePatternFromHpss(const string aChannelName, const int aTimeStart, const int aTimeEnd);
