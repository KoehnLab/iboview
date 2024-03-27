/* Copyright (c) 2015  Gerald Knizia
 * 
 * This file is part of the IboView program (see: http://www.iboview.org)
 * 
 * IboView is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IboView is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IboView (LICENSE). If not, see http://www.gnu.org/licenses/
 * 
 * Please see IboView documentation in README.txt for:
 * -- A list of included external software and their licenses. The included
 *    external software's copyright is not touched by this agreement.
 * -- Notes on re-distribution and contributions to/further development of
 *    the IboView software
 */

#include <ctime>
#include "CxTiming.h"

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace ct {

static double FindTimerInverseFrequency()
{
	LARGE_INTEGER iFreq;
	if (!QueryPerformanceFrequency(&iFreq))
		throw std::runtime_error("QueryPerformanceFrequency failed somehow.");

	return 1.0 / (double)iFreq.QuadPart;
}

double Second()
{
	static double fInvFreq = FindTimerInverseFrequency();
	LARGE_INTEGER iFreq;
	QueryPerformanceCounter(&iFreq);
	return fInvFreq * iFreq.QuadPart;
}

} // namespace ct

#else

#include <sys/time.h>

namespace ct {

double Second()
{
	timeval
		tv;
	gettimeofday( &tv, 0 );
	return tv.tv_sec + tv.tv_usec/1000000.0;
}

// get time stamp counter from current CPU. Fragile!
// Use only for debugging, if at all!
unsigned long long dbgGetTimeStampCounter()
{
	unsigned a, d;
	asm volatile("rdtsc" : "=a" (a), "=d" (d));

	return (((unsigned long long)a) | (((unsigned long long)d) << 32));
}

} // namespace ct

#endif


#include <string.h> // for memset
#include <math.h> // for fabs
#include "format.h"
#include "CxIo.h"

namespace ct {

FTimerSet::FTimerSet()
{
	m_CurrentDetailLevel = 0;
	m_TargetDetailLevel = -1;
	Reset();
}

void FTimerSet::Reset()
{
	memset(&Entries[0], 0, sizeof(Entries));
}

void FTimerSet::PrintReport(std::ostream &out, double ThrPrint)
{
	size_t nOmitted = 0;
	double fOmitted = 0.;
	for (size_t iClockId = 0; iClockId < nMaxEntries(); ++ iClockId) {
		FEntry &e = Entries[iClockId];
		if (e.iInvk != 0) {
			if (fabs((double)e.fTime) > ThrPrint) {
				out << FmtTiming(e.pDesc, e.fTime);
			} else {
				nOmitted += 1;
				fOmitted += (double)e.fTime;
			}
		}
	}
	if (nOmitted != 0) {
		std::stringstream str;
		str << fmt::format("{} entries below {:.3f} s", nOmitted, ThrPrint);
		out << FmtTiming(str.str(), fOmitted);
	}

}



}





// kate: space-indent off; tab-indent on; indent-width 4;
