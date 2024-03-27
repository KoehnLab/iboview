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

#ifndef _WIN32
    #include <unistd.h> // for getpid and pid_t. Anyone got a better idea to obtain the base path name?
#else
    #include <stdlib.h> // for malloc/free/wcstombs
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include <windows.h>
#endif

#include "CxDefs.h"
#include "CxOsInt.h"
#include "format.h"
#include <stdexcept>
#include <cstdlib> // for wcstombs
#include <cstring> // for len()
#include <algorithm> // for std::min()

namespace ct {

#ifdef _WIN32

// #ifdef UNICODE
//    #error "Sorry, this stuff assumes 8byte chars. Please compile with ASCII API bindings."
// #endif


std::string FmtLastError()
{
   // Retrieve the system error message for the last-error code
   LPTSTR lpMsgBuf;
   DWORD dw = GetLastError();

   FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
      NULL,
      dw,
      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      (LPTSTR) &lpMsgBuf,
      0, NULL );

   // convert to std::string.
   std::string
      s(lpMsgBuf, lpMsgBuf + lstrlen(lpMsgBuf));
   LocalFree(lpMsgBuf);
   return s;
}

std::string GetExePath()
{
   //assert_rt(sizeof(TCHAR) == sizeof(char)); // this may become a serious problem on non-european systems...
   //return std::string("C:\\Linux\\dev\\microscf.20180116");
   DWORD
      nBuf = 512,
      nChars;
   TCHAR
      *pwBuf = (TCHAR*)::malloc(sizeof(TCHAR) * nBuf);
   char
      *pcBuf = (char*)::malloc(sizeof(char) * nBuf);
   nChars = GetModuleFileName(0, pwBuf, nBuf);

   if (nChars == 0)
      throw std::runtime_error("Failed to obtain current path name. GetLastError() says: " + FmtLastError() );

   if (sizeof(TCHAR) == sizeof(wchar_t)) {
      std::wcstombs(pcBuf, (wchar_t*)pwBuf, std::min(size_t(nBuf), size_t(1 + std::wcslen((wchar_t*)pwBuf))));
      // ^-- get a weird compiler error with the std version. maybe it is #defined to something by default. use global one instead, without namespace
    //  size_t nw = 0;
  //    while (pwBuf[nw] != 0 && nw < nBuf) nw += 1;
//      wcstombs(pcBuf, (wchar_t*)pwBuf, std::min(nBuf, 1 + nw));
      //wcstombs(pcBuf, (wchar_t*)pwBuf, std::min(nBuf, 1 + wcslen((wchar_t*)pwBuf)));
   } else if (sizeof(TCHAR) == sizeof(char)) {
      std::memcpy(pcBuf, pwBuf, sizeof(char) * std::min(size_t(nBuf), size_t(1 + std::strlen((char*)pwBuf))));
   } else {
      throw std::runtime_error("GetExePath(): broken character format.");
   }

   // convert to std::string.
   std::string
      s(pcBuf, pcBuf + nChars);

   ::free(pcBuf);
   ::free(pwBuf);

   return s;
}


#else

// return the absolute path of the currently executed binary on
// a linux machine. Required to access our compiled algorithm files,
// because working directory is usually changed to some place else.
// if anyone got a better idea on how to do it ... [what follows now
// appears to be not some kind of obscure hack, but the actual linux
// standard way of doing this]
std::string GetExePath()
{
    pid_t
        ProcessId = getpid();
    std::string
        ExeSymLink(fmt::format("/proc/{}/exe", ProcessId)),
        ExeName(513,'\0');
    readlink(ExeSymLink.c_str(), &ExeName[0], ExeName.size() - 1);
    int
        iEnd = ExeName.find('\0');
    while ( iEnd > 0 && ExeName[iEnd] != '/' && ExeName[iEnd] != '\\' && ExeName[iEnd] != ':' )
        -- iEnd;
    ++ iEnd;
    return ExeName.substr(0, iEnd);
}

#endif

} // namespace ct
