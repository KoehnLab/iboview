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

#ifndef CT_IO_H
#define CT_IO_H

#include <string>
#include <stdexcept>
#include <sstream>

#include "CxDefs.h" // for assert
#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
   #include "CxScalarTypes.h"
#endif

#include "CxPodArray.h"
#include "format.h" // brings in namespace "fmt", including its format function (which makes std::strings)

namespace ct {

// extern char const
//    *pResultFmt,
//    *pResultFmtAnnotated,
//    *pResultFmtI,
//    *pResultFmtIAnnotated,
//    *pTimingFmt,
//    *pTimingPerFmt;

   std::string FmtResult(fmt::BasicStringRef<char> Name, double Value, fmt::BasicStringRef<char> Annotation = "");
   std::string FmtCount(fmt::BasicStringRef<char> Name, ptrdiff_t Value, fmt::BasicStringRef<char> Annotation = "");
   std::string FmtTiming(fmt::BasicStringRef<char> Name, double fTimeInSeconds, size_t nTasks = size_t(-1));
} // namespace ct


namespace ct {
   extern std::ostream
      &xerr, &xout;
   extern int
      Verbosity;
}

namespace ct {

// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro(std::ostream &out, const std::string &Name, const std::string &Version="");

// will load an entire file into a TArray<char>. Adds a 0-terminator at end.
// Returns false if failed.
// if 0 != pFileLength, *pFileLength will receive the number of loaded chars.
bool LoadFileIntoMemory(TArray<char> &pFileContent,
        std::string const &FileName, unsigned *pFileLength = 0);


// // error in some sort of input -- (e.g., syntax, brken file formats).
// class FInputError : public std::runtime_error
// {
// public:
//    typedef std::runtime_error
//       FBase;
//    explicit FInputError(std::string const &Reason);
// };
//
// ^-- moved to CxParse1.h

// // some unsupported calculation or combination of options was asked for (e.g., a tau-type meta-GGA with DFXC=1)
// class FUnsupportedError : public std::runtime_error
// {
// public:
//    typedef std::runtime_error
//       FBase;
//    explicit FUnsupportedError(std::string const &Reason);
// };
//
// ^-- moved to CxParse1.h


class FLogError : public std::runtime_error
{
public:
   typedef std::runtime_error
      FBase;
   explicit FLogError(std::string const &Reason);
};


// error in converging some sort of calculation for which there is no workaround declared
class FConvergenceError : public std::runtime_error
{
public:
   typedef std::runtime_error
      FBase;
   explicit FConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt = -1, double fFinalResidual = -1., double fFinalEnergy = -1);
};
std::string FormatConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt = -1, double fFinalResidual = -1., double fFinalEnergyChange = -1);



// a class encapsulating output streams and associated messaging.
//
// Main point about not simply using std::cout or some other global thing
// is that computations might be done on distributed machines or in some
// other non-local/non-sequential way.
class FLog
{
protected:
public:
   fmt::MemoryWriter
      w;
   typedef char
      Char;

   FLog();
   virtual ~FLog();

   enum FWriteResultFlags{
      OUTPUT_Tabulate = 0x02 // used with WriteResult(...). Marks an output meant for tabulation by an external program
   };

   void Write(fmt::BasicStringRef<Char> format) {
      w << format << '\n';
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void WriteNoNl(fmt::BasicStringRef<Char> format) {
      w << format;
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void Write(fmt::BasicStringRef<Char> format, fmt::ArgList args) {
      fmt::BasicFormatter<Char>(w).format(format, args);
      w << '\n';
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   void WriteNoNl(fmt::BasicStringRef<Char> format, fmt::ArgList args) {
      fmt::BasicFormatter<Char>(w).format(format, args);
      if (GetFlag(IOFLAG_FlushAfterWrite))
         Flush();
   }
   FMT_VARIADIC_VOID(Write, fmt::BasicStringRef<Char>)
   FMT_VARIADIC_VOID(WriteNoNl, fmt::BasicStringRef<Char>)

   virtual void Flush() = 0;

   virtual void EmitWarning(fmt::BasicStringRef<Char> Name);
   virtual void EmitError(fmt::BasicStringRef<Char> Name);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteResult(fmt::BasicStringRef<Char> Name, double fValue, unsigned Flags);
#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue);
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, unsigned Flags);
#endif // USE_GENERAL_SCALARS
   void WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value);
   void WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value, fmt::BasicStringRef<Char> Annotation);
   void WriteInfoExpf(fmt::BasicStringRef<Char> Name, double fValue);
   void WriteInfoExpf(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation);
   void WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue);
   void WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue, fmt::BasicStringRef<Char> Annotation);
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds);
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, fmt::BasicStringRef<Char> Annotation);
   // print relative timing per task.
   void WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, size_t nTasks);

   void WriteLine() { w << "\n"; };
   void WriteProgramIntro(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Version="");

   // this is our mechanism for communicating with the outside world (the other side of the log).
   // GetStatus may return something else than STATUS_Okay, for example, if an abort request was
   // placed on the other side by the user.
   // Default implementation does nothing (returns STATUS_Okay).
   enum FRunStatus {
      STATUS_Okay = 0,
      STATUS_AbortSignalled = 1
   };
   virtual FRunStatus GetStatus();
   // call GetStatus; throw FLogError if status says something differnt from STATUS_Okay.
   void CheckStatus();
   bool StatusOkay() { return GetStatus() == STATUS_Okay; };

public:
   enum FIoFlag {
      IOFLAG_FlushAfterWrite = 0x0001,
      IOFLAG_Max = IOFLAG_FlushAfterWrite
   };

   inline bool GetFlag(FIoFlag Flag) {
      assert(Flag <= IOFLAG_Max);
      return bool(m_IoFlags & Flag);
   }
   void SetFlag(FIoFlag Flag, bool NewState) {
      assert(Flag <= IOFLAG_Max);
      if (NewState)
         m_IoFlags |= Flag;
      else
         m_IoFlags &= ~Flag;
   }
   // meant for flags with more than two states; e.g., could be used for
   // local indent level.
   typedef unsigned
      FIoFlagValue;
   inline FIoFlagValue GetFlagValue(FIoFlag Flag) {
      assert(Flag <= IOFLAG_Max);
      return m_IoFlags & FIoFlagValue(Flag);
   }
   void SetFlagValue(FIoFlag Flag, FIoFlagValue NewValue) {
      assert(Flag <= IOFLAG_Max);
      m_IoFlags = (m_IoFlags & (~Flag)) | NewValue;
   }
protected:
   FIoFlagValue
      m_IoFlags;
private:
   void operator = (FLog const &); // not implemented
   FLog(FLog const &); // not implemented.
};


class FLogIoFlagOverride
{
public:
   FLogIoFlagOverride(FLog *pLog, FLog::FIoFlag Flag, bool NewState) {
      m_pLog = pLog;
      m_Flag = Flag;
      m_NewState = NewState;
      if (m_pLog) {
         m_OriginalState = m_pLog->GetFlag(m_Flag);
         m_pLog->SetFlag(m_Flag, m_NewState);
      } else {
         m_OriginalState = false;
      }
   }
   ~FLogIoFlagOverride() {
      if (m_pLog) {
         m_pLog->SetFlag(m_Flag, m_OriginalState);
         if (m_Flag == FLog::IOFLAG_FlushAfterWrite && m_NewState == false && m_OriginalState == true)
            m_pLog->Flush();
      }
   }
protected:
   FLog
      *m_pLog;
   FLog::FIoFlag
      m_Flag;
   bool
      m_OriginalState, m_NewState;
};


/// log of which the output goes to a generic std::ostream.
class FLogStdStream : public FLog
{
public:
   explicit FLogStdStream(std::ostream &TargetStream);
   void Flush(); // override
protected:
   std::ostream
      &m_TargetStream;
};



/// log which keeps the text in a memory object, without actually emitting it
/// anywhere actively.
class FMemoryLog : public FLogStdStream
{
public:
   FMemoryLog();
   ~FMemoryLog();

   enum FGetTextFlags {
      GETTEXT_RetainLog = 0x0000,
      // if set, the log buffer will be cleared out after GetText() is called.
      GETTEXT_ClearLog = 0x0001
   };

   /// assemble a single std::string containing the entire log content, and return it.
   /// Note: this does *not* clear the current text.
   std::string GetText(unsigned Flags = GETTEXT_ClearLog);
   /// flush and clera out the entire controlled log text
   void Clear();
protected:
   std::stringstream
      m_TextStream;
};


struct FFileLocatorImpl;



/// A support structure holding functionality for locating physical files (e.g.,
/// data files from the application, such as basis set libraries) relative to
/// the application directory, rather than the current working directory.
///
/// Used to allow executing programs like MicroScf to locate its auxiliary data
/// if executed via command line from other directories.
///
/// Note:
/// - The BasePath of the file search environment has to be globally set once
///   via SetBasePath().
/// - Typically, one would use the location of the executable:
///
///    FFileLocator::SetBasePath(GetExePath()); // with GetExePath() from CxOsInt.h
///
/// - Afterwards, files and additional search directories can be specified
///   relative to this base path. E.g., in AddSearchPath, '$BASEPATH' is replaced
///   by the base path of the application
/// - The implementation in principle supports additional path substitutions
///   (e.g., for HOME directories or CONFIG/PROFILE directories). But there
///   is no interface to those yet, because it is not yet used in my programs.
struct FFileLocator
{

   FFileLocator();
   ~FFileLocator();
   /// Store application base path in a global variable subsequently accessed
   /// in newly constructed FFileLocator objects.
   ///
   /// Notes:
   /// - Global to allow multiple instanciations of the object in lower level
   ///   structures of the code.
   /// - Typically this would be the application exe directory on linux, or
   ///   the application base directory on Windows.
   /// - Must be called before first invocation of FFileLocator ("." is okay
   ///   if no additional info is available).
   static void SetBasePath(std::string const &BasePath);

   static std::string JoinPath(std::string const &Path, std::string const &FileName);
   static std::string BaseName(std::string const &PathAndFileName);
   static std::string PathName(std::string const &PathAndFileName);

   enum FFindFileFlags {
      FINDFILE_RaiseErrorIfNotFound = 0x0001,
      // if set, the substitutions in m_Substitutions will not be applied
      // to the file name. These are otherwise used to register substitutions
      // like $BASEPATH or $HOME.
      FINDFILE_NoSubstitutions = 0x0002,
      FINDFILE_Default = FINDFILE_RaiseErrorIfNotFound
   };

   std::string FindFile(std::string const &FileName, unsigned Flags=FINDFILE_Default) const;
   /// note: '$BASEPATH' is replaced by the base path of the application.
   void AddSearchPath(std::string const &PathName);

   /// adds a substitution of $HOME and ~/ (at start) to the current user's
   /// home directory, if possible.
   void AddHomePath();

   /// imports environment variable sVarName as a substitution into *this.
   void ImportEnvironmentVar(std::string const &VarName);


   void PushState();
   void PopState();
protected:
   FFileLocatorImpl
      *p;
   void operator = (FFileLocator const &); // not implemented
   FFileLocator(FFileLocator const &); // not implemented
};


/// On construction, saves the current state of a FFileLocator object;
/// on destruction, restores the original state. This can be used to
/// savely add local search paths (e.g., when processing input files
/// in different directories) without compromising the structure once
/// done.
struct FFileLocatorStateGuard
{
   FFileLocatorStateGuard(FFileLocator &FileLocator_)
      : m_FileLocator(FileLocator_)
   {
      m_FileLocator.PushState();
   }
   ~FFileLocatorStateGuard()
   {
      m_FileLocator.PopState();
   }
protected:
   FFileLocator
      &m_FileLocator;
private:
   void operator = (FFileLocatorStateGuard const&); // not implemented.
   FFileLocatorStateGuard(FFileLocatorStateGuard const&); // not implemented.
};



template<class FStringLike>
struct FRepeatStrProxy {
   FStringLike const &m_str;
   size_t m_count;
   explicit FRepeatStrProxy(FStringLike const &str, size_t count) : m_str(str), m_count(count) {}
};

template<class FStringLike>
std::ostream &operator << (std::ostream &out, FRepeatStrProxy<FStringLike> const &inst) {
   for (size_t i = 0; i != inst.m_count; ++ i)
      out << inst.m_str;
   return out;
}

template<class FStringLike>
FRepeatStrProxy<FStringLike> repeat_str(FStringLike const &str, size_t count) {
   return FRepeatStrProxy<FStringLike>(str, count);
}


} // namespace ct.

#endif // CT_IO_H
