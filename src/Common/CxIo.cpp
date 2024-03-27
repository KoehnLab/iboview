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

#include "CxTypes.h" // for intrusive ptr
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cctype> // for tolower()
#include <cstdlib> // for getenv

#include <list> // used in FFileLocator
#include <map> // used in FFileLocator
#include <set> // used in FFileLocator
#include "CxParse1.h" // also used in FFileLocator...

#ifdef QT_CORE_LIB
   #include <QString>
   #include <QByteArray>
   #include <QFile>
   #include <QApplication>
   #include <QClipboard>
#endif

// #define USE_INDENT_STREAM
#ifdef USE_INDENT_STREAM
   #define INDENTSTREAM_IMPL
   #include "CxIndentStream.h"
   namespace ct {
      static fmt::FIndentStream1
         s_OutputStreamAdp(std::cout,true,1);
      std::ostream
         &xout = s_OutputStreamAdp.stream,
         &xerr = std::cerr;
   } // namespace ct
#else
   namespace ct {
      std::ostream
         &xout = std::cout,
         &xerr = std::cerr;
   } // namespace ct
#endif // USE_INDENT_STREAM

#include "CxIo.h"


namespace ct {

bool g_FileLocator_DebugPrints = false;

// char const
//    *pResultFmt = " %-32s%18.12f\n",
//    *pResultFmtAnnoted = " %-32s%18.12f  (%s)\n",
//    *pResultFmtI = " %-32s%5i\n",
//    *pResultFmtIAnnoted = " %-32s%5i  (%s)\n",
//    *pTimingFmt = " Time for %s:%40t%10.2f sec\n",
//    *pTimingPerFmt = " Time per %s:%40t%10.2f sec  (n=%i)\n";
//    *pTimingFmt = " Time for %s:%40t%10.4f sec\n",
//    *pTimingPerFmt = " Time per %s:%40t%10.4f sec  (n=%i)\n";

int
   Verbosity = 0;

void FatalError( std::string const &Message,
    char const *pFromWhere, int nLine )
{
   std::stringstream str;
   str << Message;
   if (pFromWhere) {
      str << " at " << pFromWhere << ":" << nLine;
   };
   throw std::runtime_error(str.str());
}



// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version )
{
//    out << fmt::unind();
   out << "\n*** " << Name;
   if (Version != "")
      out << " [Ver. " << Version << "]";
   out << " ***\n" << std::endl;
//    out << fmt::eind();
}


bool LoadFileIntoMemory(TArray<char> &pFileContent,
      std::string const &FileName, unsigned *pFileLength)
{
//    std::cout << "in LoadFileIntoMemory. Name: '" << FileName << "'." << std::endl;
#ifdef QT_CORE_LIB
   // Use QT functions to load the file. This gives us access to files
   // stored in resources and the application clipboard.
   QByteArray
      Text;
   if (FileName == ":/!clipboard!") {
      QClipboard
         *clipboard = QApplication::clipboard();
      QString
         Subtype = "plain",
         sText = clipboard->text(Subtype);
      Text = sText.toUtf8();
   } else {
      QFile
         StyleFile(QString(FileName.c_str()));
      if (!StyleFile.open(QFile::ReadOnly))
         return false;
      Text = StyleFile.readAll();
   }
   pFileContent.resize(Text.size() + 2);
   // 0-terminate the string.
   pFileContent[Text.size()] = 0;
   pFileContent[Text.size()+1] = 0;
   memcpy(&pFileContent[0], Text.data(), Text.size());
   if (pFileLength)
       // if requested, store actual file length at target location.
      *pFileLength = Text.size();
   return true;
#else
   // read file data with C++ standard library routines
   std::ifstream
      File(FileName.c_str());
   std::size_t
      FileLength;
   if (false == File.good())
      return false;
   File.seekg(0, std::ios::end);
   FileLength = File.tellg();
   if (0 != pFileLength)
      // if requested, store actual file length at target location.
      *pFileLength = FileLength;
   pFileContent.resize(2 + FileLength);
   // ^- cgk 2020-02-27: hm... what's with these +2 things? I guess one might
   // be for allowing interpretation as 0-terminated string (but since we don't
   // scan for other 0s in the file content, that would not be safe, either).
   // But two?
   // Update: apparently the result of this is attached to std::stringstreams
   // in several places, which require 0-termination. Still doesn't explain the
   // two 0s.
   memset(&pFileContent[0], 0, 2 + FileLength);
   File.seekg(0, std::ios::beg);
   File.read(&pFileContent[0], FileLength);
   return true;
#endif
}




std::string FormatConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt, double fFinalResidual, double fFinalEnergyChange)
{
   fmt::MemoryWriter w;
   w.write("{} failed to converge", CalcType);
   if (iFinalIt != -1)
      w.write(" within {} iterations", iFinalIt);
   w << ".";
   if (fFinalResidual != -1. || fFinalEnergyChange != -1.) {
      w << " (";
      bool bFirst = false;
      if (fFinalResidual != -1.) {
         bFirst = false;
         w.write("Final residual: {:8.2e}", fFinalResidual);
      }
      if (fFinalEnergyChange != -1.) {
         if (!bFirst) w << "; ";
         bFirst = false;
         w.write("Final energy change: {:16.8f}", fFinalEnergyChange);
      }
      w << ")";
   }
   return w.str();
}

FConvergenceError::FConvergenceError(std::string const &CalcType, ptrdiff_t iFinalIt, double fFinalResidual, double fFinalEnergy)
   : FBase(FormatConvergenceError(CalcType, iFinalIt, fFinalResidual, fFinalEnergy))
{}



FLogError::FLogError(std::string const &Reason)
   : FBase(Reason)
{}


FLog::FLog()
   : m_IoFlags(IOFLAG_FlushAfterWrite)
{
}


FLog::~FLog()
{
}


char const
   *p1ResultFmtForTabulation = "${:<32}| {:24.16e} |  ({})\n",
   // ^- used to mark quantities intended for tabulation by a script analyzing the program output
   *p1ResultFmt =            " {:<32}{:18.12f}\n",
   *p1ResultFmtAnnotated =   " {:<32}{:18.12f}  ({})\n",
   *p1CountFmt =             " {:<38}{:12}\n",
   *p1CountFmtAnnotated =    " {:<38}{:12}  ({})\n",
   *p1InfoFmt =              " {:<32}{:>18}\n",
   *p1InfoFmtAnnotated =     " {:<32}{:>18}  ({})\n",
   *p1InfoExpfFmt =          " {:<32}{:18.4e}\n",
   *p1InfoExpfFmtAnnotated = " {:<32}{:18.4e}  ({})\n",
   *p1TimingFmt =            "{:<41}{:10.2f} sec\n",        // <-- note: gets an extra leading " " before "Time for..."
   *p1TimingFmtAnnotated =   "{:<41}{:10.2f} sec  ({}{})\n";


std::string FmtResult(fmt::BasicStringRef<char> Name, double fValue, fmt::BasicStringRef<char> Annotation)
{
   if (Annotation != "")
      return fmt::format(p1ResultFmt, Name, fValue);
   else
      return fmt::format(p1ResultFmtAnnotated, Name, fValue, Annotation);
}


std::string FmtCount(fmt::BasicStringRef<char> Name, ptrdiff_t iValue, fmt::BasicStringRef<char> Annotation)
{
   if (Annotation == "")
      return fmt::format(p1CountFmt, Name, iValue);
   else
      return fmt::format(p1CountFmtAnnotated, Name, iValue, Annotation);
}


std::string FmtTiming(fmt::BasicStringRef<char> Name, double fTimeInSeconds, size_t nTasks)
{
   if (nTasks == size_t(-1))
      return fmt::format(p1TimingFmt, fmt::format(" Time for {}:", Name), fTimeInSeconds);
   else
      return fmt::format(p1TimingFmtAnnotated, fmt::format(" Time per {}:", Name), fTimeInSeconds/double(nTasks), "n=", nTasks);
}



void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue)
{
   w.write(p1ResultFmt, Name, fValue);
   Flush();
}

void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue, fmt::BasicStringRef<Char> Annotation)
{
   w.write(p1ResultFmtAnnotated, Name, fValue, Annotation);
   Flush();
}

void FLog::WriteResult(fmt::BasicStringRef<Char> Name, double fValue, unsigned Flags)
{
   if ((Flags & OUTPUT_Tabulate) != 0) {
      WriteResult(Name, fValue); // write it normally once
      w.write(p1ResultFmtForTabulation, Name, fValue, "->!Tabulate-Average");
      Flush();
   } else {
      WriteResult(Name, fValue);
   }
}

#if defined(USE_GENERAL_SCALARS) && !defined(SCALAR_IS_FLOAT_64)
void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue)
{
   w.write(p1ResultFmt, Name, fValue);
   Flush();
}

void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, fmt::BasicStringRef<Char> Annotation)
{
   w.write(p1ResultFmtAnnotated, Name, fValue, Annotation);
   Flush();
}

void FLog::WriteResult(fmt::BasicStringRef<Char> Name, FScalar fValue, unsigned Flags)
{
   if ((Flags & OUTPUT_Tabulate) != 0) {
      WriteResult(Name, fValue); // write it normally once
      w.write(p1ResultFmtForTabulation, Name, fValue, "->!Tabulate-Average");
      Flush();
   } else {
      WriteResult(Name, fValue);
   }
}
#endif // USE_GENERAL_SCALARS



void FLog::WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue)
{
   w.write(p1CountFmt, Name, iValue);
   Flush();
}

void FLog::WriteCount(fmt::BasicStringRef<Char> Name, ptrdiff_t iValue, fmt::BasicStringRef<Char> Annotation)
{
   w.write(p1CountFmtAnnotated, Name, iValue, Annotation);
   Flush();
}

void FLog::WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value)
{
   w.write(p1InfoFmt, Name, Value);
   Flush();
}

void FLog::WriteInfo(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Value, fmt::BasicStringRef<Char> Annotation)
{
   w.write(p1InfoFmtAnnotated, Name, Value, Annotation);
   Flush();
}

void FLog::WriteInfoExpf(fmt::BasicStringRef<Char> Name, double Value)
{
   w.write(p1InfoExpfFmt, Name, Value);
   Flush();
}

void FLog::WriteInfoExpf(fmt::BasicStringRef<Char> Name, double Value, fmt::BasicStringRef<Char> Annotation)
{
   if (Annotation.empty())
      return WriteInfoExpf(Name, Value);
   w.write(p1InfoExpfFmtAnnotated, Name, Value, Annotation);
   Flush();
}


void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds)
{
   w.write(p1TimingFmt, fmt::format(" Time for {}:", Name), fTimeInSeconds);
   Flush();
}

void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, fmt::BasicStringRef<Char> Annotation)
{
   if (!Annotation.empty()) {
      w.write(p1TimingFmtAnnotated, fmt::format(" Time for {}:", Name), fTimeInSeconds, Annotation, "");
      Flush();
   } else {
      return this->WriteTiming(Name, fTimeInSeconds);
   }
}

void FLog::WriteTiming(fmt::BasicStringRef<Char> Name, double fTimeInSeconds, size_t nTasks)
{
   w.write(p1TimingFmtAnnotated, fmt::format(" Time per {}:", Name), fTimeInSeconds/double(nTasks), "n=", nTasks);
   Flush();
}

void FLog::WriteProgramIntro(fmt::BasicStringRef<Char> Name, fmt::BasicStringRef<Char> Version)
{
   w << "\n*** " << Name;
   if (Version.size() != 0u)
      w << " [Ver. " << Version << "]";
   w << " ***\n\n";
   Flush();
}

FLog::FRunStatus FLog::GetStatus()
{
   // do nothing. derived classes might do more involved things.
   return STATUS_Okay;
}

void FLog::CheckStatus()
{
   // do nothing. derived classes might, however.
   FRunStatus
      Status = GetStatus();
   if (Status == STATUS_Okay)
      return;
   if (Status == STATUS_AbortSignalled)
      throw ct::FLogError("Abort signalled.");
   throw ct::FLogError("FLog signalled unrecognized error.");
}


FLogStdStream::FLogStdStream(std::ostream &TargetStream)
   : m_TargetStream(TargetStream)
{
}

void FLog::EmitWarning(fmt::BasicStringRef<Char> Name)
{
   w << " WARNING: " << Name << "\n";
   Flush();
}

void FLog::EmitError(fmt::BasicStringRef<Char> Name)
{
   w << " ERROR: " << Name << "\n";
   Flush();
}

void FLogStdStream::Flush()
{
   m_TargetStream << w.c_str();
   w.clear();
   m_TargetStream.flush();
}


FMemoryLog::FMemoryLog()
   : FLogStdStream(m_TextStream)
{}

FMemoryLog::~FMemoryLog()
{}

std::string FMemoryLog::GetText(unsigned Flags)
{
   Flush();
   std::string Text = m_TextStream.str();
   if (bool(Flags & GETTEXT_ClearLog))
      Clear();
   return Text;
}

void FMemoryLog::Clear()
{
   Flush();
   m_TextStream.str(std::string());
}



} // namespace ct


namespace ct {

static bool
   // to allow for empty base paths...
   g_FileLocator_BasePathIsSet = false;
static std::string
   // FFileLocator objects will attempt to resolve file names relative to this.
   // By default, that would be what is returned by GetExePath()
   g_FileLocator_BasePath;


void FFileLocator::SetBasePath(std::string const &BasePath)
{
   g_FileLocator_BasePath = BasePath;
   g_FileLocator_BasePathIsSet = true;
}



enum FStrSubstFlags {
   STRSUBST_MatchAtFrontOnly = 0x0001
};

struct FFileLocatorState : public FIntrusivePtrDest
{
   typedef std::list<std::string>
      FDirectoryList;
   FDirectoryList
      m_SearchPaths;
   typedef std::set<std::string>
      FDirectorySet;
   FDirectorySet
      m_SearchPaths_Unordeded;

   struct FSubstEntry {
      std::string Target;
      unsigned Flags;
      FSubstEntry() : Flags(0) {};
      FSubstEntry(std::string Target_, unsigned Flags_) : Target(Target_), Flags(Flags_) {}
   };
   typedef std::map<std::string, FSubstEntry>
      FSubstMap;
   FSubstMap
      m_Substitutions;
};
typedef TIntrusivePtr<FFileLocatorState>
   FFileLocatorStatePtr;


struct FFileLocatorImpl
{
   void AddSearchPath(std::string const &PathName);

   FFileLocatorImpl();

   std::string FindFile(std::string const &FileName, unsigned Flags) const;

   void PushState();
   void PopState();

   void AddHomePath();
   void ImportEnvironmentVar(std::string VarName);
   void AddSubstitution(std::string const &From, std::string const &To, unsigned Flags);
protected:
   typedef FFileLocatorState::FDirectoryList
      FDirectoryList;
   typedef FFileLocatorState::FDirectorySet
      FDirectorySet;
   typedef std::list<FFileLocatorStatePtr>
      FStateStack;
   typedef std::list<bool>
      FStateStackFlags;
   typedef FFileLocatorState::FSubstEntry
      FSubstEntry;
   typedef FFileLocatorState::FSubstMap
      FSubstMap;

   FStateStack
      m_StateStack;
   FStateStackFlags
      // we only copy the actual state stack contents when they are
      // modified. Before that, the current level is aliased to the
      // higher one. This list keeps track of whether or not the current
      // state is (still) aliased to the previous one.
      m_StateAliased;
   size_t
      m_StateStackDepth;

   void UnaliasState(); // if current stack level is aliased to the last one, make a copy and un-alias it.
   FFileLocatorState &State() { assert(!m_StateStack.empty()); return *m_StateStack.front(); }
   FFileLocatorState const &State() const { assert(!m_StateStack.empty()); return *m_StateStack.front(); }
   FDirectoryList &m_SearchPaths() { return State().m_SearchPaths; }
   FDirectoryList const &m_SearchPaths() const { return State().m_SearchPaths; }
   FDirectorySet &m_SearchPaths_Unordeded() { return State().m_SearchPaths_Unordeded; }
   FDirectorySet const &m_SearchPaths_Unordeded() const { return State().m_SearchPaths_Unordeded; }
   FSubstMap &m_Substitutions() { return State().m_Substitutions; }
   FSubstMap const &m_Substitutions() const { return State().m_Substitutions; }
   bool &_StateAliased() { assert(!m_StateAliased.empty()); return m_StateAliased.front(); }

   // apply substitutions and remove trailing '/' / '\' characters.
   std::string ProcessFileName(std::string const &FileName) const;
private:
   void operator = (FFileLocatorImpl const &); // not implemented
   FFileLocatorImpl(FFileLocatorImpl const &); // not implemented
};


FFileLocator::FFileLocator()
   : p(new FFileLocatorImpl)
{
}

FFileLocator::~FFileLocator()
{
   delete p;
   p = 0;
}


FFileLocatorImpl::FFileLocatorImpl()
{
   m_StateStack.push_front(FFileLocatorStatePtr(new FFileLocatorState()));
   m_StateStackDepth = 0; // no "extra" layers of state (only the current top layer, but we don't count this here).
   m_StateAliased.push_front(false); // original state is never aliased.
   if (!g_FileLocator_BasePathIsSet)
      throw std::runtime_error("FFileLocator constructor: Use SetBasePath() before instanciating FFileLocator objects.");
   std::string
      // note: overrides in the base path are used in the Debug/Release configurations of VC
      // on windows. Otherwise the relative directory stuff does not work.
      sBasePath = g_FileLocator_BasePath;
//    if (sBasePath.size() == 0)
//       sBasePath = GetExePath();
   if (sBasePath.empty())
      AddSubstitution("$BASEPATH", ".", 0);
   else
      AddSubstitution("$BASEPATH", ProcessFileName(sBasePath), 0);
}


void FFileLocatorImpl::AddSubstitution(std::string const &From, std::string const &To, unsigned Flags)
{
   if (g_FileLocator_DebugPrints)
      std::cout << fmt::format(": register file path substitution: '{}' --> '{}' (flags: {})\n", From, To, Flags);
   FSubstMap::iterator
      it = m_Substitutions().find(From);
   if (it == m_Substitutions().end() || it->second.Flags != Flags || it->second.Target != To) {
      UnaliasState();
      m_Substitutions()[From] = FSubstEntry(To, Flags);
   }
}


static std::string GetEnvironmentVar(char const *pName)
{
   assert(pName != 0);
   if (pName[0] == '$')
      // skip over leading "$" if in the variable name.
      pName += 1;
   return std::string(std::getenv(pName));
}


void FFileLocatorImpl::AddHomePath()
{
   std::string
      sHomePath;
#ifdef _WIN32
// "%HOMEDRIVE%%HOMEPATH%"
   sHomePath = FFileLocator::JoinPath(GetEnvironmentVar("HOMEDRIVE"), GetEnvironmentVar("HOMEPATH"));
#else
   sHomePath = GetEnvironmentVar("HOME");
#endif
   AddSubstitution("$HOME", sHomePath, 0);
   AddSubstitution("~/", sHomePath + "/", STRSUBST_MatchAtFrontOnly);
   // ^- the extra '/' is necessary because we remove the '/' of '~/'.
}


void FFileLocator::AddHomePath()
{
   return p->AddHomePath();
}


void FFileLocatorImpl::ImportEnvironmentVar(std::string VarName)
{
   if (VarName.empty())
      return;
   if (VarName[0] == '$')
      // remove leading "$" if supplied (e.g., turn "$HOME" into "HOME")
      VarName.erase(0, 1);
   AddSubstitution("$" + VarName, GetEnvironmentVar(VarName.c_str()), 0);
}


void FFileLocator::ImportEnvironmentVar(std::string const &VarName)
{
   return p->ImportEnvironmentVar(VarName);
}


void FFileLocatorImpl::PushState()
{
   m_StateStack.push_front(m_StateStack.front());
   // mark that current stack level is aliased to the last one
   // (the line above copies the *pointer* to the lower level state).
   m_StateAliased.push_front(true);
   m_StateStackDepth += 1;
}


void FFileLocatorImpl::PopState()
{
   if (m_StateStackDepth == 0)
      throw std::runtime_error("FFileLocatorImpl::PopState(): attempted to restore a non-existent state");
   m_StateStackDepth -= 1;
   m_StateStack.pop_front();
   m_StateAliased.pop_front();
}


void FFileLocatorImpl::UnaliasState()
{
   assert(!m_StateAliased.empty());
   if (_StateAliased()) {
      // make a copy of the (aliased) state of the current level, and replace current
      // level state pointer by a pointer to this copy.
      FFileLocatorStatePtr
         pCopyOfState(new FFileLocatorState(State()));

      m_StateStack.front() = pCopyOfState;
      _StateAliased() = false;
   }
}


void FFileLocator::PushState()
{
   p->PushState();
}


void FFileLocator::PopState()
{
   p->PopState();
}


std::string FFileLocator::JoinPath(std::string const &Path, std::string const &FileName)
{
   string_slice
      slPath(Path);
   if (slPath.empty() || slPath == "." || slPath == "./")
      return FileName;
   if (slPath.endswith("/")
#ifdef _WIN32
         || slPath.endswith("\\")
#endif
      )
      return fmt::format("{}{}", Path, FileName);
   else {
#ifdef _WIN32
      return fmt::format("{}\\{}", Path, FileName);
#else
      return fmt::format("{}/{}", Path, FileName);
#endif
   }
}


static void ReplaceSubStr(std::string &inout, std::string const &key, std::string const &value, unsigned Flags)
{
   size_t ipos = std::string::npos;
   if (bool(Flags & STRSUBST_MatchAtFrontOnly)) {
      if (StartsWith(inout, string_slice(key)))
         ipos = 0;
   } else {
      ipos = inout.find(key);
   }
   if (ipos != std::string::npos) {
      inout.replace(ipos, key.size(), value);
   }
}


std::string FFileLocatorImpl::ProcessFileName(std::string const &FileName_) const
{
   std::string
      FileName = FileName_;
   FSubstMap::const_iterator
      itSubst;
   // apply directory substitutions
   for (itSubst = m_Substitutions().begin(); itSubst != m_Substitutions().end(); ++ itSubst) {
      ReplaceSubStr(FileName, itSubst->first, itSubst->second.Target, itSubst->second.Flags);
   }

   // get rid of spaces to either side and trailing back- and front-slashes.
   string_slice
      slName(FileName);
   slName.trim();
   while (slName.first != slName.last && (slName.endswith("/") || slName.endswith("\\")))
      slName.last -= 1;
   return slName.to_str();
}


void FFileLocatorImpl::AddSearchPath(std::string const &PathName)
{
   if (PathName.empty())
      return;
   std::string
      s = ProcessFileName(PathName);
   // add the string to our search path list... unless we already have it there.
   if (m_SearchPaths_Unordeded().find(s) == m_SearchPaths_Unordeded().end()) {
      // we'll modify the state. Unlink it from higher level states now, unless already done.
      UnaliasState();
      m_SearchPaths().push_front(ProcessFileName(PathName));
      m_SearchPaths_Unordeded().insert(s);
   }
}


static bool IsFileReadable(std::string const &FileName)
{
   // well... slow, but simple.
   std::ifstream
      File(FileName.c_str());
   return File.good();
}


std::string FFileLocatorImpl::FindFile(std::string const &FileName_, unsigned Flags) const
{
   if (g_FileLocator_DebugPrints) {
      std::cout << fmt::format(": searching for {} (#paths: {})\n", FileName_, m_SearchPaths().size());
   }
   std::string
      FileName;
   if (bool(Flags & FFileLocator::FINDFILE_NoSubstitutions)) {
      FileName = FileName_;
   } else {
      FileName = ProcessFileName(FileName_);
   }

   FDirectoryList::const_iterator
      itPath;
   for (itPath = m_SearchPaths().begin(); itPath != m_SearchPaths().end(); ++ itPath) {
      std::string
         s = FFileLocator::JoinPath(*itPath, FileName);
      if (IsFileReadable(s)) {
         if (g_FileLocator_DebugPrints) {
            std::cout << fmt::format(": -> found '{}'\n", s);
         }
         return s;
      }
   }
   // couldn't find it.
   if (bool(Flags & FFileLocator::FINDFILE_RaiseErrorIfNotFound)) {
      fmt::MemoryWriter w;
      w << "! Search paths in FindFile():\n";
      for (itPath = m_SearchPaths().begin(); itPath != m_SearchPaths().end(); ++ itPath) {
         w << fmt::format("!: '{}'\n", *itPath);
      }
      w << fmt::format("FindFile: Failed to locate file '{}'", FileName);
      throw std::runtime_error(w.str());
   } else {
      return FileName;
   }
}


std::string FFileLocator::FindFile(std::string const &FileName, unsigned Flags) const
{
   return p->FindFile(FileName, Flags);
}


void FFileLocator::AddSearchPath(std::string const &PathName)
{
   return p->AddSearchPath(PathName);
}


static size_t FindPathFileSeparator(std::string const &PathAndFileName)
{
   size_t
      iSlash = PathAndFileName.rfind('/');
#ifdef _WIN32
   size_t
      iBackSlash = PathAndFileName.rfind('\\');
   if (iBackSlash != std::string::npos && (iSlash == std::string::npos || iBackSlash > iSlash))
      iSlash = iBackSlash;
#endif
   assert(iSlash == std::string::npos || iSlash < PathAndFileName.size());
   return iSlash;
}


std::string FFileLocator::BaseName(std::string const &PathAndFileName)
{
   size_t
      iSlash = FindPathFileSeparator(PathAndFileName);
   if (iSlash != std::string::npos) {
      if (iSlash + 1 >= PathAndFileName.size())
         return "";
      else
         return PathAndFileName.substr(iSlash+1);
   } else {
      // no directory separator. Return original file name.
      return PathAndFileName;
   }
}


std::string FFileLocator::PathName(std::string const &PathAndFileName)
{
   size_t
      iSlash = FindPathFileSeparator(PathAndFileName);
   if (iSlash != std::string::npos) {
      // str contains a directory separator. Return input str up to
      // this separator (excluding the separator itsef).
      // (note: That will prevent adding just "/" as a search path,
      // but I do not think this would be a great idea in any case...)
      return PathAndFileName.substr(0, iSlash);
   } else {
      // no directory separator. Return empty string.
      return std::string();
   }
}



} // namespace ct

// kate: indent-width 3;
