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

#ifndef CX_PARSE1_H
#define CX_PARSE1_H

#include <string>
#include <vector>
#include <ostream>
#include <utility>

#include "CxTypes.h"
#include "CxPodArray.h"

namespace ct {

std::string ReadFile(std::string FileName);

typedef std::string::iterator
   str_it;
typedef std::string::const_iterator
   cstr_it;

typedef std::pair<char const *, char const *>
   pchar_range;

size_t const
   // used to be called nMaxStrSplit. For standard split results, it defines
   // the number of array elements directly stored in the static storage of
   // the TPodArray type (i.e., beyond that dynamic memory allocation are
   // required and used, but for splits shorter than this, all data can be
   // kept inside the local variables on the stack)
   nSplitStaticStorage = 32;

extern char const
   *g_pDefaultPropertyList_OpenDelimCloseDef; // it's "{;}:";


typedef TArrayFix<char, nSplitStaticStorage-1>
   short_split_delim_list;
typedef TArray<char, nSplitStaticStorage-1>
   long_split_delim_list;
typedef long_split_delim_list
   split_delim_list;

// ^- NOTE:
// - originally, we only had the fixed maximum size delim lists and split
//   results which are now called short_split_delim_list and short_split_result.
//   Using these in parsing operations should be most efficient if the maximum
//   split size is known.
// - Later on, we introduced long_split_delim_list and long_split_result for
//   general size buffers.
// - Later on, the variable size POD array class TArray<> was extended to
//   optionally allow a fixed-maximum-size constant buffer to use it allocation-
//   less for small sequences. This mostly removes needs for keeping the
//   fixed size versions around.
// - 2020-04-27, I renamed the fixed-size versions to short_split_result /
//   short_split_delim_list, and made the variable size ones (with a constant-
//   size static buffer) the defaults for split_result and split_delim_list.
// - There may still be incompatibilities around in various programs depending
//   on this, if they assume split_result and long_split_result to be different,
//   for example.


struct is_path_separator_c
{
   inline bool operator() (char c) const
   {
   #ifdef _WIN32
      return c == '/' || c == '\\';
   #else
      return c == '/';
   #endif
   }
};


struct is_equal_c
{
   is_equal_c(char c_)  : m_TargetChar(c_) {}

   inline bool operator() (char c) const {
      return c == m_TargetChar;
   }
protected:
   char
      m_TargetChar;
};


static std::string const s_EmptyDummyStrToMakeDefinedIterators;


struct string_slice
{
   cstr_it
      first, last;

   string_slice() { first = s_EmptyDummyStrToMakeDefinedIterators.end(); last = first; };

   string_slice(cstr_it itFirst_, cstr_it itLast_)
      : first(itFirst_), last(itLast_)
   {
      assert(first <= last);
   }

   string_slice(std::string const &s)
      : first(s.begin()), last(s.end())
   {}

   // make this object reference the external std::string 's'. Beware that 's'
   // needs to stay alive for the lifetime of *this!
   void assign_ref(std::string const &s)
   {
      first = s.begin();
      last = s.end();
   }

//    // construct from 0-terminated string. Will explode if string is not actually
//    // null-terminated...
//    explicit string_slice(char const *p);

   cstr_it begin() const { return first; }
   cstr_it end() const { return last; }
   size_t size() const { return last - first; }
   bool empty() const { return first == last; }
   void clear() { first = last; }

   enum FToStrFlags {
      // if set, remove first level of leading/trailing quotes of returned string
      TOSTR_TrimQuotes = 0x01,
      // if set, fail if the input string is NOT quoted (implies TrimQuotes).
      // If present, also will remove one level of escapes from inner quotes
      // matching the outer quotes (e.g., for an input string
      //    "this is \"fine\""
      // it will return
      //    this is "fine"
      TOSTR_QuotesRequired = TOSTR_TrimQuotes | 0x02,
      // interpret some basic C-style/Python escape sequences in the input string;
      // At this moment:
      //    \\  -->  (single backslash)
      //    \n  ->  newline
      //    \t  ->  horizontal tab charater
      //    \'  --> '
      //    \"  --> "
      // It will not process some of the less common console-style or binary-
      // style escapes (e.g., \a (ascii bell), \b (backspace), \r/\f (carriage
      // return/form feed, \v (vertical tab)), nor will it interpret ASCII or
      // unicode symbols given as numerical escape sequences.
      TOSTR_InterpretEscapes = 0x04,
      TOSTR_EvalStringLiteral = TOSTR_QuotesRequired | TOSTR_InterpretEscapes,
      TOSTR_Default = 0
   };

   // these ones raise exceptions if invalid input is found
   std::string to_str(unsigned Flags = TOSTR_Default) const;
   ptrdiff_t to_int() const;
   size_t to_size() const;
   double to_float() const;
   bool to_bool() const;

   // same as above; the overloaded versions are for templatizing
   void convert(double *pr) const { *pr = to_float(); }
   void convert(ptrdiff_t *pr) const { *pr = to_int(); }
   void convert(size_t *pr) const { *pr = to_size(); }
   void convert(bool *pr) const { *pr = to_bool(); }
   void convert(std::string *pr, unsigned Flags = TOSTR_Default) const { *pr = to_str(Flags); }

   // thes ones return whether the conversion was successful
   bool try_convert_to_int(ptrdiff_t *pr) const;
   bool try_convert_to_size(size_t *pr) const;
   bool try_convert_to_float(double *pr) const;
   bool try_convert_to_bool(bool *pr) const;

   // same as above; the overloaded versions are for templatizing
   bool try_convert(double *pr) const { return try_convert_to_float(pr); }
   bool try_convert(ptrdiff_t *pr) const { return try_convert_to_int(pr); }
   bool try_convert(bool *pr) const { return try_convert_to_bool(pr); }
   bool try_convert(std::string *pr, unsigned Flags = TOSTR_Default) const { *pr = to_str(Flags); return true; }

   cstr_it find(char c, size_t offset = 0, unsigned flags=0) const;
   cstr_it find(char const *c, size_t offset = 0, unsigned flags=0) const;
   template<class FCharPred>
   cstr_it find_c(FCharPred const &pred, size_t offset = 0, unsigned flags=0) const;

   cstr_it rfind(char c) const;
   cstr_it rfind(char const *c) const;
   template<class FCharPred>
   cstr_it rfind_c(FCharPred const &pred) const;

   string_slice substr(size_t iFirst, size_t nCount) const {
      assert(iFirst + nCount <= size());
      return string_slice(first + iFirst, first + iFirst + nCount);
   }

   string_slice substr(cstr_it iFirst, size_t nCount) const {
      assert(iFirst + nCount <= end());
      return string_slice(iFirst, iFirst + nCount);
   }

   enum split_flags {
      // if set, whitespace at the beginning and end of each output string split will be removed.
      SPLIT_SkipWhitespace = 0x01,
      // if set, text between two parenthesis (any kind) will not be interpreted, but treated as a unit.
      SPLIT_SkipParenthesis = 0x02,
      // if set, text between two matching quotes (".." or '..') will not be interpreted, but treated as a unit.
      // This takes precedence over skip parenthesis.
      SPLIT_SkipQuotes = 0x04
   };

//    // returns out.size().
//    size_t split(TArrayFix<string_slice, nSplitStaticStorage> &out, char Delim, size_t flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;
//    // split at any of the given delimiters. store delimiters in pOutDelims if given:
//    // pOutDelims[i] is the delimiter between out[i] and out[i+1].
//    size_t split_any(TArrayFix<string_slice, nSplitStaticStorage> &out, TArrayFix<char, nSplitStaticStorage-1> *pOutDelims, char const *pDelims, size_t nDelims, size_t flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;
//
//    // convenience function for splitting a list of sub-properties into parts.
//    size_t split_list(TArrayFix<string_slice, nSplitStaticStorage> &out, char ListStart = '[', char ListEnd = ']', char ListDelim = ',') const;

   // split at any of the given delimiters. store delimiters in pOutDelims if given:
   // pOutDelims[i] is the delimiter between out[i] and out[i+1].
   // Instanciated for:
   // - split_any(split_result &, split_delim_list *, ...)   (with fixed size arrays)
   // - split_any(long_split_result &, long_split_delim_list &, ...)
   template<class TSliceList, class TDelimList>
   size_t split_any(TSliceList &out, TDelimList *pOutDelims, char const *pDelims, size_t nDelims, size_t flags) const;

   template<class TSliceList>
   size_t split_lines(TSliceList &out, size_t flags = 0) const;

   template<class TSliceList>
   size_t split(TSliceList &out, char Delim, size_t Flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const
   {
      // note: if Instanciated with long_split_result, this will invoke the
      // '<long_split_result, split_delim_list>'-version of this function (note the missing 'long_'
      // on the delimiters). However, as the output delimiter list is optional and we just pass a 0-ptr
      // (and therefore no delimiters get stored in any case) this is still fine.
      return split_any<TSliceList,split_delim_list>(out, 0, &Delim, 1, Flags);
   }

   // convenience function for splitting a list of sub-properties into parts.
   template<class TSliceList>
   size_t split_list(TSliceList &out, char ListStart='[', char ListEnd=']', char ListDelim=',') const;


   bool operator == (char const *p) const;
   bool operator != (char const *p) const { return !this->operator == (p); }

   // returns if *this is equal to *p if both are lower-cased.
   bool is_equal_except_case(char const *p) const;

   // similar to startswith(), but returns an iterator to the end of the
   // matching sequence in *this in case of a successful match.
   // Returns this->first unless the sequence in `search` matches completely.
   // NOTE:
   // - this returns `first` on failure, not `last`! `last` is a valid output
   //   if *this is a complete match with `search`.
   cstr_it match_at_start(char const *search) const;
   cstr_it match_at_start(string_slice const &search) const;
   // similar to endsswith(), but returns an iterator to the start of the
   // matching sequence in *this in case of a successful match.
   // Returns this->last unless the sequence in `search` matches completely.
   cstr_it match_at_end(char const *search) const;
   cstr_it match_at_end(string_slice const &search) const;

   bool startswith(char const *p) const;
   bool startswith(string_slice const &search) const;
   // startswith() for more complex cases, including case-insensitive matching
   // and some basic full-token matching support. Flags is a bit field of SKIP_*
   // flags (see below). Returns 'true' iff this->try_skip_literal(p, Flags)
   // would succeed.
   bool startswith(char const *p, unsigned Flags) const;
   bool startswith(string_slice const &search, unsigned Flags) const;
   bool endswith(char const *pStr) const;
   bool endswith(string_slice const &search) const;
   bool startswith(char c) const;
   bool endswith(char c) const;

   // checks if *this starts with pSearch; if not, returns `false`; if yes,
   // increments this->first to the first character behind the matching string,
   // and returns `true`. This is a convenience wrapper around try_skip_literal(),
   // which can also handle some more complex arrangements.
   //
   // Example:
   //   string_slice sl("df-rhf");
   //   bool df = false;
   //   if (sl.remove_if_at_start("df-")) {
   //      df = true;
   //   }
   //
   //   At the end, `sl` will be covering the string "rhf" only, and, `df`
   //   will be `true`. If instead, `sl` is initialized to "rhf" only at
   //   the start, it will still cover only "rhf" at the end, but `df`
   //   will be false.
   bool remove_if_at_start(char const *pSearch);
//    bool remove_if_at_start(string_slice const &search);
   // checks if *this ends with pSearch; if not, returns `false`; if yes,
   // decrements this->last to the first character of the matching string,
   // and returns `true`.
   bool remove_if_at_end(char const *pSearch);
//    bool remove_if_at_end(string_slice const &search);

   void resize(size_t n) {
      last = first + n;
   };

   // skip all characters contained in pCharsToSkip (0-terminated) at begin and end.
   // If pCharsToSkip is 0, this skips whitespace.
   void trim(char const *pCharsToSkip=0);
   // skip any characters in pCharsToSkip (if pCharsToSkip != 0) or whitespace (if pCharsToSkip==0) at end only
   void trim_right(char const *pCharsToSkip=0);
   // skip any characters in pCharsToSkip (if pCharsToSkip != 0) or whitespace (if pCharsToSkip==0) at start only
   void trim_left(char const *pCharsToSkip=0);
   // skip pair of matching quotes at begin and end of str.
   // If QuotesRequired is set, raise exception if there are no quotes.
   // returns either 0, if no quotes were found, or the opening quote character
   // of the type of quote which was removed.
   char trim_quotes(char const *pQuoteChars=0, bool QuotesRequired = true);

   // Flags for try_skip_literal and try_skip_literals
   enum skip_literal_flags {
      // if set, ignore whitespaces on the left hand side
      // of *this before searching for the target character sequence
      SKIP_SkipWhitespaceBefore = 0x0001,
      // if set, if a literal sequence was found and skipped,
      // then trim_left the resulting string before returning
      SKIP_SkipWhitespaceAfter = 0x0002,
      SKIP_SkipWhitespace = SKIP_SkipWhitespaceBefore | SKIP_SkipWhitespaceAfter,
      // compare case-insensitive when looking for the literal to skip
      SKIP_IgnoreCase = 0x0004,
      // if set, skip will not match an expression which starts with pSequence,
      // but directly continues with an alphabetic character (a-z,A-Z,_).
      //
      // For example, try_skip_literal("solve", SKIP_TokenAlpha) will match the
      // string "solve(A,x)", but will not match the string "solvent: methanol".
      // (because in "solvent", an 'n' follows the "solve")
      SKIP_TokenAlpha = 0x0008,
      SKIP_TokenDigit = 0x0010,
      SKIP_TokenUnderscore = 0x0020,
      SKIP_TokenMinus = 0x0040,
      SKIP_TokenAlphaNumeric = SKIP_TokenAlpha | SKIP_TokenDigit | SKIP_TokenUnderscore,
      SKIP_TokenIdMinus = SKIP_TokenAlpha | SKIP_TokenDigit | SKIP_TokenUnderscore | SKIP_TokenMinus,
      SKIP_Default = 0
   };

   // if controlled string_slice begins with pSequence, then
   // advance this->first to the end of the occurrence of this sequence,
   // and return the input pSequence.
   // If sequence is not found, return 0.
   char const *try_skip_literal(char const *pSequence, unsigned Flags = SKIP_Default);
   // as try_skip_literal, but will accept any of a given set of starting
   // sequences, which are provided in a double-0-terminated input string
   // (e.g., "add\0sub\0mul\0div\0\0").
   // If pSkippedSequence is != 0, it will receive a string_slice for a subset
   // of *this which describes the matched and skipped occurrence.
   char const *try_skip_literals(char const *pSequences, string_slice *pSkippedSequence = 0, unsigned Flags = SKIP_Default);

   // should this be here? not really, I guess...
   typedef std::pair<string_slice, string_slice>
      FPair;
   FPair split_path() const;
   void remove_trailing_path_separators();

   FPair rsplit1(char Delim, size_t Flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;
   FPair split1(char Delim, size_t Flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;

   pchar_range to_pchar_range() const {
      if (empty())
         return pchar_range(0,0);
      else
         return pchar_range(&*first, &*first + size());
         // ^- why &*first + size() instead of &*last? the latter would
         // formally dereference a potentially invalid iterator, and some
         // debug CRTs might take offence to that.
   }

   operator std::string () const { return to_str();  };
   operator pchar_range () const { return this->to_pchar_range();  };

   char const &operator [] (size_t i) const { assert(i < size()); return first[i]; };
protected:
   // advance iterator iPos, to next position.
   // Flags may specify SPLIT_SkipParenthesis and/or SPLIT_SkipWhitespace
   // to skip over blocks between matching parenthesis and/or whitespace
   inline void increment_it(cstr_it &iPos, unsigned Flags) const;
   void skip_over_flagged(cstr_it &iPos, unsigned Flags) const;
};

std::ostream &operator << (std::ostream &out, string_slice const &s);




typedef TArrayFix<string_slice, nSplitStaticStorage>
   short_split_result;
typedef TArray<string_slice, nSplitStaticStorage>
   long_split_result;
typedef long_split_result
   split_result;

struct FCodeLoc {
   size_t iLine;
   FCodeLoc(size_t iLine_) : iLine(iLine_) {};
};

struct FLine {
   string_slice
      Text;
   FCodeLoc
      Loc; // original location (line number) in input.
   FLine(string_slice const &Text_, FCodeLoc Loc_) : Text(Text_), Loc(Loc_) {}
};
typedef std::vector<FLine>
   FLineList;
typedef FLineList::const_iterator
   cline_it;
std::ostream &operator << (std::ostream &out, FLine const &Line);

void SplitLinesAndRemoveComments(FLineList &out, std::string const &Text);

struct FIdentifierDecl {
   enum {
      // note: none of the next allowed at start of word.
      AllowMinus = 0x0001,
      AllowColon = 0x0002,
      AllowCaret = 0x0004,
      AllowDigits = 0x0008,
      // some extra options for details...
      ForbidUnderscoreAsFirst = 0x0100
   };
   explicit FIdentifierDecl(unsigned Flags = FIdentifierDecl::AllowDigits) : m_Flags(Flags) {}
   bool IsCharAllowedAsFirst(char c) const;
   // test for additional options in Allow*.
   bool IsCharAllowed(char c) const;

   inline bool Enabled(unsigned Flag) const { return bool(m_Flags & Flag); }
protected:
   unsigned m_Flags;
};
// returns iPos when *iPos is not a valid identifier character. Note: this does *not* skip white space!
cstr_it FindIdentifierEnd(cstr_it iPos, cstr_it iEnd, FIdentifierDecl const &Decl = FIdentifierDecl());
char const *FindClosingParenthesis(char Open, char Close, char const *pFirst, char const *pLast);
cstr_it SkipWhiteSpace(cstr_it iPos, cstr_it iEnd);


struct FPropertyStr {
   string_slice
      Name, Content;
};

enum FPropertyListParseOptions {
   PROPLIST_AllowEntriesWithoutName = 0x01,
   PROPLIST_AllowOneEntryWithoutName = 0x02, // as 'PROPLIST_AllowEntriesWithoutName', but allow at most one such entry.
   PROPLIST_UnnamedList = 0x04,
   PROPLIST_AllowOmitOpenClose = PROPLIST_UnnamedList | 0x08, // note: DO NOT USE THIS if list entries could potentially start with cOpen and end with cClose! Ambiguous input then!
   PROPLIST_AllowComplexKeys = 0x10, // if set, allow for keys in (key,value) which are themselves complex parenthesized expressions
   PROPLIST_IdentifierWithOptionalArgs = 0x20 // if set, require a list name which is a identifier-class entity; but make the actual arguments (including open/close parens) optional.
};

struct FPropertyListStr : public TArray<FPropertyStr, nSplitStaticStorage>
{
   string_slice
      Name;
   // Flags: bit combination of PROPLIST_*
   explicit FPropertyListStr(string_slice const &Input, char const *pAssumedName = 0, unsigned Flags = 0, char const *pOpenDelimCloseDef = g_pDefaultPropertyList_OpenDelimCloseDef);
};


// tests if c is ' ' or '\t' (at the moment). Funky name is for avoiding conflicts
// with other functions defined for this purpose.
bool is_whitespace_cxp1(char c);


void string_slice::increment_it(cstr_it &iPos, unsigned Flags) const
{
   assert(iPos < last);
   ++ iPos;
   if (Flags != 0 && iPos < last)
      skip_over_flagged(iPos, Flags);
}

// template<class FCharPred>
// cstr_it string_slice::find_c(FCharPred const &pred, size_t offset) const
// {
//    cstr_it it(first + offset);
//    for (; it != last; ++ it) {
//       if (pred(*it))
//          return it;
//    }
//    return last;
// }


template<class FCharPred>
cstr_it string_slice::find_c(FCharPred const &pred, size_t offset, unsigned flags) const
{
   cstr_it it(first + offset);
   for (; it != last; increment_it(it, flags)) {
      if (pred(*it))
         return it;
   }
   return last;
}


template<class FCharPred>
cstr_it string_slice::rfind_c(FCharPred const &pred) const
{
   cstr_it it(last - 1);
   for ( ; it >= first; --it) {
      if (pred(*it))
         return it;
   }
   return last;
}


// returns whether this str starts with an identifier or empty string,
// directly followed by a '{'. PropertyListFlags: bit field of
// PROPLIST_*
bool looks_like_property_list(string_slice const &sl, unsigned PropertyListFlags, char const *pOpenDelimCloseDef = g_pDefaultPropertyList_OpenDelimCloseDef);


} // namespace ct


namespace ct {
   // Up next: some legacy code for basic string transformations; temporarily
   // moved here from CxIo.h/.cpp
   //
   // WARNING:
   // - these use std::isspace, which, among other things, treats newline as
   //   whitespace!
   // - And of course they are slow. Use for small stuff only, if at all!
   //
   // Notes:
   // - I wanted to get rid of boost to simplify compiling (it causes sooo many
   //   problems), ...that's why these things are here. Additionally, at least
   //   the core code needs to work in C++98 (fo the time being, anyway).
   //
   // - These used to be in CxIo.h; I moved them to CxParse1 for the time
   //   being, because these guys include the CxParse1 dependency anyway.
   //
   // - They are used in some legacy code which has not been updated.
   //   Replacing this stuff by string_slice (of CxParse1) would already
   //   go a long way in avoiding the endless copying around of std::strings.
   //
   // - The reason this is not done yet is that atm they are used in some very
   //   old and *very* dicy parsing code (libmol files and .xyz files). These
   //   are tremendously easy to subtly break, and in such a way that in 1 of
   //   100 cases, one would get not only crashes, but also possibly produce
   //   incorrect results which look almost right without ever knowing.
   //   Ultimately, this comes from the Fortran-inspired formats (of which the
   //   definitions long predate the idea that one might want to define file
   //   formats in such a way that checks of errors and consistency are possible
   //   at least in principle).
   void TrimLeft(std::string &s);
   void TrimRight(std::string &s);
   void Trim(std::string &s);
   bool StartsWith(std::string const &s, std::string const &prefix);
   bool StartsWith(std::string const &s, char const *prefix);
   bool StartsWith(std::string const &s, char prefix);
   bool EndsWith(std::string const &s, std::string const &suffix);
   bool EndsWith(std::string const &s, char const *suffix);
   void StripLineComment(std::string &s, char const *prefix);
   // Makes a string lower case. Will break with multibyte characters.
   std::string tolower(std::string const &in);
   std::string toupper(std::string const &in);
   // Returns a copy of string 'in' with all spaces removed (but not tabs or
   // newlines or other whitespace characters).
   // (note: *all* spaces, not only spaces at front and back!)
   std::string stripwhitespace(std::string const &in);

   bool IsEqual_CaseInsensitive(std::string const &a, std::string const &b);
   // these ones here are for comparing to pchar literals, without converting
   // them to std::strings first.
   bool IsEqual_CaseInsensitive(std::string const &a, char const *b);
   bool IsEqual_CaseInsensitive(char const *a, std::string const &b);
} // namespace ct


namespace ct {

// error in some sort of input -- (e.g., syntax, brken file formats).
class FInputError : public std::runtime_error
{
public:
   typedef std::runtime_error
      FBase;
   explicit FInputError(std::string const &Reason);
};


// some unsupported calculation or combination of options was asked for (e.g., a tau-type meta-GGA with DFXC=1)
class FUnsupportedError : public std::runtime_error
{
public:
   typedef std::runtime_error
      FBase;
   explicit FUnsupportedError(std::string const &Reason);
};


} // namespace ct

#endif // CX_PARSE1_H
