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

#ifndef CX_SORTED_TABLE_H
#define CX_SORTED_TABLE_H

#include <cctype> // for tolower()
#include <stddef.h> // for size_t
#include <string> // also for char_traits
#include <algorithm> // for std::sort and std::lower_bound
#include <vector>
#include <ostream> // for operator << definitions.
#include <stdexcept>
#include <iostream> // FIXME: remove this.

namespace ct {

// don't want to add dependency to CxAssert.h, but don't want
// the <assert> one, either.
static void _sorted_table_assert(bool condition) {
   if (!condition)
      throw std::runtime_error("CxSortedTable.inl: internal error.");
}

class EKeyNotFound : public std::runtime_error {
public:
   EKeyNotFound(std::string const &s) : FBase(s) {}
private:
   typedef std::runtime_error
      FBase;
};


static void RaiseNameNotFoundError(char const *pName)
{
   throw EKeyNotFound("Failed to look up symbol '"+std::string(pName)+"'");
}


enum TShortNameKeyFlags {
   // if set, convert all input names to lower case for making keys.
   NAMEKEY_CaseIndependent = 0x0001
};

/// A directly (by value) stored short string of fixed maximum length, and
/// associated comparison functions (for use as map key)
template<size_t MaxNameLength, unsigned KeyFlags>
struct TShortNameKey
{
   typedef char const
      // type of input argument we use to construct the key from. Here a char const *.
      // This is the type the static lookup table gets as actual arguments.
      *FArgType;
   int _compare(TShortNameKey const &other) const {
      return std::char_traits<char>::compare(&m_Chars[0], &other.m_Chars[0], MaxNameLength);
   }
   // ^-- WARNING: this does *not* just return -1, 0, +1... but possibly large
   //     positive/negative values!

   bool operator < (TShortNameKey const &other) const { return _compare(other) < 0; }
   bool operator == (TShortNameKey const &other) const { return _compare(other) == 0; }
   bool operator != (TShortNameKey const &other) const { return _compare(other) > 0; }

   TShortNameKey() {} // <- needs to have one of these to qualify as aggregate type.
   explicit TShortNameKey(char const *pName) {
      size_t
         i = 0;
      bool
         NameEnded = false;
      for (; i < MaxNameLength; ++ i) {
         char
            pi = 0;
         if (!NameEnded) {
            pi = pName[i];
            if (bool(KeyFlags & NAMEKEY_CaseIndependent))
               pi = static_cast<char>(std::tolower(pi));
         }
         if (pi == 0) {
            // ^- the names are 0-terminated. Fill the rest of
            // the this->m_Chars entries with zeros.
            NameEnded = true;
         }
         m_Chars[i] = pi;
      }
//       std::cout << "    | TShortNameKey: '" << pName << "'  !NameEnded = " << (!NameEnded) << "  i = " << i << "  " << std::endl;
//       if (!NameEnded && pName[i] != 0)
      if (!NameEnded) {
//          throw std::runtime_error("do I ever get here?");
         if (pName[i] != 0)
            RaiseNameNotFoundError(pName);
      }
      // ^-- cgk: 2021-02-02: hm.. g++ 9.3.1's static code analysis says:
      //     CxSortedTable.h:92:23: warning: array subscript 4 is outside array bounds of 'const char [4]' [-Warray-bounds]
      //     However, I do not see how a place where the invalid array element
      //     is deferenced can ever be reached (I think !NameEnded cannot be true
      //     if a valid too-short input str is provided). Just to make sure (normally
      //     g++ is spot on with such messages), I also checked explicitly experimentally,
      //     and the offending place is indeed not reached. So I'll just ignore this for now.
   }

   size_t MaxLength() const { return MaxNameLength; }
   void Print(std::ostream &out) const {
      for (size_t i = 0; i < MaxNameLength; ++ i)
         if (m_Chars[i] != 0)
            out << m_Chars[i];
   }
protected:
   char
      m_Chars[MaxNameLength];
};


template<size_t MaxNameLength, unsigned KeyFlags>
std::ostream &operator << (std::ostream &out, TShortNameKey<MaxNameLength, KeyFlags> const &NameKey)
{
   NameKey.Print(out);
   return out;
}



template<class FKeyType, class FKeyArg, class FMappedType>
struct TLookupTableEntry {
   TLookupTableEntry() {}

   TLookupTableEntry(FKeyArg const &KeyArg_, FMappedType const &Value_)
      : m_Key(KeyArg_), m_Value(Value_)
   {}

   FKeyType
      m_Key;
   FMappedType
      m_Value;
};


template<class FKeyType, class FKeyArg, class FMappedType>
std::ostream &operator << (std::ostream &out, TLookupTableEntry<FKeyType, FKeyArg, FMappedType> const &LookupEntry)
{
   out << LookupEntry.m_Key << ": " << LookupEntry.m_Value;
   return out;
}


template<class TLookupTableEntry, class FKeyCompare, bool SortValuesToo>
struct TLookupTableEntryCompare {
   FKeyCompare m_KeyCompare;
   explicit TLookupTableEntryCompare(FKeyCompare const &KeyCompare_) : m_KeyCompare(KeyCompare_) {}
   bool operator () (TLookupTableEntry const &A, TLookupTableEntry const &B) const {
      if (SortValuesToo) {
         if (m_KeyCompare(A.m_Key, B.m_Key)) return true;
         if (m_KeyCompare(B.m_Key, A.m_Key)) return false;
         return (A.m_Value < B.m_Value);
      } else {
         return m_KeyCompare(A.m_Key, B.m_Key);
      }
   }
};



enum FLookupTableFlags {
   // if set, eliminate duplicate value entries by std::swap instead of direct
   // assign. That is preferable for value types being complex objects (e.g.,
   // std::vectors, maps, etc), but worse for simple direct data types (e.g.,
   // ints, doubles, etc)
   TABLEFLAGS_PreferValueSwapOverAssign = 0x0001,
   // if set, entries with duplicate keys will be retained instead
   // of deleted (i.e., this will be the binary search table equivalent
   // of a std::multimap rather than a regular std::map)
   TABLEFLAGS_AllowEqualKeys = 0x0002,
   // if set, sort by (key,value) instead of just by key.
   // This only makes sense if AllowEqualKeys is set, because otherwise there
   // will be at most one value entry for a given key.
   TABLEFLAGS_SortValues = TABLEFLAGS_AllowEqualKeys | 0x0004
};


/// An associative container-like object (i.e., a map/multimap-like object
/// describing the mapping of "key"  elements to "value" objects), which
/// internally represents its data sets as a flat sorted array. The data
/// structure provides the following formal properties:
///
/// - Element lookup is amortized O(log N)
///   (and it is *extremely* efficient* at that! This holds as long as the
///   elements have reasonably non-equivalent keys)
///
/// - Element insertion is amortized O(log N)
///   (for each element, if adding N elements to a previously empty container)
///
/// - Element removal is amortized O(log N)
///   (for each element, provided a sizeable fraction of the elements do get
///   removed)
///
/// ...however, unlike in standard tree-based associative containers (e.g.
/// RbTrees, which most std::map implementations use), in the present structure
/// these different kinds of operations can not be efficiently mixed.
/// In fact, most changes of the container content invalidate the "sorted"
/// property of the internal array, and require a O(N log N) data preparation
/// step before the next lookup operation can be performed.
///
/// As a consequence, the data structure is only efficient in situations where the "filling phase(s)"
/// and the "lookup phases" can be cleanly separated. However, *if* that is the
/// case --- in particular if a table is filled only once, and then used many
/// times afterwards --- this data structure can be significantly more efficient
/// than any other associative container.
///
/// Template arguments:
///
/// - FKeyType: data type to use in internal data structures to represent the
///   keys of the association
///
/// - FKeyArg: either identical to FKeyType, or another data type from
///   which FKeyType objects can be efficiently constructed. FKeyType must have
///   a constructor taking FKeyArg arguments. All "key" arguments in TLookupTable actually
///   employ FKeyArg types as arguments for communication (e.g., operator [], find,
///   contains, etc.).
///
///   For example, above we define the TShortNameKey template class to use as
///   internal map keys for representing fixed maximum-size strings. However,
///   the as actual interface objects of TLookupTable one can use other string-
///   like objects which are easier to handle (e.g. char const * arguments, by
///   combining FKeyType = TShortNameKey<...> with FKeyArg = char const *).
///
/// - FMappedType: ...what in other languages would be called the "value type":
///   the types of objects the key objects are mapped to.
///
/// - TableFlags: bit field of TABLEFLAGS_* objects controlling details of the
///   behavior. In particular, setting TABLEFLAGS_AllowEqualKeys turns this from
///   a std::map-like object into a std::multimap-like object, and if setting
///   TABLEFLAGS_SortValues additionally to TABLEFLAGS_AllowEqualKeys, then
///   table entries with equivalent keys, are sorted by value.
///
/// - FKeyCompare: std::less<FKeyType>-like predicate imposing a strict weak
///   order on FKeyType. Must be copyable (all the std:: algorithms copy the
///   predicates around, e.g., std::stable_sort or std::lower_bound)
///
template<class FKeyType, class FKeyArg, class FMappedType, unsigned TableFlags, class FKeyCompare = std::less<FKeyType> >
struct TLookupTable
{
   typedef TLookupTableEntry<FKeyType, FKeyArg, FMappedType>
      FValueType;
   typedef TLookupTableEntryCompare<FValueType, FKeyCompare, bool(TableFlags & TABLEFLAGS_SortValues)>
      FValueCompare;

   TLookupTable() : m_Compare(FKeyCompare()) { _CommonInit(); }
   explicit TLookupTable(FKeyCompare const &Compare_) : m_Compare(Compare_) { _CommonInit(); }
   virtual ~TLookupTable() {}

   FMappedType &at(FKeyArg const &KeyArg)
   {
      iterator it = _FindIter(KeyArg);
      if (it == m_Entries.end())
         RaiseKeyError(KeyArg);
      return it->m_Value;
   }

   FMappedType const &at(FKeyArg const &KeyArg) const {
      return const_cast<TLookupTable*>(this)->at(KeyArg);
   }

   bool &contains(FKeyArg const &KeyArg) const {
      return _FindIter(KeyArg) != m_Entries.end();
   }


   void insert_or_assign(FKeyArg const &KeyArg, FMappedType const &Value) {
      // note: potential duplicate entries are cleaned up in _PrepareForSearch. If
      // multiple entries with the same key are added, only the last one is
      // retained.
      m_Entries.push_back(FValueType(KeyArg, Value));
      // mark that we added elements which could potentially collide with
      // earlier elements (in which case the earlier elements must be removed if
      // keys are meant to be unique)
      if (_NeedUniqueKeys())
         m_DuplicateCheckRequired = true;
      _MarkAsNotReady();
   }

   void reserve(size_t nEntries_) { m_Entries.reserve(nEntries_); }

   FMappedType &operator[] (FKeyArg const &KeyArg) { return at(KeyArg); }
   FMappedType const &operator[] (FKeyArg const &KeyArg) const { return at(KeyArg); }

   void Print(std::ostream &out) const {
      out << "{";
      for (size_t i = 0; i != m_Entries.size(); ++ i) {
         if (i != 0)
            out << ", ";
         out << m_Entries[i];
      }
      out << "}";
   };
protected:
   bool
      m_Ready,
      // set if elements were added after last _PrepareForSearch, and this table
      // is meant to contain unique keys. This is maintained to diagnose a
      // potential error condition in which element addition and element removal
      // are mixed.
      m_DuplicateCheckRequired;
   FValueCompare
      m_Compare;
   typedef std::vector<FValueType>
      FEntryTable;
   FEntryTable
      m_Entries;
public:
   void swap(TLookupTable &other)
   {
      using std::swap;
      swap(m_Ready, other.m_Ready);
      swap(m_Compare, other.m_Compare);
      m_Entries.swap(other.m_Entries);
   }
public:
   // some definitions for std container compatibility.
   typedef FKeyType
      key_type;
   typedef FValueType
      value_type;
   typedef FMappedType
      mapped_type;
   typedef typename FEntryTable::iterator
      iterator;
   typedef typename FEntryTable::const_iterator
      const_iterator;
   bool empty() const { return m_Entries.empty(); }
   size_t size() const { return m_Entries.size(); }
   void clear() { _MarkAsNotReady(); m_Entries.clear(); }
   const_iterator begin() const { _PrepareForSearch(); return m_Entries.begin(); }
   const_iterator end() const { return m_Entries.end(); }
   iterator begin() { _PrepareForSearch(); return m_Entries.begin(); }
   // ^-- WARNING: *do not* modify the keys inside the table!
   // (you may modify they values inside, unless TABLEFLAGS_SortValues is set)
   iterator end() { return m_Entries.end(); }

   std::pair<iterator, iterator> equal_range(FKeyArg const &KeyArg) {
      _PrepareForSearch();
      FValueType
         SearchKey(KeyArg, FMappedType());
      return std::equal_range(m_Entries.begin(), m_Entries.end(), SearchKey, m_Compare);
   }
   std::pair<const_iterator, const_iterator> equal_range(FKeyArg const &KeyArg) const {
      return const_cast<TLookupTable*>(this)->equal_range(KeyArg);
   }
   iterator find(FKeyArg const &KeyArg) {
      return _FindIter(KeyArg);
   }
   const_iterator find(FKeyArg const &KeyArg) const {
      return _FindIter(KeyArg);
   }

   void erase(iterator it) {
      // as in std::map, `it` must point to a valid and dereferenceable element
      _sorted_table_assert(!empty() && begin() <= it && it < end());
      iterator
         itLastValidEntry = begin() + (size() - 1);

      // deleting the last element?
      if (it == itLastValidEntry) {
         // this case is special: in this case we do not need to re-sort the data.
         // we just remove the entry (happens below)
      } else {
         using std::swap;
         if (m_DuplicateCheckRequired)
            throw std::runtime_error("TLookupTable:: detected unsafe combinaton of element add/erase. Please ensure table is made consistent before deleting elements away from the end.");
         // ^-- we're moving elements around, and this may break the
         // "in case of multiple equivalent (key,value) pairs, keep the
         // last element pair which was added" property we otherwise.
         // would have
         _MarkAsNotReady();
         // swap current element with element at the end of the table
         swap(*it, *itLastValidEntry);
      }
      // delete last element of the table.
      m_Entries.pop_back();
   }

   void erase(iterator first, iterator last) {
      _sorted_table_assert(begin() <= first && last <= end());
      size_t
         nRemove = last - first,
         nLeftAfterRange = end() - last;
      if (nRemove == 0)
         return;
      if (nLeftAfterRange <= nRemove * 10) {
         // compared to the number of elements to erase, the number
         // of elements which need to be copied around in a direct
         // std::vector::erase is not excessive. So we do this:
         // Note also that in this case the order of the remaining elements
         // stays intact (if it presently is), so there is no need
         // to mark the table as NotReady.
         m_Entries.erase(first, last);
      } else {
         // entries are somewhere in the middle of the array.
         // In this case do a sequence of O(1) erase operations,
         // which, however, will break the order of the array so
         // that is needs to be re-sorted afterwards.
         // (the single-element erase takes care of that).
         for (iterator it = last - 1; ; --it) {
            erase(it);
            if (it == first) break;
         }
      }
   }

   // delete all table entries which have keys which are equivalent to
   // FKeyType(KeyArg) under the sorting predicate. Returns number of
   // entries deleted (may be 0 if no element in the table matches).
   size_t erase(FKeyArg const &KeyArg) {
      std::pair<iterator, iterator>
         range = this->equal_range(KeyArg);
      this->erase(range.first, range.second);
      return range.second - range.first;
   }

protected:
   void _CommonInit() {
      m_Ready = false;
      m_DuplicateCheckRequired = false;
   }

   // mark that the table was changed, and is currently not ready for lookup
   // operations.
   void _MarkAsNotReady() {
      m_Ready = false;
   }

   bool _NeedUniqueKeys() const {
      return !bool(TableFlags & TABLEFLAGS_AllowEqualKeys);
   }

   void _PrepareForSearch() const {
      // ^- I think it's fair to hack it const-accessible because it does not
      // technically modify the "logical" content of *this, even though it does
      // modify the internal data structures.
      return const_cast<TLookupTable*>(this)->_PrepareForSearchImpl();
   }

   void _PrepareForSearchImpl()
   {
      if (m_Ready)
         // data is already consistent---wasn't modified since last sorting
         return;

      std::stable_sort(m_Entries.begin(), m_Entries.end(), m_Compare);

      if (_NeedUniqueKeys() && m_DuplicateCheckRequired) {
         _PrepareForSearchImpl_RemoveDuplicatesAfterSort();
      }
      m_Ready = true;
   }

   void _PrepareForSearchImpl_RemoveDuplicatesAfterSort() {
      // remove duplicates: we performed a stable sort, so entries added later
      // come later in the array. So we only retain the value of the last
      // FMappedType which was added within an equal key range.
      size_t
         iOut = 0,
         iBeg = 0,
         iEnd = 0;
      while (iBeg < m_Entries.size()) {
         // find range [iBeg,iEnd) of entries with equal keys.
         while (iEnd < m_Entries.size() && m_Entries[iEnd].m_Key == m_Entries[iBeg].m_Key)
            iEnd += 1;
         iBeg = iEnd;
         _sorted_table_assert(iEnd > 0);
         if (iEnd - 1 != iOut) {
            using std::swap;
            if (bool(TableFlags & TABLEFLAGS_PreferValueSwapOverAssign)) {
               swap(m_Entries[iOut].m_Value, m_Entries[iEnd - 1].m_Value);
               // ^- for complex value types that is better... for direct PODs
               // (in particular, ints, doubles, etc.), it is worse.
            } else {
               m_Entries[iOut].m_Value = m_Entries[iEnd - 1].m_Value; // or use std::swap?
            }
         }
         iOut += 1;
      }
      m_Entries.resize(iOut);
      m_DuplicateCheckRequired = false;
   }

   iterator _FindIter(FKeyArg const &KeyArg)
   {
      _PrepareForSearch();
      _sorted_table_assert(m_Ready);
      FValueType
         SearchKey(KeyArg, FMappedType());
      typename FEntryTable::iterator
         it = std::lower_bound(m_Entries.begin(), m_Entries.end(), SearchKey, m_Compare);
      if (it != m_Entries.end() && it->m_Key != SearchKey.m_Key)
         it = m_Entries.end();
      return it;
   }

   const_iterator _FindIter(FKeyArg const &KeyArg) const {
      return const_cast<TLookupTable*>(this)->_FindIter(KeyArg);
   }

   virtual void RaiseKeyError(FKeyArg const &KeyArg) {
      std::stringstream str;
      str << "TLookupTable: key '" << KeyArg << " not found.";
      throw EKeyNotFound(str.str());
   }
};

template<class FKeyType, class FKeyArg, class FMappedType, unsigned TableFlags, class FKeyCompare>
std::ostream &operator << (std::ostream &out, TLookupTable<FKeyType, FKeyArg, FMappedType, TableFlags, FKeyCompare> const &LookupTable)
{
   LookupTable.Print(out);
   return out;
}


// note: meant to be used with "using std::swap; swap(a,b) etc." via implicit
// namespace references derived from argument-dependent lookup (see comments on swap in CxPodArray.h)
template<class FKeyType, class FKeyArg, class FMappedType, unsigned TableFlags, class FKeyCompare>
void swap(TLookupTable<FKeyType, FKeyArg, FMappedType, TableFlags, FKeyCompare> &A,
          TLookupTable<FKeyType, FKeyArg, FMappedType, TableFlags, FKeyCompare> &B)
{
   A.swap(B);
}


} // namespace ct




#endif // CX_SORTED_TABLE_H
