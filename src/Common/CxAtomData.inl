#include <string> // also for char_traits
#include <stdexcept>
#include <sstream>
#include <ctype.h> // for tolower()
#include <stddef.h> // for size_t
#include <assert.h>

#ifndef NO_CX_SORTED_TABLE
   #include "CxSortedTable.h"
   #include <memory> // for unique_ptr
   // ^-- hm... it's literally the only thing in this file which is not
   //     available in C++98. And auto_ptr, which would work just fine
   //     for my current case, generates annoying deprecation warnings
   //     (although I agree that auto_ptr is a piece of higher insanity which
   //     should never have been designed in the first place, let alone included
   //     in the C++ standard. I really wonder what they were thinking...)
   //
   //     But well.. have to move on at some point. C++11 it is, then.
#endif

namespace ct {

EDataRequestError::EDataRequestError(std::string const &s)
   : FBase("CxAtomData.inl: Data request failed. " + s)
{
}


#ifdef INCLUDE_ABANDONED
// // don't want to add dependency to CxAssert.h, but don't want
// // the <assert> one, either.
// static void assert1(bool condition) {
//    if (!condition)
//       throw std::runtime_error("CxAtomData.inl: internal error.");
// }
#endif // INCLUDE_ABANDONED


void CheckElement(int iElement)
{
   if (iElement < 0 || size_t(iElement) >= g_nMaxElement) {
//       IvNotify(NOTIFY_Error, IvFmt("Element '%1' not recognized. Returning mass 0.", iElement));
      std::stringstream str;
      str << "Element '" << iElement << " not recognized.";
      throw EDataRequestError(str.str());
   }
}

double GetAtomicMass(int iElement, FAtomicMassType Type)
{
   if (iElement == 0)
      return 0.; // dummy atom
   CheckElement(iElement);
   assert(sizeof(g_ElementAverageMasses) == sizeof(g_ElementMostCommonIsotopeMasses));
   if (Type == ATMASS_MostCommonIsotope)
      return g_ElementMostCommonIsotopeMasses[size_t(iElement)];
   else if (Type == ATMASS_StandardAtomicWeight)
      return g_ElementAverageMasses[iElement];
   else if (Type == ATMASS_ChargeInsteadOfMass)
      return double(iElement);
   else
      throw EDataRequestError("element mass type not recognized.");
}

// case insensitive compare for simple ASCII strings. There is no standard C function
// for doing this (there are nonstandard functions in POSIX (strcasecmp) and
// win32 (stricmp), but neither one is portable)
static int my_strcicmp(char const *a, char const *b)
{
    for (;; ++ a, ++ b) {
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}

// conversion element name <-> atomic number
char const *ElementNameFromNumber(int iElement)
{
//    assert(( AtomicNumber >= 1 ) && ( AtomicNumber <= s_nElementNames ));
   CheckElement(iElement);
   return g_ElementNames[iElement];
}


static void RaiseElementNotFoundError(char const *pName)
{
   throw EDataRequestError("Failed to look up atomic number for element symbol '"+std::string(pName)+"'");
}

#ifndef NO_CX_SORTED_TABLE

   static size_t const
      s_MaxElementNameLength = 4;
      // ^- looks like char_traits::compare might be more efficient with
      // 4 than with 3. And >= 3 we need for the old Uux names.
   typedef TShortNameKey<s_MaxElementNameLength, NAMEKEY_CaseIndependent>
      FElementNameKey;
   struct FElementLookupTable : protected TLookupTable<FElementNameKey, char const *, int, 0>
   {
      FElementLookupTable();
      int iElement(char const *pElementName) const { 
         try {
            return this->at(pElementName); 
         } catch (EKeyNotFound const &e) {
            RaiseElementNotFoundError(pElementName);
            return 0;
         }
      }
   };


#ifdef INCLUDE_ABANDONED
//    template<class T>
//    struct TMiniUniquePtr {
//       TMiniUniquePtr() : m_pObject(0) {}
//       TMiniUniquePtr(T *p) : m_pObject(0) { reset(p); }
//       ~TMiniUniquePtr() { delete m_pObject; m_pObject = 0; } // <-- note: delete 0 is okay.
// 
//       void reset(T *pNewObject) {
//          T *pOldObject = m_pObject;
//          m_pObject = pNewObject;
//          delete pOldObject;
//       }
//       
//       T *get() { return m_pObject; }
//       T const *get() const { return m_pObject; }
//       operator bool () const { return m_pObject != 0; }
//       
//       T &operator * () { return *m_pObject; }
//       T const &operator * () const { return *m_pObject; }
//       T *operator -> () { return m_pObject; }
//       T const *operator -> () const { return m_pObject; }
//    protected:
//       T *m_pObject;
//    private:
//       TMiniUniquePtr(TMiniUniquePtr const &); // not implemented
//       void operator = (TMiniUniquePtr const &); // not implemented
//    };
#endif // INCLUDE_ABANDONED


   FElementLookupTable::FElementLookupTable()
   {
      reserve(200);

      assert(0 == my_strcicmp(g_ElementNames[0], "X") && 0 == my_strcicmp(g_ElementNames[1], "H"));
      // ^-- make sure array indices are equivalent to element indices (i.e.,
      // that first entry is 'X', not 'H'.
      for (size_t i = 0; i < sizeof(g_ElementNames)/sizeof(g_ElementNames[0]); ++i)
         insert_or_assign(g_ElementNames[i], i);

      // add some additional alternative names for elements. These are, in
      // particular, special names for deuterium and tritium, and older provisorial
      // names for higher synthetic elements. We might still encounter those
      // guys in input data (in particular, in basis set/ecp library files...).

      insert_or_assign("D", 1); // deuterium... from an electronic structure perspective, it's just a H.
      insert_or_assign("T", 1); // tritium
                  
      insert_or_assign("Ha", 105); // "Hahnium"... temporary name for what is now called Dubnium.
      insert_or_assign("Uut", 113); // 113 is now named 'Nihonium'
      insert_or_assign("Uuq", 114); // 114 is now named 'Flerovium'
      insert_or_assign("Uup", 115); // 114 is now named 'Moscovium'
      insert_or_assign("Uuh", 116); // 116 is now named 'Livermorium'
      insert_or_assign("Uus", 117); // 117 is now named 'Tennessine'
      insert_or_assign("Uuo", 118); // 118 is now named 'Oganesson'
      (void)my_strcicmp; // suppress unused warning
   }

   static std::unique_ptr<FElementLookupTable>
//    static TMiniUniquePtr<FElementLookupTable>
      s_pElementLookupTable;

   int ElementNumberFromName(char const *pName)
   {
      if (!s_pElementLookupTable)
         s_pElementLookupTable.reset(new FElementLookupTable);
      return s_pElementLookupTable->iElement(pName);
   }

#else
   int ElementNumberFromName(char const *pName)
   {
      if (my_strcicmp(pName, "X") == 0)
         return 0;
      // deuterium and tritium... we can't deal with those in a
      // more dignified manner. The current way may allow using element
      // labels to encode additional information (e.g., D/T or specific
      // basis sets) to element numbers, of both are stored separately.
      if (my_strcicmp(pName, "D") == 0 || my_strcicmp(pName, "T") == 0)
         return 1;
      if (my_strcicmp(pName, "Ha") == 0) // "Hahnium"... temporary name for what is now called Dubnium.
         return 105;
      if (my_strcicmp(pName, "Uut") == 0) // 113 is now named 'Nihonium'
         return 113;
      if (my_strcicmp(pName, "Uuq") == 0) // 114 is now named 'Flerovium'
         return 114;
      if (my_strcicmp(pName, "Uup") == 0) // 114 is now named 'Moscovium'
         return 115;
      if (my_strcicmp(pName, "Uuh") == 0) // 116 is now named 'Livermorium'
         return 116;
      if (my_strcicmp(pName, "Uus") == 0) // 117 is now named 'Tennessine'
         return 117;
      if (my_strcicmp(pName, "Uuo") == 0) // 118 is now named 'Oganesson'
         return 118;

      for (int i = 0; i < g_nMaxElement; ++i)
         if (0 == my_strcicmp(pName, g_ElementNames[i]))
            return i;
      RaiseElementNotFoundError(pName);
   }
#endif // NO_CX_SORTED_TABLE


double GetCovalentRadius(int iElement)
{
   CheckElement(iElement);
   return g_CovalentRadii[size_t(iElement)];
}


double GetVdwRadius(int iElement, unsigned Flags)
{
   CheckElement(iElement);
   double fVdwRadius = g_VdwRadii[size_t(iElement)];
   if (fVdwRadius == -1. && (Flags & DATAREQUEST_AssertDataExists))
      throw EDataRequestError( "VdW-radius for element '"+std::string(ElementNameFromNumber(iElement))+"' not tabulated.");
   return fVdwRadius;
}


// vdW radii for H - Cm, from 10.1002/chem.201602949 supp info, converted to
// a_bohr. these ones are approx 0.001 e^{-}/a_bohr^3 iso-density values for
// free atoms note: in some grid tests, I recomputed them from MINAO data and
// average free-atom occupancies. The results I obtained there were near
// identical to these ones:
// (2.87,2.49,4.16,4.16,3.89,3.6 ,3.36,3.21... by my ad-hoc computation, and
//  2.91,2.53,4.16,4.14,3.87,3.59,3.38,3.23... from the reference).
// So these do not appear to be terribly sensitive to the
// calculation parameters... (note: made by ~/dev/migrid/radial/vdw_table.py)
static double const s_VdwRadiiIsoDensity[97] = {
   -1.00,
   2.91,                                            2.53, 
   4.16, 4.14,        3.87, 3.59, 3.38, 3.23, 3.08, 2.95,
   4.25, 4.54,        4.52, 4.38, 4.21, 4.04, 3.89, 3.72,
   4.42, 5.10, 4.97, 4.86, 4.76,
   4.40, 4.57, 4.27, 4.20, 4.14, 4.10, 4.20, 4.40, 4.42, 4.37, 4.23, 4.14,
   4.01, 4.54, 5.27, 5.18, 5.06, 4.74, 4.61, 4.55, 4.48, 4.40, 4.06, 4.25,
   4.50, 4.65, 4.69, 4.65, 4.57, 4.50, 4.38, 4.71, 5.54, 5.37, 5.33, 5.40,
   5.37, 5.35, 5.29, 5.29, 5.23, 5.22, 5.20, 5.16, 5.14, 5.12, 5.23, 5.10,
   4.99, 4.88, 4.78, 4.71, 4.61, 4.40, 4.35, 4.27, 4.33, 4.57, 4.71, 4.72,
   4.72, 4.67, 4.59, 4.88, 5.52, 5.54, 5.46, 5.39, 5.35, 5.29, 5.25, 5.22,
   5.22
};


double GetVdwRadius_IsoDensity(int iElement) {
   size_t
      nRadiiInTable = sizeof(s_VdwRadiiIsoDensity)/sizeof(s_VdwRadiiIsoDensity[0]);
   if (iElement <= 0)
      return s_VdwRadiiIsoDensity[1]; // replace by H
   else if (size_t(iElement) >= nRadiiInTable)
      // replace by last element we actually have
      return s_VdwRadiiIsoDensity[nRadiiInTable-1];
   else
      return s_VdwRadiiIsoDensity[size_t(iElement)]; // <- note: NOT -1! There is a dummy item at [0] such that table index and element correspond 1:1.
}


} // namespace ct

// kate: syntax c++;
