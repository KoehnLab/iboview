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

#ifndef _CX_TYPES_H
#define _CX_TYPES_H

#ifndef _for_each
    #define _for_each(it,con) for ( (it) = (con).begin(); (it) != (con).end(); ++(it) )
#endif

#include <cstddef> // for size_t and ptrdiff_t

#define __STDC_CONSTANT_MACROS
// ^- ask stdint.h to include fixed-size literal macros (e.g., UINT64_C).

#ifdef USE_BOOST_STDINT_H
    #include <boost/cstdint.hpp>
    using boost::uint64_t;
    using boost::uint32_t;
    using boost::uint16_t;
    using boost::uint8_t;
    using boost::int64_t;
    using boost::int32_t;
    using boost::int16_t;
    using boost::int8_t;
#else
//     #include <cstdint>
//     // ^- not quite compatible with C++98...
//     using std::uint64_t;
//     using std::uint32_t;
//     using std::uint16_t;
//     using std::uint8_t;
//     using std::int64_t;
//     using std::int32_t;
//     using std::int16_t;
//     using std::int8_t;
    #include <stdint.h>
    // ^- technically a C99 header, but since VC 2015 also supports it now,
    //    I guess we can just go with it.
#endif


using std::size_t;
using std::ptrdiff_t;

typedef unsigned int
    uint;
typedef unsigned char
    uchar;

#include "CxDefs.h"


#ifdef USE_BOOST_INTRUSIVE_PTR
    namespace ct {
        struct FIntrusivePtrDest;
    }

    void intrusive_ptr_add_ref(ct::FIntrusivePtrDest const *pExpr);
    void intrusive_ptr_release(ct::FIntrusivePtrDest const *pExpr);


    namespace ct {
        /// A base class for reference counted objects. Classes derived from this can
        /// be used as target for boost::intrusive_ptr.
        struct FIntrusivePtrDest
        {
            FIntrusivePtrDest() : m_RefCount(0) {};
            inline virtual ~FIntrusivePtrDest() = 0;

            // making copies of the targeted objects does not affect how many references
            // are kept to them. E.g., in this:
            //
            //   pOne = FIntrusivePtr(new ...)
            //   FObj Two = *pOne;
            //
            // there is no reference to object `Two`, regardless how many there were in
            // `pOne`. Similarly, assignment the other way around, i.e., *pOne = Two,
            // would not change the references in either `*pOne` or `Two`. So just like
            // in the regular construtor, in copy construction we need to set the ref
            // count to zero, and in assignment, we need to leave the ref count alone at
            // whatever it is before making the copy.
            FIntrusivePtrDest1(FIntrusivePtrDest1 const &other) : m_RefCount(0) {}
            void operator = (FIntrusivePtrDest1 const &other) {} // <-- must leave ref count alone on copy.

            mutable int m_RefCount;
            friend void ::intrusive_ptr_add_ref(FIntrusivePtrDest const *Expr);
            friend void ::intrusive_ptr_release(FIntrusivePtrDest const *Expr);
        };

        inline FIntrusivePtrDest::~FIntrusivePtrDest()
        {
        }
    } // namespace ct

    inline void intrusive_ptr_add_ref(ct::FIntrusivePtrDest const *pExpr) {
        pExpr->m_RefCount += 1;
    }

    inline void intrusive_ptr_release(ct::FIntrusivePtrDest const *pExpr) {
        assert(pExpr->m_RefCount > 0);
        pExpr->m_RefCount -= 1;
        if (pExpr->m_RefCount == 0)
            delete pExpr;
    }

    #include <boost/intrusive_ptr.hpp>
    // ^- that's for gcc 4.7., which otherwise refuses to instantiate intrusive_ptr,
    //    since due to changes in two-phase name lookup, it cannot anymore instantiate *any*
    //    templates depending on (non-argument dependent) code which was declared only
    //    later in the translation unit.
    //    Yes, you heard right: In the case of intrusive_ptr, it means that all
    //    reference counting functions must be declared *BEFORE* intrusive_ptr.hpp
    //    is #included *ANYWHERE*. Sounds like a fun can of worms. For the
    //    beginning: don't include #intrusive_ptr.hpp directly, that might break
    //    random code.
#else
    #include "CxIntrusivePtr.h"

    namespace ct {
        typedef FIntrusivePtrDest1
            FIntrusivePtrDest;
    }

#endif


namespace ct {
    // note: the implementation of that is presently in CxAssertFail.cpp. Should be moved.
    // And I wonder if *here* is really the right place for it. It's a low level function,
    // but not actually related to any types. Maybe PodArray would be a better place?
    // Or make a separate CxBasicAlgos file set or something? I guess its small and benign
    // enough that we could just leave it in CxAssertFail...
    void ArgSort1(size_t *pOrd, double const *pVals, size_t nValSt, size_t nVals, bool Reverse = false);
}



#endif // _CX_TYPES_H

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4; indent-width 4; mixedindent off; indent-mode normal;
