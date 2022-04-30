#pragma once

#include <stdint.h>
#include "succinct_config.hpp"

typedef unsigned __int128 uint128_t;

#define SUCCINCT_USE_INTRINSICS 1

#if SUCCINCT_USE_INTRINSICS

#include <xmmintrin.h>

#if defined(__GNUC__) || defined(__clang__)
#    define __INTRIN_INLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#    define __INTRIN_INLINE inline __forceinline
#else
#    define __INTRIN_INLINE inline
#endif

#endif

#if SUCCINCT_USE_POPCNT
#    if !SUCCINCT_USE_INTRINSICS
#        error "Intrinsics support needed for popcnt"
#    endif
#include <smmintrin.h>
#endif


namespace succinct {
    namespace intrinsics {


//#if SUCCINCT_USE_INTRINSICS

        __INTRIN_INLINE uint64_t byteswap64(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
            return __builtin_bswap64(value);
#elif defined(_MSC_VER)
            return _byteswap_uint64(value);
#else
#     error Unsupported platform
#endif
        }

        __INTRIN_INLINE bool bsf64(unsigned long *const index, const uint64_t mask) {
#if defined(__GNUC__) || defined(__clang__)
            if (mask) {
                *index = (unsigned long) __builtin_ctzll(mask);
                return true;
            } else {
                return false;
            }
#elif defined(_MSC_VER)
            return _BitScanForward64(index, mask) != 0;
#else
#     error Unsupported platform
#endif
        }

        __INTRIN_INLINE bool bsr64(unsigned long *const index, const uint64_t mask) {
#if defined(__GNUC__) || defined(__clang__)
            if (mask) {
                *index = (unsigned long) (63 - __builtin_clzll(mask));
                return true;
            } else {
                return false;
            }
#elif defined(_MSC_VER)
            return _BitScanReverse64(index, mask) != 0;
#else
#     error Unsupported platform
#endif
        }

        inline int clz_u128(uint128_t u) {
            uint64_t hi = u >> 64;
            uint64_t lo = (uint64_t) u;
            int retval[3] = {
                    __builtin_clzll(hi),
                    __builtin_clzll(lo) + 64,
                    128
            };
            int idx = !hi + ((!lo) & (!hi));
            return retval[idx];
        }

        __INTRIN_INLINE bool bsr128(unsigned long *const index, const uint128_t mask) {
#if defined(__GNUC__) || defined(__clang__)
            if (mask) {
                *index = (unsigned long) (127 - clz_u128(mask));
                return true;
            } else {
                return false;
            }
#elif defined(_MSC_VER)
            return _BitScanReverse64(index, mask) != 0;
#else
#     error Unsupported platform
#endif
        }


        template<typename T>
        __INTRIN_INLINE void prefetch(T const *ptr) {
            _mm_prefetch((const char *) ptr, _MM_HINT_T0);
        }

//#else /* SUCCINCT_USE_INTRINSICS */

//        template <typename T>
//        inline void prefetch(T const* /* ptr */)
//        {
//            /* do nothing */
//        }

//#endif /* SUCCINCT_USE_INTRINSICS */

#if SUCCINCT_USE_POPCNT

        __INTRIN_INLINE uint64_t popcount(uint64_t x)
        {
            return uint64_t(_mm_popcnt_u64(x));
        }

#endif /* SUCCINCT_USE_POPCNT */

    }
}
