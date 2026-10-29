#pragma once

#include "archive.hpp"
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cmath>
#include <succinct/bit_vector.hpp>
#include <iostream>
#include <cassert>
#include <sys/time.h>
#include <chrono>

#include <immintrin.h>

#include "global_parameters.hpp"
#include "sux/support/common.hpp"

#define __SUX_SELECT__

#define BIT(n)  (1ull << (n))
#define MASK(n) (BIT(n) - 1ull)
#define GET_BIT(bits, pos) ((bits) & BIT(pos))
#define SET_BIT(bits, pos) ((bits) |= BIT(pos))
#define CLEAR_BIT(bits, pos) ((bits) &= ~BIT(pos))

#define EXPECT(cond) \
do { \
  if (!(cond)) { \
    printf("%s:%d: assertion failed: %s\n", __FILE__, __LINE__, (#cond)); \
    fflush(stdout); \
    exit(-1); \
  } \
} while (0)


// gcc attributues
#define __ALWAYS_INLINE __attribute__((always_inline))
#define __NOINLINE __attribute__((noinline))
#define __HOT __attribute__((hot))
#ifdef __PROFILE__
# define __NOINLINE_IF_PROFILE __NOINLINE  // used for profiling specific functions
#else
# define __NOINLINE_IF_PROFILE
#endif



namespace coco_paper {

inline ds2i::global_parameters params;

static constexpr uint8_t terminator_ = 0;
static constexpr int mb_bits = 1024*1024*8;
static constexpr int mb_bytes = 1024*1024;

inline auto usec() -> uint64_t {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

using uint128_t = unsigned __int128;
using code_t = uint64_t;
#ifdef COCO_HAS_UINT256
constexpr int width_in_bits(uint256_t x) {
    for(int word=3;word>=0;--word) {
        const uint64_t value=uint64_t(x>>(word*64));
        if(value) return word*64+64-__builtin_clzll(value);
    }
    return 1;
}
constexpr int log2_ceil(uint256_t x) {
    assert(x);
    return x == 1 ? 0 : width_in_bits(x-1);
}
#endif

// passing 0 is undefined behavior
constexpr auto clz128(uint128_t x) -> int {
  uint64_t high = static_cast<uint64_t>(x >> 64);
  uint64_t low = static_cast<uint64_t>(x);
  if (high) {
    return __builtin_clzll(high);
  } else {
    return 64 + __builtin_clzll(low);
  }
}

// passing 0 is undefined behavior
constexpr auto ctz128(uint128_t x) -> int {
  uint64_t high = static_cast<uint64_t>(x >> 64);
  uint64_t low = static_cast<uint64_t>(x);
  if (low) {
    return __builtin_ctzll(low);
  } else {
    return 64 + __builtin_ctzll(high);
  }
}

constexpr auto width_in_bits(uint32_t x) -> int {
  if (x == 0) {
    return 1;
  }
  return 32 - __builtin_clz(x);
}

constexpr auto width_in_bits(uint64_t x) -> int {
  if (x == 0) {
    return 1;
  }
  return 64 - __builtin_clzll(x);
}

constexpr auto width_in_bits(uint128_t x) -> int {
  if (x == 0) {
    return 1;
  }
  return 128 - clz128(x);
}

constexpr auto log2_ceil(uint32_t x) -> int {
  assert(x > 0);
  if ((x & (x-1)) == 0) {
    return __builtin_ctz(x);
  }
  return 32 - __builtin_clz(x);
}

constexpr auto log2_ceil(uint64_t x) -> int {
  assert(x > 0);
  if ((x & (x-1)) == 0) {
    return __builtin_ctzll(x);
  }
  return 64 - __builtin_clzll(x);
}

constexpr auto log2_ceil(uint128_t x) -> int {
  assert(x > 0);
  if ((x & (x-1)) == 0) {
    return ctz128(x);
  }
  return 128 - clz128(x);
}

// 64-bit select 1; refer to https://graphics.stanford.edu/~seander/bithacks.html
inline auto selectll(uint64_t bits, int rank) -> int {
#ifdef __SUX_SELECT__
  return rank == 0 ? 0 : sux::select64(bits, rank - 1);
#else
  int pos;
  uint64_t c2, c4, c8, c16;  // popcounts of each 2-bit, 4-bit, 8-bit and 16-bit chunks
  uint32_t mid;  // popcount of half of the current chunk; used for binary search

  c2 = bits - ((bits >> 1) & ~0UL/3);
  c4 = (c2 & ~0UL/5) + ((c2 >> 2) & ~0UL/5);
  c8 = (c4 + (c4 >> 4)) & ~0UL/0x11;
  c16 = (c8 + (c8 >> 8)) & ~0UL/0x101;

  // binary search, starting with chunk size 64
  pos = 0;
  mid = ((c16 >> 16) + c16) & 0xff;
  if (mid < rank) {
    pos += 32;
    rank -= mid;
  }

  mid = (c16 >> pos) & 0xff;
  if (mid < rank) {
    pos += 16;
    rank -= mid;
  }

  mid = (c8 >> pos) & 0xf;
  if (mid < rank) {
    pos += 8;
    rank -= mid;
  }

  mid = (c4 >> pos) & 0x7;
  if (mid < rank) {
    pos += 4;
    rank -= mid;
  }

  mid = (c2 >> pos) & 0x3;
  if (mid < rank) {
    pos += 2;
    rank -= mid;
  }

  mid = (bits >> pos) & 0x1;
  if (mid < rank) {
    pos++;
  }

  return pos;
#endif
}

// 64-bit select 0; refer to https://graphics.stanford.edu/~seander/bithacks.html
inline auto select0ll(uint64_t bits, int rank) -> int {
  return selectll(~bits, rank);
}

inline auto rank00ll(uint64_t bits, bool prev) -> int {
  uint64_t v = bits | (bits << 1) | prev;
  return __builtin_popcountll(~v);
}

inline auto select00ll(uint64_t bits, int rank, bool prev) -> int {
  uint64_t v = bits | (bits << 1) | prev;
  return select0ll(v, rank);
}

inline auto rank(const uint64_t *bits, size_t begin, size_t end) -> size_t {
  uint64_t block = bits[begin / 64] & ~MASK(begin % 64);
  size_t ret = 0;
  int first_block_idx = begin / 64, last_block_idx = end / 64;
  for (int i = first_block_idx; i < last_block_idx; i++) {
    ret += __builtin_popcountll(block);
    block = bits[i + 1];
  }
  return ret + __builtin_popcountll(block & MASK(end % 64));
};

/**
 * bit copying
 */
// copy `len` bits starting from the `src_start`-th bit of `src` to the `dst_start`-th bit of `dst`
// returns the number of bits copied (which is always `len`)
inline auto copy_bits(uint64_t *dst, size_t dst_start, const uint64_t *src, size_t src_start, size_t len) -> size_t {
  const size_t total=len;
  while(len) {
    const size_t d=dst_start%64,s=src_start%64,n=std::min(len,64-d);
    uint64_t value=src[src_start/64]>>s;
    if(s+n>64) value|=src[src_start/64+1]<<(64-s);
    const uint64_t mask=n==64?~uint64_t(0):(uint64_t(1)<<n)-1;
    const uint64_t low=d?(dst[dst_start/64]&((uint64_t(1)<<d)-1)):0;
    dst[dst_start/64]=low|((value&mask)<<d);
    dst_start+=n;src_start+=n;len-=n;
  }
  return total;
}

// similar to `copy_bits`, except that bits are traversed and copied in reverse order
inline auto copy_bits_backward(uint64_t *dst, size_t dst_end, const uint64_t *src, size_t src_end, size_t len) -> size_t {
  assert(dst_end >= len);
  assert(src_end >= len);

  int offset1 = dst_end % 64, offset2 = src_end % 64;
  int remainder1 = (dst_end - len) % 64, remainder2 = (src_end - len) % 64;
  int last_block1 = dst_end / 64, first_block1 = (dst_end - len) / 64;
  int last_block2 = src_end / 64, first_block2 = (src_end - len) / 64;

  if (offset1 == offset2) {
    if (offset1 == 0) {  // prevent array out of bound
      uint64_t high_bits = src[last_block2 - 1];
      for (int i = 1; i < last_block1 - first_block1; i++) {  // TODO
        dst[last_block1 - i] = high_bits;
        high_bits = src[last_block2 - i - 1];
      }
      dst[first_block1] = high_bits & ~MASK(remainder1);
    } else {
      uint64_t high_bits = dst[last_block2] & ~MASK(offset1);
      for (int i = 0; i < last_block1 - first_block1; i++) {
        dst[last_block1 - i] = (high_bits | (src[last_block2 - i] & MASK(offset2)));
        high_bits = src[last_block2 - i - 1] & ~MASK(offset1);
      }
      dst[first_block1] = (high_bits | (src[first_block2] & MASK(offset2))) & ~MASK(remainder1);
    }
  } else if (offset1 < offset2) {
    int rshift = offset2 - offset1;

    if (offset1 == 0) {  // prevent array out of bound
      uint64_t high_bits = (src[last_block2] << (64 - rshift));
      for (int i = 1; i < last_block2 - first_block2; i++) {  // TODO
        dst[last_block1 - i] = (high_bits | (src[last_block2 - i] >> rshift));
        high_bits = ((src[last_block2 - i] << (64 - rshift)) | (src[last_block2 - i - 1] >> rshift));
      }
      if (last_block2 - first_block2 == last_block1 - first_block1) {
        dst[first_block1] = (high_bits | ((src[first_block2] & ~MASK(remainder2)) >> rshift));
      } else {  // one more destination block
        assert(last_block2 - first_block2 == last_block1 - first_block1 - 1);
        dst[first_block1 + 1] = (high_bits | (src[first_block2] >> rshift));
        dst[first_block1] = (src[first_block2] & ~MASK(remainder2)) << (64 - rshift);
      }
    } else {
      uint64_t high_bits = dst[last_block1] & ~MASK(offset1);
      for (int i = 0; i < last_block2 - first_block2; i++) {
        dst[last_block1 - i] = (high_bits | ((src[last_block2 - i] & MASK(offset2)) >> rshift));
        high_bits = ((src[last_block2 - i] << (64 - rshift)) | (src[last_block2 - i - 1] >> offset2) << offset1);
      }
      if (last_block2 - first_block2 == last_block1 - first_block1) {
        dst[first_block1] = (high_bits | ((src[first_block2] & MASK(offset2) & ~MASK(remainder2)) >> rshift));
      } else {  // one more destination block
        assert(last_block2 - first_block2 == last_block1 - first_block1 - 1);
        dst[first_block1 + 1] = (high_bits | ((src[first_block2] & MASK(offset2)) >> rshift));
        dst[first_block1] = (src[first_block2] & ~MASK(remainder2)) << (64 - rshift);
      }
    }
  } else {  // offset1 > offset2
    int lshift = offset1 - offset2;

    if (offset2 == 0) {  // prevent array out of bound
      uint64_t high_bits = dst[last_block1] & ~MASK(offset1);
      for (int i = 0; i < last_block1 - first_block1; i++) {
        dst[last_block1 - i] = (high_bits | (src[last_block2 - i - 1] >> (64 - lshift)));
        high_bits = src[last_block2 - i - 1] << lshift;
      }
      if (last_block2 - first_block2 == last_block1 - first_block1) {
        dst[first_block1] = high_bits & ~MASK(remainder1);
      } else {  // one more source block
        assert(last_block2 - first_block2 == last_block1 - first_block1 + 1);
        dst[first_block1] = (high_bits | (src[first_block2] >> (64 - lshift))) & ~MASK(remainder1);
      }
    } else {
      uint64_t high_bits = dst[last_block1] & ~MASK(offset1);
      for (int i = 0; i < last_block1 - first_block1; i++) {
        uint64_t low_bits = ((src[last_block2 - i] & MASK(offset2)) << lshift) | (src[last_block2 - i - 1] >> (64 - lshift));
        dst[last_block1 - i] = (high_bits | low_bits);
        high_bits = (src[last_block2 - i - 1] & ~MASK(offset2)) << lshift;
      }
      if (last_block2 - first_block2 == last_block1 - first_block1) {
        uint64_t low_bits = ((src[first_block2] & MASK(offset2)) << lshift);
        dst[first_block1] = (high_bits | low_bits) & ~MASK(remainder1);
      } else {  // one more source block
        assert(last_block2 - first_block2 == last_block1 - first_block1 + 1);
        uint64_t low_bits = ((src[first_block2 + 1] & MASK(offset2)) << lshift) | (src[first_block2] >> (64 - lshift));
        dst[first_block1] = (high_bits | low_bits) & ~MASK(remainder1);
      }
    }
  }

  return len;
}

}  // namespace coco_paper
