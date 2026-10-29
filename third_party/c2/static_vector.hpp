#pragma once

#include "utils.hpp"

#include <type_traits>
#include <cstring>
#include <cassert>


namespace coco_paper {

template <typename Value>
class StaticVector {
 public:
  using value_type = Value;

  StaticVector() = default;

  ~StaticVector() {
    free(values_);
  }

  auto size() const -> size_t {
    return size_;
  }

  auto data() const -> const value_type * {
    return values_;
  }

  auto data() -> value_type * {
    return values_;
  }

  auto is_empty() const -> bool {
    return size_ == 0;
  }

  auto operator[](size_t pos) const -> const value_type & {
    return at(pos);
  }

  auto operator[](size_t pos) -> value_type & {
    return at(pos);
  }

  auto at(size_t pos) const -> const value_type & {
    assert(pos < size_);
    return values_[pos];
  }

  auto at(size_t pos) -> value_type & {
    assert(pos < size_);
    return values_[pos];
  }

  void reserve(size_t size) {
    if (size > capacity_) {
      capacity_ = size;
      values_ = reinterpret_cast<value_type *>(realloc(values_, sizeof(value_type) * capacity_));
    }
  }

  void emplace_back(const value_type &value) {
    if (size_ == capacity_) {
      capacity_ = std::max<uint64_t>(size_ * 3 / 2, 8ull);
      values_ = reinterpret_cast<value_type *>(realloc(values_, sizeof(value_type) * capacity_));
    }
    values_[size_++] = value;
  }

  void emplace_back(value_type &&value) {
    if (size_ == capacity_) {
      capacity_ = std::max<uint64_t>(size_ * 3 / 2, 8ull);
      values_ = reinterpret_cast<value_type *>(realloc(values_, sizeof(value_type) * capacity_));
    }
    values_[size_++] = value;
  }

  auto find(const value_type &value, size_t begin, size_t end) const -> size_t {
    assert(end <= size());

    if constexpr (std::is_same<value_type, uint8_t>::value) {  // special optimization for uint8_t: SIMD search
      __m128i val_vec = _mm_set1_epi8(value);

      while (end - begin >= 16) {  // search in batch of 16
        __m128i src_vec = _mm_loadu_si128(reinterpret_cast<const __m128i *>(values_ + begin));
        __m128i res_vec = _mm_cmpeq_epi8(src_vec, val_vec);
        uint32_t mask = _mm_movemask_epi8(res_vec);
        if (mask) {
          size_t idx = begin + __builtin_ctz(mask);
          return idx;
        }
        begin += 16;
      }
    }
    while (begin < end) {  // if Value is uint8_t, this loop deals with the remainder
      if (values_[begin] == value) {
        break;
      }
      begin++;
    }
    return begin;
  }

  void clear() {
    free(values_);
    values_ = nullptr;
  }

  void shrink_to_fit() {
    if (capacity_ > size_) {
      values_ = reinterpret_cast<value_type *>(realloc(values_, sizeof(value_type) * size_));
      capacity_ = size_;
    }
  }

  auto size_in_bytes() const -> size_t {
    return sizeof(StaticVector) + sizeof(value_type) * capacity_;
  }

  auto size_in_bits() const -> size_t {
    return size_in_bytes() * 8;
  }

  size_t capacity_{0};
  size_t size_{0};
  value_type *values_{nullptr};
};

}  // namespace coco_paper
