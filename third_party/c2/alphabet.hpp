#pragma once

#include "utils.hpp"


namespace coco_paper {

struct Alphabet {
  uint32_t rank_{0};
  uint8_t  subrank_[4]{0};
  uint64_t bits_[4]{0};

  Alphabet() = default;

  Alphabet(const Alphabet &other) {
    load(other.bits_, 0, 256);
    build_index();
  }

  auto operator=(const Alphabet &other) -> Alphabet & {
    load(other.bits_, 0, 256);
    build_index();
    return *this;
  }

  auto operator+=(const Alphabet &other) -> Alphabet & {
    bits_[0] |= other.bits_[0];
    bits_[1] |= other.bits_[1];
    bits_[2] |= other.bits_[2];
    bits_[3] |= other.bits_[3];
    build_index();
    return *this;
  }

  auto operator+(const Alphabet &other) const -> Alphabet {
    Alphabet ret;
    ret.bits_[0] = bits_[0] | other.bits_[0];
    ret.bits_[1] = bits_[1] | other.bits_[1];
    ret.bits_[2] = bits_[2] | other.bits_[2];
    ret.bits_[3] = bits_[3] | other.bits_[3];
    ret.build_index();
    return ret;
  }

  void clear() {
    std::memset(this, 0, sizeof(Alphabet));
  }

  void load(const uint64_t *bitmap, size_t start, size_t size) {
    copy_bits(bits_, 0, bitmap, start, size);
  }

  void serialize(uint64_t *bitmap, size_t start, size_t size) {
    copy_bits(bitmap, start, bits_, 0, size);
  }

  void build_index() {
    rank_ = 0;
    for (int i = 0; i < 4; i++) {
      subrank_[i] = rank_;
      rank_ += __builtin_popcountll(bits_[i]);
    }
  }

  auto get(uint8_t pos) const -> bool {
    return GET_BIT(bits_[pos/64], pos % 64);
  }

  void set1(uint8_t pos) {
    SET_BIT(bits_[pos/64], pos % 64);
  }

  // equivalent to rank
  auto encode(uint8_t label) const -> uint8_t {
    assert(get(label));
    return subrank_[label/64] + __builtin_popcountll(bits_[label/64] & MASK(label % 64));
  }

  // equivalent to select
  auto decode(uint8_t label) const -> uint8_t {
    assert(label < rank_);
    if (label <= subrank_[2]) {
      if (label <= subrank_[1]) {
        return 64*0 + selectll(bits_[0], label);
      } else {
        return 64*1 + selectll(bits_[1], label - subrank_[1]);
      }
    } else {
      if (label <= subrank_[3]) {
        return 64*2 + selectll(bits_[2], label - subrank_[2]);
      } else {
        return 64*3 + selectll(bits_[3], label - subrank_[3]);
      }
    }
  }

  auto alphabet_size() const -> uint32_t {
    return rank_;
  }
};

}  // namespace coco_paper
