#pragma once

#include "utils.hpp"

// Optional upstream MARISA comparison disabled in this integration.
#ifdef __COMPARE_MARISA__
# include "louds_marisa.hpp"
# include <queue>
#endif


namespace coco_paper {

// louds-sparse tree topology; used by patricia
class LS4Patricia {
 private:
  static constexpr int sample_rate_ = 256;
  static constexpr int spill_threshold_ = 128;

  struct BitVector {
    struct Block {  // packed has-child, louds and is-link bitvectors
      uint32_t rank0_{0};    // cumulative block ranks of bv<0> (i.e. has-child)
      uint32_t rank1_{0};    // cumulative block ranks of bv<1> (i.e. louds)
      uint32_t rank2_{0};    // cumulative block ranks of bv<2> (i.e. is-link)
      uint32_t select0_{0};  // index of the last block such that blocks_[select0_].rank1_ + 1 < rank0_; used for parent navigation
      uint32_t select1_{0};  // index of the last block such that blocks_[select1_].rank0_ - 1 < rank1_; used for child navigation
      uint8_t subrank0_[4]{0};
      uint8_t subrank1_[4]{0};
      uint8_t subrank2_[4]{0};
      uint64_t bits0_[4]{0};  // has-child
      uint64_t bits1_[4]{0};  // louds
      uint64_t bits2_[4]{0};  // is-link

      template <int bvnum>
      auto build_index() -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);

        uint32_t rank = 0;
        if constexpr (bvnum == 0) {
          for (int i = 0; i < 4; i++) {
            subrank0_[i] = rank;
            rank += __builtin_popcountll(bits0_[i]);
          }
          return rank;
        } else if constexpr (bvnum == 1) {
          for (int i = 0; i < 4; i++) {
            subrank1_[i] = rank;
            rank += __builtin_popcountll(bits1_[i]);
          }
          return rank;
        } else {
          for (int i = 0; i < 4; i++) {
            subrank2_[i] = rank;
            rank += __builtin_popcountll(bits2_[i]);
          }
          return rank;
        }
      }

      template <int bvnum>
      auto get(uint32_t pos) const -> bool {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
        assert(pos < 256);
        if constexpr (bvnum == 0) {
          return GET_BIT(bits0_[pos/64], pos % 64);
        } else if constexpr (bvnum == 1) {
          return GET_BIT(bits1_[pos/64], pos % 64);
        } else {
          return GET_BIT(bits2_[pos/64], pos % 64);
        }
      }

      template <int bvnum>
      auto rank1(uint32_t size) const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
        if constexpr (bvnum == 0) {
          return subrank0_[size/64] + __builtin_popcountll(bits0_[size/64] & MASK(size%64));
        } else if constexpr (bvnum == 1) {
          return subrank1_[size/64] + __builtin_popcountll(bits1_[size/64] & MASK(size%64));
        } else {
          return subrank2_[size/64] + __builtin_popcountll(bits2_[size/64] & MASK(size%64));
        }
      }

      template <int bvnum>
      auto rank0(uint32_t size) const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
        return size - rank1<bvnum>(size);
      }

      template <int bvnum>
      auto select1(uint32_t rank) const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1);
        if constexpr (bvnum == 0) {
          assert(rank <= rank1<0>());
          if (rank <= subrank0_[2]) {
            if (rank <= subrank0_[1]) {
              return 64*0 + selectll(bits0_[0], rank);
            } else {
              return 64*1 + selectll(bits0_[1], rank - subrank0_[1]);
            }
          } else {
            if (rank <= subrank0_[3]) {
              return 64*2 + selectll(bits0_[2], rank - subrank0_[2]);
            } else {
              return 64*3 + selectll(bits0_[3], rank - subrank0_[3]);
            }
          }
        } else {
          assert(rank <= rank1<1>());
          if (rank <= subrank1_[2]) {
            if (rank <= subrank1_[1]) {
              return 64*0 + selectll(bits1_[0], rank);
            } else {
              return 64*1 + selectll(bits1_[1], rank - subrank1_[1]);
            }
          } else {
            if (rank <= subrank1_[3]) {
              return 64*2 + selectll(bits1_[2], rank - subrank1_[2]);
            } else {
              return 64*3 + selectll(bits1_[3], rank - subrank1_[3]);
            }
          }
        }
      }

      template <int bvnum>
      auto select0(uint32_t rank) const -> uint32_t {
        static_assert(bvnum == 0);
        assert(rank <= rank0<0>());
        if (rank <= 64*2 - subrank0_[2]) {
          if (rank <= 64*1 - subrank0_[1]) {
            return 64*0 + select0ll(bits0_[0], rank);
          } else {
            return 64*1 + select0ll(bits0_[1], rank - (64*1 - subrank0_[1]));
          }
        } else {
          if (rank <= 64*3 - subrank0_[3]) {
            return 64*2 + select0ll(bits0_[2], rank - (64*2 - subrank0_[2]));
          } else {
            return 64*3 + select0ll(bits0_[3], rank - (64*3 - subrank0_[3]));
          }
        }
      }

      template <int bvnum>
      auto rank1() const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
        if constexpr (bvnum == 0) {
          return subrank0_[3] + __builtin_popcountll(bits0_[3]);
        } else if constexpr (bvnum == 1) {
          return subrank1_[3] + __builtin_popcountll(bits1_[3]);
        } else {
          return subrank2_[3] + __builtin_popcountll(bits2_[3]);
        }
      }

      template <int bvnum>
      auto rank0() const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
        return 256 - rank1<bvnum>();
      }
    };  // 128 bytes

    Block *blocks_{nullptr};

    uint32_t *sample_{nullptr};  // sample for has-child.select0

    uint32_t *spill_{nullptr};

    uint32_t capacity_{0};

    uint32_t size_{0};

    uint32_t rank_[3]{0};

    uint32_t spill_size_{0};

    BitVector() = default;

    ~BitVector() {
      clear();
    }

    auto size() const -> uint32_t {
      return size_;
    }

    template <int bvnum>
    auto rank1() const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
      return rank_[bvnum];
    }

    template <int bvnum>
    auto rank0() const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
      return size_ - rank1<bvnum>();
    }

    auto is_empty() const -> bool {
      return size_ == 0;
    }

    void clear() {
      free(blocks_);
      free(sample_);
      free(spill_);
      blocks_ = nullptr;
      sample_ = nullptr;
      spill_ = nullptr;
    }

    void reserve(uint32_t size) {
      if (size > capacity_) {
        capacity_ = (size + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    void load_bits(const uint64_t *bits0, const uint64_t *bits1, const uint64_t *bits2, size_t start, uint32_t size) {
      if (size_ + size > capacity_) {
        capacity_ = std::max<uint32_t>(((size_ + size) * 2 + 255) / 256 * 256, 256 * 8);
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }

      uint32_t leftover = (size_ + 255) / 256 * 256 - size_;
      if (size <= leftover) {
        copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, size);
        copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, size);
        copy_bits(blocks_[size_/256].bits2_, size_ % 256, bits2, start, size);
        size_ += size;
        return;
      }

      size_t end = start + size;
      copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, leftover);
      copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, leftover);
      copy_bits(blocks_[size_/256].bits2_, size_ % 256, bits2, start, leftover);
      size_ += leftover;
      start += leftover;
      while (start + 256 <= end) {
        copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, 256);
        copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, 256);
        copy_bits(blocks_[size_/256].bits2_, size_ % 256, bits2, start, 256);
        size_ += 256;
        start += 256;
      }
      copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, end - start);
      copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, end - start);
      copy_bits(blocks_[size_/256].bits2_, size_ % 256, bits2, start, end - start);
      size_ += end - start;
    }

    auto size_in_bytes() const -> size_t {
      return sizeof(BitVector) + sizeof(Block) * capacity_/256;
    }

    auto size_in_bits() const -> size_t {
      return size_in_bytes() * 8;
    }

    void shrink_to_fit() {
      if (capacity_ - size_ >= 256) {
        capacity_ = (size_ + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    void build(bool sample) {
      if (blocks_ == nullptr) {
        return;
      }

      shrink_to_fit();

      // clear trailing bits
      uint32_t remainder = size_ % 256;
      for (uint32_t i = (remainder + 63) / 64; i < 4; i++) {
        blocks_[size_ / 256].bits0_[i] = 0;
        blocks_[size_ / 256].bits1_[i] = 0;
        blocks_[size_ / 256].bits2_[i] = 0;
      }
      for (uint32_t i = 0; i < 4; i++) {
        blocks_[capacity_ / 256].subrank0_[i] = 0;
        blocks_[capacity_ / 256].subrank1_[i] = 0;
        blocks_[capacity_ / 256].subrank2_[i] = 0;
        blocks_[capacity_ / 256].bits0_[i] = 0;
        blocks_[capacity_ / 256].bits1_[i] = 0;
        blocks_[capacity_ / 256].bits2_[i] = 0;
      }

      // build rank index
      rank_[0] = rank_[1] = rank_[2] = 0;
      for (uint32_t i = 0; i < capacity_ / 256; i++) {
        blocks_[i].rank0_ = rank_[0];
        blocks_[i].rank1_ = rank_[1];
        blocks_[i].rank2_ = rank_[2];
        rank_[0] += blocks_[i].template build_index<0>();
        rank_[1] += blocks_[i].template build_index<1>();
        rank_[2] += blocks_[i].template build_index<2>();
      }
      blocks_[capacity_ / 256].rank0_ = rank_[0];
      blocks_[capacity_ / 256].rank1_ = rank_[1];
      blocks_[capacity_ / 256].rank2_ = rank_[2];

      // build select index
      init_select();
      if (sample) {
        init_sample();
        init_spill();
      } else {
        sample_ = nullptr;
        spill_ = nullptr;
      }
    }

    // get the `pos`-th bit of bv<bvnum>
    template <int bvnum>
    auto get(uint32_t pos) const -> bool {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
      assert(pos < size_);
      return blocks_[pos/256].template get<bvnum>(pos % 256);
    }

    // compute the rank of the first `size` bits of bv<bvnum>
    template <int bvnum>
    auto rank1(uint32_t size) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
      assert(size <= size_);

      const auto &block = blocks_[size/256];
      if constexpr (bvnum == 0) {
        return block.rank0_ + block.template rank1<0>(size % 256);
      } else if constexpr (bvnum == 1) {
        return block.rank1_ + block.template rank1<1>(size % 256);
      } else {
        return block.rank2_ + block.template rank1<2>(size % 256);
      }
    }

    template <int bvnum>
    auto rank0(uint32_t size) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);
      assert(size <= size_);
      return size - rank1<bvnum>(size);
    }

    // select the `rank`-th 1 bit from bv<bvnum>; used in build phase only
    template <int bvnum>
    auto select1(uint32_t rank) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);
      assert(rank <= rank_[bvnum]);

      uint32_t left = 0, right = capacity_/256;
      if constexpr (bvnum == 0) {
        while (left + 8 < right) {
          uint32_t mid = (left + right + 1) / 2;
          if (blocks_[mid].rank0_ < rank) {
            left = mid;
          } else {
            right = mid - 1;
          }
        }
        while (blocks_[left+1].rank0_ < rank) {
          left++;
        }
        return left*256 + blocks_[left].template select1<0>(rank - blocks_[left].rank0_);
      } else {
        while (left + 8 < right) {
          uint32_t mid = (left + right + 1) / 2;
          if (blocks_[mid].rank1_ < rank) {
            left = mid;
          } else {
            right = mid - 1;
          }
        }
        while (blocks_[left+1].rank1_ < rank) {
          left++;
        }
        return left*256 + blocks_[left].template select1<1>(rank - blocks_[left].rank1_);
      }
    }

    template <int bvnum>
    auto select0(uint32_t rank) const -> uint32_t {
      static_assert(bvnum == 0);
      assert(rank <= size_ - rank_[bvnum]);
      assert(sample_ != nullptr);

      int sample_idx = rank / sample_rate_;
      uint32_t sample = sample_[sample_idx];
      if (GET_BIT(sample, 31)) {  // spill
        uint32_t spill_idx = sample & MASK(31);
        assert(spill_idx + rank % sample_rate_ < spill_size_);
        return spill_[spill_idx + rank % sample_rate_];
      }

      int left = sample & MASK(24), right = left + (sample >> 24);
      while (right - left > 8) {
        int mid = (left + right + 1) / 2;
        if (mid*256 - blocks_[mid].rank0_ < rank) {
          left = mid;
        } else {
          right = mid - 1;
        }
      }
      while ((left+1)*256 - blocks_[left+1].rank0_ < rank) {
        left++;
      }
      uint32_t remainder = rank - (left*256 - blocks_[left].rank0_);
      uint32_t ret = left*256 + blocks_[left].template select0<0>(remainder);
      return ret;
    }

    // return the position of the closest 1 bit before position `pos` (inclusive)
    // if not found, return 0
    template <int bvnum>
    auto prev1(uint32_t pos) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);

      if (pos >= size_) {
        return size_;
      }
      uint32_t block_idx = pos / 64;
      uint32_t remainder = pos % 64;
      uint32_t dist = 0;
      if constexpr (bvnum == 0) {
        const uint64_t &block = blocks_[block_idx / 4].bits0_[block_idx % 4];
        if (block << (63 - remainder)) {
          return pos - __builtin_clzll(block << (63 - remainder));
        }
        dist += remainder + 1;
        while (pos - dist >= 0) {
          block_idx = (pos - dist) / 64;
          const uint64_t &block = blocks_[block_idx / 4].bits0_[block_idx % 4];
          if (block) {
            return pos - dist - __builtin_clzll(block);
          }
          dist += 64;
        }
      } else {
        const uint64_t &block = blocks_[block_idx / 4].bits1_[block_idx % 4];
        if (block << (63 - remainder)) {
          return pos - __builtin_clzll(block << (63 - remainder));
        }
        dist += remainder + 1;
        while (pos - dist >= 0) {
          block_idx = (pos - dist) / 64;
          const uint64_t &block = blocks_[block_idx / 4].bits1_[block_idx % 4];
          if (block) {
            return pos - dist - __builtin_clzll(block);
          }
          dist += 64;
        }
      }
    }

    // return the position of the closest 1 bit after position `pos` (inclusive)
    // if not found, return `size_`
    template <int bvnum>
    auto next1(uint32_t pos) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1 || bvnum == 2);

      if (pos >= size_) {
        return size_;
      }
      uint32_t block_idx = pos / 64;
      uint32_t remainder = pos % 64;
      uint32_t dist = 0;
      if constexpr (bvnum == 0) {
        const uint64_t &block = blocks_[block_idx / 4].bits0_[block_idx % 4];
        if (block >> remainder) {
          return pos + __builtin_ctzll(block >> remainder);
        }
        dist += 64 - remainder;
        while (pos + dist < size_) {
          block_idx = (pos + dist) / 64;
          const uint64_t &block = blocks_[block_idx / 4].bits0_[block_idx % 4];
          if (block) {
            return pos + dist + __builtin_ctzll(block);
          }
          dist += 64;
        }
      } else if constexpr (bvnum == 1) {
        const uint64_t &block = blocks_[block_idx / 4].bits1_[block_idx % 4];
        if (block >> remainder) {
          return pos + __builtin_ctzll(block >> remainder);
        }
        dist += 64 - remainder;
        while (pos + dist < size_) {
          block_idx = (pos + dist) / 64;
          const uint64_t &block = blocks_[block_idx / 4].bits1_[block_idx % 4];
          if (block) {
            return pos + dist + __builtin_ctzll(block);
          }
          dist += 64;
        }
      } else {
        const uint64_t &block = blocks_[block_idx / 4].bits2_[block_idx % 4];
        if (block >> remainder) {
          return pos + __builtin_ctzll(block >> remainder);
        }
        dist += 64 - remainder;
        while (pos + dist < size_) {
          block_idx = (pos + dist) / 64;
          const uint64_t &block = blocks_[block_idx / 4].bits2_[block_idx % 4];
          if (block) {
            return pos + dist + __builtin_ctzll(block);
          }
          dist += 64;
        }
      }
      return size_;
    }

    template <int bvnum>
    auto next0(uint32_t pos) const -> uint32_t {
      static_assert(bvnum == 0);
      uint32_t idx = pos / 64;
      uint64_t elt = ~blocks_[idx/4].bits0_[idx%4] & ~MASK(pos%64);
      while (elt == 0) {  // terminates at sentinel block since it's all 0
        idx++;
        elt = ~blocks_[idx/4].bits0_[idx%4];
      }
      uint32_t ret = idx*64 + __builtin_ctzll(elt);
      return ret;
    }

    // bvnum == 0: returns bv<0>.select1(bv<1>.rank1(pos + 1) - 1); used for parent navigation
    // bvnum == 1: returns bv<1>.select1(bv<0>.rank1(pos + 1) + 1); used for child navigation
    template <int bvnum>
    auto rs(uint32_t pos) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);
      assert(pos < size_);

      Block &block = blocks_[(pos+1)/256];
      if constexpr (bvnum == 0) {
        uint32_t rank = block.rank1_ + block.template rank1<1>((pos+1)%256) - 1;
        assert(rank != static_cast<uint32_t>(-1));
        uint32_t left = block.select0_, right = blocks_[(pos+1)/256 + 1].select0_ + 1;
        assert(rank == 0 || blocks_[left].rank0_ < rank);
        assert(blocks_[right].rank0_ >= rank);
        while (left + 8 < right) {
          uint32_t mid = (left + right + 1) / 2;
          if (blocks_[mid].rank0_ < rank) {
            left = mid;
          } else {
            right = mid - 1;
          }
        }
        while (blocks_[left+1].rank0_ < rank) {
          left++;
        }
        return left*256 + blocks_[left].template select1<0>(rank - blocks_[left].rank0_);
      } else {
        uint32_t rank = block.rank0_ + block.template rank1<0>((pos+1)%256) + 1;
        uint32_t left = block.select1_, right = blocks_[(pos+1)/256 + 1].select1_ + 1;
        assert(blocks_[left].rank1_ < rank);
        assert(blocks_[right].rank1_ >= rank);
        while (left + 8 < right) {
          uint32_t mid = (left + right + 1) / 2;
          if (blocks_[mid].rank1_ < rank) {
            left = mid;
          } else {
            right = mid - 1;
          }
        }
        while (blocks_[left+1].rank1_ < rank) {
          left++;
        }
        return left*256 + blocks_[left].template select1<1>(rank - blocks_[left].rank1_);
      }
    }

    void prefetch_block(uint32_t pos) const {
      assert(pos < size());
      _mm_prefetch((const char*)&blocks_[pos/256], _MM_HINT_T0);
      _mm_prefetch((const char*)&blocks_[pos/256] + 64, _MM_HINT_T0);
    }
   private:
    void init_select() {
      int num_blocks = capacity_ / 256;

      // build index for bv<0>.select1(bv<1>.rank1(pos) - 1)
      uint32_t i = 0, j = 0;
      while (i <= num_blocks) {
        while (j < num_blocks && blocks_[j+1].rank0_ + 1 < blocks_[i].rank1_) {
          j++;
        }
        blocks_[i].select0_ = j;
        i++;
      }

      // build index for bv<1>.select1(bv<0>.rank1(pos) + 1)
      i = j = 0;
      while (i <= num_blocks) {
        while (j < num_blocks && blocks_[j+1].rank1_ < blocks_[i].rank0_ + 1) {
          j++;
        }
        blocks_[i].select1_ = j;
        i++;
      }
    }

    void init_sample() {
      uint32_t num_blocks = capacity_ / 256;
      uint32_t num_samples = (size_ - rank_[0] + sample_rate_ - 1) / sample_rate_;

      sample_ = reinterpret_cast<uint32_t *>(realloc(sample_, (num_samples + 1) * sizeof(uint32_t)));
      sample_[0] = 0;
      for (uint32_t i = 1; i < num_samples; i++) {
        uint32_t j = sample_[i-1];
        assert(j*256 - blocks_[j].rank0_ < i*sample_rate_);
        while ((j+1)*256 - blocks_[j+1].rank0_ < i*sample_rate_) {
          j++;
        }
        sample_[i] = j;
      }
      sample_[num_samples] = num_blocks;
    }

    void init_spill() {
      int num_samples = (size_ - rank_[0] + sample_rate_ - 1) / sample_rate_;
      int num_spills = 0;

      for (int i = 0; i < num_samples; i++) {
        uint32_t dist = sample_[i+1] - sample_[i];
        if (dist >= spill_threshold_) {
          sample_[i] |= BIT(31);
          num_spills += std::min(static_cast<uint32_t>(sample_rate_), size_ - rank_[0] - i*sample_rate_) + 1;
        } else {
          sample_[i] |= (dist << 24);
        }
      }

      spill_ = reinterpret_cast<uint32_t *>(realloc(spill_, num_spills * sizeof(uint32_t)));
      spill_size_ = num_spills;
      int spilled = 0;
      for (int i = 0; i < num_samples; i++) {
        if (!GET_BIT(sample_[i], 31)) {
          continue;
        }
        uint32_t block_idx = sample_[i] & MASK(31);
        sample_[i] = BIT(31) | spilled;

        uint32_t remainder = i*sample_rate_ - (block_idx*256 - blocks_[block_idx].rank0_);
        uint32_t pos = block_idx*256 + blocks_[block_idx].template select0<0>(remainder);
        uint32_t end = std::min(static_cast<uint32_t>(sample_rate_), size_ - rank_[0] - i*sample_rate_);
        for (uint32_t j = 0; j <= end; j++) {
          spill_[spilled++] = pos;
          pos = next0<0>(pos + 1);
          assert(spilled <= spill_size_);
        }
      }
      assert(spilled == spill_size_);
    }
  };
  using bitvec_t = BitVector;

 public:
  LS4Patricia() = default;

  void add_node(const uint64_t has_child[4], const uint64_t is_link[4], size_t size) {
    assert(size <= 256);
    const uint64_t louds[4]{1};
    bv_.load_bits(has_child, louds, is_link, 0, size);
  }

  void build(bool sample = false) {
    bv_.build(sample);
  }

  auto size() const -> uint32_t {
    return bv_.size();
  }

  auto node_start(uint32_t pos) const -> uint32_t {
    return bv_.template prev1<1>(pos);
  }

  auto node_end(uint32_t pos) const -> uint32_t {
    return bv_.template next1<1>(pos + 1);
  }

  auto node_degree(uint32_t pos) const -> uint32_t {
    return node_end(pos) - node_start(pos);
  }

  auto node_id(uint32_t pos) const -> uint32_t {
    assert(louds(pos));
    return bv_.template rank1<1>(pos);
  }

  auto leaf_id(uint32_t pos) const -> uint32_t {
    assert(!has_child(pos));
    return bv_.template rank0<0>(pos);
  }

  auto link_id(uint32_t pos) const -> uint32_t {
    assert(is_link(pos));
    return bv_.template rank1<2>(pos);
  }

  auto label_id(uint32_t pos) const -> uint32_t {
    assert(!is_link(pos));
    return bv_.template rank0<2>(pos);
  }

  auto child_pos(uint32_t pos) const -> uint32_t {
    return bv_.template rs<1>(pos);
  }

  auto parent_pos(uint32_t pos) const -> uint32_t {
    return bv_.template rs<0>(pos);
  }

  auto select_leaf(uint32_t leaf_id) const -> uint32_t {
    return bv_.template select0<0>(leaf_id + 1);
  }

  auto has_child(uint32_t pos) const -> bool {
    return bv_.template get<0>(pos);
  }

  auto has_parent(uint32_t pos) const -> bool {
    return bv_.template rank1<1>(pos + 1) > 1;
  }

  auto louds(uint32_t pos) const -> bool {
    return bv_.template get<1>(pos);
  }

  auto is_link(uint32_t pos) const -> bool {
    return bv_.template get<2>(pos);
  }

  auto next_link(uint32_t pos) const -> uint32_t {
    return bv_.template next1<2>(pos);
  }

  auto num_nodes() const -> uint32_t {
    return bv_.template rank1<1>();
  }

  auto num_leaves() const -> uint32_t {
    return bv_.template rank0<0>();
  }

  auto num_children() const -> uint32_t {
    return bv_.template rank1<0>();
  }

  auto num_links() const -> uint32_t {
    return bv_.template rank1<2>();
  }

  auto size_in_bytes() const -> size_t {
    return bv_.size_in_bytes() + sizeof(uint32_t);
  }

  auto size_in_bits() const -> size_t {
    return size_in_bytes() * 8;
  }

  void clear() {
    bv_.clear();
  }

#ifdef __COMPARE_MARISA__
  void to_louds_marisa(std::unique_ptr<LoudsMarisa> &out) const {
    out = std::make_unique<LoudsMarisa>();
    std::queue<uint32_t> queue;
    queue.push(0);
    while (!queue.empty()) {
      auto pos = queue.front();
      queue.pop();
      if (pos == -1) {
        out->add_node(0);
      } else {
        uint32_t deg = node_degree(pos);
        out->add_node(deg);
        for (uint32_t i = 0; i < deg; i++) {
          if (has_child(pos + i)) {
            queue.push(child_pos(pos + i));
          } else {
            queue.push(-1);
          }
        }
      }
    }
    for (uint32_t i = 0; i < size(); i++) {
      out->add_bit(!has_child(i), is_link(i));
    }
    out->build();
  }

  void prefetch_block(uint32_t pos) const {
    bv_.prefetch_block(pos);
  }
#endif

 private:
  bitvec_t bv_;
 public:
  void serialize(std::ostream& out) const {
    archive::bits(out,size(),[&](uint32_t i){return has_child(i);});
    archive::bits(out,size(),[&](uint32_t i){return louds(i);});
    archive::bits(out,size(),[&](uint32_t i){return is_link(i);});
  }
  void load(std::istream& in) {
    std::vector<bool> children,louds_bits;
    archive::bits(in,[&](bool bit){children.push_back(bit);});
    archive::bits(in,[&](bool bit){louds_bits.push_back(bit);});
    size_t i=0;
    archive::bits(in,[&](bool bit){
      if(i>=children.size() || i>=louds_bits.size()) throw std::runtime_error("Invalid Patricia topology");
      uint64_t c=children[i],l=louds_bits[i],k=bit;
      bv_.load_bits(&c,&l,&k,0,1);++i;
    });
    if(i!=children.size() || i!=louds_bits.size()) throw std::runtime_error("Invalid Patricia topology");
    build(true);
  }

};

}  // namespace coco_paper
