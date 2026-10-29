#pragma once

#include "utils.hpp"

#include <vector>


namespace coco_paper {

class LoudsSparseCC {
 private:
  struct BitVector {
    struct Block {  // packed has-child and louds bitvectors
      uint32_t rank0_{0};    // cumulative block ranks of bv<0> (i.e. has-child)
      uint32_t rank1_{0};    // cumulative block ranks of bv<1> (i.e. louds)
      uint32_t select0_{0};  // index of the last block such that blocks_[select0_].rank1_ + 1 < rank0_; used for parent navigation
      uint32_t select1_{0};  // index of the last block such that blocks_[select1_].rank0_ - 1 < rank1_; used for child navigation
      uint8_t subrank0_[4]{0};
      uint8_t subrank1_[4]{0};
      uint64_t bits0_[4]{0};  // has-child
      uint64_t bits1_[4]{0};  // louds

      template <int bvnum>
      auto build_index() -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1);

        uint32_t rank = 0;
        if constexpr (bvnum == 0) {
          for (int i = 0; i < 4; i++) {
            subrank0_[i] = rank;
            rank += __builtin_popcountll(bits0_[i]);
          }
          return rank;
        } else {
          for (int i = 0; i < 4; i++) {
            subrank1_[i] = rank;
            rank += __builtin_popcountll(bits1_[i]);
          }
          return rank;
        }
      }

      template <int bvnum>
      auto get(uint32_t pos) const -> bool {
        static_assert(bvnum == 0 || bvnum == 1);
        assert(pos < 256);
        if constexpr (bvnum == 0) {
          return GET_BIT(bits0_[pos/64], pos % 64);
        } else {
          return GET_BIT(bits1_[pos/64], pos % 64);
        }
      }

      template <int bvnum>
      auto rank1(uint32_t size) const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1);
        if constexpr (bvnum == 0) {
          return subrank0_[size/64] + __builtin_popcountll(bits0_[size/64] & MASK(size%64));
        } else {
          return subrank1_[size/64] + __builtin_popcountll(bits1_[size/64] & MASK(size%64));
        }
      }

      template <int bvnum>
      auto rank0(uint32_t size) const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1);
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
      auto rank1() const -> uint32_t {
        static_assert(bvnum == 0 || bvnum == 1);
        if constexpr (bvnum == 0) {
          return subrank0_[3] + __builtin_popcountll(bits0_[3]);
        } else {
          return subrank1_[3] + __builtin_popcountll(bits1_[3]);
        }
      }
    };  // 88 bytes

    Block *blocks_{nullptr};

    uint32_t capacity_{0};

    uint32_t size_{0};

    uint32_t rank_[2]{0};

    BitVector() = default;

    ~BitVector() {
      clear();
    }

    auto size() const -> uint32_t {
      return size_;
    }

    template <int bvnum>
    auto rank1() const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);
      return rank_[bvnum];
    }

    auto is_empty() const -> bool {
      return size_ == 0;
    }

    void clear() {
      free(blocks_);
      blocks_ = nullptr;
    }

    void reserve(uint32_t size) {
      if (size > capacity_) {
        capacity_ = (size + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    void push_back(bool has_child, bool louds) {
      if (size_ >= capacity_) {
        capacity_ = std::max<uint32_t>((size_ * 2 + 255) / 256 * 256, 256 * 8);
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }

      if (has_child) {
        SET_BIT(blocks_[size_/256].bits0_[(size_%256)/64], size_ % 64);
      } else {
        CLEAR_BIT(blocks_[size_/256].bits0_[(size_%256)/64], size_ % 64);
      }

      if (louds) {
        SET_BIT(blocks_[size_/256].bits1_[(size_%256)/64], size_ % 64);
      } else {
        CLEAR_BIT(blocks_[size_/256].bits1_[(size_%256)/64], size_ % 64);
      }

      size_++;
    }

    void load_bits(const uint64_t *bits0, const uint64_t *bits1, size_t start, uint32_t size) {
      if (size_ + size > capacity_) {
        capacity_ = std::max<uint32_t>(((size_ + size) * 2 + 255) / 256 * 256, 256 * 8);
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }

      uint32_t leftover = (size_ + 255) / 256 * 256 - size_;
      if (size <= leftover) {
        copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, size);
        copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, size);
        size_ += size;
        return;
      }

      size_t end = start + size;
      copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, leftover);
      copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, leftover);
      size_ += leftover;
      start += leftover;
      while (start + 256 <= end) {
        copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, 256);
        copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, 256);
        size_ += 256;
        start += 256;
      }
      copy_bits(blocks_[size_/256].bits0_, size_ % 256, bits0, start, end - start);
      copy_bits(blocks_[size_/256].bits1_, size_ % 256, bits1, start, end - start);
      size_ += end - start;
    }

    void shrink_to_fit() {
      if (capacity_ - size_ >= 256) {
        capacity_ = (size_ + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    auto size_in_bytes() const -> size_t {
      return sizeof(BitVector) + sizeof(Block) * capacity_/256;
    }

    auto size_in_bits() const -> size_t {
      return size_in_bytes() * 8;
    }

    void build() {
      if (blocks_ == nullptr) {
        return;
      }

      shrink_to_fit();

      // clear trailing bits
      uint32_t remainder = size_ % 256;
      blocks_[size_/256].bits0_[remainder/64] &= MASK(remainder%64);
      blocks_[size_/256].bits1_[remainder/64] &= MASK(remainder%64);
      for (uint32_t i = (remainder + 63) / 64; i < 4; i++) {
        blocks_[size_ / 256].bits0_[i] = 0;
        blocks_[size_ / 256].bits1_[i] = 0;
      }
      for (uint32_t i = 0; i < 4; i++) {
        blocks_[capacity_ / 256].subrank0_[i] = 0;
        blocks_[capacity_ / 256].subrank1_[i] = 0;
        blocks_[capacity_ / 256].bits0_[i] = 0;
        blocks_[capacity_ / 256].bits1_[i] = 0;
      }

      // build rank index
      rank_[0] = rank_[1] = 0;
      for (uint32_t i = 0; i < capacity_ / 256; i++) {
        blocks_[i].rank0_ = rank_[0];
        blocks_[i].rank1_ = rank_[1];
        rank_[0] += blocks_[i].template build_index<0>();
        rank_[1] += blocks_[i].template build_index<1>();
      }
      blocks_[capacity_ / 256].rank0_ = rank_[0];
      blocks_[capacity_ / 256].rank1_ = rank_[1];

      // build select index
      init_select();
    }

    // get the `pos`-th bit of bv<bvnum>
    template <int bvnum>
    auto get(uint32_t pos) const -> bool {
      static_assert(bvnum == 0 || bvnum == 1);
      assert(pos < size_);
      return blocks_[pos/256].template get<bvnum>(pos % 256);
    }

    // compute the rank of the first `size` bits of bv<bvnum>
    template <int bvnum>
    auto rank1(uint32_t size) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);
      assert(size <= size_);

      const auto &block = blocks_[size/256];
      if constexpr (bvnum == 0) {
        return block.rank0_ + block.template rank1<0>(size % 256);
      } else {
        return block.rank1_ + block.template rank1<1>(size % 256);
      }
    }

    template <int bvnum>
    auto rank0(uint32_t size) const -> uint32_t {
      static_assert(bvnum == 0 || bvnum == 1);
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
      static_assert(bvnum == 0 || bvnum == 1);

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
      } else {
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
      }
      return size_;
    }

    // bvnum == 0: returns bv<0>.select1(bv<1>.rank1(pos + 1) - 1); used for child navigation
    // bvnum == 1: returns bv<1>.select1(bv<0>.rank1(pos + 1) + 1); used for parent navigation
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
  };
  using bitvec_t = BitVector;

 public:
  void reserve(uint32_t size) {
    bv_.reserve(size);
  }

  void add_node(uint64_t has_child[4], uint32_t size) {
    assert(size <= 256);
    const uint64_t louds[4]{1};
    bv_.load_bits(has_child, louds, 0, size);
  }

  void push_back(bool has_child, bool louds) {
    bv_.push_back(has_child, louds);
  }

  void build() {
    bv_.build();
  }

  auto node_start(uint32_t pos) const -> uint32_t {
    return bv_.template prev1<1>(pos);
  }

  auto node_end(uint32_t pos) const -> uint32_t {
    return bv_.template next1<1>(pos + 1);
  }

  auto prev_child(uint32_t pos) const -> uint32_t {
    return bv_.template prev1<0>(pos);
  }

  auto next_child(uint32_t pos) const -> uint32_t {
    return bv_.template next1<0>(pos);
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

  auto select_node(uint32_t node_id) const -> uint32_t {
    assert(node_id < num_nodes());
    return bv_.template select1<1>(node_id + 1);
  }

  auto child_pos(uint32_t pos) const -> uint32_t {
    return bv_.template rs<1>(pos);
  }

  auto has_child_rank1(uint32_t pos) const -> uint32_t {
    return bv_.template rank1<0>(pos + 1);
  }

  auto parent_pos(uint32_t pos) const -> uint32_t {
    return bv_.template rs<0>(pos);
  }

  auto has_child(uint32_t pos) const -> bool {
    return bv_.template get<0>(pos);
  }

  auto louds(uint32_t pos) const -> bool {
    return bv_.template get<1>(pos);
  }

  auto num_nodes() const -> uint32_t {
    return bv_.template rank1<1>();
  }

  auto num_children() const -> uint32_t {
    return bv_.template rank1<0>();
  }

  auto size() const -> uint32_t {
    return bv_.size();
  }

  auto size_in_bytes() const -> size_t {
    return bv_.size_in_bytes();
  }

  auto size_in_bits() const -> size_t {
    return bv_.size_in_bits();
  }

  auto get_level_boundaries() const -> std::vector<uint32_t> {
    std::vector<uint32_t> ret;
    uint32_t num_children = 0;
    uint32_t pos = 0;
    ret.emplace_back(pos);
    while (pos < size()) {
      pos = node_end(bv_.template select1<1>(num_children + 1));
      ret.emplace_back(pos);
      num_children = bv_.template rank1<0>(pos);
    }
    return ret;
  }

  void clear() {
    bv_.clear();
  }

 private:
  bitvec_t bv_;
 public:
  void serialize(std::ostream& out) const {
    archive::bits(out,size(),[&](uint32_t i){return has_child(i);});
    archive::bits(out,size(),[&](uint32_t i){return louds(i);});
  }
  void load(std::istream& in) {
    std::vector<bool> child;
    archive::bits(in,[&](bool bit){child.push_back(bit);});
    size_t i=0;
    archive::bits(in,[&](bool bit){if(i>=child.size()) throw std::runtime_error("Invalid LOUDS-Sparse");push_back(child[i++],bit);});
    if(i!=child.size()) throw std::runtime_error("Invalid LOUDS-Sparse");
    build();
  }

};

}  // namespace coco_paper
