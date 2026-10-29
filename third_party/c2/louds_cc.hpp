#pragma once

#include "utils.hpp"


#define __COMPARE_COCO__
#ifdef __COMPARE_COCO__
# include "louds_sux.hpp"
#endif


namespace coco_paper {

class LoudsCC {
 private:
  static constexpr int spill_threshold_ = 64;
  static_assert(spill_threshold_ <= 128);

  struct BitVector {
    struct Block {
      uint32_t rank1_;       // rank1 before block
      uint32_t rank00_;      // rank00 before block; used for value navigation
      /** select index for computing select0(rank1(pos))
       * msb: indicates if the select index is spilled or not
       * not spilled - bits 0-23: maximum block index with rank0 < this->rank1_; bits 24-30: distance to the upper bound
       * spilled - bits 0-30: spill index (only half the bits are 0, so no overflow)
       */
      uint32_t select0_;
      uint8_t subrank_[4];    // accumulative rank1 before each uint64_t element
      uint8_t subrank00_[4];  // accumulative rank00 before each uint64_t element (block-relative)
      /** carries_ packs 4 carry bits into one byte:
       *  bit 0: carry into bits_[0] (= MSB of previous block's bits_[3])
       *  bit i (i=1,2,3): carry into bits_[i] (= MSB of bits_[i-1])
       *  used by rank00(size) to avoid loading bits_[word-1] entirely
       */
      uint8_t carries_;
      uint8_t _pad[3];
      uint64_t bits_[4];      // 256 actual bits

      auto build_index() -> uint32_t {
        uint32_t subrank = 0;
        for (int i = 0; i < 4; i++) {
          subrank_[i] = subrank;
          subrank += __builtin_popcountll(bits_[i]);
        }
        return subrank;
      }

      auto build_subrank00(bool prev) -> void {
        carries_ = prev;  // bit 0: carry into word 0
        uint32_t acc = 0;
        for (int i = 0; i < 4; i++) {
          subrank00_[i] = acc;
          acc += rank00ll(bits_[i], prev);
          prev = (bits_[i] >> 63);
          if (i < 3) carries_ |= (prev << (i + 1));  // bits 1-3: carries into words 1-3
        }
      }

      auto get(uint32_t pos) const -> bool {
        assert(pos < 256);
        return GET_BIT(bits_[pos/64], pos%64);
      }

      auto rank1(uint32_t size) const -> uint32_t {
        return subrank_[size/64] + __builtin_popcountll(bits_[size/64] & MASK(size%64));
      }

      auto rank0(uint32_t size) const -> uint32_t {
        return size - rank1(size);
      }

      auto rank00(bool prev) const -> uint32_t {
        return rank00ll(bits_[0], prev) +
               rank00ll(bits_[1], bits_[0] >> 63) +
               rank00ll(bits_[2], bits_[1] >> 63) +
               rank00ll(bits_[3], bits_[2] >> 63);
      }

      auto rank00(uint32_t size, bool prev) const -> uint32_t {
        assert(size < 256);
        uint32_t word = size / 64;
        bool carry = (word > 0) ? (bool)(bits_[word - 1] >> 63) : prev;
        return subrank00_[word] + rank00ll(bits_[word] | ~MASK(size % 64), carry);
      }

      // fast path: all data from this block only, no cross-block or cross-word loads
      auto rank00(uint32_t size) const -> uint32_t {
        assert(size < 256);
        uint32_t word = size / 64;
        bool carry = (carries_ >> word) & 1;
        return subrank00_[word] + rank00ll(bits_[word] | ~MASK(size % 64), carry);
      }

      auto select0(uint32_t rank) const -> uint32_t {
        if (rank <= 64*2 - subrank_[2]) {
          if (rank <= 64*1 - subrank_[1]) {
            return 64*0 + select0ll(bits_[0], rank);
          } else {
            return 64*1 + select0ll(bits_[1], rank - (64*1 - subrank_[1]));
          }
        } else {
          if (rank <= 64*3 - subrank_[3]) {
            return 64*2 + select0ll(bits_[2], rank - (64*2 - subrank_[2]));
          } else {
            return 64*3 + select0ll(bits_[3], rank - (64*3 - subrank_[3]));
          }
        }
      }
    };  // 56 bytes

    Block *blocks_{nullptr};

    // spill list when 0s are sparse, stores precomputed select0 results for EVERY possible input
    uint32_t *spill0_{nullptr};

    uint32_t spill0_size_{0};

    uint32_t capacity_{0};

    uint32_t size_{0};

    uint32_t rank1_{0};

    uint32_t rank00_{0};

    BitVector() = default;

    BitVector(uint32_t size) {
      reserve(size);
    }

    ~BitVector() {
      clear();
    }

    auto size() const -> uint32_t {
      return size_;
    }

    auto rank1() const -> uint32_t {
      return rank1_;
    }

    auto rank0() const -> uint32_t {
      return size_ - rank1_;
    }

    auto rank00() const -> uint32_t {
      return rank00_;
    }

    auto is_empty() const -> bool {
      return size_ == 0;
    }

    void reserve(uint32_t size) {
      if (size > capacity_) {
        capacity_ = (size + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    void append1() {
      if (size_ == capacity_) {
        capacity_ = std::max<uint32_t>((size_ * 2 + 255) / 256 * 256, 256 * 8);
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
      SET_BIT(blocks_[size_/256].bits_[(size_%256)/64], size_ % 64);
      size_++;
    }

    void append0() {
      if (size_ == capacity_) {
        capacity_ = std::max<uint32_t>((size_ * 2 + 255) / 256 * 256, 256 * 8);
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
      CLEAR_BIT(blocks_[size_/256].bits_[(size_%256)/64], size_ % 64);
      size_++;
    }

    void shrink_to_fit() {
      if (size_ + 256 <= capacity_) {
        capacity_ = (size_ + 255) / 256 * 256;
        blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
      }
    }

    void clear() {
      free(blocks_);
      free(spill0_);
      blocks_ = nullptr;
      spill0_ = nullptr;
    }

    auto size_in_bytes() const -> size_t {
      return sizeof(Block) * (capacity_/256 + 1) + sizeof(uint32_t) * spill0_size_ + sizeof(BitVector);
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
      blocks_[size_/256].bits_[remainder/64] &= MASK(remainder%64);
      for (uint32_t i = (remainder + 63) / 64; i < 4; i++) {
        blocks_[size_/256].bits_[i] = 0;
      }
      for (uint32_t i = 0; i < 4; i++) {
        blocks_[capacity_/256].subrank_[i] = 0;
        blocks_[capacity_/256].subrank00_[i] = 0;
      }

      rank1_ = rank00_ = 0;
      bool prev = true;
      for (uint32_t i = 0; i < capacity_ / 256; i++) {
        uint32_t block_rank1 = blocks_[i].build_index();
        blocks_[i].rank1_ = rank1_;
        rank1_ += block_rank1;

        blocks_[i].rank00_ = rank00_;
        blocks_[i].build_subrank00(prev);
        rank00_ += blocks_[i].rank00(prev);
        prev = blocks_[i].bits_[3] >> 63;
      }
      if (size_ < capacity_) {
        rank00_ -= capacity_ - size_ - get(size_ - 1);
      }
      blocks_[capacity_/256].rank1_ = rank1_;
      blocks_[capacity_/256].rank00_ = rank00_;

      init_select0();
      init_spill0();
    }

    // get the `pos`-th bit
    auto get(uint32_t pos) const -> bool {
      assert(pos < size_);
      bool ret = blocks_[pos/256].get(pos%256);
      return ret;
    }

    // return the number of 1 bits in the first `size` bits
    auto rank1(uint32_t size) const -> uint32_t {
      assert(size <= size_);
      const auto &block = blocks_[size / 256];
      uint32_t ret = block.rank1_ + block.rank1(size % 256);
      return ret;
    }

    // return the number of 0 bits in the first `size` bits
    auto rank0(uint32_t size) const -> uint32_t {
      return size - rank1(size);
    }

    // return the number of 00 patterns in the first `size` bits
    auto rank00(uint32_t size) const -> uint32_t {
      assert(size <= size_);
      const auto &block = blocks_[size / 256];
      return block.rank00_ + block.rank00(size % 256);  // uses precomputed carry_
    }

    // return the position of the closest 0 bit after position `pos` (inclusive)
    // if not found, return `size_`
    auto next0(uint32_t pos) const -> uint32_t {
      // TODO: maybe search by block using rank, then by word?
      uint32_t idx = pos / 64;
      uint64_t elt = ~blocks_[idx/4].bits_[idx%4] & ~MASK(pos%64);
      while (elt == 0) {
        idx++;
        elt = ~blocks_[idx/4].bits_[idx%4];
      }
      uint32_t ret = idx*64 + __builtin_ctzll(elt);
      return ret;
    }

    // return select0(rank1(pos))
    auto r1s0(uint32_t pos) const -> uint32_t {
      assert(pos < size_);
      Block &block = blocks_[pos / 256];
      uint32_t rank = rank1(pos);
      if (GET_BIT(block.select0_, 31)) {
        uint32_t spill_idx = block.select0_ & MASK(31);
        assert(spill_idx + (rank - block.rank1_) < spill0_size_);
        return spill0_[spill_idx + (rank - block.rank1_)];
      }

      int left = block.select0_ & MASK(24), right = left + (block.select0_ >> 24) + 1;
      assert(left*256 - blocks_[left].rank1_ < rank);
      assert(right*256 - blocks_[right].rank1_ >= rank);
      while (right - left > 8) {
        int mid = (left + right + 1) / 2;
        if (mid*256 - blocks_[mid].rank1_ < rank) {
          left = mid;
        } else {
          right = mid - 1;
        }
      }
      while ((left+1)*256 - blocks_[left+1].rank1_ < rank) {
        left++;
      }
      uint32_t remainder = rank - (left*256 - blocks_[left].rank1_);
      uint32_t ret = left*256 + blocks_[left].select0(remainder);
      return ret;
    }

    // used for debugging only
    auto select0(uint32_t rank) const -> uint32_t {
      assert(rank <= rank0());

      int left = 0;
      int right = size_ / 256;

      while (left + 8 < right) {  // binary search when the gap is large
        int mid = (left + right + 1) / 2;
        if (mid*256 - blocks_[mid].rank1_ < rank) {
          left = mid;
        } else {
          right = mid - 1;
        }
      }
      int block_idx = left;  // switch to linear search
      while ((block_idx+1)*256 - blocks_[block_idx+1].rank1_ < rank) {
        block_idx++;
      }

      assert((rank == 0 || block_idx*256 - blocks_[block_idx].rank1_ < rank));
      assert((block_idx+1)*256 - blocks_[block_idx+1].rank1_ >= rank);

      auto &block = blocks_[block_idx];
      uint32_t remainder = rank - (block_idx*256 - block.rank1_);
      return block_idx*256 + block.select0(remainder);
    }

   private:
    void init_select0() {
      int num_blocks = capacity_ / 256;

      // build select index
      int i = 0, j = 0;
      while (i < num_blocks) {
        while ((j+1)*256 - blocks_[j+1].rank1_ < blocks_[i].rank1_) {
          j++;
        }
        // j is the maximum block index satisfying blocks_[j].rank0_ < blocks_[i].rank1_
        blocks_[i].select0_ = j;
        i++;
      }
      blocks_[num_blocks].select0_ = capacity_/256;  // sentinel
    }

    void init_spill0() {
      int num_blocks = capacity_ / 256;
      int num_spills = 0;
      for (int i = 0; i < num_blocks; i++) {
        uint32_t dist = blocks_[i+1].select0_ - blocks_[i].select0_;
        if (dist >= spill_threshold_) {
          blocks_[i].select0_ |= BIT(31);
          // also precompute blocks_[i].rank1_ just in case block starts with 0 bits
          num_spills += blocks_[i+1].rank1_ - blocks_[i].rank1_ + 1;
        } else {
          blocks_[i].select0_ |= (dist << 24); // # bits fits in 32 bits, so # blocks fits in 24 bits
        }
      }

      spill0_ = reinterpret_cast<uint32_t *>(realloc(spill0_, num_spills * sizeof(uint32_t)));
      spill0_size_ = num_spills;
      int spilled = 0;
      for (int i = 0; i < num_blocks; i++) {
        if (!GET_BIT(blocks_[i].select0_, 31)) {
          continue;
        }
        uint32_t block_idx = blocks_[i].select0_ & MASK(31);
        blocks_[i].select0_ = BIT(31) | spilled;

        uint32_t remainder = blocks_[i].rank1_ - (block_idx*256 - blocks_[block_idx].rank1_);
        uint32_t pos = block_idx*256 + blocks_[block_idx].select0(remainder);
        for (uint32_t j = 0; j <= blocks_[i+1].rank1_ - blocks_[i].rank1_; j++) {
          spill0_[spilled++] = pos;
          pos = next0(pos + 1);
          assert(spilled <= spill0_size_);
        }
      }
      assert(spilled == spill0_size_);
    }
  };
  using bitvec_t = BitVector;

 public:
  LoudsCC() = default;

  LoudsCC(uint32_t size) : bv_(size == 0 ? 0 : (size*2 - 1)) {}

  ~LoudsCC() = default;

  void add_node(uint32_t degree) {
    for (int i = 0; i < degree; i++) {
      bv_.append1();
    }
    bv_.append0();
  }

  void build() {
    bv_.build();
  }

  auto node_id(uint32_t pos) const -> uint32_t {
    assert(pos < bv_.size());
    return bv_.rank0(pos);
  }

  auto child_id(uint32_t pos) const -> uint32_t {
    assert(pos < bv_.size());
    return bv_.rank1(pos + 1);
  }

  auto leaf_id(uint32_t pos) const -> uint32_t {
    assert(pos < bv_.size());
    return bv_.rank00(pos);
  }

  auto internal_id(uint32_t pos) const -> uint32_t {
    assert(pos < bv_.size());
    return bv_.rank0(pos) - bv_.rank00(pos);
  }

  auto has_parent(uint32_t pos) const -> uint32_t {
    assert(pos < bv_.size());
    return bv_.rank0(pos) > 0;
  }

  auto num_nodes() const -> uint32_t {
    return bv_.rank0();
  }

  auto num_leaves() const -> uint32_t {
    return bv_.rank00();
  }

  auto num_internals() const -> uint32_t {
    return num_nodes() - num_leaves();
  }

  auto is_leaf(uint32_t pos) const -> bool {
    assert(pos == 0 || bv_.get(pos - 1) == 0);
    assert(pos < bv_.size());
    return bv_.get(pos) == 0;
  }

  // return # of children
  auto node_degree(uint32_t pos) const -> uint32_t {
    assert(pos == 0 || bv_.get(pos - 1) == 0);
    uint32_t next0 = bv_.next0(pos);
    return next0 - pos;
  }

  auto child_pos(uint32_t pos) const -> uint32_t {
    assert(bv_.get(pos));
    auto ret = bv_.r1s0(pos + 1) + 1;
    return ret;
  }

  void shrink_to_fit() {
    bv_.shrink_to_fit();
  }

  void reserve(uint32_t num_nodes) {
    bv_.reserve(num_nodes*2 - 1);
  }

  auto size() const -> uint32_t {
    return bv_.size();
  }

  auto size_in_bytes() const -> uint32_t {
    return bv_.size_in_bytes();
  }

  auto size_in_bits() const -> uint32_t {
    return size_in_bytes() * 8;
  }

  void clear() {
    bv_.clear();
  }

#ifdef __COMPARE_COCO__
  void to_louds_sux(std::unique_ptr<LoudsSux<>> &out) {
    out = std::make_unique<LoudsSux<>>(num_nodes());
    uint32_t pos = 0;
    while (pos < bv_.size()) {
      uint32_t deg = node_degree(pos);
      out->add_node(deg);
      pos += deg + 1;
    }
    out->build();
  }

  auto get(uint32_t pos) const -> bool {
    return bv_.get(pos);
  }
#endif

 private:
  bitvec_t bv_;
};

}  // namespace coco_paper
