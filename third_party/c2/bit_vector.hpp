#pragma once

#include "utils.hpp"


namespace coco_paper {

struct BitVector {
  struct Block {
    uint32_t rank_{0};
    uint8_t subrank_[4]{0};
    uint64_t bits_[4]{0};

    auto build_index() -> uint32_t {
      uint32_t rank = 0;
      for (int i = 0; i < 4; i++) {
        subrank_[i] = rank;
        rank += __builtin_popcountll(bits_[i]);
      }
      return rank;
    }

    auto get(uint32_t pos) const -> bool {
      assert(pos < 256);
      return GET_BIT(bits_[pos/64], pos % 64);
    }

    auto rank1(uint32_t size) const -> uint32_t {
      return subrank_[size/64] + __builtin_popcountll(bits_[size/64] & MASK(size%64));
    }

    auto rank0(uint32_t size) const -> uint32_t {
      return size - rank1(size);
    }

    auto select1(uint32_t rank) const -> uint32_t {
      assert(rank <= rank1());
      if (rank <= subrank_[2]) {
        if (rank <= subrank_[1]) {
          return 64*0 + selectll(bits_[0], rank);
        } else {
          return 64*1 + selectll(bits_[1], rank - subrank_[1]);
        }
      } else {
        if (rank <= subrank_[3]) {
          return 64*2 + selectll(bits_[2], rank - subrank_[2]);
        } else {
          return 64*3 + selectll(bits_[3], rank - subrank_[3]);
        }
      }
    }

    auto select0(uint32_t rank) const -> uint32_t {
      assert(rank <= rank0());
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

    auto rank1() const -> uint32_t {
      return subrank_[3] + __builtin_popcountll(bits_[3]);
    }

    auto rank0() const -> uint32_t {
      return 256 - rank1();
    }
  };  // 40 bytes

  struct Select1Index {
    static constexpr uint32_t sample_rate_ = 256;
    static constexpr uint32_t spill_threshold_ = 128;

    Select1Index() = default;

    Select1Index(const BitVector &bv) : sample_{nullptr}, spill_{nullptr} {
      build(bv);
    }

    ~Select1Index() {
      free(sample_);
      free(spill_);
    }

    void build(const BitVector &bv) {
      bv_ = &bv;
      if (bv_ == nullptr) {
        clear();
      }
      init_sample();
      init_spill();
    }

    void clear() {
      free(sample_);
      free(spill_);
      bv_ = nullptr;
      sample_ = nullptr;
      spill_ = nullptr;
      spill_size_ = 0;
    }

    auto select1(uint32_t rank) const -> uint32_t {
      assert(bv_ != nullptr);
      assert(sample_ != nullptr);
      assert(rank <= bv_->rank_);

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
        if (bv_->blocks_[mid].rank_ < rank) {
          left = mid;
        } else {
          right = mid - 1;
        }
      }
      while (bv_->blocks_[left+1].rank_ < rank) {
        left++;
      }
      uint32_t remainder = rank - bv_->blocks_[left].rank_;
      uint32_t ret = left*256 + bv_->blocks_[left].select1(remainder);
      return ret;
    }

   private:
    void init_sample() {
      uint32_t num_blocks = bv_->capacity_ / 256;
      uint32_t num_samples = (bv_->rank_ + sample_rate_ - 1) / sample_rate_;

      sample_ = reinterpret_cast<uint32_t *>(realloc(sample_, (num_samples + 1) * sizeof(uint32_t)));
      sample_[0] = 0;
      for (uint32_t i = 1; i < num_samples; i++) {
        uint32_t j = sample_[i-1];
        assert(bv_->blocks_[j].rank_ < i*sample_rate_);
        while (bv_->blocks_[j+1].rank_ < i*sample_rate_) {
          j++;
        }
        sample_[i] = j;
      }
      sample_[num_samples] = num_blocks;
    }

    void init_spill() {
      int num_samples = (bv_->rank_ + sample_rate_ - 1) / sample_rate_;
      int num_spills = 0;

      for (int i = 0; i < num_samples; i++) {
        uint32_t dist = sample_[i+1] - sample_[i];
        if (dist >= spill_threshold_) {
          sample_[i] |= BIT(31);
          num_spills += std::min(static_cast<uint32_t>(sample_rate_), bv_->rank_ - i*sample_rate_) + 1;
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

        uint32_t remainder = i*sample_rate_ - bv_->blocks_[block_idx].rank_;
        uint32_t pos = block_idx*256 + bv_->blocks_[block_idx].select1(remainder);
        uint32_t end = std::min(static_cast<uint32_t>(sample_rate_), bv_->rank_ - i*sample_rate_);
        for (uint32_t j = 0; j <= end; j++) {
          spill_[spilled++] = pos;
          pos = bv_->next1(pos + 1);
          assert(spilled <= spill_size_);
        }
      }
      assert(spilled == spill_size_);
    }

    const BitVector *bv_{nullptr};

    uint32_t *sample_{nullptr};

    uint32_t *spill_{nullptr};

    uint32_t spill_size_{0};
  };

  struct Select0Index {
    static constexpr uint32_t sample_rate_ = 256;
    static constexpr uint32_t spill_threshold_ = 128;

    Select0Index() = default;

    Select0Index(const BitVector &bv) : sample_{nullptr}, spill_{nullptr} {
      build(bv);
    }

    ~Select0Index() {
      free(sample_);
      free(spill_);
    }

    void build(const BitVector &bv) {
      bv_ = &bv;
      if (bv_ == nullptr) {
        clear();
      }
      init_sample();
      init_spill();
    }

    void clear() {
      free(sample_);
      free(spill_);
      bv_ = nullptr;
      sample_ = nullptr;
      spill_ = nullptr;
      spill_size_ = 0;
    }

    auto select0(uint32_t rank) const -> uint32_t {
      assert(bv_ != nullptr);
      assert(sample_ != nullptr);
      assert(rank <= bv_->size_ - bv_->rank_);

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
        if (mid*256 - bv_->blocks_[mid].rank_ < rank) {
          left = mid;
        } else {
          right = mid - 1;
        }
      }
      while ((left+1)*256 - bv_->blocks_[left+1].rank_ < rank) {
        left++;
      }
      uint32_t remainder = rank - (left*256 - bv_->blocks_[left].rank_);
      uint32_t ret = left*256 + bv_->blocks_[left].select0(remainder);
      return ret;
    }

   private:
    void init_sample() {
      uint32_t num_blocks = bv_->capacity_ / 256;
      uint32_t num_samples = (bv_->size_ - bv_->rank_ + sample_rate_ - 1) / sample_rate_;

      sample_ = reinterpret_cast<uint32_t *>(realloc(sample_, (num_samples + 1) * sizeof(uint32_t)));
      sample_[0] = 0;
      for (uint32_t i = 1; i < num_samples; i++) {
        uint32_t j = sample_[i-1];
        assert(j*256 - bv_->blocks_[j].rank_ < i*sample_rate_);
        while ((j+1)*256 - bv_->blocks_[j+1].rank_ < i*sample_rate_) {
          j++;
        }
        sample_[i] = j;
      }
      sample_[num_samples] = num_blocks;
    }

    void init_spill() {
      int num_samples = (bv_->size_ - bv_->rank_ + sample_rate_ - 1) / sample_rate_;
      int num_spills = 0;

      for (int i = 0; i < num_samples; i++) {
        uint32_t dist = sample_[i+1] - sample_[i];
        if (dist >= spill_threshold_) {
          sample_[i] |= BIT(31);
          num_spills += std::min(static_cast<uint32_t>(sample_rate_), bv_->size_ - bv_->rank_ - i*sample_rate_) + 1;
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

        uint32_t remainder = i*sample_rate_ - (block_idx*256 - bv_->blocks_[block_idx].rank_);
        uint32_t pos = block_idx*256 + bv_->blocks_[block_idx].select0(remainder);
        uint32_t end = std::min(static_cast<uint32_t>(sample_rate_), bv_->size_ - bv_->rank_ - i*sample_rate_);
        for (uint32_t j = 0; j <= end; j++) {
          spill_[spilled++] = pos;
          pos = bv_->next0(pos + 1);
          assert(spilled <= spill_size_);
        }
      }
      assert(spilled == spill_size_);
    }

    const BitVector *bv_{nullptr};

    uint32_t *sample_{nullptr};

    uint32_t *spill_{nullptr};

    uint32_t spill_size_{0};
  };

  BitVector() = default;

  ~BitVector() {
    clear();
  }

  auto size() const -> uint32_t {
    return size_;
  }

  auto rank1() const -> uint32_t {
    return rank_;
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

  void load_bits(const uint64_t *bits, size_t start, uint32_t size) {
    if (size_ + size > capacity_) {
      capacity_ = std::max<uint32_t>(((size_ + size) * 2 + 255) / 256 * 256, 256 * 8);
      blocks_ = reinterpret_cast<Block *>(realloc(blocks_, sizeof(Block) * (capacity_ / 256 + 1)));
    }

    uint32_t leftover = (size_ + 255) / 256 * 256 - size_;
    if (size <= leftover) {
      copy_bits(blocks_[size_/256].bits_, size_ % 256, bits, start, size);
      size_ += size;
      return;
    }

    size_t end = start + size;
    copy_bits(blocks_[size_/256].bits_, size_ % 256, bits, start, leftover);
    size_ += leftover;
    start += leftover;
    while (start + 256 <= end) {
      copy_bits(blocks_[size_/256].bits_, size_ % 256, bits, start, 256);
      size_ += 256;
      start += 256;
    }
    copy_bits(blocks_[size_/256].bits_, size_ % 256, bits, start, end - start);
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

  void build() {
    if (blocks_ == nullptr) {
      return;
    }

    shrink_to_fit();

    // clear trailing bits
    uint32_t remainder = size_ % 256;
    for (uint32_t i = (remainder + 63) / 64; i < 4; i++) {
      blocks_[size_ / 256].bits_[i] = 0;
    }
    for (uint32_t i = 0; i < 4; i++) {
      blocks_[capacity_ / 256].subrank_[i] = 0;
      blocks_[capacity_ / 256].bits_[i] = 0;
    }

    // build rank index
    rank_ = 0;
    for (uint32_t i = 0; i < capacity_ / 256; i++) {
      blocks_[i].rank_ = rank_;
      rank_ += blocks_[i].build_index();
    }
    blocks_[capacity_ / 256].rank_ = rank_;
  }

  // get the `pos`-th bit of bv<bvnum>
  auto get(uint32_t pos) const -> bool {
    assert(pos < size_);
    return blocks_[pos/256].get(pos % 256);
  }

  auto rank1(uint32_t size) const -> uint32_t {
    assert(size <= size_);
    const auto &block = blocks_[size/256];
    return block.rank_ + block.rank1(size % 256);
  }

  auto rank0(uint32_t size) const -> uint32_t {
    assert(size <= size_);
    return size - rank1(size);
  }

  auto next1(uint32_t pos) const -> uint32_t {
    if (pos >= size_) {
      return size_;
    }
    uint32_t block_idx = pos / 64;
    uint32_t remainder = pos % 64;
    uint32_t dist = 0;
    const uint64_t &block = blocks_[block_idx / 4].bits_[block_idx % 4];
    if (block >> remainder) {
      return pos + __builtin_ctzll(block >> remainder);
    }
    dist += 64 - remainder;
    while (pos + dist < size_) {
      block_idx = (pos + dist) / 64;
      const uint64_t &block = blocks_[block_idx / 4].bits_[block_idx % 4];
      if (block) {
        return pos + dist + __builtin_ctzll(block);
      }
      dist += 64;
    }
    return size_;
  }

  auto next0(uint32_t pos) const -> uint32_t {
    uint32_t idx = pos / 64;
    uint64_t elt = ~blocks_[idx/4].bits_[idx%4] & ~MASK(pos%64);
    while (elt == 0) {  // terminates at sentinel block since it's all 0
      idx++;
      elt = ~blocks_[idx/4].bits_[idx%4];
    }
    uint32_t ret = idx*64 + __builtin_ctzll(elt);
    return ret;
  }

 private:
  Block *blocks_{nullptr};

  uint32_t capacity_{0};

  uint32_t size_{0};

  uint32_t rank_{0};

  friend struct Select1Index;
  friend struct Select0Index;
 public:
  void serialize(std::ostream& out) const {
    archive::bits(out,size(),[&](uint32_t i){return get(i);});
  }
  void load(std::istream& in) {
    archive::bits(in,[&](bool bit){if(bit) append1();else append0();});
    build();
  }

};

}  // namespace coco_paper
