// Taken and modified from https://github.com/aboffa/CoCo-trie;
// Used for ablation study

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <sux/bits/SimpleSelectZero.hpp>
#include <sux/bits/SimpleSelectZeroHalf.hpp>

#pragma once


namespace coco_paper {

template<typename rank_support1 = sdsl::rank_support_v<1>,
        typename rank_support00 = sdsl::rank_support_v<0, 2>,
        typename select0_type = sux::bits::SimpleSelectZero<>>
struct LoudsSux {
  sdsl::bit_vector bv;
  // pointer used while building
  size_t build_p = 0;
  rank_support1 rs1;
  rank_support00 rs00;
  select0_type ss0;

  LoudsSux() = default;

  LoudsSux(size_t n) {
    // n is the number of nodes
    bv = sdsl::bit_vector(2 * n - 1, 0);
    build_p = 0;
  }

  ~LoudsSux() {}

  LoudsSux(sdsl::bit_vector &_bv) {
    size_t _bv_size = _bv.size();
    bv.swap(_bv);
    assert(bv.size() == _bv_size);
    rs1 = rank_support1(&bv);
    rs00 = rank_support00(&bv);
    if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZero<>>)
      ss0 = select0_type(bv.data(), bv.size(), 2);
    else if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZeroHalf<>>)
      ss0 = select0_type(bv.data(), bv.size());
    else
      ss0 = select0_type(&bv);
  }

  void reserve(uint32_t n) {
    bv = sdsl::bit_vector(n, 0);
    build_p = 0;
  }

  void add_node(size_t num_child) {
    for (auto i = 0; i < num_child; ++i) {
      bv[build_p++] = 1;
    }
    bv[build_p++] = 0;
  }

  void build() {
    // init rank select support
    rs1 = rank_support1(&bv);
    rs00 = rank_support00(&bv);
    if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZero<>>)
      ss0 = select0_type(bv.data(), bv.size(), 2);
    else if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZeroHalf<>>)
      ss0 = select0_type(bv.data(), bv.size());
    else
      ss0 = select0_type(&bv);
    assert(bv.size() == this->build_p);
  }

  // from index of nodes to position in the bv
  size_t node_select(size_t i) const {
    if constexpr (std::is_same_v<select0_type, sdsl::select_support_mcl<0>>)
      return const_cast<LoudsSux *>(this)->ss0.select(i) + 1;
    else
      return const_cast<LoudsSux *>(this)->ss0.selectZero(i - 1) + 1;
  }

  size_t size_in_bytes() const {
    size_t to_return = 0;
    if constexpr (std::is_same_v<select0_type, sdsl::select_support_mcl<0>>)
      to_return += sdsl::size_in_bytes(ss0);
    else
      to_return += ss0.bitCount() / CHAR_BIT;
    to_return += sdsl::size_in_bytes(rs1);
    to_return += sdsl::size_in_bytes(rs00);
    to_return += sdsl::size_in_bytes(bv);
    to_return += sizeof(LoudsSux);
    return to_return;
  }

  inline size_t size_in_bits() const {
    size_t to_return = size_in_bytes();
    return to_return * CHAR_BIT;
  }

  size_t nextZero(const uint64_t bv_idx) const {
    uint32_t idx = bv_idx / 64;
    uint64_t elt = ~bv.data()[idx] & ~MASK(bv_idx%64);
    while (elt == 0) {
      idx++;
      elt = ~bv.data()[idx];
    }
    uint32_t ret = idx*64 + __builtin_ctzll(elt);
    return ret;
  }

  // following are the performance critical functions of cache optimized LOUDS
  auto node_id(uint32_t pos) const -> uint32_t {
    return pos - rs1.rank(pos);
  }

  auto leaf_id(uint32_t pos) const -> uint32_t {
    return rs00.rank(pos);
  }

  auto internal_id(uint32_t pos) const -> uint32_t {
    return node_id(pos) - leaf_id(pos);
  }

  auto is_leaf(uint32_t pos) const -> bool {
    return bv[pos] == 0;
  }

  // return # of children
  auto node_degree(uint32_t pos) const -> uint32_t {
    uint32_t next0 = nextZero(pos);
    return next0 - pos;
  }

  auto child_pos(uint32_t pos) const -> uint32_t {
    uint32_t rank = rs1.rank(pos + 1);
    uint32_t ret = node_select(rank);
    return ret;
  }

  auto get(uint32_t pos) const -> bool {
    return bv[pos];
  }

  auto has_parent(uint32_t pos) const -> uint32_t {
    return pos - rs1.rank(pos) > 0;
  }

  auto num_nodes() const -> uint32_t {
    return (bv.size() + 1) / 2;
  }

  auto num_leaves() const -> uint32_t {
    return rs00.rank(bv.size());
  }

  auto num_internals() const -> uint32_t {
    return num_nodes() - num_leaves();
  }

  auto size() const -> uint32_t {
    return build_p;
  }
 public:
  void serialize(std::ostream& out) const { bv.serialize(out); }
  void load(std::istream& in) { bv.load(in); if(!in) throw std::runtime_error("Invalid LOUDS"); build_p=bv.size();build(); }

};

}  // namespace coco_paper
