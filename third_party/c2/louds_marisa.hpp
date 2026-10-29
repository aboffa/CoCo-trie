// Taken and modified from https://github.com/s-yata/marisa-trie
// Used for ablation study
#pragma once

#include "utils.hpp"
#include "../baseline_marisa/marisa/lib/marisa/grimoire/vector.h"


namespace coco_paper {

class LoudsMarisa {
 public:
  using bitvec_t = marisa::grimoire::vector::BitVector;

  LoudsMarisa() = default;

  void add_node(uint32_t degree) {
    for (uint32_t i = 0; i < degree; i++) {
      louds_.push_back(true);
    }
    louds_.push_back(false);
  }

  void add_bit(bool term, bool link) {
    terminal_flags_.push_back(term);
    link_flags_.push_back(link);
  }

  void build() {
    louds_.build(true, true);
    terminal_flags_.build(false, false);
    link_flags_.build(false, false);

    uint32_t ones = 0, zeros = 0;
    for (uint32_t i = 0; i < size(); i++) {
      ones += louds(i);
      zeros += !louds(i);
    }
    printf("ones: %d, zeros: %d\n", ones, zeros);
  }

  auto size() const -> uint32_t {
    return louds_.size();
  }

  auto louds(uint32_t pos) const -> bool {
    return louds_[pos];
  }

  auto is_term(uint32_t pos) const -> bool {
    return terminal_flags_[pos];
  }

  auto is_link(uint32_t pos) const -> bool {
    return link_flags_[pos];
  }

  auto leaf_id(uint32_t pos) const -> uint32_t {
    return terminal_flags_.rank1(pos);
  }

  auto link_id(uint32_t pos) const -> uint32_t {
    return link_flags_.rank1(pos);
  }

  auto child_pos(uint32_t pos) const -> uint32_t {
    return louds_.select0(louds_.rank1(pos)) + 1;
  }

  auto parent_pos(uint32_t pos) const -> uint32_t {
    return louds_.select1(louds_.rank0(pos));
  }

  bitvec_t louds_;
  bitvec_t terminal_flags_;
  bitvec_t link_flags_;
};

}  // namespace coco_paper
