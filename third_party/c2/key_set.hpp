#pragma once

#include "utils.hpp"



namespace coco_paper {

template <typename Key>
struct KeySet {
  using key_type = Key;

  struct Fragment {
    const key_type *key_{nullptr};
    uint32_t offset_{0};
    uint32_t length_{0};
    uint32_t id_{0};

    Fragment() = default;

    ~Fragment() = default;

    Fragment(const Fragment &other) = default;

    Fragment(Fragment &&other) = default;

    auto operator=(const Fragment &other) -> Fragment & = default;

    auto operator=(Fragment &&other) -> Fragment & = default;

    Fragment(const key_type *key, uint32_t offset, uint32_t length, uint32_t id)
        : key_(key), offset_(offset), length_(length), id_(id) {}

    Fragment(const key_type *key, uint32_t id)
        : key_(key), offset_(0), length_(key->size()), id_(id) {}

    auto get_label(uint32_t idx, bool reverse = false) const -> uint8_t {
      assert(idx < length_);
      return reverse ? (*key_)[offset_ + length_ - idx - 1] : (*key_)[offset_ + idx];
    }

    auto size() const -> uint32_t {
      return length_;
    }

    auto lt(const Fragment &rhs, bool reverse) const -> bool {
      int len = std::min(length_, rhs.length_);
      for (int i = 0; i < len; i++) {
        if (get_label(i, reverse) < rhs.get_label(i, reverse)) {
          return true;
        } else if (get_label(i, reverse) > rhs.get_label(i, reverse)) {
          return false;
        }
      }
      return length_ < rhs.length_;
    }

    auto operator==(const Fragment &rhs) const -> bool {
      if (length_ != rhs.length_) {
        return false;
      }
      return !memcmp(key_->c_str() + offset_, rhs.key_->c_str() + rhs.offset_, length_);
    }

    auto materialize(bool reverse = false) const -> key_type {
      key_type ret = key_->substr(offset_, length_);
      if (reverse) {
        std::reverse(ret.begin(), ret.end());
      }
      return ret;
    }

    template <typename T>
    void append_to(T &t, bool terminator = true, bool reverse = false) const {
      if (!reverse) {
        auto begin = key_->begin() + offset_, end = key_->begin() + (offset_ + length_);
        t.insert(t.end(), begin, end);
      } else {
        auto begin = key_->rbegin() + (key_->size() - offset_ - length_), end = key_->rbegin() + (key_->size() - offset_);
        t.insert(t.end(), begin, end);
      }
      if (terminator) {
        t.push_back(terminator_);
      }
    }

    auto substr_range(uint32_t idx, uint32_t length, bool reverse) const -> std::pair<uint32_t, uint32_t> {
      assert(idx + length <= length_);
      if (reverse) {
        return std::make_pair(offset_ + length_ - idx - length, length);
      } else {
        return std::make_pair(offset_ + idx, length);
      }
    }

    auto c_str() const -> const uint8_t * {
      return key_->c_str() + offset_;
    }
  };

  std::vector<Fragment> fragments_;
  size_t space_cost_{0};
  bool reverse_{false};

  KeySet() = default;

  KeySet(const KeySet &other) = default;

  KeySet(KeySet &&other) = default;

  auto operator=(const KeySet &other) -> KeySet & = default;

  auto operator=(KeySet &&other) -> KeySet & = default;

  KeySet(bool reverse) : reverse_(reverse) {}

  void emplace_back(const key_type *key) {
    fragments_.emplace_back(key, fragments_.size());
    space_cost_ += key->size();
  }

  void emplace_back(const key_type *key, uint32_t offset, uint32_t length) {
    assert(offset + length <= key->size());
    fragments_.emplace_back(key, offset, length, fragments_.size());
    space_cost_ += length;
  }

  void push_back(const Fragment &frag, bool reset_id = false) {
    fragments_.push_back(frag);
    space_cost_ += frag.length_;
    if (reset_id) {
      fragments_.back().id_ = fragments_.size() - 1;
    }
  }

  auto operator[](int idx) const -> const Fragment & {
    return get(idx);
  }

  auto get(int idx) const -> const Fragment & {
    assert(idx >=0 && idx < size());
    return fragments_[idx];
  }

  auto get_label(int idx, int depth) const -> uint8_t {
    assert(idx >= 0 && idx < size());
    return fragments_[idx].get_label(depth, reverse_);
  }

  auto materialize(int idx) const -> key_type {
    assert(idx >= 0 && idx < size());
    return fragments_[idx].materialize(reverse_);
  }

  auto substr_range(uint32_t idx, uint32_t offset, uint32_t length) const -> std::pair<uint32_t, uint32_t> {
    assert(idx >= 0 && idx < size());
    return fragments_[idx].substr_range(offset, length, reverse_);
  }

  auto front() const -> const Fragment & {
    return get(0);
  }

  auto back() const -> const Fragment & {
    return get(size() - 1);
  }

  auto empty() const -> bool {
    return fragments_.empty();
  }

  auto size() const -> size_t {
    return fragments_.size();
  }

  auto space_cost() const -> size_t {
    return space_cost_;
  }

  void sort() {
    std::sort(fragments_.begin(), fragments_.end(),
              [&](const Fragment &f1, const Fragment &f2) -> bool {
                return f1.lt(f2, reverse_);
              });
  }

  // must be sorted; return size after deduplication
  auto deduped_size() const -> size_t {
    size_t ret = 0;
    const Fragment *next = nullptr;
    for (size_t i = fragments_.size(); i > 0; i--) {
      const KeySet<key_type>::Fragment &cur = fragments_[i - 1];
      if (next == nullptr || cur != *next) {
        ret += cur.size() + 1;  // key + terminator
      }
    }
    return ret;
  }

  // must be sorted; returns {total size of lcp, size after sorting and deduplication}
  auto lcp_size() const -> std::pair<size_t, size_t> {
    size_t total_lcp = 0, sorted_size = 0;
    const Fragment *next = nullptr;
    for (size_t i = fragments_.size(); i > 0; i--) {
      const KeySet<key_type>::Fragment &cur = fragments_[i - 1];
      if (next != nullptr) {
        uint32_t len = std::min(cur.size(), next->size());
        uint32_t match = 0;
        while (match < len && cur.get_label(match, reverse_) == next->get_label(match, reverse_)) {
          match++;
        }
        total_lcp += match;
        if (match == cur.size()) {
          continue;
        }
      }
      sorted_size += cur.size() + 1;  // key + terminator
      next = &cur;
    }
    return std::make_pair(total_lcp, sorted_size);
  }

  void reverse() {
    reverse_ = !reverse_;
  }

  void set_reverse(bool rev) {
    reverse_ = rev;
  }

  auto make_sample(size_t sample_size, bool sort_rev = false) const -> KeySet<key_type> {
    if (sample_size*2 >= space_cost()) {
      return KeySet<key_type>(*this);
    }
    KeySet<key_type> sample;
    size_t sampled_len = (size() + 63) / 64 + 1;
    auto sampled = new uint64_t[sampled_len];
    memset(sampled, 0, sampled_len*sizeof(uint64_t));
    size_t cur_size = 0;
    while (cur_size < sample_size) {
      size_t i = std::rand() % size();
      if (!GET_BIT(sampled[i/64], i%64)) {
        SET_BIT(sampled[i/64], i%64);
        sample.push_back(fragments_[i], true);
        cur_size += fragments_[i].size();
      }
    }
    if (sort_rev) {
      sample.reverse();
      sample.sort();
    }
    return sample;
  }
};

}
