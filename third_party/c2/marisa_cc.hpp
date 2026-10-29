#pragma once

#include "utils.hpp"
#include "static_vector.hpp"
#include "ls4patricia.hpp"
#include "key_set.hpp"
#include <sdsl/int_vector.hpp>

#include <algorithm>
#include <type_traits>
#include <vector>
#include <queue>

#include "strpool.hpp"


// #define __NO_BRANCHING_LABEL__
#define __ENABLE_CACHE__


namespace coco_paper {

template <typename Key, bool reverse = false>
class MarisaCC : public StringPool<Key> {
 public:
  using key_type = Key;
  using label_vec = StaticVector<uint8_t>;
  using topo_t = LS4Patricia;
  using strpool_t = StringPool<key_type>;

  static constexpr bool reverse_ = reverse;
  static constexpr uint32_t link_cutoff_ = 3;  // unary paths must be at least this long to be considered for recursive compression
  static_assert(link_cutoff_ >= 2);

#ifdef __ENABLE_CACHE__
  // static constexpr uint32_t cache_ratio_ = 128;  // MARISA_HUGE_CACHE
  // static constexpr uint32_t cache_ratio_ = 256;  // MARISA_LARGE_CACHE
  static constexpr uint32_t cache_ratio_ = 512;  // MARISA_DEFAULT_CACHE — matches baseline MarisaWrapper
#endif

  void print_space_cost_breakdown() const {
    size_t topo = 0, link = 0, data = 0;
    space_cost_breakdown(topo, link, data);
    printf("topology: %lf MB, link: %lf MB, data: %lf MB\n", (double)topo/mb_bits, (double)link/mb_bits, (double)data/mb_bits);
  }

  struct UnaryPathStats {
    uint64_t num_single{0};       // path_len == 1: plain branching label, no unary extension
    uint64_t num_short{0};        // 1 < path_len < link_cutoff_: unary path too short to compress
    uint64_t total_short_len{0};
    uint64_t num_links{0};        // path_len >= link_cutoff_: compressed to string pool
    uint64_t total_link_len{0};
    uint64_t max_link_len{0};
    std::vector<uint64_t> len_hist;  // len_hist[i] = # links with length (link_cutoff_ + i)
  };

  void print_unary_path_stats() const {
    if constexpr (reverse_) return;
    uint64_t total = stats_.num_single + stats_.num_short + stats_.num_links;
    if (total == 0) return;
    printf("--- Unary Path Analysis (Marisa, cutoff=%u) ---\n", link_cutoff_);
    printf("  branches: %llu  |  single-char: %llu (%.1f%%)  short-unary: %llu (%.1f%%)  links: %llu (%.1f%%)\n",
           total,
           stats_.num_single, 100.0 * stats_.num_single / total,
           stats_.num_short,  100.0 * stats_.num_short  / total,
           stats_.num_links,  100.0 * stats_.num_links  / total);
    if (stats_.num_links > 0) {
      printf("  link len: avg=%.2f  max=%llu  total_chars=%llu\n",
             (double)stats_.total_link_len / stats_.num_links,
             stats_.max_link_len, stats_.total_link_len);
      printf("  link length distribution (len:count):");
      for (size_t i = 0; i < stats_.len_hist.size(); i++) {
        if (stats_.len_hist[i] > 0)
          printf("  %u:%llu", (uint32_t)(link_cutoff_ + i), stats_.len_hist[i]);
      }
      printf("\n");
    }
    if (stats_.num_short > 0) {
      printf("  short unary avg len: %.2f\n", (double)stats_.total_short_len / stats_.num_short);
    }
  }

 private:
  struct Range {
    uint32_t begin_{0};
    uint32_t end_{0};
    uint32_t depth_{0};

    Range() = default;

    Range(uint32_t begin, uint32_t end, uint32_t depth) : begin_(begin), end_(end), depth_(depth) {}
  };

#ifdef __ENABLE_CACHE__
  struct Cache {
    uint32_t parent_{0};
    uint32_t child_{0};
    union {
      uint32_t link_;
      float weight_;
    } union_;

    Cache() { set_weight(-1); }

    auto link() const -> uint32_t {
      return union_.link_;
    }

    auto weight() const -> float {
      return union_.weight_;
    }

    void set_link(uint32_t link) {
      union_.link_ = link;
    }

    void set_weight(float weight) {
      union_.weight_ = weight;
    }
  };
#endif

 public:
  MarisaCC() = default;

  ~MarisaCC() {
    delete next_;
  }

  template <typename Iterator, bool rev = reverse, typename = std::enable_if_t<!rev>>
  void build(Iterator begin, Iterator end, bool sorted = false, int max_recursion = 0, int mask = 0) {
    KeySet<key_type> key_set;
    while (begin != end) {
      key_set.emplace_back(&(*begin));
      ++begin;
    }
    assert(!key_set.empty());
    if (!sorted) {
      key_set.sort();
    }
    build(key_set, nullptr, max_recursion, mask);
  }

  auto size() const -> uint32_t override {
    if constexpr (!reverse_) {
      return topo_.num_leaves();
    } else {
      return links_.size();
    }
  }

  auto size_in_bytes() const -> size_t override {
    size_t ret = topo_.size_in_bytes() + labels_.size_in_bytes() + sdsl::size_in_bytes(links_);
    if (next_ != nullptr) {
      ret += next_->size_in_bytes();
    }
  #ifdef __ENABLE_CACHE__
    ret += cache_.size() * sizeof(Cache);
  #endif
    return ret;
  }

  auto size_in_bits() const -> size_t override {
    return size_in_bytes() * 8;
  }

  auto trie_size_in_bits() const -> size_t {
    size_t ret;
    ret = topo_.size_in_bits() + labels_.size_in_bytes() * 8 + sdsl::size_in_bytes(links_) * 8;
  #ifdef __ENABLE_CACHE__
    ret += cache_.size() * sizeof(Cache);
  #endif
    return ret;
  }

  void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const override {
    topo += topo_.size_in_bits();
  #ifdef __ENABLE_CACHE__
    topo += cache_.size() * sizeof(Cache);
  #endif
    link += sdsl::size_in_bytes(links_) * 8;
    data += labels_.size_in_bytes() * 8;
    next_->space_cost_breakdown(topo, link, data);
  }

#ifdef __NO_BRANCHING_LABEL__
  // returns leaf ID (-1 if not found)
  template <bool rev = reverse, typename = std::enable_if_t<!rev>>
  auto lookup(const key_type &key) const -> uint32_t {
    uint32_t len = key.size(), matched_len = 0;
    uint32_t pos = 0;

    while (matched_len < len) {
    #ifdef __ENABLE_CACHE__
      auto ret = search_cache(pos, key, matched_len);
      if (ret == -1) {  // mismatch
        return -1;
      } else if (ret == 1) {  // branch terminates
        return matched_len == len ? topo_.leaf_id(pos) : -1;
      } else if (ret == 2) {
        continue;
      }  // else ret == 0, i.e. not cached
    #endif
      bool found = false;
      do {  // search for label
        // _mm_prefetch(&labels_[pos], _MM_HINT_T0);
        if (topo_.is_link(pos)) {
          uint32_t link_len = next_->match(key, matched_len, topo_.link_id(pos), labels_[pos]);
          if (link_len != -1) {
            matched_len += link_len;
            found = true;
            break;
          }
        } else if (labels_[pos] == key[matched_len]) {
          matched_len++;
          found = true;
          break;
        }
        pos++;
      } while (!topo_.louds(pos));

      if (!found) {
        return -1;
      }
      if (!topo_.has_child(pos)) {  // branch terminates
        return matched_len == len ? topo_.leaf_id(pos) : -1;
      }
      pos = topo_.child_pos(pos);
    }

    if (labels_[pos] == terminator_) {  // prefix key
      return topo_.leaf_id(pos);
    }  // else early termination
    return -1;
  }
#else
  // returns leaf ID (-1 if not found)
  template <bool rev = reverse, typename = std::enable_if_t<!rev>>
  auto lookup(const key_type &key) const -> uint32_t {
    uint32_t len = key.size(), matched_len = 0;
    uint32_t pos = 0;

    while (matched_len < len) {
    #ifdef __ENABLE_CACHE__
      auto ret = search_cache(pos, key, matched_len);
      if (ret == -1) {  // mismatch
        return -1;
      } else if (ret == 1) {  // branch terminates
        return matched_len == len ? topo_.leaf_id(pos) : -1;
      } else if (ret == 2) {
        continue;
      }  // else ret == 0, i.e. not cached
    #endif
      uint32_t end = topo_.node_end(pos);
      pos = labels_.find(key[matched_len], pos, end);
      if (pos == end) {  // mismatch
        return -1;
      }
      assert(labels_[pos] == key[matched_len]);
      matched_len++;

      if (topo_.is_link(pos)) {
        uint32_t link_len = next_->match(key, matched_len, topo_.link_id(pos));
        if (link_len == -1) {
          return -1;
        }
        matched_len += link_len;
      }
      if (!topo_.has_child(pos)) {  // branch terminates
        return matched_len == len ? topo_.leaf_id(pos) : -1;
      }
      pos = topo_.child_pos(pos);
    }

    if (labels_[pos] == terminator_) {  // prefix key
      return topo_.leaf_id(pos);
    }  // else early termination
    return -1;
  }
#endif

  // Returns true if any stored key starts with q.
  template <bool rev = reverse, typename = std::enable_if_t<!rev>>
  auto contains_prefix(const key_type &q) const -> bool {
    uint32_t len = q.size(), matched_len = 0;
    uint32_t pos = 0;

    while (matched_len < len) {
      uint32_t end = topo_.node_end(pos);
      pos = labels_.find(q[matched_len], pos, end);
      if (pos == end) return false;
      matched_len++;

      if (topo_.is_link(pos)) {
        uint32_t lr = topo_.link_id(pos);
        uint32_t link_len = next_->match(q, matched_len, lr);
        if (link_len == (uint32_t)-1) {
          return next_->starts_with(q, matched_len, lr);
        }
        matched_len += link_len;
        if (matched_len >= len) return true;
      }
      if (!topo_.has_child(pos)) return matched_len == len;
      pos = topo_.child_pos(pos);
    }
    return true;  // query exhausted at internal node
  }

  auto leftmost_leaf_m(uint32_t pos) const -> int32_t {
    while (topo_.has_child(pos))
      pos = topo_.child_pos(pos);
    return (int32_t)topo_.leaf_id(pos);
  }

  template <bool rev = reverse, typename = std::enable_if_t<!rev>>
  auto successor(const key_type &key) const -> int32_t {
    std::vector<uint32_t> stack;
    uint32_t pos = 0, depth = 0;
    while (true) {
      uint32_t end = topo_.node_end(pos);
      if (depth >= key.size()) return leftmost_leaf_m(pos);
      uint8_t target = (uint8_t)key[depth];
      uint32_t p = pos;
      if (labels_[p] == terminator_) p++;
      while (p < end && labels_[p] < target) p++;
      if (p >= end) break;
      if (labels_[p] > target) return leftmost_leaf_m(p);
      // exact label match at p
      if (topo_.is_link(p))
        return (int32_t)topo_.leaf_id(p);  // link: treat as exact match
      if (!topo_.has_child(p)) {
        uint32_t lid = topo_.leaf_id(p);
        if (depth + 1 == key.size()) return (int32_t)lid;
        if (p + 1 < end) return leftmost_leaf_m(p + 1);
        break;
      }
      stack.push_back(p);
      pos = topo_.child_pos(p);
      depth++;
    }
    while (!stack.empty()) {
      uint32_t ap = stack.back(); stack.pop_back();
      if (ap + 1 < topo_.node_end(ap)) return leftmost_leaf_m(ap + 1);
    }
    return -1;
  }

  // returns matched length (-1 on mismatch)
  auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t override {
    if constexpr (!reverse_) {
      return -1;
    }
    assert(key_id < size());
    uint32_t link = links_[key_id];
    return match_link(key, begin, link);
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id, uint8_t partial_link) const -> uint32_t override {
    if constexpr (!reverse_) {
      return -1;
    }
    assert(key_id < size());
    uint32_t link = (links_[key_id] << 8) | partial_link;
    return match_link(key, begin, link);
  }

  auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool override {
    if constexpr (!reverse_) {
      return false;
    }
    assert(key_id < size());
    uint32_t link = links_[key_id];
    return starts_with_link(key, begin, link);
  }

  auto starts_with_link(const key_type &key, uint32_t begin, size_t link) const -> bool {
    uint32_t pos = link;
    uint32_t matched_len = begin;
    while (true) {
      if (matched_len >= (uint32_t)key.size()) return true;  // query exhausted
      if (topo_.is_link(pos)) {
        if (!next_->starts_with(key, matched_len, topo_.link_id(pos))) return false;
        // advance matched_len by the length of the nested pool entry
        uint32_t len = next_->match(key, matched_len, topo_.link_id(pos), labels_[pos]);
        if (len == (uint32_t)-1) return true;  // query was a prefix of nested entry
        matched_len += len;
      } else {
        uint8_t label = labels_[pos];
        if (label != terminator_) {
          if (matched_len >= (uint32_t)key.size()) return true;
          if (label != (uint8_t)key[matched_len]) return false;
          matched_len++;
        }
      }
      if (!topo_.has_parent(pos)) break;
      pos = topo_.parent_pos(pos);
    }
    return matched_len >= (uint32_t)key.size();
  }

  auto match_link(const key_type &key, uint32_t begin, size_t link) const -> uint32_t {
    uint32_t pos = link;
    uint32_t matched_len = begin;
    while (true) {
      if (topo_.is_link(pos)) {
        uint32_t link_len = next_->match(key, matched_len, topo_.link_id(pos), labels_[pos]);
        if (link_len == -1) {
          return -1;
        }
        matched_len += link_len;
      } else {
        uint8_t label = labels_[pos];
        if (label != terminator_) {
          if (label != key[matched_len]) {
            return -1;
          }
          matched_len++;
        }
      }
      if (!topo_.has_parent(pos)) {
        break;
      }
      pos = topo_.parent_pos(pos);
    }
    return matched_len - begin;
  }

#ifdef __COMPARE_MARISA__
  void to_louds_marisa(std::unique_ptr<LoudsMarisa> &out) const {
    topo_.to_louds_marisa(out);
  }

  auto get_topo() const -> const topo_t * {
    return &topo_;
  }
#endif

 private:
  void build_current_trie(const KeySet<key_type> &key_set, KeySet<key_type> &next_keys,
                          std::vector<uint8_t> *partial_links = nullptr) {
    if constexpr (reverse_) {
      links_.resize(key_set.size());
    }
    if (partial_links != nullptr) {
      partial_links->resize(key_set.size());
    }
  #ifdef __ENABLE_CACHE__
    if constexpr (!reverse_) {
      reserve_cache(key_set.size());
    }
  #endif

    uint32_t key_id = 0;
    uint32_t branch_id = 0;
    std::queue<Range> queue;
    queue.push(Range(0, key_set.size(), 0));
    while (!queue.empty()) {
      auto range = queue.front();
      queue.pop();
      assert(range.begin_ < range.end_);

      uint64_t has_child[4]{0}, is_link[4]{0};  // each range corresponds to a node
      uint32_t node_start = branch_id, num_branches = 0;

      uint32_t begin = range.begin_, end = range.begin_;  // group common fragments

      while (end < range.end_ && range.depth_ == key_set[end].length_) {  // skip empty suffixes
        if constexpr (reverse_) {
          if (partial_links == nullptr) {
            links_[key_set[end].id_] = branch_id;
          } else {
            links_[key_set[end].id_] = branch_id >> 8;
            (*partial_links)[key_set[end].id_] = branch_id & MASK(8);
          }
        }
        end++;
      }
      if (end > begin) {
        num_branches++;
        branch_id++;
        key_id++;
        labels_.emplace_back(terminator_);
        begin = end;
      }

      while (end < range.end_) {
        while (end < range.end_) {  // horizontal expansion
          if (key_set.get_label(end, range.depth_) != key_set.get_label(begin, range.depth_)) {
            break;
          }
          end++;
        }
        uint8_t branch_label = key_set.get_label(begin, range.depth_);
        assert(end > begin);

      #ifdef __ENABLE_CACHE__
        if constexpr (!reverse_) {
          cache_branch(node_start, num_branches, end - begin, branch_label);
        }
      #endif

        uint32_t depth;
        if (end == begin + 1) {  // single key; skip suffix
          depth = key_set[begin].length_;
        } else {
          depth = range.depth_ + 1;
          while (depth < key_set[begin].length_) {  // vertical extension
            if (key_set.get_label(begin, depth) != key_set.get_label(end - 1, depth)) {
              break;
            }
            depth++;
          }
        }

        if constexpr (!reverse_) {
          uint32_t path_len = depth - range.depth_;
          if (path_len < 2) {
            stats_.num_single++;
          } else if (path_len < link_cutoff_) {
            stats_.num_short++;
            stats_.total_short_len += path_len;
          } else {
            stats_.num_links++;
            stats_.total_link_len += path_len;
            if (path_len > stats_.max_link_len) stats_.max_link_len = path_len;
            uint32_t bucket = path_len - link_cutoff_;
            if (bucket >= (uint32_t)stats_.len_hist.size())
              stats_.len_hist.resize(bucket + 1, 0);
            stats_.len_hist[bucket]++;
          }
        }

        if (depth - range.depth_ >= link_cutoff_) {  // link
          if constexpr (!reverse_) {
          #ifdef __NO_BRANCHING_LABEL__
            auto [pos, len] = key_set.substr_range(begin, range.depth_, depth - range.depth_);
            labels_.emplace_back(terminator_);
            next_keys.emplace_back(key_set[begin].key_, pos, len);
          #else
            auto [pos, len] = key_set.substr_range(begin, range.depth_ + 1, depth - range.depth_ - 1);
            labels_.emplace_back(branch_label);  // store branching label in place for fast lookup
            next_keys.emplace_back(key_set[begin].key_, pos, len);
          #endif
          } else {
            auto [pos, len] = key_set.substr_range(begin, range.depth_, depth - range.depth_);
            labels_.emplace_back(terminator_);
            next_keys.emplace_back(key_set[begin].key_, pos, len);
          }
          SET_BIT(is_link[num_branches / 64], num_branches % 64);
        } else {  // label
          labels_.emplace_back(branch_label);
          depth = range.depth_ + 1;
        }

        bool empty = true;
        for (uint32_t i = begin; i < end; i++) {
          if (key_set[i].length_ > depth) {
            empty = false;
          } else if constexpr (reverse_) {
            if (partial_links == nullptr) {
              links_[key_set[i].id_] = branch_id;
            } else {
              links_[key_set[i].id_] = branch_id >> 8;
              (*partial_links)[key_set[i].id_] = branch_id & MASK(8);
            }
          }
        }
        if (!empty) {
          SET_BIT(has_child[num_branches / 64], num_branches % 64);
          queue.push(Range(begin, end, depth));
        } else {
          key_id++;
        }
        num_branches++;
        branch_id++;
        begin = end;
      }
      topo_.add_node(has_child, is_link, num_branches);
    }
    topo_.build(false);
    labels_.shrink_to_fit();
    sdsl::util::bit_compress(links_);

  }

  void build_next(const KeySet<key_type> &next_keys, int max_recursion, int mask = 0) {
    auto build_next_with_partial_links = [&]() {
      std::vector<uint8_t> next_partial_links;
      next_ = strpool_t::build_optimal(next_keys, &next_partial_links, trie_size_in_bits(), max_recursion, mask);
      uint32_t pos = -1;
      for (auto partial_link : next_partial_links) {
        pos = topo_.next_link(pos + 1);
        labels_[pos] = partial_link;
      }
    };
  #ifdef __NO_BRANCHING_LABEL__
    build_next_with_partial_links();
  #else
    if constexpr (reverse_) {
      build_next_with_partial_links();
    } else {
      next_ = strpool_t::build_optimal(next_keys, nullptr, trie_size_in_bits(), max_recursion, mask);
    }
  #endif
  #ifdef __ENABLE_CACHE__
    if constexpr (!reverse_) {
      fill_cache();
    }
  #endif
  }

  void build(const KeySet<key_type> &key_set, std::vector<uint8_t> *partial_links = nullptr,
             int max_recursion = 0, int mask = 0) override {
    KeySet<key_type> next_keys;
    build_current_trie(key_set, next_keys, partial_links);
    build_next(next_keys, max_recursion, mask);
  }

#ifdef __ENABLE_CACHE__
  void reserve_cache(uint32_t num_keys) {
    uint32_t cache_size = 256;
    while (cache_size < num_keys / cache_ratio_) {
      cache_size *= 2;
    }
    cache_.resize(cache_size);
    cache_mask_ = cache_size - 1;
  }

  void cache_branch(uint32_t parent, uint32_t child_id, float weight, uint8_t label) {
    auto cache_id = get_cache_id(parent, label);
    auto &cache = cache_[cache_id];
    if (weight > cache.weight()) {
      cache.parent_ = parent;
      cache.child_ = child_id;
      cache.set_weight(weight);
    }
  }

  auto get_cache_id(uint32_t pos, uint8_t label) const -> uint32_t {
    return (pos ^ (pos << 5) ^ label) & cache_mask_;
  }

  auto restore_label(uint32_t cache_id) const -> uint8_t {
    return restore_label(cache_id, cache_[cache_id].parent_);
  }

  auto restore_label(uint32_t cache_id, uint32_t pos) const -> uint8_t {
    return (cache_id ^ pos ^ (pos << 5)) & cache_mask_;
  }

  void fill_cache() {
    uint32_t cache_size = cache_.size();
    for (uint32_t i = 0; i < cache_size; i++) {
      auto &cache = cache_[i];
      auto parent = cache.parent_, branch_pos= cache.parent_ + cache.child_;
      if (cache.weight() < 0) {
        cache.parent_ = -1;
        cache.child_ = -1;
      } else {
        cache.child_ = (topo_.has_child(branch_pos) ? topo_.child_pos(branch_pos) : 0);
        if (!topo_.is_link(branch_pos)) {
          assert(restore_label(i) == labels_[branch_pos]);
          cache.set_link(-1);
        } else {
        #ifdef __NO_BRANCHING_LABEL__
          assert(topo_.num_links() < (1 << 24));
          cache.set_link((topo_.link_id(branch_pos) << 8) | labels_[branch_pos]);
        #else
          cache.set_link(topo_.link_id(branch_pos));
        #endif
        }
      }
    }
  }

  auto search_cache(uint32_t &pos, const key_type &key, uint32_t &matched_len) const -> uint32_t {
    auto cache_id = get_cache_id(pos, key[matched_len]);
    const auto &cache = cache_[cache_id];
    if (cache.parent_ != pos) {  // not cached
      return 0;
    }

    if (cache.link() == -1) {
      matched_len++;
    } else {
    #ifdef __NO_BRANCHING_LABEL__
      uint32_t link_len = next_->match(key, matched_len, cache.link() >> 8, cache.link() & MASK(8));
    #else
      matched_len++;
      uint32_t link_len = next_->match(key, matched_len, cache.link());
    #endif
      if (link_len == -1) {  // mismatch
        return -1;
      }
      matched_len += link_len;
    }

    if (cache.child_ == 0) {  // no child
      return 1;
    } else {
      pos = cache.child_;
      return 2;
    }
  }
#endif

  topo_t topo_;
  label_vec labels_;
  sdsl::int_vector<> links_;
  strpool_t *next_{nullptr};
  UnaryPathStats stats_;

#ifdef __ENABLE_CACHE__
  std::vector<Cache> cache_;
  uint32_t cache_mask_{0};
#endif

  template <typename K> friend class StringPool;
 public:
  void serialize(std::ostream& out) const override {
    archive::write(out,uint8_t(3));topo_.serialize(out);
    archive::write(out,uint64_t(labels_.size()));
    for(size_t i=0;i<labels_.size();++i) archive::write(out,labels_[i]);
    links_.serialize(out);next_->serialize(out);
#ifdef __ENABLE_CACHE__
    archive::vector(out,cache_);archive::write(out,cache_mask_);
#endif
  }
  void load_data(std::istream& in,unsigned depth) {
    topo_.load(in);uint64_t n;archive::read(in,n);
    if(n>UINT32_MAX) throw std::runtime_error("Too many labels");
    for(uint64_t i=0;i<n;++i) { uint8_t ch;archive::read(in,ch);labels_.emplace_back(ch); }
    links_.load(in);next_=strpool_t::load(in,depth+1);
#ifdef __ENABLE_CACHE__
    archive::vector(in,cache_);archive::read(in,cache_mask_);
#endif
  }

};

}  // namespace coco_paper

namespace coco_paper {
template<class Key> StringPool<Key>* StringPool<Key>::load(std::istream& in,unsigned depth) {
  if(depth>16) throw std::runtime_error("Too many tail recursion levels");
  uint8_t type;archive::read(in,type);
  if(type==0) { auto p=std::make_unique<SortedStringPool<Key>>();p->load_data(in);return p.release(); }
  if(type==2) { auto p=std::make_unique<FsstStringPool<Key>>();p->load_data(in);return p.release(); }
  if(type==3) { auto p=std::make_unique<MarisaCC<Key,true>>();p->load_data(in,depth);return p.release(); }
  throw std::runtime_error("Unknown tail container");
}
}
