#pragma once

#define __COMPARE_FST__

#include "utils.hpp"
#include "static_vector.hpp"
#include "louds_sparse_cc.hpp"
#include "bit_vector.hpp"
#include "key_set.hpp"
#include "marisa_cc.hpp"

#include <limits>
#include <vector>
#include <queue>


namespace coco_paper {

// louds-sparse trie used for building CoCo-trie
template <typename Key>
class FstCC {
 public:
  using key_type = Key;
  using strpool_t = StringPool<key_type>;
  using label_vec = StaticVector<uint8_t>;
  using topo_t = LoudsSparseCC;
  using bitvec_t = BitVector;

  static constexpr uint32_t link_cutoff_ = 4;  // suffixes of length above this value will be moved to string pool

  void print_space_cost_breakdown() const {
    size_t topo = topo_.size_in_bits();
    size_t link = is_link_.size_in_bits();
    size_t data = labels_.size_in_bits();
    next_->space_cost_breakdown(topo, link, data);
    printf("topology: %lf MB, link: %lf MB, data: %lf MB\n", (double)topo/mb_bits, (double)link/mb_bits, (double)data/mb_bits);
  }

  // FST distinguishes two kinds of unary paths:
  //   shared: multiple keys share a common prefix extension → kept in main trie char-by-char
  //   terminal: only one key remains (unique suffix) → may be compressed if len >= link_cutoff_
  struct UnaryPathStats {
    uint64_t num_shared_steps{0};   // individual char-steps through shared unary extensions
    uint64_t num_inplace_term{0};   // short terminal suffixes (len < link_cutoff_), kept in trie
    uint64_t total_inplace_len{0};
    uint64_t num_links{0};          // long terminal suffixes (len >= link_cutoff_), compressed
    uint64_t total_link_len{0};
    uint64_t max_link_len{0};
    std::vector<uint64_t> len_hist;  // len_hist[i] = # links with length (link_cutoff_ + i)
  };

  void print_unary_path_stats() const {
    printf("--- Unary Path Analysis (FST, cutoff=%u) ---\n", link_cutoff_);
    printf("  shared unary steps in main trie: %llu\n", stats_.num_shared_steps);
    printf("  terminal suffixes: inplace=%llu (avg len=%.2f)  links=%llu (%.1f%% of terminal)\n",
           stats_.num_inplace_term,
           stats_.num_inplace_term > 0 ? (double)stats_.total_inplace_len / stats_.num_inplace_term : 0.0,
           stats_.num_links,
           (stats_.num_links + stats_.num_inplace_term) > 0
               ? 100.0 * stats_.num_links / (stats_.num_links + stats_.num_inplace_term) : 0.0);
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
  }

 public:
  // helper class for walking down the trie and traversing macro-node keys; used by the optimizer
  struct walker {
    using trie_t = FstCC;

    key_type key_;
    const trie_t *trie_;
    uint32_t pos_;
    uint32_t level_;  // number of walked levels

    walker(const trie_t *trie, uint32_t pos) : trie_(trie), pos_(pos), level_(1) {
      uint8_t label = trie_->get_label(pos_);
      if (label != terminator_) {
        key_.push_back(label);
      }
    }

    walker(const walker &other) : key_(other.key_), trie_(other.trie_), pos_(other.pos_), level_(other.level_) {}

    // is the key legitimately terminated?
    auto valid() const -> bool {
      return !trie_->topo_.has_child(pos_);
    }

    // is the current key a prefix key?
    auto prefix_key() const -> bool {
      return trie_->get_label(pos_) == terminator_;
    }

    auto key() const -> const key_type & {
      return key_;
    }

    // move to the leftmost label in node
    void move_to_front() {
      uint32_t front = trie_->topo_.node_start(pos_);
      if (trie_->get_label(front) == terminator_) {
        pos_ = front;
        key_.pop_back();
      } else {
        pos_ = front;
        key_.back() = trie_->get_label(pos_);
      }
    }

    // move to the rightmost label in node
    void move_to_back() {
      uint32_t back = trie_->topo_.node_end(pos_) - 1;
      if (trie_->get_label(pos_) == terminator_) {
        pos_ = back;
        key_.push_back(trie_->get_label(pos_));
      } else {
        pos_ = back;
        key_.back() = trie_->get_label(pos_);
      }
    }

    // does NOT regress if current level is already greater than `max_level`
    void get_min_key(uint32_t max_level = std::numeric_limits<uint32_t>::max()) {
      while (level_ < max_level && trie_->topo_.has_child(pos_)) {
        pos_ = trie_->topo_.child_pos(pos_);  // keep taking the leftmost branch
        uint8_t label = trie_->get_label(pos_);
        if (label != terminator_) {
          key_.push_back(label);
        }
        level_++;
      }
    }

    // does NOT regress if current level is already greater than `max_level`
    void get_max_key(uint32_t max_level = std::numeric_limits<uint32_t>::max()) {
      while (level_ < max_level && trie_->topo_.has_child(pos_)) {
        pos_ = trie_->topo_.child_pos(pos_);
        pos_ = trie_->topo_.node_end(pos_) - 1;  // keep taking the rightmost branch
        uint8_t label = trie_->get_label(pos_);
        if (label != terminator_) {
          key_.push_back(label);
        }
        level_++;
      }
    }

    // make sure to call `get_min_key(max_level)` before calling this
    // return true on success and false if there is no more key
    auto next(uint32_t max_level = std::numeric_limits<uint32_t>::max()) -> bool {
      while (level_ > 0) {
        if (trie_->get_label(pos_) == terminator_) {  // terminator is never the last label
          pos_++;
          key_.push_back(trie_->get_label(pos_));
          get_min_key(max_level);
          return true;
        } else if (pos_ + 1 < trie_->topo_.size() && !trie_->topo_.louds(pos_ + 1)) {  // not the last label in node
          pos_++;
          key_.back() = trie_->get_label(pos_);
          get_min_key(max_level);
          return true;
        }
        // last label in node; regress to parent
        key_.pop_back();
        pos_ = trie_->topo_.parent_pos(pos_);
        level_--;
      }
      return false;
    }

    // move to the leftmost next-level node in subtrie
    // return true if found and false otherwise
    // calling any of the `move_down` functions after false is returned is undefined behavior
    auto move_down_one_level_left() -> bool {
      uint32_t next_level = level_ + 1;

      get_min_key(next_level);
      while (level_ < next_level) {
        assert(!trie_->topo_.has_child(pos_));  // key terminates before `next_level`
        while (true) {
          if (trie_->get_label(pos_) != terminator_) {
            key_.pop_back();
          }
          uint32_t end = trie_->topo_.node_end(pos_);
          uint32_t next = trie_->topo_.next_child(pos_ + 1);
          if (next < end) {  // trace next branch
            pos_ = next;
            key_.push_back(trie_->get_label(pos_));
            break;
          }
          // regress to parent
          pos_ = trie_->topo_.parent_pos(pos_);
          level_--;
          if (level_ == 0) {
            return false;
          }
        }
        get_min_key(next_level);
      }
      return true;
    }

    // move to the rightmost next-level node in subtrie
    // return true if found and false otherwise
    // calling any of the `move_down` functions after false is returned is undefined behavior
    auto move_down_one_level_right() -> bool {
      uint32_t next_level = level_ + 1;

      get_max_key(next_level);
      while (level_ < next_level) {
        assert(!trie_->topo_.has_child(pos_));  // key terminates before `next_level`
        while (true) {
          if (trie_->get_label(pos_) != terminator_) {
            key_.pop_back();
          }
          uint32_t start = trie_->topo_.node_start(pos_);
          uint32_t prev = trie_->topo_.prev_child(pos_ - 1);
          if (prev >= start) {  // trace previous branch
            pos_ = prev;
            key_.push_back(trie_->get_label(pos_));
            break;
          }
          // regress to parent
          pos_ = trie_->topo_.parent_pos(pos_);
          level_--;
          if (level_ == 0) {
            return false;
          }
        }
        get_max_key(next_level);
      }
      return true;
    }
  };

 private:
  class TempStringPool : public StringPool<key_type> {
   public:
    using strpool_t = typename succinct::tries::compressed_string_pool<uint8_t>;

    TempStringPool() = default;

    ~TempStringPool() = default;

    void build(const KeySet<key_type> &key_set, std::vector<uint8_t> *partial_links,
               int max_recursion = 0, int mask = 0) override {
      keys_ = key_set;
    }

    void build(KeySet<key_type> &&key_set) {
      keys_ = key_set;
    }

    auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t override {
      return -1;  // unused
    }

    auto match(const key_type &key, uint32_t begin, uint32_t key_id,
               uint8_t partial_link) const -> uint32_t override {
      return -1;  // unused
    }

    auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool override {
      return false;  // unused during build
    }

    auto size() const -> uint32_t override {
      return keys_.size();
    }

    auto size_in_bytes() const -> size_t override {
      return 0;  // unused
    }

    auto size_in_bits() const -> size_t override {
      return 0; // unused
    }

    void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const override {}  // unused
   private:
    KeySet<key_type> keys_;

    template <typename K> friend class FstCC;
    template <typename K, typename T, typename C> friend class CoCoCC;
  };

 private:
  struct Range {
    uint32_t begin_{0};
    uint32_t end_{0};
    uint32_t depth_{0};
    uint32_t lcp_{0};

    Range() = default;

    Range(uint32_t begin, uint32_t end, uint32_t depth, uint32_t lcp)
        : begin_(begin), end_(end), depth_(depth), lcp_(lcp) {}
  };

 public:
  FstCC() = default;

  ~FstCC() {
    delete next_;
  }

  template <typename Iterator>
  void build(Iterator begin, Iterator end, bool sorted = false,
             int max_recursion = 0, int mask = 0) {
    KeySet<key_type> key_set;
    while (begin != end) {
      key_set.emplace_back(&(*begin));
      ++begin;
    }
    assert(!key_set.empty());
    if (!sorted) {
      key_set.sort();
    }
    build(key_set, false, max_recursion, mask);
  }

  void clear() {
    labels_.clear();
    topo_.clear();
  }

  auto get_label(uint32_t idx) const -> uint8_t {
    assert(idx < labels_.size());
    return labels_.at(idx);
  }

  auto lookup(const key_type &key) const -> uint32_t {
    uint16_t len = key.size(), matched_len = 0;
    uint32_t pos = 0;

    while (matched_len < len) {
      uint32_t end = topo_.node_end(pos);
      pos = labels_.find(key[matched_len], pos, end);
      if (pos == end) {  // mismatch
        return -1;
      }
      assert(get_label(pos) == key[matched_len]);
      matched_len++;

      if (!topo_.has_child(pos)) {
        auto leaf_id = topo_.leaf_id(pos);
        if (!is_link_.get(leaf_id)) {
          return matched_len == len ? leaf_id : -1;
        } else {
          uint32_t link = is_link_.rank1(leaf_id);
          return next_->match(key, matched_len, link) == len - matched_len ? leaf_id : -1;
        }
      }
      pos = topo_.child_pos(pos);
    }

    if (get_label(pos) == terminator_) {  // prefix key
      return topo_.leaf_id(pos);
    }
    return -1;
  }

  // Returns true if any key in the trie starts with q (i.e., q is a prefix of some stored key).
  auto contains_prefix(const key_type &q) const -> bool {
    uint32_t pos = 0;
    uint32_t matched_len = 0, len = q.size();

    while (matched_len < len) {
      uint32_t end = topo_.node_end(pos);
      pos = labels_.find(q[matched_len], pos, end);
      if (pos == end) return false;
      matched_len++;

      if (!topo_.has_child(pos)) {
        uint32_t lid = topo_.leaf_id(pos);
        if (!is_link_.get(lid)) {
          return matched_len == len;  // non-link leaf: exact match required
        }
        uint32_t lr = is_link_.rank1(lid);
        return next_->starts_with(q, matched_len, lr);
      }
      pos = topo_.child_pos(pos);
    }
    return true;  // consumed all of q: subtree below guarantees at least one key
  }

  // return the cutoffs between levels
  auto get_level_boundaries() const -> std::vector<uint32_t> {
    return topo_.get_level_boundaries();
  }

  auto size_in_bytes() const -> size_t {
    return labels_.size_in_bytes() + topo_.size_in_bytes() + next_->size_in_bytes() + is_link_.size_in_bytes();
  }

  auto size_in_bits() const -> size_t {
    return size_in_bytes() * 8;
  }

  auto trie_size_in_bits() const -> size_t {
    return (labels_.size_in_bytes() + topo_.size_in_bytes() + is_link_.size_in_bytes()) * 8;
  }

  auto get_topo() const -> const topo_t* { return &topo_; }

 private:
  auto leftmost_leaf(uint32_t pos) const -> int32_t {
    while (topo_.has_child(pos))
      pos = topo_.child_pos(pos);
    return (int32_t)topo_.leaf_id(pos);
  }

 public:
  auto successor(const key_type &key) const -> int32_t {
    std::vector<uint32_t> stack;
    uint32_t pos = 0, depth = 0;
    while (true) {
      uint32_t end = topo_.node_end(pos);
      if (depth >= key.size()) return leftmost_leaf(pos);
      uint8_t target = (uint8_t)key[depth];
      uint32_t p = pos;
      if (get_label(p) == terminator_) p++;
      while (p < end && get_label(p) < target) p++;
      if (p >= end) break;
      if (get_label(p) > target) return leftmost_leaf(p);
      // exact label match at p
      if (!topo_.has_child(p)) {
        uint32_t lid = topo_.leaf_id(p);
        if (!is_link_.get(lid)) {
          if (depth + 1 == key.size()) return (int32_t)lid;
          if (p + 1 < end) return leftmost_leaf(p + 1);
          break;
        }
        return (int32_t)lid;  // link: treat as exact match
      }
      stack.push_back(p);
      pos = topo_.child_pos(p);
      depth++;
    }
    while (!stack.empty()) {
      uint32_t ap = stack.back(); stack.pop_back();
      if (ap + 1 < topo_.node_end(ap)) return leftmost_leaf(ap + 1);
    }
    return -1;
  }

 // -------------------------------------------------------------------------
  // Iterator-based range query
  // -------------------------------------------------------------------------
  struct RangeIter {
    const FstCC *trie_  = nullptr;
    uint32_t pos_       = 0;
    int32_t  rank_      = -1;
    uint32_t level_     = 0;
    key_type key_;
    key_type end_key_;   // exclusive upper bound
    bool     valid_     = false;

    bool valid()             const { return valid_; }
    int32_t      rank()      const { return rank_; }
    const key_type& key()    const { return key_; }

    // Advance to the next key in sorted order; returns valid().
    bool next() {
      while (level_ > 0) {
        uint8_t lbl = trie_->get_label(pos_);
        if (lbl == terminator_) {
          // terminator is always the first edge in its node and never the last
          pos_++;
          uint8_t nl = trie_->get_label(pos_);
          if (nl != terminator_) key_.push_back(nl);
          go_leftmost();
        } else if (pos_ + 1 < trie_->topo_.size() && !trie_->topo_.louds(pos_ + 1)) {
          pos_++;                             // next sibling in same node
          key_.back() = trie_->get_label(pos_);
          go_leftmost();
        } else {
          key_.pop_back();                    // last edge — backtrack to parent
          pos_ = trie_->topo_.parent_pos(pos_);
          level_--;
          continue;
        }
        rank_  = (int32_t)trie_->topo_.leaf_id(pos_);
        valid_ = is_below_end(pos_);
        return valid_;
      }
      rank_  = -1;
      valid_ = false;
      return false;
    }

   private:
    void go_leftmost() {
      while (trie_->topo_.has_child(pos_)) {
        pos_ = trie_->topo_.child_pos(pos_);
        level_++;
        uint8_t l = trie_->get_label(pos_);
        if (l != terminator_) key_.push_back(l);
      }
    }

    // Stop condition: full key at p is strictly less than end_key_.
    // For link nodes key_ is a truncated main-trie prefix; if key_ happens to be a
    // proper prefix of end_key_ we consult the string pool to avoid a false positive.
    // This pool call only fires when the iterator reaches end_key_'s own position
    // (at most once per range_count_iter call), so it is O(1) amortised.
    bool is_below_end(uint32_t p) const {
      if (key_ >= end_key_) return false;
      uint32_t lid = trie_->topo_.leaf_id(p);
      if (!trie_->is_link_.get(lid)) return true;
      // Link: key_ == full_key[0:key_.size()]. Check whether key_ is a proper prefix
      // of end_key_ (the only case where key_ < end_key_ could be a false positive).
      size_t k = key_.size(), e = end_key_.size();
      if (k >= e) return true;   // key_ < end_key_ and no prefix ambiguity
      for (size_t i = 0; i < k; i++) {
        if (key_[i] != (uint8_t)end_key_[i]) return true;  // differ before k → clear
      }
      // key_ == end_key_[0:k]: consult pool to decide direction of the suffix.
      // match() returns pool_len if end_key_[k:k+pool_len] fully matches the pool
      // entry, or (uint32_t)-1 on any mismatch / key too short.
      // matched == e-k  ↔  full_key == end_key_  →  NOT valid (exclusive bound).
      // matched != e-k  ↔  full_key  < end_key_  →  valid (sorted-order invariant).
      uint32_t lr = trie_->is_link_.rank1(lid);
      uint32_t matched = trie_->next_->match(end_key_, (uint32_t)k, lr);
      return matched != (uint32_t)(e - k);
    }

    friend class FstCC;
  };

  // Return an iterator over keys in [key, r). One trie traversal for key;
  // stop condition uses string comparison + O(1)-amortised pool check.
  auto lower_bound(const key_type& key, const key_type& r) const -> RangeIter {
    struct Frame { uint32_t pos; bool pushed; };
    std::vector<Frame> bt;

    RangeIter it;
    it.trie_    = this;
    it.end_key_ = r;

    uint32_t pos = 0, depth = 0;

    auto land = [&](uint32_t p, bool go_left) {
      uint8_t lbl = get_label(p);
      it.pos_ = p;
      if (lbl != terminator_) it.key_.push_back(lbl);
      it.level_++;
      if (go_left) it.go_leftmost();
      it.rank_  = (int32_t)topo_.leaf_id(it.pos_);
      it.valid_ = it.is_below_end(it.pos_);
    };

    while (true) {
      uint32_t end = topo_.node_end(pos);

      if (depth >= (uint32_t)key.size()) { land(pos, true); return it; }

      uint8_t target = (uint8_t)key[depth];
      uint32_t p = pos;
      if (get_label(p) == terminator_) p++;
      while (p < end && get_label(p) < target) p++;

      if (p >= end) break;

      uint8_t lbl = get_label(p);
      if (lbl > target) { land(p, true); return it; }

      // exact label match at p
      if (!topo_.has_child(p)) {
        uint32_t lid = topo_.leaf_id(p);
        if (!is_link_.get(lid)) {
          if (depth + 1 == (uint32_t)key.size()) { land(p, false); return it; }
          if (p + 1 < end)                        { land(p + 1, true); return it; }
          break;
        }
        land(p, false); return it;    // link leaf: treat as inclusive match
      }

      bool pushed = (lbl != terminator_);
      if (pushed) it.key_.push_back(lbl);
      it.level_++;
      bt.push_back({p, pushed});
      pos   = topo_.child_pos(p);
      depth++;
    }

    while (!bt.empty()) {
      auto [ap, pushed] = bt.back(); bt.pop_back();
      if (pushed) it.key_.pop_back();
      it.level_--;
      if (ap + 1 < topo_.node_end(ap)) { land(ap + 1, true); return it; }
    }
    return it;   // rank_ = -1, no successor
  }

 private:
  void build(const KeySet<key_type> &key_set, bool temp = false,
             int max_recursion = 0, int mask = 0) {
    KeySet<key_type> suffixes;

    auto lcp = [&](uint32_t begin, uint32_t end, uint32_t depth) -> uint32_t {
      assert(end > begin);
      assert(depth <= key_set[begin].length_ && depth <= key_set[end - 1].length_);
      if (end == begin + 1) {
        return key_set[begin].length_ - depth;
      }
      uint32_t len = std::min(key_set[begin].length_, key_set[end - 1].length_) - depth;
      uint32_t ret = 0;
      while (ret < len) {
        if (key_set.get_label(begin, depth + ret) != key_set.get_label(end - 1, depth + ret)) {
          break;
        }
        ret++;
      }
      return ret;
    };

    auto is_same_key = [&](const Range &range) -> bool {
      return range.lcp_ == key_set[range.end_ - 1].length_ - range.depth_;
    };

    std::queue<Range> queue;
    queue.push(Range(0, key_set.size(), 0, lcp(0, key_set.size(), 0)));
    while (!queue.empty()) {
      auto range = queue.front();
      queue.pop();
      assert(range.begin_ < range.end_);

      uint64_t has_child[4]{0};    // each range corresponds to a node
      uint32_t num_branches = 0;

      if (range.lcp_ > 0) {
        labels_.emplace_back(key_set.get_label(range.begin_, range.depth_));
        if (!is_same_key(range)) {  // not a suffix: shared unary path, kept char-by-char in main trie
          stats_.num_shared_steps++;
          SET_BIT(has_child[0], 0);
          topo_.add_node(has_child, 1);
          queue.push(Range(range.begin_, range.end_, range.depth_ + 1, range.lcp_ - 1));
        } else if (range.lcp_ < link_cutoff_) {  // in place
          if (range.lcp_ == 1) {  // last label of short terminal suffix
            stats_.num_inplace_term++;
            stats_.total_inplace_len += 1;
            topo_.add_node(has_child, 1);
            is_link_.append0();
          } else {  // not last label: continuation of short terminal suffix
            stats_.total_inplace_len++;
            SET_BIT(has_child[0], 0);
            topo_.add_node(has_child, 1);
            queue.push(Range(range.begin_, range.end_, range.depth_ + 1, range.lcp_ - 1));
          }
        } else {  // link: long terminal suffix, compressed to string pool
          uint32_t link_len = range.lcp_ - 1;  // first char stored in main trie, rest in pool
          stats_.num_links++;
          stats_.total_link_len += link_len;
          if (link_len > stats_.max_link_len) stats_.max_link_len = link_len;
          uint32_t bucket = range.lcp_ - link_cutoff_;  // bucket by full lcp (including first char)
          if (bucket >= (uint32_t)stats_.len_hist.size())
            stats_.len_hist.resize(bucket + 1, 0);
          stats_.len_hist[bucket]++;
          topo_.add_node(has_child, 1);
          is_link_.append1();
          suffixes.emplace_back(key_set[range.begin_].key_, range.depth_ + 1, range.lcp_ - 1);
        }
        continue;
      }

      // group common fragments
      uint32_t begin = range.begin_, end = range.begin_;

      while (end < range.end_ && range.depth_ == key_set[end].length_) {  // skip empty suffixes
        end++;
      }
      if (end > begin) {
        num_branches++;
        is_link_.append0();
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
        assert(end > begin);

        labels_.emplace_back(key_set.get_label(begin, range.depth_));

        uint32_t depth = range.depth_ + 1;
        if (key_set[end - 1].length_ > depth) {  // subtree not empty
          SET_BIT(has_child[num_branches / 64], num_branches % 64);
          queue.push(Range(begin, end, depth, lcp(begin, end, depth)));
        } else {
          is_link_.append0();
        }
        num_branches++;
        begin = end;
      }
      topo_.add_node(has_child, num_branches);
    }
    topo_.build();
    is_link_.build();
    labels_.shrink_to_fit();

    if (!temp) {
      next_ = strpool_t::build_optimal(suffixes, nullptr, trie_size_in_bits(), max_recursion, mask);
    } else {
      auto temp = new TempStringPool();
      temp->build(std::move(suffixes));
      next_ = temp;
    }
  }

  label_vec labels_;

  topo_t topo_;

  bitvec_t is_link_;

  strpool_t *next_{nullptr};

  UnaryPathStats stats_;

  friend class walker;
  template <typename K, typename C> friend class CoCoOptimizer;
  template <typename K, typename T, typename C> friend class CoCoCC;
};

}  // namespace coco_paper
