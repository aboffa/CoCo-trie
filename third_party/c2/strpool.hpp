#pragma once

#include "utils.hpp"
#include "key_set.hpp"
#include "compressed_string_pool.hpp"
#include "succinct/mapper.hpp"
#include "../fsst/fsst.h"

#include <sdsl/int_vector.hpp>
#include <vector>


namespace coco_paper {

template <typename Key, bool reverse>
class MarisaCC;

template <typename Key>
class FsstStringPool;

template <typename Key>
class SortedStringPool;

template <typename Key>
class RepairStringPool;

template <typename Key>
class StringPool {
 public:
  using key_type = Key;
  using sorted_t = SortedStringPool<key_type>;
  using repair_t = RepairStringPool<key_type>;
  using fsst_t = FsstStringPool<key_type>;
  using trie_t = MarisaCC<key_type, true>;

  enum class Type {
    SORTED = 0, REPAIR = 1, FSST = 2, TRIE = 3,
  };

  static constexpr int SORTED_FLAG   = BIT(static_cast<int>(Type::SORTED));
  static constexpr int REPAIR_FLAG   = BIT(static_cast<int>(Type::REPAIR));
  static constexpr int FSST_FLAG     = BIT(static_cast<int>(Type::FSST));

  // sample size used for estimation
  static constexpr size_t sample_size_ = (1 << 22);  // 4MB
  // recurse only if improvement is above this threshold
  static constexpr size_t recurse_threshold_ = 10;

  StringPool() = default;

  virtual ~StringPool() = default;
  virtual void serialize(std::ostream&) const { throw std::runtime_error("Unsupported tail container"); }
  static StringPool* load(std::istream&, unsigned depth = 0);


 private:
  static auto build_tail(const KeySet<key_type> &keys, const KeySet<key_type> &sorted_rev_keys,
                         std::vector<uint8_t> *partial_links, int mask) -> StringPool * {
    StringPool *ret;
    if (mask & FSST_FLAG) {
      ret = new fsst_t();
      ret->build(sorted_rev_keys, partial_links);
    } else if (mask & REPAIR_FLAG) {
      ret = new repair_t();
      ret->build(keys, partial_links);
    } else {
      ret = new sorted_t();
      ret->build(sorted_rev_keys, partial_links);
    }
    return ret;
  }

  static auto build_recursive(const KeySet<key_type> &keys, const KeySet<key_type> &sorted_rev_keys,
                              std::vector<uint8_t> *partial_links, size_t total_trie_size, size_t prev_tail_estimate,
                              int max_recursion, int mask) -> StringPool * {
    if (max_recursion == 0) {
      return build_tail(keys, sorted_rev_keys, partial_links, mask);
    }

    KeySet<key_type> next_keys;
    auto next_trie = new MarisaCC<key_type, true>();
    next_trie->build_current_trie(sorted_rev_keys, next_keys, partial_links);
    auto sorted_rev_next_keys = next_keys;
    sorted_rev_next_keys.set_reverse(true);
    sorted_rev_next_keys.sort();
    size_t tail_estimate = (mask & FSST_FLAG) ? fsst_t::estimate_space_cost(next_keys, false) :
                           (mask & REPAIR_FLAG) ? repair_t::estimate_space_cost(next_keys, false) :
                           sorted_t::estimate_space_cost(sorted_rev_next_keys, false);
    size_t trie_size = next_trie->trie_size_in_bits();
    size_t prev_estimate = total_trie_size + prev_tail_estimate;
    size_t cur_estimate = total_trie_size + trie_size + tail_estimate;
    if ((100 + recurse_threshold_)*cur_estimate >= 100*prev_estimate) {
      delete next_trie;
      return build_tail(keys, sorted_rev_keys, partial_links, mask);
    }

    std::vector<uint8_t> next_partial_links;
    next_trie->next_ = build_recursive(next_keys, sorted_rev_next_keys, &next_partial_links, total_trie_size +
                                       next_trie->size_in_bits(), tail_estimate, max_recursion - 1, mask);
    uint32_t pos = -1;
    for (auto partial_link : next_partial_links) {
      pos = next_trie->topo_.next_link(pos + 1);
      next_trie->labels_[pos] = partial_link;
    }
    return next_trie;
  }

 public:
  static auto build_optimal(const KeySet<key_type> &keys, std::vector<uint8_t> *partial_links = nullptr,
                            size_t first_trie_size = 0, int max_recursion = 0, int mask = 0) -> StringPool * {
    if (keys.size() == 0) {
      auto *empty = new sorted_t();
      empty->build(keys, partial_links);
      return empty;
    }
    if (mask == 0) {
      mask = FSST_FLAG;
    }
    StringPool *ret;
    auto sorted_rev_keys = keys;
    sorted_rev_keys.set_reverse(true);
    sorted_rev_keys.sort();
    size_t tail_estimate = (mask & FSST_FLAG) ? fsst_t::estimate_space_cost(keys, false) :
                           (mask & REPAIR_FLAG) ? repair_t::estimate_space_cost(keys, false) :
                           sorted_t::estimate_space_cost(sorted_rev_keys, false);
    return build_recursive(keys, sorted_rev_keys, partial_links, first_trie_size, tail_estimate, max_recursion, mask);
  }

  virtual void build(const KeySet<key_type> &keys, std::vector<uint8_t> *partial_links = nullptr,
                     int max_recursion = 0, int mask = 0) = 0;

  virtual auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t = 0;

  virtual auto match(const key_type &key, uint32_t begin, uint32_t key_id,
                     uint8_t partial_link) const -> uint32_t = 0;

  // Returns true if the pool entry at key_id starts with key[begin:].
  virtual auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool = 0;

  virtual auto size() const -> uint32_t = 0;

  virtual auto size_in_bytes() const -> size_t = 0;

  virtual auto size_in_bits() const -> size_t = 0;

  virtual void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const = 0;
};


template <typename Key>
class SortedStringPool : public StringPool<Key> {
 public:
  using key_type = Key;

  SortedStringPool() = default;

  ~SortedStringPool() = default;

  // partial link: 8 low bits of the actual link
  void build(const KeySet<key_type> &sorted_rev_keys, std::vector<uint8_t> *partial_links = nullptr,
             int max_recursion = 0, int mask = 0) override {
    links_.resize(sorted_rev_keys.size());
    if (partial_links != nullptr) {
      partial_links->resize(sorted_rev_keys.size());
    }

    const typename KeySet<key_type>::Fragment *next = nullptr;
    for (size_t i = sorted_rev_keys.size(); i > 0; i--) {
      auto &cur = sorted_rev_keys[i - 1];
      if (next != nullptr) {
        uint32_t len = std::min(cur.size(), next->size());
        uint32_t match = 0;
        while (match < len && cur.get_label(match, true) == next->get_label(match, true)) {
          match++;
        }
        if (match == cur.size()) {  // deduplicate prefix key
          size_t link = labels_.size() - match - 1;
          if (partial_links == nullptr) {
            links_[cur.id_] = link;
          } else {
            links_[cur.id_] = link >> 8;
            (*partial_links)[cur.id_] = link & MASK(8);
          }
          continue;
        }
      }
      if (partial_links == nullptr) {
        links_[cur.id_] = labels_.size();
      } else {
        links_[cur.id_] = labels_.size() >> 8;
        (*partial_links)[cur.id_] = labels_.size() & MASK(8);
      }
      cur.append_to(labels_, true);
      next = &cur;
    }
    sdsl::util::bit_compress(links_);
    labels_.shrink_to_fit();
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t override {
    assert(key_id < size());
    size_t link = links_[key_id];
    return match_link(key, begin, link);
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id, uint8_t partial_link) const -> uint32_t override {
    assert(key_id < size());
    size_t link = (links_[key_id] << 8) | partial_link;
    return match_link(key, begin, link);
  }

  auto match_link(const key_type &key, uint32_t begin, size_t link) const -> uint32_t {
    uint32_t matched_len = 0;
    while (labels_[link + matched_len] != terminator_) {
      if (begin + matched_len >= key.size() || labels_[link + matched_len] != uint8_t(key[begin + matched_len])) {
        return -1;
      }
      matched_len++;
    }
    return matched_len;
  }

  auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool override {
    assert(key_id < size());
    size_t link = links_[key_id];
    uint32_t i = 0;
    while (labels_[link + i] != terminator_) {
      if (begin + i >= key.size()) return true;
      if (labels_[link + i] != (uint8_t)key[begin + i]) return false;
      i++;
    }
    return begin + i == key.size();
  }

  auto size() const -> uint32_t override {
    return links_.size();
  }

  auto size_in_bytes() const -> size_t override {
    return sdsl::size_in_bytes(links_) + labels_.size() * sizeof(uint8_t) + sizeof(labels_);
  }

  auto size_in_bits() const -> size_t override {
    return size_in_bytes() * 8;
  }

  void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const {
    link += sdsl::size_in_bytes(links_) * 8;
    data += labels_.size() * sizeof(uint8_t) * 8;
  }

  static auto estimate_space_cost(const KeySet<key_type> &sorted_rev_keys, bool partial_links = false) -> size_t {
    auto [lcp_size, sorted_size] = sorted_rev_keys.lcp_size();
    return estimate_space_cost(sorted_size, sorted_rev_keys.size(), partial_links);
  }

  static auto estimate_space_cost(size_t sorted_size, size_t num_keys, bool partial_links = false) -> size_t {
    if (sorted_size == 0) {
      return 0;
    }
    int link_width = partial_links ? std::max(56 - __builtin_clzll(sorted_size), 1) :
                     std::max(64 - __builtin_clzll(sorted_size), 1);
    return sorted_size * 8 + link_width * num_keys;
  }
  void serialize(std::ostream& out) const override {
    archive::write(out,uint8_t(0));archive::vector(out,labels_);links_.serialize(out);
  }
  void load_data(std::istream& in) { archive::vector(in,labels_);links_.load(in); }
 private:
  std::vector<uint8_t> labels_;
  sdsl::int_vector<> links_;
};

template <typename Key>
class RepairStringPool : public StringPool<Key> {
 public:
  using key_type = Key;
  using strpool_t = typename succinct::tries::compressed_string_pool<uint8_t>;

  RepairStringPool() = default;

  ~RepairStringPool() = default;

  // partial links: (4 low bits of this link) | (4 low bits of next link)
  void build(const KeySet<key_type> &key_set, std::vector<uint8_t> *partial_links = nullptr,
             int max_recursion = 0, int mask = 0) override {
    std::vector<uint8_t> keys;
    for (const auto &frag : key_set.fragments_) {
      frag.append_to(keys);
    }
    build(keys, partial_links);
  }

  void build(const std::vector<uint8_t> &keys, std::vector<uint8_t> *partial_links = nullptr) {
    if (!keys.empty()) {
      strpool_t strpool(keys);
      strpool_.swap(strpool);
      if (partial_links == nullptr) {
        links_.swap(strpool_.m_positions);
      } else {
        auto enu = typename succinct::elias_fano::select_enumerator(strpool_.m_positions, 1);

        size_t n = (strpool_.m_positions.select(strpool_.size()) >> 4) + 1;
        size_t m = strpool_.size() + 1;
        typename succinct::elias_fano::elias_fano_builder builder(n, m);
        uint8_t prev = 0;
        builder.push_back(0);
        partial_links->reserve(m - 1);
        for (size_t i = 0; i < m - 1; i++) {
          size_t link = enu.next();
          builder.push_back(link >> 4);
          partial_links->push_back(prev | link & MASK(4));
          prev = (link & MASK(4)) << 4;
        }
        strpool_.release_positions();
        typename succinct::elias_fano(&builder, false).swap(links_);
      }
    } else {
      strpool_t strpool;
      strpool_.swap(strpool);
    }
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t override {
    assert(key_id < size());
    auto [link, end] = links_.select_range(key_id);
    return match_link(key, begin, link, end);
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id, uint8_t partial_link) const -> uint32_t override {
    assert(key_id < size());
    auto [link, end] = links_.select_range(key_id);
    link = (link << 4) | (partial_link >> 4);
    end = (end << 4) | (partial_link & MASK(4));
    return match_link(key, begin, link, end);
  }

  auto match_link(const key_type &key, uint32_t begin, size_t link, size_t end) const -> uint32_t {
    auto enu = strpool_.get_string_enumerator(link, end);
    uint32_t matched_len = 0;
    uint8_t label;
    while ((label = enu.next()) != terminator_) {
      if (begin + matched_len >= key.size() || key[begin + matched_len] != label) {
        return -1;
      }
      matched_len++;
    }
    return matched_len;
  }

  auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool override {
    assert(key_id < size());
    auto [link, end] = links_.select_range(key_id);
    auto enu = strpool_.get_string_enumerator(link, end);
    uint32_t i = 0;
    uint8_t label;
    while ((label = enu.next()) != terminator_) {
      if (begin + i >= key.size()) return true;
      if ((uint8_t)key[begin + i] != label) return false;
      i++;
    }
    return begin + i == key.size();
  }

  auto size() const -> uint32_t override {
    return strpool_.size();
  }

  auto size_in_bytes() const -> size_t override {
    return succinct::mapper::size_of(const_cast<strpool_t &>(strpool_)) +
           succinct::mapper::size_of(const_cast<succinct::elias_fano &>(links_));
  }

  auto size_in_bits() const -> size_t override {
    return size_in_bytes() * 8;
  }

  void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const {
    link += succinct::mapper::size_of(const_cast<succinct::elias_fano &>(links_)) * 8;
    data += succinct::mapper::size_of(const_cast<strpool_t &>(strpool_)) * 8;
  }

  static auto estimate_space_cost(const KeySet<key_type> &keys, bool partial_links = false) -> size_t {
    auto sample = keys.make_sample(StringPool<key_type>::sample_size_, false);
    RepairStringPool<key_type> temp;
    if (partial_links) {
      std::vector<uint8_t> links;
      temp.build(sample, &links);
    } else {
      temp.build(sample);
    }
    return keys.space_cost() * 8 * temp.size_in_bits() / (sample.space_cost() * 8);
  }
 private:
  strpool_t strpool_;
  succinct::elias_fano links_;
};

template <typename Key>
class FsstStringPool : public StringPool<Key> {
 public:
  using key_type = Key;

  FsstStringPool() = default;

  ~FsstStringPool() = default;

  void build(const KeySet<key_type> &sorted_rev_keys, std::vector<uint8_t> *partial_links = nullptr,
             int max_recursion = 0, int mask = 0) override {
    links_.resize(sorted_rev_keys.size());
    if (partial_links != nullptr) {
      partial_links->resize(sorted_rev_keys.size());
    }

    std::vector<const uint8_t *> in_str;
    std::vector<uint8_t *> out_str;
    std::vector<size_t> in_len, out_len;
    std::vector<uint8_t> in_buf;  // zero terminated
    std::vector<uint32_t> mapping(sorted_rev_keys.size());  // map key IDs to positions in `in_str` and `out_str`
    const typename KeySet<key_type>::Fragment *next = nullptr;
    in_str.reserve(sorted_rev_keys.size());
    in_len.reserve(sorted_rev_keys.size());
    in_buf.reserve(sorted_rev_keys.space_cost() + sorted_rev_keys.size());  // key + terminator
    for (size_t i = sorted_rev_keys.size(); i > 0; i--) {
      auto &cur = sorted_rev_keys[i - 1];
      assert(cur.id_ < mapping.size());
      if (next != nullptr && cur == *next) {  // duplicate; overlap previous key
        mapping[cur.id_] = in_str.size() - 1;
        continue;
      }
      // new key
      in_str.push_back(in_buf.data() + in_buf.size());
      in_len.push_back(cur.size() + 1);
      cur.append_to(in_buf, true);
      mapping[cur.id_] = in_str.size() - 1;
      next = &cur;
    }

    auto encoder = fsst_create(in_str.size(), in_len.data(), in_str.data(), true);
    codes_.resize(in_buf.size() * 2 + 16);
    out_str.resize(in_str.size());
    out_len.resize(in_str.size());
    fsst_compress(encoder, in_str.size(), in_len.data(), in_str.data(), codes_.size(), codes_.data(), out_len.data(), out_str.data());
    size_t compressed_len = (in_str.size() == 0) ? 0 : (out_str[in_str.size() - 1]  + out_len[in_str.size() - 1] - codes_.data());
    codes_.resize(compressed_len + 8);

    for (uint32_t i = 0; i < sorted_rev_keys.size(); i++) {
      auto id = sorted_rev_keys[i].id_;
      size_t link = out_str[mapping[id]] - codes_.data();
      assert(link < compressed_len);
      if (partial_links == nullptr) {
        links_[id] = link;
      } else {
        links_[id] = (link >> 8);
        (*partial_links)[id] = (link & MASK(8));
      }
    }
    sdsl::util::bit_compress(links_);
    codes_.shrink_to_fit();

    uint8_t buf[sizeof(fsst_decoder_t)];
    fsst_export(encoder, buf);
    fsst_destroy(encoder);
    fsst_import(&decoder_, buf);
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id) const -> uint32_t override {
    assert(key_id < size());
    auto link = links_[key_id];
    return match_link(key, begin, link);
  }

  auto match(const key_type &key, uint32_t begin, uint32_t key_id, uint8_t partial_link) const -> uint32_t override {
    assert(key_id < size());
    auto link = (links_[key_id] << 8) | partial_link;
    return match_link(key, begin, link);
  }

  auto starts_with(const key_type &key, uint32_t begin, uint32_t key_id) const -> bool override {
    assert(key_id < size());
    #define FSSTSTRPOOL_BATCH_SIZE (32)
    #define FSSTSTRPOOL_BUFFER_SIZE (FSSTSTRPOOL_BATCH_SIZE * 8)
    uint8_t out_buf[FSSTSTRPOOL_BUFFER_SIZE];
    size_t link = links_[key_id];
    size_t pool_end = link + strlen(reinterpret_cast<const char *>(codes_.data() + link));
    uint32_t pos = begin;
    while (link < pool_end) {
      if (pos >= (uint32_t)key.size()) return true;
      size_t decode_size = std::min((size_t)FSSTSTRPOOL_BATCH_SIZE, pool_end - link);
      if (decode_size > 1 && codes_[link + decode_size - 1] == FSST_ESC) decode_size--;
      size_t out_size = fsst_decompress(&decoder_, decode_size, codes_.data() + link, FSSTSTRPOOL_BUFFER_SIZE, out_buf);
      uint32_t remaining = (uint32_t)key.size() - pos;
      if (out_size >= remaining) {
        return memcmp(key.c_str() + pos, out_buf, remaining) == 0;
      }
      if (pos > key.size() || out_size > key.size() - pos || memcmp(key.c_str() + pos, out_buf, out_size)) return false;
      pos += (uint32_t)out_size;
      link += decode_size;
    }
    return pos == (uint32_t)key.size();
    #undef FSSTSTRPOOL_BATCH_SIZE
    #undef FSSTSTRPOOL_BUFFER_SIZE
  }

  auto match_link(const key_type &key, uint32_t begin, size_t link) const -> uint32_t {
  #define FSSTSTRPOOL_BATCH_SIZE (32)  // decode this many input bytes at a time
  #define FSSTSTRPOOL_BUFFER_SIZE (FSSTSTRPOOL_BATCH_SIZE * 8)  // each input byte decodes to at most 8 bytes
    uint8_t out_buf[FSSTSTRPOOL_BUFFER_SIZE];
    size_t out_size, end = link + strlen(reinterpret_cast<const char *>(codes_.data() + link));  // zero terminated
    uint32_t pos = begin;
    while (end - link > FSSTSTRPOOL_BATCH_SIZE) {
      // make sure FSST_ESC is not the last code
      uint32_t decode_size = FSSTSTRPOOL_BATCH_SIZE - (codes_[link + FSSTSTRPOOL_BATCH_SIZE - 1] == FSST_ESC);
      out_size = fsst_decompress(&decoder_, decode_size, codes_.data() + link, FSSTSTRPOOL_BUFFER_SIZE, out_buf);
      if (pos > key.size() || out_size > key.size() - pos || memcmp(key.c_str() + pos, out_buf, out_size)) {
        return -1;
      }
      link += decode_size;
      pos += out_size;
    }
    if (end > link) {
      out_size = fsst_decompress(&decoder_, end - link, codes_.data() + link, FSSTSTRPOOL_BUFFER_SIZE, out_buf);
      if (pos > key.size() || out_size > key.size() - pos || memcmp(key.c_str() + pos, out_buf, out_size)) {
        return -1;
      }
      pos += out_size;
    }
    return pos - begin;
  #undef FSSTSTRPOOL_BATCH_SIZE
  #undef FSSTSTRPOOL_BUFFER_SIZE
  }

  auto size() const -> uint32_t override {
    return links_.size();
  }

  auto size_in_bytes() const -> size_t override {
    return sdsl::size_in_bytes(links_) + codes_.size() * sizeof(uint8_t) + sizeof(fsst_decoder_t);
  }

  auto size_in_bits() const -> size_t override {
    return size_in_bytes() * 8;
  }

  void space_cost_breakdown(size_t &topo, size_t &link, size_t &data) const override {
    link += (sdsl::size_in_bytes(links_) + sizeof(fsst_decoder_t)) * 8;
    data += codes_.size() * sizeof(uint8_t) * 8;
  }

  static auto estimate_space_cost(const KeySet<key_type> &keys, bool partial_links = false) {
    if (keys.space_cost() == 0) return size_t(0);
    auto sample = keys.make_sample(StringPool<key_type>::sample_size_, true);
    FsstStringPool<key_type> temp;
    if (partial_links) {
      std::vector<uint8_t> links;
      temp.build(sample, &links);
    } else {
      temp.build(sample);
    }
    return keys.space_cost() * 8 * temp.size_in_bits() / (sample.space_cost() * 8);
  }
  void serialize(std::ostream& out) const override {
    archive::write(out,uint8_t(2));links_.serialize(out);archive::vector(out,codes_);
    archive::write(out,decoder_.version);archive::write(out,decoder_.zeroTerminated);
    archive::write(out,decoder_.len);archive::write(out,decoder_.symbol);
  }
  void load_data(std::istream& in) {
    links_.load(in);archive::vector(in,codes_);
    archive::read(in,decoder_.version);archive::read(in,decoder_.zeroTerminated);
    archive::read(in,decoder_.len);archive::read(in,decoder_.symbol);
  }
 private:
  sdsl::int_vector<> links_;
  std::vector<uint8_t> codes_;
  fsst_decoder_t decoder_{};
};

}  // namespace coco_paper
