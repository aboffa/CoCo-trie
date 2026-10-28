// This file is part of CoCo-trie, licensed under GPL-3.0.
#pragma once

#include <array>
#include <deque>
#include "CoCo-trie_v1.hpp"

struct coco_sorted_unique_t {};
inline constexpr coco_sorted_unique_t coco_sorted_unique{};

// Build from prefix intervals in a sorted dictionary. No character-node tree
// is materialized: only the active LCP stack and two-byte decisions are kept.
template<uint8_t MIN_L = 1, typename Code = uint128_t,
         uint8_t MaxLevels = MAX_L_THRS, uint8_t Relaxation = 0, class Bitmap = Code>
class sorted_trie_builder {
    static_assert(MIN_L == 1, "The streaming builder currently supports MIN_L=1");
    static_assert(MaxLevels > 0 && MaxLevels <= MAX_L_THRS_amap);
    static_assert(BASE_COST == 0, "Frontier summaries assume zero-cost leaves");
    using alphamap = basic_alphamap<Bitmap>;
    using Legacy = Trie_lw<MIN_L, Code, MaxLevels, Relaxation, Bitmap>;
    using Encoders = CoCo_v1<MIN_L, Code, MaxLevels, Relaxation, sdsl::rank_support_v<1>, sdsl::rank_support_v<0, 2>, sux::bits::SimpleSelectZero<>, Bitmap>;
    static constexpr size_t levels = MAX_L_THRS_amap;

    struct Frame {
        size_t first_key = 0, decision = 0;
        size_t height = 0, max_child_levels = 0, safe_levels = 0, best = 0;
        bool terminal = false;
        std::array<size_t, levels> internal{}, leaves{}, costs{};
        std::array<Bitmap, levels + 1> maps{};
    };
    const std::vector<std::string>& keys_;
    std::vector<size_t> offsets_;
    std::vector<uint16_t> decisions_;
    std::array<size_t, 258> safe_limits_{};

    size_t common_prefix(size_t i) const {
        return keys_[i].size() - (offsets_[i + 1] - offsets_[i]);
    }

    size_t safe_limit(size_t radix) {
        if (safe_limits_[radix]) return safe_limits_[radix];
        Code maximum = ~Code(0);
#ifdef COCO_BITINT_CAP_128
        if constexpr (sizeof(Code) > 16) maximum = (Code(1) << 128) - 1;
#endif
        // Preserve the existing builder's conservative overflow rule.
        size_t limit = levels;
        for (size_t length = 1; length <= levels; ++length) {
            Code acc = radix;
            bool valid = true;
            for (size_t i = 1; i <= length; ++i) {
                Code next = acc * Code(radix);
                if (next >= maximum / Code(radix)) { valid = false; break; }
                acc = next;
            }
            if (!valid) { limit = length - 1; break; }
        }
        safe_limits_[radix] = limit;
        return limit;
    }

    void prepare() {
        offsets_.resize(keys_.size() + 1);
        offsets_[0] = 1; // Decision zero belongs to the root.
        for (size_t i = 0; i < keys_.size(); ++i) {
            size_t common = 0;
            if (i) {
                if (!(keys_[i - 1] < keys_[i]))
                    throw std::invalid_argument("sorted builder requires strictly increasing, deduplicated keys");
                while (common < std::min(keys_[i - 1].size(), keys_[i].size()) &&
                       keys_[i - 1][common] == keys_[i][common]) ++common;
            }
            const size_t added = keys_[i].size() - common;
            if (added > std::numeric_limits<size_t>::max() - offsets_[i])
                throw std::length_error("too many prefixes");
            offsets_[i + 1] = offsets_[i] + added;
        }
        decisions_.resize(offsets_.back());
        prefix_count = decisions_.size();
    }

    void optimize(size_t fixed_length) {
        std::deque<Frame> stack(1);
        size_t depth = 0;
        // One reusable candidate workspace, rather than candidate vectors on
        // every character node. Scoring is shared with the original builder.
        typename Legacy::TrieNode_lw scratch;
        scratch.alphamaps.resize(levels + 1);
        scratch.u_vec.resize(levels); scratch.u_vec_real.resize(levels);
        scratch.n_vec.resize(levels); scratch.node_type_vec.resize(levels);
        std::vector<size_t> counts(levels), costs(levels), scores(levels);
        std::vector<Code> low(levels), high(levels), low_real(levels), high_real(levels);

        auto finish = [&](Frame& f, size_t d, size_t end) {
            if (!f.height) return;
            size_t limit = std::min(levels, f.max_child_levels + 1);
            for (size_t i = 1; i <= limit; ++i) {
                f.maps[i] |= f.maps[i - 1];
                const size_t radix = alphamap(f.maps[i]).rankmax() + 1;
                if (i > safe_limit(radix)) { limit = i - 1; break; }
            }
            f.safe_levels = limit;
            const size_t global_limit = std::min({limit, size_t(MaxLevels), sizeof(Code) * CHAR_BIT / log_sigma});
            const size_t first_child_key = f.first_key + size_t(f.terminal);
            const auto left = std::string_view(keys_[first_child_key]).substr(d);
            const auto right = std::string_view(keys_[end - 1]).substr(d);
            size_t ended = 0;
            for (size_t i = 0; i < limit; ++i) {
                ended += f.leaves[i];
                counts[i] = ended + f.internal[i];
                costs[i] = f.costs[i];
                scores[i] = size_t(-1);
                scratch.node_type_vec[i] = index_types;
                scratch.alphamaps[i] = alphamap(f.maps[i]);
                if (i < global_limit) {
                    const Code a = i < left.size() ? Code(alphabet.offset(left[i]) + 1) : Code(0);
                    const Code b = i < right.size() ? Code(alphabet.offset(right[i]) + 1) : Code(0);
                    low[i] = (i ? low[i - 1] * Code(alphabet.size) : Code(0)) + a;
                    high[i] = (i ? high[i - 1] * Code(alphabet.size) : Code(0)) + b;
                }
                const auto& map = scratch.alphamaps[i];
                if (i && f.maps[i] == f.maps[i - 1]) {
                    const Code radix = map.rankmax() + 1;
                    low_real[i] = low_real[i - 1] * radix +
                        (i < left.size() ? Code(map.rank(alphabet.offset(left[i])) + 1) : Code(0));
                    high_real[i] = high_real[i - 1] * radix +
                        (i < right.size() ? Code(map.rank(alphabet.offset(right[i])) + 1) : Code(0));
                } else {
                    low_real[i] = Legacy::enc_real(left, i, map, alphabet);
                    high_real[i] = Legacy::enc_real(right, i, map, alphabet);
                }
            }
            Legacy::evaluate_space_cost(&scratch, global_limit, limit, counts, costs,
                                        low, high, low_real, high_real, scores, alphabet, log_sigma);
            size_t selected = 0;
            // The legacy <= comparator selects the last equal-cost candidate.
            for (size_t i = 1; i < limit; ++i) if (scores[i] <= scores[selected]) selected = i;
            if (fixed_length) selected = std::min(fixed_length - 1, limit - 1);
            else if constexpr (Relaxation > 0) {
                const auto best = scores[selected];
                for (size_t i = selected + 1; i < limit; ++i)
                    if ((scores[i] - best) / double(best) * 100.0 < Relaxation) selected = i;
            }
            f.best = scores[selected];
            decisions_[f.decision] = uint16_t(selected | (size_t(scratch.node_type_vec[selected]) << 6));
            max_l_idx = std::max(max_l_idx, selected);
        };

        auto close = [&](size_t end) {
            Frame& child = stack[depth];
            finish(child, depth, end);
            Frame& parent = stack[depth - 1];
            const auto ch = keys_[child.first_key][depth - 1];
            parent.maps[0] |= Bitmap(1) << (alphabet.offset(ch));
            parent.height = std::max(parent.height, std::min(levels, child.height + 1));
            parent.max_child_levels = std::max(parent.max_child_levels, child.safe_levels);
            parent.costs[0] += child.best;
            if (!child.height) ++parent.leaves[0];
            else {
                ++parent.internal[0];
                if (child.terminal) ++parent.leaves[1];
                for (size_t i = 0; i < std::min(child.height, levels - 1); ++i) {
                    parent.internal[i + 1] += child.internal[i];
                    parent.leaves[i + 1] += child.leaves[i];
                    parent.costs[i + 1] += child.costs[i];
                }
                // Accumulate map changes; finish() makes them cumulative.
                for (size_t i = 0; i <= std::min(child.safe_levels, levels - 1); ++i)
                    parent.maps[i + 1] |= child.maps[i];
            }
            --depth;
        };

        for (size_t i = 0; i < keys_.size(); ++i) {
            const size_t common = common_prefix(i);
            while (depth > common) close(i);
            while (depth < keys_[i].size()) {
                ++depth;
                if (depth == stack.size()) stack.emplace_back();
                else stack[depth] = Frame{};
                stack[depth].first_key = i;
                stack[depth].decision = offsets_[i] + depth - common - 1;
            }
            stack[depth].terminal = true;
        }
        while (depth) close(keys_.size());
        finish(stack[0], 0, keys_.size());
        root_cost = stack[0].best;
        stack_prefixes = stack.size();
    }

    void emit() {
        struct Interval { size_t lo, hi, depth; };
        std::deque<Interval> queue{{0, keys_.size(), 0}};
        std::vector<bool> topology{true, false};
        std::vector<size_t> pointers;
        std::vector<Code> codes;
        succinct::bit_vector_builder encoding;
        bits_for_L = log_universe(std::max<uint64_t>(1, max_l_idx));
        while (!queue.empty()) {
            const auto node = queue.front(); queue.pop_front();
            const bool terminal = keys_[node.lo].size() == node.depth;
            const size_t first = node.lo + size_t(terminal);
            if (first == node.hi) { topology.push_back(false); continue; }
            const size_t decision = node.depth ? offsets_[node.lo] + node.depth - common_prefix(node.lo) - 1 : 0;
            const size_t index = decisions_[decision] & 63;
            const auto type = node_type(decisions_[decision] >> 6);
            const bool remapped = type >= elias_fano_amap;
            alphamap map(0);
            if (remapped) {
                for (size_t i = first; i < node.hi; ++i)
                    for (char ch : std::string_view(keys_[i]).substr(node.depth, index + 1))
                        map.template set_bit<1>(alphabet.offset(ch));
            }
            codes.clear();
            for (size_t begin = first; begin < node.hi;) {
                const auto label = std::string_view(keys_[begin]).substr(node.depth, index + 1);
                size_t end = begin + 1;
                while (end < node.hi && std::string_view(keys_[end]).substr(node.depth, index + 1) == label) ++end;
                codes.push_back(remapped ? Legacy::enc_real(label, index, map, alphabet) : Legacy::enc(label, index, alphabet));
                queue.push_back({begin, end, node.depth + label.size()});
                topology.push_back(true);
                begin = end;
            }
            topology.push_back(false);
            if (!node.depth) root_children = codes.size();
            pointers.push_back(encoding.size());
            encoding.append_bits(uint64_t(index), bits_for_L);
            encoding.append_bits(uint64_t(type), NUM_BIT_TYPE);
            encoding.append_bits(uint64_t(terminal), 1);
            if (remapped) encoding.append_bits(map.bitmap, alphabet.bitmap_bits(alphamap::width));
            encoding.append_bits(codes.front(), remapped ? bits_first_code<MIN_L, Code>(map.rankmax(), index) : log_sigma * (index + 1));
            const Code first_code = codes.front(), universe = codes.back() - first_code;
            for (size_t i = 1; i < codes.size(); ++i) codes[i - 1] = codes[i] - first_code - 1;
            codes.pop_back();
            switch (type) {
                case elias_fano: case elias_fano_amap: Encoders::write_elias_fano(universe, codes, encoding); break;
                case bitvector: case bitvector_amap: Encoders::write_bitvector(universe, codes, encoding); break;
                case packed: case packed_amap: Encoders::write_packed(universe, codes, encoding); break;
                case all_ones: case all_ones_amap: break;
                default: throw std::logic_error("invalid sorted-builder encoding decision");
            }
        }
        topology_bits = sdsl::bit_vector(topology.size(), 0);
        for (size_t i = 0; i < topology.size(); ++i) topology_bits[i] = topology[i];
        encoding_pointers = std::make_unique<sdsl::int_vector<>>(pointers.size());
        for (size_t i = 0; i < pointers.size(); ++i) (*encoding_pointers)[i] = pointers[i];
        sdsl::util::bit_compress(*encoding_pointers);
        labels = std::make_unique<succinct::bit_vector>(&encoding);
    }

public:
    alphabet_config alphabet;
    bool root_terminal = false;
    size_t log_sigma, max_l_idx = 0, bits_for_L = 1, root_children = 0;
    size_t prefix_count = 0, root_cost = 0, stack_prefixes = 0;
    sdsl::bit_vector topology_bits;
    std::unique_ptr<sdsl::int_vector<>> encoding_pointers;
    std::unique_ptr<succinct::bit_vector> labels;

    explicit sorted_trie_builder(const std::vector<std::string>& keys, size_t fixed_length = 0)
        : keys_(keys), alphabet(alphabet_config::from_dataset(keys, alphamap::width)), log_sigma(log_universe(alphabet.size)) {
        prepare();
        root_terminal = !keys.empty() && keys.front().empty();
        if (!keys.empty() && !keys.back().empty()) {
            optimize(fixed_length);
            emit();
        }
        std::vector<uint16_t>().swap(decisions_);
        std::vector<size_t>().swap(offsets_);
    }
};
