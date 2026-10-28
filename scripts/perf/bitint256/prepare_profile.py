#!/usr/bin/env python3
"""Add an untimed lookup trace to an already prepared experimental overlay."""
import argparse
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('overlay', type=Path)
a = p.parse_args()
path = a.overlay / 'include/CoCo-trie_v2.hpp'
s = path.read_text()
start = s.index('    size_t look_up(const std::string &to_search) const {')
opening = s.index('{', start)
end, depth = opening + 1, 1
while depth:
    depth += (s[end] == '{') - (s[end] == '}')
    end += 1
body = s[start:end]
body = body.replace('look_up(const std::string &to_search)', 'look_up_profiled(const std::string &to_search, lookup_profile& profile)')
body = body.replace('            bool is_end_of_world = it.take(1);', '''            bool is_end_of_world = it.take(1);
            ++profile.visited_nodes;
            ++profile.node_types[nt];''')
body = body.replace('            if (!alphabet_.contains(substr))', '''            profile.characters += substr.size();
            profile.jump_capacity += l + MIN_L;
            if (!alphabet_.contains(substr))''')
body = body.replace('                to_search_code = utrie_t::enc_real(substr, l, am, alphabet_);', '''                to_search_code = utrie_t::enc_real(substr, l, am, alphabet_);
                profile.work.push_back({substr, l, am.bitmap, true, to_search_code});''')
body = body.replace('                to_search_code = utrie_t::enc(substr, l, alphabet_);', '''                to_search_code = utrie_t::enc(substr, l, alphabet_);
                profile.work.push_back({substr, l, 0, false, to_search_code});''')
addition = '''
    struct lookup_work {
        std::string_view text;
        size_t level_index;
        uint128_t bitmap;
        bool remapped;
        code_type expected;
    };
    struct lookup_profile {
        size_t visited_nodes = 0, characters = 0, jump_capacity = 0;
        size_t node_types[8] = {};
        std::vector<lookup_work> work;
    };
    code_type encode_work(const lookup_work& work) const {
        if (work.remapped) {
            alphamap am(work.bitmap);
            return utrie_t::enc_real(work.text, work.level_index, am, alphabet_);
        }
        return utrie_t::enc(work.text, work.level_index, alphabet_);
    }
'''
s = s[:end] + addition + '\n' + body + '\n' + s[end:]
path.write_text(s)
