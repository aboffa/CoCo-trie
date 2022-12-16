# This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
# Copyright (c) 2022 Antonio Boffa.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
    #
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


import sys

def return_only_ascii(str):
    return ''.join([x for x in str if x.isascii()])

def lcp(s1, s2):
    n = min(len(s1), len(s2))
    for i in range(n):
        if s1[i] != s2[i]:
            return i
    return n

def distinguishing_prefix(l):
    assert len(l) >= 2

    def trunc(s1, s2, s3):
        m = max(lcp(s1, s2), lcp(s2, s3))
        m = min(len(s2), m + 1)
        return s2[:m]

    result = [trunc('', l[0], l[1])]
    for i in range(1, len(l) - 1):
        result.append(trunc(l[i-1], l[i], l[i+1]))

    result.append(trunc(l[-2], l[-1], ''))
    result.sort()
    return result


if __name__ == "__main__":
    for arg in sys.argv:
        with open(arg) as f:
            original_lines = [x.strip() for x in f.readlines()]

        for i, s in enumerate(original_lines):
            original_lines[i] = return_only_ascii(s)

        original_lines.sort()
        original_lines = list(set(original_lines))
        original_lines.sort()
        dpa = distinguishing_prefix(original_lines)
        dpa = list(set(dpa))
        dpa.sort()

        sum = 0
        for i, s in enumerate(dpa):
            sum += len(s)
            if sum > 10**10:
                dpa = dpa[:i]

        with open(arg + '_distinguishing_new', mode='w') as f:
            f.write('\n'.join(dpa))