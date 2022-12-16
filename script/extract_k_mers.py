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


K = 31
with open('dna', 'r') as file:
    data = file.read().replace('\n', '')
    res = set(data[i:i+K] for i in range(len(data)-K+1))
    res = sorted(list(res))

    with open('dna-'+str(K)+'-mer.txt', mode='w') as f:
        f.write('\n'.join(res))
