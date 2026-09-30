// AfAME - parallel-tempering search and algebraic certification of
// absolutely maximally entangled (AME) states.
// Copyright (C) 2026 Zakaria Dahbi
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
using namespace std;

// Independent rank calculation: enumerate the image of every row combination.
int cut_rank(const array<uint8_t,8>& graph, int mask) {
    uint8_t row[4]{};
    int k = 0;
    for (int v = 0; v < 8; ++v)
        if (mask & (1 << v)) row[k++] = graph[v] & (255 ^ mask);
    bool seen[256]{};
    int distinct = 0;
    for (int combination = 0; combination < (1 << k); ++combination) {
        int value = 0;
        for (int j = 0; j < k; ++j)
            if (combination & (1 << j)) value ^= row[j];
        if (!seen[value]) { seen[value] = true; ++distinct; }
    }
    int rank = 0;
    while ((1 << rank) < distinct) ++rank;
    return rank;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    ifstream in(argv[1]);
    array<uint8_t,8> graph{};
    int min_cost = 1000000, seven = 0, cases = 0, winners = 0;
    while (true) {
        int value;
        if (!(in >> value)) break;
        graph[0] = value;
        for (int v = 1; v < 7; ++v) { in >> value; graph[v] = value; }
        ++seven;
        for (int neighbors = 0; neighbors < 128; ++neighbors) {
            graph[7] = neighbors;
            for (int v = 0; v < 7; ++v) {
                graph[v] &= 127;
                if (neighbors & (1 << v)) graph[v] |= 128;
            }
            int cost = 0;
            for (int mask = 1; mask < 256; ++mask) {
                int k = __builtin_popcount((unsigned)mask);
                if (k > 4) continue;
                int deficiency = k - cut_rank(graph, mask);
                cost += deficiency * deficiency;
            }
            if (cost < min_cost) { min_cost = cost; winners = 1; }
            else if (cost == min_cost) ++winners;
            ++cases;
        }
    }
    cout << "atlas7=" << seven << " extensions=" << cases
         << " min_cost=" << min_cost << " winning_extensions=" << winners << '\n';
}
