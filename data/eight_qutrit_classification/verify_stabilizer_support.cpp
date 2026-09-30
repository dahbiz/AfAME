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
// Independent exact check: enumerate stabilizer supports, without elimination.
// For x in GF(3)^8, a graph stabilizer has support {j : (x_j,(P*x)_j) != (0,0)}.
// The number supported within S is 3^(|S|-rank(P[S,Sbar])).
// This program parses the original extracted database independently of Python.
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> out;
    std::stringstream stream(line);
    std::string part;
    while (std::getline(stream, part, delimiter)) out.push_back(part);
    return out;
}
int bits(unsigned int mask) {
    int count = 0;
    while (mask) { count += mask & 1U; mask >>= 1U; }
    return count;
}
int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: verify_by_stabilizer database.tsv output.csv");
        std::ifstream input(argv[1]);
        std::ofstream output(argv[2]);
        require(bool(input) && bool(output), "Cannot open input/output");
        output << "class,mask,size,rank,deficit\n";
        std::string line;
        int index = 0, minimum = 99999, indecomposable = 0;
        std::vector<int> minimizers;
        std::map<int, int> histogram;
        while (std::getline(input, line)) {
            ++index;
            const auto parts = split(line, '\t');
            require(parts.size() == 6 && parts[0] == "8", "Malformed database row");
            require(parts[1] == "I" || parts[1] == "D", "Invalid code type");
            indecomposable += parts[1] == "I";
            std::array<std::array<int, 8>, 8> matrix{};
            if (parts[5] != "-") for (const auto& edge : split(parts[5], ',')) {
                const auto star = edge.find('*');
                const int weight = star == std::string::npos ? 1 : std::stoi(edge.substr(0, star));
                const auto vertices = split(edge.substr(star == std::string::npos ? 0 : star + 1), '-');
                require(vertices.size() == 2, "Malformed edge");
                const int i = std::stoi(vertices[0]), j = std::stoi(vertices[1]);
                require(0 <= i && i < j && j < 8 && 1 <= weight && weight <= 2, "Invalid edge");
                require(matrix[i][j] == 0, "Duplicate edge");
                matrix[i][j] = matrix[j][i] = weight;
            }
            std::array<int, 256> support_counts{};
            std::array<int, 9> weights{};
            for (int label = 0; label < 6561; ++label) {
                std::array<int, 8> x{};
                int quotient = label;
                for (auto& value : x) { value = quotient % 3; quotient /= 3; }
                unsigned int support = 0;
                for (int i = 0; i < 8; ++i) {
                    int y = 0;
                    for (int j = 0; j < 8; ++j) y += matrix[i][j] * x[j];
                    if (x[i] || y % 3) support |= 1U << i;
                }
                ++support_counts[support];
                ++weights[bits(support)];
            }
            const auto expected = split(parts[3].substr(1, parts[3].size()-2), ',');
            require(expected.size() == weights.size(), "Weight enumerator length");
            for (int w = 0; w <= 8; ++w)
                require(weights[w] == std::stoi(expected[w]), "Published weight enumerator mismatch at class " + std::to_string(index));
            int distance = 1;
            while (distance <= 8 && weights[distance] == 0) ++distance;
            require(distance == std::stoi(parts[2]), "Distance mismatch");
            // Subset zeta transform: counts of stabilizers contained in each S.
            for (int bit = 0; bit < 8; ++bit)
                for (int mask = 0; mask < 256; ++mask)
                    if (mask & (1 << bit)) support_counts[mask] += support_counts[mask ^ (1 << bit)];
            int cost = 0, cuts = 0;
            for (int mask = 1; mask < 256; ++mask) {
                const int k = bits(mask);
                if (k > 4) continue;
                int value = support_counts[mask], deficit = 0;
                require(value > 0, "Identity stabilizer missing");
                while (value > 1 && value % 3 == 0) { value /= 3; ++deficit; }
                require(value == 1 && deficit <= k, "Support count is not a valid power of 3");
                cost += deficit * deficit;
                ++cuts;
                output << index << ',' << mask << ',' << k << ',' << k-deficit << ',' << deficit << '\n';
            }
            require(cuts == 162, "Cut count mismatch");
            ++histogram[cost];
            if (cost < minimum) { minimum = cost; minimizers.clear(); }
            if (cost == minimum) minimizers.push_back(index);
        }
        require(index == 817 && indecomposable == 659, "Incomplete database");
        std::cout << "classes=" << index << " stabilizers_per_class=6561 cut_checks=" << index*162
                  << " minimum=" << minimum << "\nminimizing_classes=";
        for (int id : minimizers) std::cout << id << ' ';
        std::cout << "\nAll published weight enumerators and distances matched.\n";
        require(bool(output), "Output failure");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
