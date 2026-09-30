#!/usr/bin/env bash
# AfAME - parallel-tempering search and algebraic certification of
# absolutely maximally entangled (AME) states.
# Copyright (C) 2026 Zakaria Dahbi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
if [[ ! -x "$package_dir/AfAME" || "$package_dir/AfAME.cpp" -nt "$package_dir/AfAME" ]]; then
    c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -pthread \
        "$package_dir/AfAME.cpp" -o "$package_dir/AfAME"
fi
exec "$package_dir/AfAME" "$@"
