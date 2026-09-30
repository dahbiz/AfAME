#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
if [[ ! -x "$package_dir/AfAME" || "$package_dir/AfAME.cpp" -nt "$package_dir/AfAME" ]]; then
    c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -pthread \
        "$package_dir/AfAME.cpp" -o "$package_dir/AfAME"
fi
exec "$package_dir/AfAME" "$@"
