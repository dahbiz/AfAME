#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
if [[ ! -x "$package_dir/AME7DIA_PT" || "$package_dir/AME7DIA_PT.cpp" -nt "$package_dir/AME7DIA_PT" ]]; then
    c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -pthread \
        "$package_dir/AME7DIA_PT.cpp" -o "$package_dir/AME7DIA_PT"
fi
exec "$package_dir/AME7DIA_PT" "$@"
