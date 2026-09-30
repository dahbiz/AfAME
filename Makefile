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
CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON ?= python3

.PHONY: all test sanitize clean
all: AfAME
AfAME: AfAME.cpp
	$(CXX) $(CXXFLAGS) -pthread $< -o $@
tests/test_core: tests/test_core.cpp AfAME.cpp
	$(CXX) $(CXXFLAGS) -pthread $< -o $@
test: AfAME tests/test_core
	./tests/test_core
	$(PYTHON) tests/test_integration.py ./AfAME
sanitize:
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -pthread tests/test_core.cpp -o tests/test_core_sanitized
	./tests/test_core_sanitized
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -pthread AfAME.cpp -o tests/AfAME_sanitized
	$(PYTHON) tests/test_integration.py ./tests/AfAME_sanitized
clean:
	rm -f AfAME tests/test_core tests/test_core_sanitized tests/AfAME_sanitized
