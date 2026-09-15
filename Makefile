CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON ?= python3

.PHONY: all test sanitize clean
all: AME7DIA_PT
AME7DIA_PT: AME7DIA_PT.cpp
	$(CXX) $(CXXFLAGS) -pthread $< -o $@
tests/test_core: tests/test_core.cpp AME7DIA_PT.cpp
	$(CXX) $(CXXFLAGS) -pthread $< -o $@
test: AME7DIA_PT tests/test_core
	./tests/test_core
	$(PYTHON) tests/test_integration.py ./AME7DIA_PT
sanitize:
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -pthread tests/test_core.cpp -o tests/test_core_sanitized
	./tests/test_core_sanitized
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -pthread AME7DIA_PT.cpp -o tests/AME7DIA_PT_sanitized
	$(PYTHON) tests/test_integration.py ./tests/AME7DIA_PT_sanitized
clean:
	rm -f AME7DIA_PT tests/test_core tests/test_core_sanitized tests/AME7DIA_PT_sanitized
