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
