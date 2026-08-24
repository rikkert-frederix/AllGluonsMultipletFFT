PYTHON = python3
BENCHMARK_OPTIONS =
NATIVE_FFLAGS = -Ofast -march=native -std=f2018 -Wall -Wextra

.PHONY: all benchmark benchmark-native clean test

all: benchmark

benchmark:
	$(PYTHON) run_benchmark.py $(BENCHMARK_OPTIONS) --output build/results.md

benchmark-native:
	$(PYTHON) run_benchmark.py $(BENCHMARK_OPTIONS) \
		--fflags '$(NATIVE_FFLAGS)' --build-dir build-native \
		--output build-native/results.md

test:
	$(PYTHON) -m py_compile run_benchmark.py
	$(PYTHON) -m unittest discover -s tests -v
	$(PYTHON) run_benchmark.py --min-gluons 4 --max-gluons 6 \
		--points 1 --target-seconds 0.001 --batches 1 \
		--initialization-runs 1 --output build/test-results.md > /dev/null
	@echo 'ten-variant proxy-weighted sampled-helicity benchmark regression: PASS'

clean:
	rm -rf build build-debug build-native __pycache__
