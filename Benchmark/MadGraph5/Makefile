PYTHON ?= python3
OPTIONS ?=

.PHONY: all benchmark smoke clean

all: benchmark

benchmark:
	$(PYTHON) ../run_madgraph_benchmark.py --skip-reference-check $(OPTIONS)

smoke:
	$(PYTHON) ../run_madgraph_benchmark.py --min-gluons 4 --max-gluons 6 \
		--target-seconds 0.01 --batches 1 --skip-reference-check \
		--output ../build/results-madgraph-smoke.md $(OPTIONS)

clean:
	$(RM) -r ../build/madgraph ../build/results-madgraph-fixed.md \
		../build/results-madgraph-smoke.md
