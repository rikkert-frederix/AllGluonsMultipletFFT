LATEXMK ?= latexmk

.PHONY: all clean

all: gluon_colour_sums.pdf

gluon_colour_sums.pdf: gluon_colour_sums.tex references.bib data/colour_sum_times.dat tables/evaluation_results.tex tables/memory_results.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error gluon_colour_sums.tex

clean:
	$(LATEXMK) -C gluon_colour_sums.tex
