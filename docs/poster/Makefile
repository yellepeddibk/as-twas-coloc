# Makefile for ISU Honors Poster
# Usage: make poster

TEX = poster.tex
PDF = poster.pdf
LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error

.PHONY: poster clean figures

poster: $(PDF)

$(PDF): $(TEX) beamerthemeISU.sty figures
	$(LATEXMK) $(TEX)

figures:
	cd ../.. && python scripts/build_poster_figures.py

clean:
	latexmk -C
	rm -f *.nav *.snm *.vrb *.bbl *.blg *.fls *.fdb_latexmk
