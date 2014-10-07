TEX	= latex
BIBTEX  = bibtex
DVIPS	= dvips
PSPDF	= ps2pdf

SOURCES = paper.tex introduction.tex algorithm.tex characterization.tex conclusion.tex references.bib

all : paper.pdf

paper.pdf : $(SOURCES)
	$(TEX) paper
	$(BIBTEX) paper
	$(TEX) paper
	$(TEX) paper
	$(DVIPS) paper
	$(PSPDF) paper.ps
clean :
	rm -f *aux *log *dvi *bbl *blg *toc *.ps *.pdf
