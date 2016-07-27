TEX	= pdflatex
BIBTEX  = bibtex

SOURCES = paper.tex introduction.tex method.tex algorithm.tex limitations.tex characterization.tex conclusion.tex parameters.tex references.bib

all : paper.pdf

paper.pdf : $(SOURCES)
	$(TEX) paper
	$(BIBTEX) paper
	$(TEX) paper
	$(TEX) paper
clean :
	rm -f *aux *out *log *dvi *bbl *blg *toc paper.pdf
