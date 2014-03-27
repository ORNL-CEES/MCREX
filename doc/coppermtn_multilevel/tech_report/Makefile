TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: 
	$(PDFTEX) multilevel_mc
	$(BIB) multilevel_mc
	$(PDFTEX) multilevel_mc
	$(PDFTEX) multilevel_mc
	$(PDFTEX) multilevel_mc

clean: 
	-rm -f *~
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.dvi
	-rm -f *.rel
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.out
	-rm -f *.toc *.lof *.lot
	-rm -f *.bak
	-rm -f *.pdftex_t
	-rm -f *.pdf
	-rm -f TAGS
