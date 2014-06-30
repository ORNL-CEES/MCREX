TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all:
	$(PDFTEX) mrb_mc.tex
	$(PDFTEX) mrb_mc.tex

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
