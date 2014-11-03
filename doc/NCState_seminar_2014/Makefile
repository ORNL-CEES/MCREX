TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: ncstate_seminar.tex
	$(PDFTEX) ncstate_seminar
	$(PDFTEX) ncstate_seminar
	$(PDFTEX) ncstate_seminar

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
	-rm -f *.nav
	-rm -f *.snm
	-rm -f TAGS

realclean: clean
	-rm -f ncstate_seminar.pdf
