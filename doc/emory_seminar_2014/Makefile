TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: emory_talk.tex
	$(PDFTEX) emory_talk
	$(PDFTEX) emory_talk
	$(PDFTEX) emory_talk

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
	-rm -f emory_talk.pdf
