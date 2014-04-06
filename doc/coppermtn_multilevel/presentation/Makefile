TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: heat_eq_setup.pdf
	$(PDFTEX) copper_mountain_14.tex
	$(PDFTEX) copper_mountain_14.tex

heat_eq_setup.pdf: heat_eq_setup.fig
	fig2dev -L pdftex -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdf
	fig2dev -L pdftex_t -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdftex_t

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
