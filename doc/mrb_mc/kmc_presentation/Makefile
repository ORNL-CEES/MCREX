TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: orthogonal_residual.pdf
	$(PDFTEX) krylov_mc.tex
	$(PDFTEX) krylov_mc.tex

orthogonal_residual.pdf: orthogonal_residual.fig
	fig2dev -L pdftex -p orthogonal_residual.pdf -m 0.4 orthogonal_residual.fig orthogonal_residual.pdf
	fig2dev -L pdftex_t -p orthogonal_residual.pdf -m 0.4 orthogonal_residual.fig orthogonal_residual.pdftex_t

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
