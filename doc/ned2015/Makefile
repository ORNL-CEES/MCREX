TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: stencil_graph.pdf
	$(PDFTEX) slattery_wilson_evans
	$(BIB) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans

stencil_graph.pdf: stencil_graph.fig
	fig2dev -L pdftex -p stencil_graph.pdf -m 0.4 stencil_graph.fig stencil_graph.pdf
	fig2dev -L pdftex_t -p stencil_graph.pdf -m 0.4 stencil_graph.fig stencil_graph.pdftex_t

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
	-rm -f TAGS
	-rm -f *.pdf
