TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all:
	$(PDFTEX) slattery_wilson_evans
	$(BIB) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans

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
