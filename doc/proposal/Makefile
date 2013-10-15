##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

PSTEX  := latex
TEX    := pdflatex
BIB    := bibtex
DVI    := dvips
PDF    := ps2pdf13
SED    := sed
PSBOOK := psbook

XFIG   := $(wildcard *fig)
PDFFIG := $(patsubst %.fig, %.pdf, $(XFIG))
TEXSRC := $(wildcard *.tex)
STYSRC := $(wildcard *.sty)

SUPSRC := $(wildcard supp/*.tex)
BIOSRC := $(wildcard bios/*.tex)
CNPSRC := $(wildcard cnp/*.tex)

TEXSRC += $(SUPSRC) $(BIOSRC) $(CNPSRC)

##---------------------------------------------------------------------------##

.PHONY: clean pdf ps realclean all

all: proposal.pdf

white: white.pdf
emory: emory.pdf

white.pdf: white.tex $(PDFFIG)
	$(TEX) white
	$(TEX) white

proposal.pdf: proposal.tex $(PDFFIG) $(TEXSRC) $(STYSRC)
	$(TEX) proposal
	$(BIB) proposal
	$(TEX) proposal
	$(TEX) proposal
	$(TEX) proposal

emory.pdf: emory.tex $(PDFFIG) $(TEXSRC) $(STYSRC)
	$(TEX) emory
	$(BIB) emory
	$(TEX) emory
	$(TEX) emory
	$(TEX) emory

#set.pdf: set.fig 
#	fig2dev -L pdftex   -p set.pdf -m 0.4 set.fig set.pdf
#	fig2dev -L pdftex_t -p set.pdf -m 0.4 set.fig set.pdftex_t

#set.pdf: set.fig 
#	fig2dev -L pdftex   -p set.pdf -m 0.3 set.fig set.pdf
#	fig2dev -L pdftex_t -p set.pdf -m 0.3 set.fig set.pdftex_t

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
	-rm -f set.pdf set.pdftex_t
	-rm -f TAGS
	-rm -f *.brf

realclean: clean
	-rm -f white.pdf
	-rm -f proposal.pdf
	-rm -f emory.pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
