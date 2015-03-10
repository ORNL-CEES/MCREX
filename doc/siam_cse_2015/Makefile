TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: ddnu_example.pdf domain_to_domain.pdf transport_loop.pdf binary_comm_tree.pdf master_comm_tree.pdf
	$(PDFTEX) siam_cse_2015
	$(PDFTEX) siam_cse_2015
	$(PDFTEX) siam_cse_2015

ddnu_example.pdf: ddnu_example.fig
	fig2dev -L pdftex -p ddnu_example.pdf -m 0.4 ddnu_example.fig ddnu_example.pdf
	fig2dev -L pdftex_t -p ddnu_example.pdf -m 0.4 ddnu_example.fig ddnu_example.pdftex_t

domain_to_domain.pdf: domain_to_domain.fig
	fig2dev -L pdftex -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdf
	fig2dev -L pdftex_t -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdftex_t

transport_loop.pdf: transport_loop.fig
	fig2dev -L pdftex -p transport_loop.pdf -m 0.4 transport_loop.fig transport_loop.pdf
	fig2dev -L pdftex_t -p transport_loop.pdf -m 0.4 transport_loop.fig transport_loop.pdftex_t

binary_comm_tree.pdf: binary_comm_tree.fig
	fig2dev -L pdftex -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdf
	fig2dev -L pdftex_t -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdftex_t

master_comm_tree.pdf: master_comm_tree.fig
	fig2dev -L pdftex -p master_comm_tree.pdf -m 0.4 master_comm_tree.fig master_comm_tree.pdf
	fig2dev -L pdftex_t -p master_comm_tree.pdf -m 0.4 master_comm_tree.fig master_comm_tree.pdftex_t

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
	-rm -f siam_cse_2015.pdf
