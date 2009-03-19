# $id: $

EPS_FILES = 
#EPS/adaptiveintegration_alberta2d.eps EPS/adaptiveintegration_ug3d.eps

DOCSOURCE = pdelab-howto-article.tex pdelab-howto-beamer.tex pdelab-howto-main.tex

MAINDOCS = pdelab-howto-article.dvi pdelab-howto-article.pdf \
	pdelab-howto-article.ps pdelab-howto-beamer.dvi pdelab-howto-beamer.pdf \
	pdelab-howto-beamer.ps

EXTRA_DIST = pdelab-howto-article.pdf pdelab-howto-beamer.pdf Makefile.dist.in Makefile.dist.am

if BUILD_DOCS
if CONVERT
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/pdelab-howto

docdir=$(datadir)/doc/dune-pdelab-howto

doc_DATA = pdelab-howto-article.pdf pdelab-howto-beamer.pdf

EXTRAINSTALL = pdelab-howto-article.pdf

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps *.nav *.snm *.vrb
CLEANFILES += $(DOCS)

if CONVERT
CLEANFILES += $(EPS_FILES)

#%.eps: $(srcdir)/%.jpg
#	convert $< $@
pdelab-howto-article.tex: $(EPS_FILES) lit.bib
pdelab-howto-beamer.tex: $(EPS_FILES) lit.bib
endif

dist-hook:
	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS:=$(top_srcdir)

pdelab-howto-article.dvi: pdelab-howto-main.tex $(EPS_FILES) lit.bib
pdelab-howto-beamer.dvi: pdelab-howto-main.tex $(EPS_FILES) lit.bib

art: pdelab-howto-article.pdf
	xpdf -z 200 pdelab-howto-article.pdf
beamer: pdelab-howto-beamer.pdf
	xpdf -z 200 pdelab-howto-beamer.pdf

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
