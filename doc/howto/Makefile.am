# $id: $

# setting like in dune-web
CURDIR=doc/pdelab-howto
BASEDIR=../..

# install the html pages
docdir=$(datadir)/doc/dune-pdelab-howto
if BUILD_DOCS
  DOCFILES = pdelab-howto-article.dvi pdelab-howto-article.pdf \
	pdelab-howto-article.ps pdelab-howto-beamer.dvi pdelab-howto-beamer.pdf \
	pdelab-howto-beamer.ps
  EXTRA_DIST = $(DOCFILES)
  EXTRAINSTALL = $(DOCFILES)
endif

if BUILD_DOCS
if CONVERT
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif

TEX_FILES = functionspaces.tex istl-backend.tex pdelab-howto-beamer.tex \
	summary.tex instationaryproblems.tex  operators.tex \
	pdelab-howto-main.tex systems.tex introduction.tex \
	pdelab-howto-article.tex  stationaryproblems.tex   weightedresidual.tex \
        worklow.tex

pdelab-howto-article.dvi: pdelab-howto-main.tex lit.bib $(TEX_FILES)
pdelab-howto-beamer.dvi: pdelab-howto-main.tex lit.bib $(TEX_FILES)

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps *.nav *.snm *.vrb
CLEANFILES += $(DOCSFILES)

EXTRA_TEXINPUTS = $(top_srcdir)

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
