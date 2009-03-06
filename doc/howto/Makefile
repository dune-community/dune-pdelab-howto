SRCP = pdelab-howto-beamer
SRCA = pdelab-howto-article

beamer: 
	latex $(SRCP)
	bibtex $(SRCP)
	latex $(SRCP)
	dvipdf $(SRCP)
	xpdf -z 200 $(SRCP).pdf

art:
	latex $(SRCA)
	bibtex $(SRCA)
	latex $(SRCA)
	dvipdf $(SRCA)
	xpdf $(SRCA).pdf

clean:
	rm *.bbl *.blg *.dvi *.aux *.idx *.log *.lof *.toc *.lot *.ind *.ilg *.snm *.nav *.out
