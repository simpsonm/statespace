# LaTeX Makefile

## Define the file name (.tex)

DOC=dlmasis
APP=appendices
DOCB=dlmasisblind
APPB=appendicesblind
FULL=dlmasisfull
FULLB=dlmasisblindfull

all:	$(DOC).pdf $(APP).pdf $(DOCB).pdf $(APPB).pdf \
 $(FULL).pdf $(FULLB).pdf

.PHONY: clean

clean:
	\rm *.aux *.blg *.out *.bbl *.log $(DOC).pdf $(DOCB).pdf \
 $(APP).pdf $(APPB).pdf $(FULL).pdf $(FULLB).pdf

$(DOC).pdf: 	$(DOC).tex
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	bibtex $(DOC)
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(DOC).tex}"

$(APP).pdf: 	$(APP).tex
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	bibtex $(APP)
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"

$(DOCB).pdf: 	$(DOC).pdf
	cp $(DOC).pdf temp.pdf
	pdflatex "\def\blind{1} \input{$(DOC).tex}"
	pdflatex "\def\blind{1} \input{$(DOC).tex}"
	bibtex $(DOC)
	pdflatex "\def\blind{1} \input{$(DOC).tex}"
	pdflatex "\def\blind{1} \input{$(DOC).tex}"
	cp $(DOC).pdf $(DOCB).pdf
	cp temp.pdf $(DOC).pdf
	rm temp.pdf

$(APPB).pdf: 	$(APP).pdf
	cp $(APP).pdf temp.pdf
	pdflatex "\def\blind{1} \input{$(APP).tex}"
	pdflatex "\def\blind{1} \input{$(APP).tex}"
	bibtex $(APP)
	pdflatex "\def\blind{1} \input{$(APP).tex}"
	pdflatex "\def\blind{1} \input{$(APP).tex}"
	cp $(APP).pdf $(APPB).pdf
	cp temp.pdf $(APP).pdf
	rm temp.pdf

$(FULL).pdf: $(DOC).pdf $(APP).pdf
	pdfunite $(DOC).pdf $(APP).pdf $(FULL).pdf

$(FULLB).pdf: $(DOCB).pdf $(APPB).pdf
	pdfunite $(DOCB).pdf $(APPB).pdf $(FULLB).pdf
