# Makefile for LaTeX
PDFC = xelatex
BIBC = bibtex
READER = mupdf
 
TARGET = report.pdf
OUTGROWTH = *.out *.log *.aux *.toc *.lof *.lot *.lol *.bbl *.blg
#OUTGROWTH = *.out *.log *.aux
SOURCE = report.tex
SOURCES = *.tex
BSOURCE = report
  
STRIP = detex
COUNTER = wc

# Here is the compilation of documents
# with bibtex.
#
$(TARGET): $(SOURCE)
	$(PDFC)    $(SOURCE)
	$(BIBC)    $(BSOURCE)
	$(PDFC)    $(SOURCE)
	$(PDFC)    $(SOURCE)

#$(TARGET): $(SOURCES)
#	$(PDFC)    $(SOURCE)  && $(PDFC)    $(SOURCE)

read: report.pdf
	$(READER) $(TARGET)

count: report.tex
	$(STRIP)  $(SOURCES) | uniq | $(COUNTER) -m

clean:
	-rm $(OUTGROWTH)
	 
purge:
	-rm $(OUTGROWTH) $(TARGET)

love:
	echo "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nOh... Je t'aime !\n\nJe t'aime bien.\n\n\n"
	 
.PHONY:read count clean cleanall
