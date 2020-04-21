danotes:
	pdflatex $@
	bibtex $@
	pdflatex $@
	pdflatex $@

refs:
	python ~/python/bibclean.py references_raw.bib -o references.bib --fix-thesis
