#!/bin/sh

# pdf
Rscript rmd2md.R 
pandoc cea-vignette.md --latex-engine /Library/TeX/texbin/pdflatex --bibliography cea-references.bib --variable graphics=yes -o cea-vignette.pdf
pandoc cea-vignette.md --bibliography cea-references.bib -o cea-vignette.html
