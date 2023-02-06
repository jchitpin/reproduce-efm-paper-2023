#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/figure-01"

echo "Compiling Figure 1"
cd $DIR_FIGURE
pdflatex figure-01.tex
rm figure-01.aux figure-01.log
echo ""
echo "Figure 1 located under $DIR_FIGURE."

