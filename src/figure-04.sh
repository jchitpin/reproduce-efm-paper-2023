#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/figure-04"

echo "Compiling Figure 4"
cd $DIR_FIGURE
pdflatex figure-04.tex
rm figure-04.aux figure-04.log
echo ""
echo "Figure 4 located under $DIR_FIGURE."

