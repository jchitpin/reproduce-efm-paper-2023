#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/figure-03"

echo "Compiling Figure 3"
cd $DIR_FIGURE
pdflatex figure-03.tex
rm figure-03.aux figure-03.log
echo ""
echo "Figure 3 located under $DIR_FIGURE."

