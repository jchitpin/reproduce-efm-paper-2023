#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/figure-02"

echo "Compiling Figure 2"
cd $DIR_FIGURE
pdflatex figure-02.tex
rm figure-02.aux figure-02.log
echo ""
echo "Figure 2 located under $DIR_FIGURE."

