#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/supplementary-04"

echo "Compiling Supplementary Figure 4"
cd $DIR_FIGURE
pdflatex supplementary-04.tex
rm supplementary-04.aux supplementary-04.log
echo ""
echo "Supplementary Figure 4 under $DIR_FIGURE."

