#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/supplementary-01"

echo "Compiling Supplementary Figure 1 (wildtype/disease)"
cd $DIR_FIGURE
pdflatex supplementary-01-wildtype.tex
pdflatex supplementary-01-disease.tex
rm supplementary-01-wildtype.aux supplementary-01-wildtype.log
rm supplementary-01-disease.aux supplementary-01-disease.log
echo ""
echo "Supplementary Figure 1 (wildtype/disease) under $DIR_FIGURE."

