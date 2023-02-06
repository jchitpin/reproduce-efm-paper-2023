#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/supplementary-02"

echo "Compiling Supplementary Figure 2 (wildtype/disease)"
cd $DIR_FIGURE
pdflatex supplementary-02-wildtype.tex
pdflatex supplementary-02-disease.tex
rm supplementary-02-wildtype.aux supplementary-02-wildtype.log
rm supplementary-02-disease.aux supplementary-02-disease.log
echo ""
echo "Supplementary Figure 2 (wildtype/disease) under $DIR_FIGURE."

