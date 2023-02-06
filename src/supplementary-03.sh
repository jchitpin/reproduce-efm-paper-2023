#!/usr/bin/env

DIR_FIGURE="$(pwd)/../figures/supplementary-03"
DIR_DATA="$(pwd)/../data/"

echo "Compiling Supplementary Figure 3"
cd $DIR_DATA

# Manually modify table to pretty print scientific notation
# Do this rather than modifying the matrix2table julia function
# Changing line 2 with the following and saving to a new file:
sed '2s/.*/\\begin{tabularx}{\\textwidth}{X *{4}{S[table-format=-1.3e-2]}lX}/' table-markov-sphingo.tex > table-markov-sphingo-edited.tex

cd $DIR_FIGURE
pdflatex supplementary-03.tex
rm supplementary-03.aux supplementary-03.log
echo ""
echo "Supplementary Figure 3 under $DIR_FIGURE."

