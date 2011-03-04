# !/bin/bash

latex note.tex
latex note.tex

dvipdf note.dvi

gv --scale=2 note.pdf
#gv note.pdf
