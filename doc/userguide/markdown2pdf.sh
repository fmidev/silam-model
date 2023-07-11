#!/bin/bash

title="SILAMv5 User Manual" 
subtitle=`LC_ALL=en_US.UTF-8 date +"%a %d %b, %Y "` 
version="5.8" #x.x.x

ouf="SILAMv${version}-UserGuide.pdf"
inf="silam_user-guide.md"

#fonts:
mainfont="DejaVu Serif" #"Liberation Sans"
sansfont="Liberation Sans Narrow"
monofont="Liberation Mono"
fontsize="14pt"

pandoc -s -N --template=./templates/mytemplate.tex --variable mainfont="$mainfont" --variable sansfont="$sansfont" --variable monofont="$monofont" --variable fontsize="$fontsize" --variable version="${version}" --variable title="${title}" --variable subtitle="${subtitle}" --variable geometry:margin=0.5in --toc --pdf-engine=xelatex -f markdown-implicit_figures -s -o ${ouf} ${inf}
