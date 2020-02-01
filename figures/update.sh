#!/bin/bash

in=$1

cmd="R -e \"rmarkdown::render('~/lab/repos/composite-mutations/pdf/"$in".Rmd', output_file = '~/lab/repos/composite-mutations/pdf/"$in".pdf')\""
echo $cmd
eval $cmd
