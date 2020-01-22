#!/bin/bash

in=$1

cmd="R -e \"rmarkdown::render('~/lab/repos/compound-mutations/pdf/"$in".Rmd', output_file = '~/lab/repos/compound-mutations/pdf/"$in".pdf')\""
echo $cmd
eval $cmd
