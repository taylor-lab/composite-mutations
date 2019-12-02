#!/bin/bash

in=$1
out=$2

cmd="R -e \"rmarkdown::render('~/lab/repos/compound-mutations/do/"$in"', output_file = '~/lab/repos/compound-mutations/do/"$out"')\""
echo $cmd
eval $cmd
