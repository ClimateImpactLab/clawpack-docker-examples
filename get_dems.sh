#!/usr/bin/env bash

IFS=''

while read line;
do
    if [[ $line != "#"* ]]; then
        if [[ ! -f examples/data/$(basename $line) ]]; then
            echo "downloading $(basename $line)"
            curl $line > examples/data/$(basename $line);
        fi;
    fi;
done < "highres_dems.txt"
