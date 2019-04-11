#!/usr/bin/env bash

cd clawpack
git submodule init
git submodule update

for submod in amrclaw classic clawutil geoclaw pyclaw riemann visclaw; do
    cd $submod
    git submodule init
    git submodule update
    cd ..
done

cd ..

bash get_dems.sh

docker pull climateimpactlab/geoclaw-example:latest
docker build \
    -t climateimpactlab/geoclaw-example:latest \
    --cache-from climateimpactlab/geoclaw-example:latest \
    .
docker push climateimpactlab/geoclaw-example:latest

docker run -i --rm \
    --name "geoclaw-local" \
    -v $(pwd)/examples:/home/examples \
    -v $(pwd)/clawpack:/clawpack \
    -t climateimpactlab/geoclaw-example:latest  \
    /bin/bash
