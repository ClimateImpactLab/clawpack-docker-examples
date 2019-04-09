#!/usr/bin/env bash

docker pull climateimpactlab/geoclaw-example:latest
docker build \
    -t climateimpactlab/geoclaw-example:latest \
    --cache-from climateimpactlab/geoclaw-example:latest \
    .

docker push climateimpactlab/geoclaw-example:latest

docker run -i --rm \
    --name "geoclaw-local" \
    -v $(pwd)/examples:/home/examples \
    -t climateimpactlab/geoclaw-example:latest  \
    /bin/bash
