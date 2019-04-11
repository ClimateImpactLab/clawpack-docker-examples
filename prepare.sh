#!/bin/bash

cd /clawpack

pip install -e .

cd /home

# run extra commands
$@
