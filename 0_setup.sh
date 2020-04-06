#!/bin/bash

set -o nounset -o errexit

## Basic setup of results directory structure.
## Assumes that $GOATSDIR is defined and exists.

cd $RAREDIR
mkdir data_v8
mkdir features_v8
mkdir figures
mkdir paper_figures
mkdir preprocessing_v8
