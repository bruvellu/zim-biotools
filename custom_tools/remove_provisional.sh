#!/bin/bash

ROOT=`dirname $1`

grep -rl "@provisional" $ROOT | xargs rm
