#!/bin/bash

# Name: Remove provisional
# Description: Removes notes with tag @provisional from current directory
# Command: ~/.local/share/zim/plugins/zim-biotools/custom_tools/remove_provisional.sh %s

ROOT=`dirname $1`

grep -rl "@provisional" $ROOT | xargs rm
