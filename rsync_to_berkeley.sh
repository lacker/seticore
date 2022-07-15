#!/bin/bash -e
CMD="rsync -av --append-verify $1 rsync://blpd18.ssl.berkeley.edu/datax/"
echo $CMD
$CMD
