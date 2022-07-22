#!/bin/bash -e
# Usage:
#   ./install.sh <directory>
# Defaults to ~/bin
# Put the most recently compiled code in the current user's ~/bin directory


if [ "$1" == "" ]; then
    INSTALL_DIR=/home/$USER/bin
else
    INSTALL_DIR=$1
fi

BIN=./build/seticore

VERSION=`$BIN --help 2>&1 | grep version | sed 's:.* ::'`

INSTALLED_BIN=${INSTALL_DIR}/seticore
VERSIONED_BIN=${INSTALLED_BIN}-${VERSION}

if test -e "$VERSIONED_BIN"; then
    echo $VERSIONED_BIN already exists - delete it or bump the version
    exit 1
fi

echo installing seticore $VERSION
cp $BIN $VERSIONED_BIN
rm -f $INSTALLED_BIN
ln -s $VERSIONED_BIN $INSTALLED_BIN
echo version $VERSION now available at $INSTALLED_BIN
