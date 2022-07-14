#!/bin/bash -e
# Deploy the most recently compiled code to meerkat

DEPLOY_DIR=/home/lacker/bin

BIN=./build/seticore

VERSION=`$BIN --help 2>&1 | grep version | sed 's:.* ::'`

DEPLOYED_BIN=${DEPLOY_DIR}/seticore
VERSIONED_BIN=${DEPLOYED_BIN}-${VERSION}

if test -e "$VERSIONED_BIN"; then
    echo $VERSIONED_BIN already exists - delete it or bump the version
    exit 1
fi

echo deploying seticore $VERSION
cp $BIN $VERSIONED_BIN
rm -f $DEPLOYED_BIN
ln -s $VERSIONED_BIN $DEPLOYED_BIN
echo version $VERSION now available at $DEPLOYED_BIN
