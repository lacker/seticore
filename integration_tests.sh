#!/bin/bash -e

# Runs the integration tests. These require a bunch of specific local data but
# they catch more bugs than just the unit tests.

echo integration testing...

./build/pipeline_integration_test

./test_extract.sh
