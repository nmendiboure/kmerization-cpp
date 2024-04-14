#!/usr/bin/env sh

# This script is used to generate the proto files.
protoc --cpp_out=. matrix.proto
protoc --cpp_out=. index.proto