#!/bin/bash

cd src/
g++ -pg -O3 main.cpp Locus.cpp Individual.cpp model.cpp matrix.cpp readBinData.cpp -o ../mplink
cd ../