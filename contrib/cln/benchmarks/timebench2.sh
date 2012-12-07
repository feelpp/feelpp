#!/bin/sh

./timebench2a -n 100 -r 100000 > /dev/null
./timebench2a -n 1000 -r 10000 > /dev/null
./timebench2a -n 10000 -r 100 > /dev/null
./timebench2a -n 100000 -r 100 > /dev/null
./timebench2ap -r 100
