#!/bin/sh
set -eu

for i in core coefficientformpdes heat electric fluid solid heatfluid thermoelectric; do
    dest="$i"
    if [ "$i" = "coefficientformpdes" ]; then
        dest=cfpdes
    fi
    cat > python3-feelpp-toolboxes-$i.install << EOF
usr/lib/python3*/*-packages/feelpp/toolboxes/$dest/*.py
usr/lib/python3*/*-packages/feelpp/toolboxes/$dest/*.so
EOF
done
