#!/bin/sh
set -eu

for i in core coefficientformpdes heat electric fluid solid heatfluid thermoelectric; do
    depends='${shlibs:Depends}, ${misc:Depends}, ${python3:Depends}, libfeelpp-toolboxes1-'"$i"' (= ${binary:Version}), python3-feelpp-toolboxes-core'
    case "$i" in
        core)
            depends='${shlibs:Depends}, ${misc:Depends}, ${python3:Depends}, libfeelpp-toolboxes1-core (= ${binary:Version}), python3-feelpp'
            ;;
        heatfluid)
            depends='${shlibs:Depends}, ${misc:Depends}, ${python3:Depends}, libfeelpp-toolboxes1-heatfluid (= ${binary:Version}), python3-feelpp-toolboxes-heat, python3-feelpp-toolboxes-fluid'
            ;;
        thermoelectric)
            depends='${shlibs:Depends}, ${misc:Depends}, ${python3:Depends}, libfeelpp-toolboxes1-thermoelectric (= ${binary:Version}), python3-feelpp-toolboxes-heat, python3-feelpp-toolboxes-electric'
            ;;
    esac
    cat <<EOF
Package: python3-feelpp-toolboxes-$i
Section: python
Architecture: amd64 i386 ia64 powerpc sparc
Depends: $depends
Recommends: feelpp-toolboxes-data
Description: Feel++ toolboxes $i Python bindings
 Provides Feel++ $i toolbox via Python bindings

EOF
done
