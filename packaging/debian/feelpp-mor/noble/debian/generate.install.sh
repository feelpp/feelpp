#!/bin/bash

for _tb in core coefficientformpdes heat electric fluid solid hdg fsi heatfluid thermoelectric
do
    if [ "$_tb" == "core"  ] ; then
        # libraries
        cat >  libfeelpp-toolboxes1-${_tb}.install <<EOF
usr/lib/*/libfeelpp_model*.so.*
EOF

        # dev files
        cat > libfeelpp-toolboxes1-${_tb}-dev.install <<EOF
usr/lib/*/libfeelpp_model*.so
usr/include/feelpp/toolboxes/feel/feelmodels/modelcore/*.hpp
usr/include/feelpp/toolboxes/feel/feelmodels/modelmesh/*.hpp
usr/include/feelpp/toolboxes/feel/feelmodels/modelvf/*.hpp
usr/include/feelpp/toolboxes/feel/feelmodels/modelmaterials/*.hpp
usr/share/feelpp/toolboxes/cmake
EOF


    elif [ "$_tb" == "coefficientformpdes"  ] ; then
        # libraries
        cat > libfeelpp-toolboxes1-${_tb}.install <<EOF
usr/lib/*/libfeelpp_toolbox_coefficientformpde*.so.*
EOF

        # dev files
        cat > libfeelpp-toolboxes1-${_tb}-dev.install <<EOF
usr/lib/*/libfeelpp_toolbox_coefficientformpde*.so
usr/include/feelpp/toolboxes/feel/feelmodels/${_tb}/*.hpp
EOF

        # executables
        cat > feelpp-toolboxes-${_tb}.install <<EOF
usr/bin/feelpp_toolbox_${_tb}
EOF
#        usr/share/man/man1/feelpp_toolbox_${_tb}.*
#        usr/share/doc/feelpp/toolboxes/feelpp_toolbox_${_tb}.*

    else
        # libraries
        echo "usr/lib/*/libfeelpp_toolbox_${_tb}_*.so.*" > libfeelpp-toolboxes1-${_tb}.install

        # dev files
        cat > libfeelpp-toolboxes1-${_tb}-dev.install <<EOF
usr/lib/*/libfeelpp_toolbox_${_tb}_*.so
usr/include/feelpp/toolboxes/feel/feelmodels/${_tb}/*.hpp
EOF

        # executables
        if [ "$_tb" == "fsi"  ] || [ "$_tb" == "hdg"  ] ; then
            cat > feelpp-toolboxes-${_tb}.install <<EOF
usr/bin/feelpp_toolbox_${_tb}_*
EOF
        else
        cat > feelpp-toolboxes-${_tb}.install <<EOF
usr/bin/feelpp_toolbox_${_tb}
usr/share/man/man1/feelpp_toolbox_${_tb}.*
usr/share/doc/feelpp/toolboxes/feelpp_toolbox_${_tb}.*
EOF
        fi
    fi
done

