apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install gcc g++ gfortran wget cpio && \
    cd /tmp && \
    wget -q http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/12414/l_mkl_2018.1.163.tgz && \
    tar -xzf l_mkl_2018.1.163.tgz && \
    cd l_mkl_2018.1.163 && \
    sed -i 's/ACCEPT_EULA=decline/ACCEPT_EULA=accept/g' silent.cfg && \
    sed -i 's/ARCH_SELECTED=ALL/ARCH_SELECTED=INTEL64/g' silent.cfg && \
    #  sed -i 's/COMPONENTS=DEFAULTS/COMPONENTS=;intel-comp-l-all-vars__noarch;intel-openmp-l-all__x86_64;intel-openmp-l-ps-libs__x86_64;intel-openmp-l-ps-libs-jp__x86_64;intel-tbb-libs__noarch;intel-mkl-common__noarch;intel-mkl-sta-common__noarch;intel-mkl__x86_64;intel-mkl-rt__x86_64;intel-mkl-ps-rt-jp__x86_64;intel-mkl-doc__noarch;intel-mkl-ps-doc__noarch;intel-mkl-ps-doc-jp__noarch;intel-mkl-gnu__x86_64;intel-mkl-gnu-rt__x86_64;intel-mkl-ps-common__noarch;intel-mkl-ps-common-jp__noarch;intel-mkl-ps-common-64bit__x86_64;intel-mkl-common-c__noarch;intel-mkl-common-c-64bit__x86_64;intel-mkl-ps-common-c__noarch;intel-mkl-doc-c__noarch;intel-mkl-ps-doc-c-jp__noarch;intel-mkl-ps-ss-tbb__x86_64;intel-mkl-ps-ss-tbb-rt__x86_64;intel-mkl-gnu-c__x86_64;intel-mkl-ps-common-f__noarch;intel-mkl-ps-common-f-64bit__x86_64;intel-mkl-ps-doc-f__noarch;intel-mkl-ps-doc-f-jp__noarch;intel-mkl-ps-gnu-f-rt__x86_64;intel-mkl-ps-gnu-f__x86_64;intel-mkl-ps-f95-common__noarch;intel-mkl-ps-f__x86_64;intel-mkl-psxe__noarch;intel-psxe-common__noarch;intel-psxe-common-doc__noarch;intel-compxe-pset/g' silent.cfg && \
    sed -i 's/COMPONENTS=DEFAULTS/COMPONENTS=;intel-comp-l-all-vars__noarch;intel-comp-nomcu-vars__noarch;intel-openmp__x86_64;intel-tbb-libs__x86_64;intel-mkl-common__noarch;intel-mkl-installer-license__noarch;intel-mkl-core__x86_64;intel-mkl-core-rt__x86_64;intel-mkl-doc__noarch;intel-mkl-doc-ps__noarch;intel-mkl-gnu__x86_64;intel-mkl-gnu-rt__x86_64;intel-mkl-common-ps__noarch;intel-mkl-core-ps__x86_64;intel-mkl-common-c__noarch;intel-mkl-core-c__x86_64;intel-mkl-common-c-ps__noarch;intel-mkl-tbb__x86_64;intel-mkl-tbb-rt__x86_64;intel-mkl-gnu-c__x86_64;intel-mkl-common-f__noarch;intel-mkl-core-f__x86_64;intel-mkl-gnu-f-rt__x86_64;intel-mkl-gnu-f__x86_64;intel-mkl-f95-common__noarch;intel-mkl-f__x86_64;intel-mkl-psxe__noarch;intel-psxe-common__noarch;intel-psxe-common-doc__noarch;intel-compxe-pset/g' silent.cfg && \
    ./install.sh -s silent.cfg && \
    cd .. && rm -rf * && \
    rm -rf /opt/intel/.*.log /opt/intel/compilers_and_libraries_2018.1.163/licensing && \
    echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && \
    ldconfig && \
      echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc
