#!/bin/bash
sudo wget -c https://www.scientificlinux.org/documentation/gpg/RPM-GPG-KEY-cern -O RPM-GPG-KEY-cern.asc
sudo rpm --import RPM-GPG-KEY-cern.asc
wget -c http://people.centos.org/tru/devtools-1.1/devtools-1.1.repo -O /etc/yum.repos.d/devtools-1.1.repo
wget -c http://linuxsoft.cern.ch/cern/devtoolset/slc6-devtoolset.repo -O /etc/yum.repos.d/devtools-2.repo
wget -c http://dl.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
wget -c http://rpms.famillecollet.com/enterprise/remi-release-6.rpm
wget -c ftp://ftp.pbone.net/mirror/ftp5.gwdg.de/pub/opensuse/repositories/home:/cathay4t:/misc-rhel6/CentOS_CentOS-6/x86_64/cmake-2.8.9-3.1.x86_64.rpm
sudo rpm -U cmake-2.8.9-3.1.x86_64.rpm
sudo rpm -Uvh remi-release-6*.rpm epel-release-6*.rpm
sudo yum -y update
sudo yum -y install devtoolset-1.1 devtoolset-2 openmpi openmpi-devel git bison bison-devel flex lapack lapack-devel automake autoconf
export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
scl enable devtoolset-2 ./coreInstall.sh
