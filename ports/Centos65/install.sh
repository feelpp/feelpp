#!/bin/bash
wget http://people.centos.org/tru/devtools-1.1/devtools-1.1.repo -O /etc/yum.repos.d/devtools-1.1.repo
wget http://dl.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
wget http://rpms.famillecollet.com/enterprise/remi-release-6.rpm
sudo rpm -Uvh remi-release-6*.rpm epel-release-6*.rpm
sudo yum -y update
sudo yum -y install devtoolset-1.1 cmake28 cmake28-gui openmpi openmpi-devel git bison bison-devel flex lapack lapack-devel automake autoconf
sudo ln -s /usr/bin/cmake28 /usr/bin/cmake
export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
scl enable devtoolset-1.1 ./coreInstall.sh
