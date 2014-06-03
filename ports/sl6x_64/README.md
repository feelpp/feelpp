Using LinuxBrew on "Scientific" Linux or RedHat 6.5

**THAT DOES NOT WORK AT THIS POINT !**

The main problem is the lack of a compiler that support c++11.

# First step - Preparing system
```
sudo yum update -y
sudo yum groupinstall -y "Development Tools"
sudo yum install -y \
        autoconf automake19 libtool gettext \
        git scons cmake flex bison \
        libcurl-devel curl \
        ncurses-devel ruby bzip2-devel expat-devel
```

# Second step - Installing Linux brew

Do it as user !
```
git clone https://github.com/Homebrew/linuxbrew.git ~/.linuxbrew
```

and cc that in /etc/bashrc

```
# Until LinuxBrew is fixed, the following is required.
# See: https://github.com/Homebrew/linuxbrew/issues/47
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/lib64/pkgconfig:/usr/lib64/pkgconfig:/usr/lib/pkgconfig:/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/lib64/pkgconfig:/usr/share/pkgconfig:$PKG_CONFIG_PATH
## Setup linux brew
export LINUXBREWHOME=$HOME/.linuxbrew
export PATH=$LINUXBREWHOME/bin:$PATH
export MANPATH=$LINUXBREWHOME/man:$MANPATH
export PKG_CONFIG_PATH=$LINUXBREWHOME/lib64/pkgconfig:$LINUXBREWHOME/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$LINUXBREWHOME/lib64:$LINUXBREWHOME/lib:$LD_LIBRARY_PATH
```

source it
```
source ~/.bashrc
```
or close/open the terminal
and test the brew install:
```
which brew
/home/ubuntu/.linuxbrew/bin/brew
echo $PKG_CONFIG_PATH
/home/ubuntu/.linuxbrew/lib64/pkgconfig:/home/ubuntu/.linuxbrew/lib/pkgconfig:/usr/local/lib/pkgconfig:/usr/local/lib64/pkgconfig:/usr/lib64/pkgconfig:/usr/lib/pkgconfig:/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/lib64/pkgconfig:/usr/share/pkgconfig:
```

# Third step - Install feel++

su -c "yum install zlib-devel"
su -c "yum install libtool"
su -c "yum -y install glibc-devel.i686 glibc-devel"

```
brew install gcc
ln -s .linuxbrew/Cellar/gcc....64 .linuxbrew/bin/gcc64
for i in `ls ~/.linuxbrew/Cellar/gcc/4.8.3/lib64/`; do ln -s ~/.linuxbrew/Cellar/gcc/4.8.3/lib64/$i ~/.linuxbrew/lib/; done
brew install --with-openblas --whith-shared --withou-check scalapack
brew install gmsh
brew install homebrew/versions/llvm34 --with-clang --without-libffi
brew install feelpp
```

