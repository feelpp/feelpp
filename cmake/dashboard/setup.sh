#! /bin/sh

mkdir -p $1
cd $1
if [ ! -d $1/feelpp/.git ]; then
    #git clone https://code.google.com/p/feelpp feelpp
    git clone https://github.com/feelpp/feelpp.git feelpp
fi

echo "Feel++ is now setup in $1 and ready for CDash dashboard"
if [ ! -f $1/feelpp/cmake/dashboard/feelpp.site.`hostname -s`.cmake ]; then
    echo "make sure that the file $1/feelpp/cmake/dashboard/feelpp.site.`hostname -s`.cmake exists and setup correctly"
fi
echo "edit crontab using 'crontab -e' and add the following line"
echo "0 1 * * * $1/feelpp/cmake/dashboard/linux.sh $1/feelpp/ Nightly"
echo
echo "here is the content of crontab:"
crontab -l
