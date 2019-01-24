tar -xvf libpng-1.6.2.tar.gz
cd libpng-1.6.2
./configure --prefix=`pwd`
make
make install
LIBPNGDIR=`pwd`
cd ..
unzip blatSrc35.zip
cd blatSrc
cp $LIBPNGDIR/png.h lib/
cp $LIBPNGDIR/pngconf.h lib/
cp $LIBPNGDIR/pnglibconf.h lib/
echo $MACHTYPE
MACHTYPE=x86_64
export MACHTYPE
mkdir -p ~/bin/$MACHTYPE
make

