= Third party library dependencies and installation instruction =


== Overview of library dependencies ==

* BOOST: Graph classes, various iterator classes, ublas
* BOOST numeric bindings: bindings to lapack, etc (using syev, etc)
* opt++: Numeric optimization rutines
* newmat (v11): Linear algebra (Eigenvalue decomposition)
* NTL: Number theoretic library (provides float type with extended exponent range)

* bison++: grammar parser
* flex++: lexer 

Note that bison++ and flex++ are not needed to compile the code as such. They are only needed to generate the code for the newickParser.

Some of the below can probably be obtained as packaged binaries. In
longer terms the number of third parties dependencies should be
reduced. It is also possible to optionally include some libraries
using pre-processing directives, say opt++, which are used only in
isolated parts of the code.

== Installation instructions and web-pointers ==


=== BOOST ===

see http://sourceforge.net/projects/boost/files/boost/1.48.0/boost_1_48_0.tar.bz2/download

 # fetch 
 wget http://sourceforge.net/projects/boost/files/boost/1.48.0/boost_1_48_0.tar.gz/download -O boost_1_48_0.tar.gz
 tar xvzf boost_1_48_0.tar.gz
 cd boost_1_48_0
 
 # compile and install
 ./bootstrap.sh
 ./bjam
 sudo ./bjam install
 

=== BOOST numeric bindings ===

see http://mathema.tician.de/node/391

 # fetch
 git clone http://git.tiker.net/trees/boost-numeric-bindings.git
 cd boost-numeric-bindings/

 # compile and install
 ./configure
 sudo make install
 # fix permissions
 sudo find /usr/local/include/boost-numeric-bindings/ -type d -exec chmod a+x {} \;
 sudo mv /usr/local/include/boost-numeric-bindings/boost/numeric/bindings /usr/local/include/boost/numeric/


=== opt++ and newmat ===

====Overview====
opt++ comes bundles with newmat v11. We replace the bundled version
with the latest version and try to clean things up in a minimalistic
fashion.

====Homepages:====
* opt++: https://software.sandia.gov/opt++/
* newmt: http://www.robertnz.net/index.html

Registration is needed for download of opt++:
 https://software.sandia.gov/opt++/opt++_download.html

====compile and install ====
  # fetch
  tar xvzf tar xvzf optpp-2.4.tar.gz 
  cd optpp-2.4

  # compile and install
  configure --with-pic --includedir=/usr/local/include/optpp
  make
  sudo make install
  sudo chmod a+x /usr/local/include/optpp
  sudo cp include/OPT++_config.h /usr/local/include/optpp/

  # now the newmat include headers live in usr/local/include/optpp. We thus perform a
  # little hack to allow more standard calling of the include files
  # from code
  sudo ln -s /usr/local/include/optpp/ /usr/local/include/newmat

====Future improvements====

It would be nice to intall newmat indendently of opt++, which would
also allow newer versions to be used.

 
=== NTL ===
see www.shoup.net/ntl/

 # fetch
 wget http://www.shoup.net/ntl/ntl-6.0.0.tar.gz
 tar xvzf ntl-6.0.0.tar.gz
 cd ntl-6.0.0/src

 # Compile and install
 ./configure
 make # this takes 5 min...
 sudo make install


=== bison++ and flex++ === 
see http://code.google.com/p/flexpp-bisonpp/

 # fetch
 wget http://flexpp-bisonpp.googlecode.com/files/bisonpp-1.21-45.tar.gz
 tar xvzf bisonpp-1.21-45.tar.gz 
 cd bison++-1.21/

 # compile and install
 ./configure
 make
 sudo make install
 cd ..

 # fetch
 wget http://flexpp-bisonpp.googlecode.com/files/flexpp-2.3.8-45.tar.gz
 tar xvzf flexpp-2.3.8-45.tar.gz
 cd flex++-2.3.8-45/

 # compile and install
 ./configure 
 make
 sudo make install
