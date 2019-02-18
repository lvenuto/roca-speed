# roca-speed
Fast, multiprocess implementation in Python of the RoCa vulnerability of RSA keys

## Installation instructions
### Ubuntu 18.04

Install Python 2.7 and pip
```
sudo apt-get install python-minimal && sudo apt-get install python-pip
```

We need Cython because it's a requirement for python-flint, grab it with pip
```
sudo pip install Cython
```

[libgmp](https://gmplib.org/) should be already installed in Ubuntu 16.04, but if it isn't:

```
sudo apt-get install libgmp-dev
```

We also need the multiprecision library for floats [mpfr](https://www.mpfr.org/), it's a requirement for FLINT C. It should be already installed in ubuntu 18.04, but not in ubuntu 16.04
```
sudo apt-get install libmpfr6
```

This script uses the [FLINT C library](http://www.flintlib.org/) to factor polynomials, it's way faster and easier to use than Sympy, it's in the ubuntu repo:
```
sudo apt-get install libflint-dev
```

We also need [arb](https://github.com/fredrik-johansson) because it's a requirement for python-flint (TODO: remove this dependancy, arb is not needed for our stuff..)
Additionally the version in the repository is too old, and we need to compile from source, ARGH!
```
wget https://github.com/fredrik-johansson/arb/archive/2.16.0.tar.gz
tar -xvf 2.16.0.tar.gz
cd arb-2.16.0/
./configure
make
sudo make install
```

We need numpy because the python-flint installer uses it... TODO: remove this dependency
```
sudo pip install numpy
```

Now install python-flint
```
sudo pip install python-flint
```

Now we need to install fplll(fastest lll implementation)
```
sudo apt-get install libfplll-dev
```

Now we install the python wrapper for lll, from the package manager and not pip, because the version in the pip repo has some [issues](https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=905434)
```
sudo apt-get install python-fpylll
```

Install gmpy2 (used to test primality) and the library associated with it:
```
sudo apt-get install libmpc-dev
sudo pip install gmpy2
```

And finally clone this repo and run some tests:
```
git clone https://github.com/lvenuto/roca-speed
cd roca-speed
python2 test.py
```

### Mac os X

The easiest way to install the requirements is through [Brew](https://brew.sh/)
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
This script runs on Python 2.7, install it if you haven't already
```
brew install python2
```
We need Cython because it's a requirement for python-flint 
```
brew install cython
```
We need the multiprecision library [gmp](https://gmplib.org/), because it's a requirement for FLINT C
```
brew install gmp
```
We also need the multiprecision library for floats [mpfr](https://www.mpfr.org/), it's a requirement for FLINT C 

This script uses the [FLINT C library](http://www.flintlib.org/) to factor polynomials, it's way faster and easier to use than Sympy, but it's not in the brew repo, so we need to compile it.
```
brew install wget
wget http://www.flintlib.org/flint-2.5.2.tar.gz
tar -xvf flint-2.5.2.tar.gz
cd flint-2.5.2
./configure
./make
#Takes 10 minutes
./make install
```

We also need [arb](https://github.com/fredrik-johansson) because it's a requirement for python-flint

```
wget https://github.com/fredrik-johansson/arb/archive/2.15.1.tar.gz
cd arb-2.15.1
./configure
make
make install
```

Now with pip we can install python-flint
```
pip install python-flint
```

Now we need [fpylll](https://github.com/fplll/fpylll) and [fplll](https://github.com/fplll/fplll) for the lll reduction
First we need to create the configure script, so either get it from: [Download](https://www.dropbox.com/s/ohprvleybgvgk3n/configure?dl=0)
or 
```
brew install autoconf
brew install automake
brew install libtool
```
Then clone the repo containing fplll (fastest lll implementation)
and there it is, another round of compile-action!
```
git clone https://github.com/fplll/fplll
cd fplll
./autogen.sh #necessary only if you haven't downloaded the configure file
./configure
./make
./make install
```
Now pip comes in handy again:
```
pip install cysignals
pip install fpylll
```
The installation of fpylll fails because:
```
  build/src/fpylll/fplll/integer_matrix.cpp:691:10: fatal error: 'fplll/pruner.h' file not found
  #include "fplll/pruner.h"
           ^~~~~~~~~~~~~~~~
  1 error generated.
  error: command 'clang' failed with exit status 1
```
`pruner.h` has been moved, the simplest way is to create a symbolic link:

```
cd /usr/local/include
ln -s pruner/pruner.h pruner.h
```


