# roca-speed
Fast, multiprocess implementation in Python of the RoCa vulnerability of RSA keys

## Installation instructions
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

And finally clone this repo and run some tests:
```
git clone https://github.com/lvenuto/roca-speed
cd roca-speed
./test.py
```

### Linux

