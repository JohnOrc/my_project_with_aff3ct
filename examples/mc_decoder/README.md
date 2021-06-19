# How to compile this example

Make sure to have done the instructions from the `README.md` file at the root of this repository before doing this.

Copy the cmake configuration files from the AFF3CT build

	$ mkdir cmake && mkdir cmake/Modules
	$ cp ../../lib/aff3ct/build/lib/cmake/aff3ct-*/* cmake/Modules

Compile the code on Linux/MacOS/MinGW:

	$ mkdir build
	$ cd build
	$ cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-funroll-loops -march=native"
	$ make

The source code of this mini project is in `src/main.cpp`.
The compiled binary is in `build/bin/my_project`.

Then cd directory `../build`. Set run.sh permission

    $ chmod +755 ./run.sh

Run and debug

    $ ./run.sh