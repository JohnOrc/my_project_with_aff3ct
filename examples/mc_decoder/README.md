# How to compile this example

Make sure to have done the instructions from the `README.md` file at the root of this repository before doing this.

Copy the cmake configuration files from the AFF3CT build

	$ mkdir cmake && mkdir cmake/Modules
	$ cp ../../lib/aff3ct/build/lib/cmake/aff3ct-*/* cmake/Modules

Compile the code on Linux/MacOS/MinGW:

	$ mkdir build
	$ cd build
	$ cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-funroll-loops -march=native"
	$ make

The source code of this mini project is in `src/main.cpp`.
The compiled binary is in `build/bin/my_project`.

Then cd directory `../build`. Set run.sh permission

    $ chmod +755 ./run.sh

Run and debug

    $ ./run.sh

## Decoder_polar_SCL_mcfast_sys: variables

**path**

	vector L, delete_path()

**metrics**

	vector L

**l**
	
	vector L, N. llr values

**s**

	vector L, N. partial sums

**metrics_vec**

	vector[0] for Rep	2 * L
	vector[1] for Rate1	4 * L
	vector[2] for SPC	

**dup_count**

	vector L, 0. erase_bad_path()

**bit_flips**

	vector 4 * L

**is_even**

	vector L

**best_path**

	_store(), select_best_path()

**n_active_paths**

	init_buffers(), duplicate_tree(), delete_path()

**n_array_ref**

	init_buffers(), up_ref_array_idx(), duplicate_tree()

**path_2_array**

	vector L, m. init_buffers()

**sorter**

	N

**best_idx**

	vector L

**l_tmp**

	vector N

* * *

## Decoder_polar_SCL_mcfast_sys: functions

### Decoder_polar_SCL_mcfast_sys: constructor

### ~Decoder_polar_SCL_mcfast_sys

### clone

### set_frozen_bits

### get_frozen_bits

### init_buffers

### _decode

### _decode_siho

### _decode_siho_cw

### recursive_decode

### _store

### _store_cw

### update_paths_r0

### update_paths_r1

### flip_bits_r1

### update_paths_rep

### update_paths_spc

### flip_bits_spc

### delete_path

### select_best_path

### up_ref_array_idx

### erase_bad_paths

### duplicate_tree


# add glog
> git@github.com:google/glog.git

use debug mode 

	$ cmake .. -DCMAKE_BUILD_TYPE="Debug"
	