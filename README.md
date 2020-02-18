[![DOI](https://zenodo.org/badge/140095242.svg)](https://zenodo.org/badge/latestdoi/140095242)

# DSP library
This project contains several different functions to apply DSP algorithms
for optical communications. These functions are suitable both for coherent and
non-coherent (PAM, DMT) optical communications.

Most of the functions are independent from each other. Therefore, functions
in this code can be easily used and combined with other DSP functions.

## Usage
* The list of functions, together with a short description, is
  in the [Contents.m](Contents.m) file.
* The description of input and output parameters is in the header of each
  function.
  
### Signals 
In general, input (and output) signals have time in the **first** dimension
(e.g. column vectors), and the second dimension is used to manage multiple
signals at a time (e.g. different polarizations, different parameters, ...). 

This convention has been chosen since MATLAB stores a matrix by keeping
the columns in contiguous parts of memory, therefore this convention is
*faster* than the other way around (row vectors).

### Parameters 
Most of the functions use a parameters structure for input parameters. 
The description of the parameters used in the function is usually found in
the header, and the default parameters for the coherent-dsp functions 
is in the [get_default_cohdsp_params.m](get_default_cohdsp_params.m) file.

## References
- [D. Pilori, "Advanced Digital Signal Processing Techniques for High-Speed Optical Links," Ph.D. Thesis, Mar 22, 2019.](https://hdl.handle.net/11583/2729814)

## License
This code is released under [MIT License](https://opensource.org/licenses/MIT).
