# Transparent lazy data access for matlab

By default, [Matlab][1] matrices must be fully loaded into memory. This can make allocating and working with
huge matrices a pain, especially if you only _really_ need access to a small portion of the matrix at a time.
[`memmapfile`][2] allows the data for a matrix to be stored on disk, but you can't access the matrix transparently
in functions that don't expect a [`memmapfile`][2] object without reading in the whole matrix. `MappedTensor` is
a matlab class that looks like a simple matlab tensor, with all the data stored on disk.

A few extra niceties over [`memmapfile`][2] are included, such as built-in per-slice access; fast addition,
subtraction, multiplication and division by scalars; fast negation; permution; complex support.

Tensor data is automatically allocated on disk in a temporary file, which is removed when all referencing
objects are cleared. Existing binary files can also be accessed. `MappedTensor` is a handle class, which means
that assigning an existing mapped tensor to another variable _will not_ make a copy, but both variables will point
to the same data. Changing the data in one variable will change both variables.

MappedTensor internally uses `mex` functions, which need to be compiled the first time MappedTensor is used. If
compilation fails then slower, non-mex versions will be used.

## Download and install

Download [MappedTensor][3]  
Unzip the `@MappedTensor` directory to somewhere on the [Matlab][1] path. The *@* ampersand symbol is important,
as it signals to [Matlab][2] that this is a class directory.

## Creating a MappedTensor object

    mtVariable = MappedTensor(vnTensorSize)
    mtVariable = MappedTensor(nDim1, nDim2, nDim3, ...)
    mtVariable = MappedTensor(strExistingFilename, ...)
    mtVariable = MappedTensor(..., 'Class', strClassName)

`vnTensorSize`, or [`nDim1 nDim2 nDim3 ...]` defines the desired size of the variable. By default, a new binary
temporary file will be generated, and deleted when the `mtVariable` is destroyed. `strExistingFilename` can be
used to map an existing file on disk, but the full size (and class) of the file must be known and specified in
advance. This file will not be removed when all handle references are destroyed.

By default the tensor will have class `double`. This can be specified as an argument to `MappedTensor`. Supported
classes: `char`, `int8`, `uint8`, `logical`, `int16`, `uint16`, `int32`, `uint32`, `single`, `int64`, `uint64`, `double`.

## Usage examples

    size(mtVariable)
    mtVariable(:) = rand(100, 100, 100);
    mfData = mtVariable(:, :, 34, 2);
    mtVariable(12) = 1+6i;

**Note**: `mtVariable = rand(100, 100, 100);` would over-write the mapped tensor with a standard matlab tensor!
To assign to the entire tensor you must use colon referencing: `mtVariable(:) = ...`

It's not clear why you would do this anyway, because the right hand side of the assignment would already allocate
enough space for the full tensor... which is presumably what you're trying to avoid.

Permute is supported. Complex numbers are supported (a definite improvement over [`memmapfile`][2]). `transpose`
(`A.'`) and `ctranspose` (`A'`) are both supported. Transposition just swaps the first two dimensions, leaving
the trailing dimensions unpermuted.

Unary plus (`+A`) and minus (`-A`) are supported. Binary plus (`A+B`), minus (`A-B`), times (`A*B`, `A.*B`) as
long as one of `A` or `B` is a scalar. Division (`A/B`, `A./B`, `B/A`, `B./A`) is supported, as long as `B` is a scalar.

Save and load is minimally supported — data is _not_ saved, but on load a new mapped tensor will be generated and
filled with zeros. Both save and load generate warnings.

Dot referencing (`A.something`) is not supported.

`sum(mtVar <, nDimension>)` is implemented internally, to avoid having to read the entire tensor into memory.

## Convenience methods

`SliceFunction`: Execute a function on the entire tensor, by slicing it along a specified dimension, and store the
results back in the tensor.

Usage: [`<mtnewvar>] = SliceFunction(mtVar, fhFunctionHandle, nSliceDim <, vnSliceSize,> ...)`

`mtVar` is a MappedTensor. This tensor will be sliced up along dimensions `nSliceDim`, with each slice passed
individually to `fhFunctionHandle`, along with any trailing argments (`...`). If no return argument is supplied, the
results will be stored back in `mtVar`. If a return argument is supplied, a new `MappedTensor` will be created to
contain the results. The optional argument `vnSliceSize` can be used to call a function that returns a different sized
output than the size of a single slice of `mtVar`. In that case, a new tensor `mtNewVar` will be generated, and it
will have the size `vnSliceSize`, with the dimension `nSliceDim` having the same length as in the original tensor `mtVar`.

"Slice assign" operations can be performed by passing in a function than takes no input arguments for `fhFunctionHandle`.

For example:

    	mtVar(:) = abs(fft2(mtVar(:, :, :)));

is equivalent to

    	SliceFunction(mtVar, @(x)(abs(fft2(x)), 3);

Each slice of the third dimension of `mtVar`, taken in turn, is passed to `fft2` and the result stored back into the
same slice of `mtVar`.

    	mtVar2 = SliceFunction(mtVar, @(x)fft2(x), 3);

This will return the result in a new, complex `MappedTensor` with temporary storage.

    	mtVar2 = SliceFunction(mtVar, @(x)sum(x), 3, [1 10 1]);

This will create a new `MappedTensor` with size [`1 10 N]`, where `N` is the length along dimension 3 of `mtVar`.

    	SliceFunction(mtVar, @()(randn(10, 10)), 3);

This will assign random numbers to each slice of `mtVar` independently.

    	SliceFunction(mtVar, @(x, n)(x .* vfFactor(n)), 3);

The second argument to the function is passed the index of the current slice. This line will multiply each slice in
mtVar by a scalar corresponding to that slice index.

## Publications

This work was published in [Frontiers in Neuroinformatics][4]: DR Muir and BM Kampa. 2015. [_FocusStack and StimServer:
A new open source MATLAB toolchain for visual stimulation and analysis of two-photon calcium neuronal imaging data_][5],
**Frontiers in Neuroinformatics** 8 _85_. DOI: [10.3389/fninf.2014.00085](http://dx.doi.org/10.3389/fninf.2014.00085).
Please cite our publication in lieu of thanks, if you use this code.

[1]: http://www.mathworks.com
[2]: http://www.mathworks.com/help/techdoc/ref/memmapfile.html
[3]: /resources/code/MappedTensor.zip
[4]: http://www.frontiersin.org/neuroinformatics
[5]: http://dx.doi.org/10.3389/fninf.2014.00085
