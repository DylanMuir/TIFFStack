# README #

This is the source repository for FocusStack, a Matlab toolbox for analysis of two-photon calcium imaging data. See the [FocusStack wiki](https://bitbucket.org/DylanMuir/twophotonanalysis/wiki) for more information.

### How do I get set up? ###

* Check out the repository to your local machine, using your git client of choice
* Add the base path of the repository to your Matlab path. You do not need to add sub-directories to the Matlab path
* FocusStack relies on several low-level accelerated disk access routines. These will be compiled automatically when first run, using "mex". This implies that "mex" must be able to compile C source files on your setup.

### How to contribute ###

* Firstly, thanks for your willingness to help! Documentation improvements, especially on the [FocusStack wiki](https://bitbucket.org/DylanMuir/twophotonanalysis/wiki), are especially useful. Bug fixes and feature contributions are an amazing contribution.
* If you would like to contribute a patch or bug fix, write to Dylan Muir (<dylan.muir@unibas.ch>) to request write access to the repository.
* Create a feature branch to encapsulate your code changes.
* After completing, documenting and testing your code, create a pull request for your branch to be integrated into the master branch.
* You may also of course fork the repository, work on your fork (under a feature branch, of course), and then create a pull request as above.

### Who do I talk to? ###

* In general the toolbox is supplied "as-is", with only moderate support
* Please use the [Issues tracker](https://bitbucket.org/DylanMuir/twophotonanalysis/issues) to report bugs and feature suggestions
* For questions or feedback, contact Dylan Muir (<dylan.muir@unibas.ch>)