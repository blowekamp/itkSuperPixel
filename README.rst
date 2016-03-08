.. image:: https://circleci.com/gh/blowekamp/itkSuperPixel/tree/master.svg?style=svg
    :target: https://circleci.com/gh/blowekamp/itkSuperPixel/tree/master
.. image:: https://travis-ci.org/blowekamp/itkSuperPixel.svg?branch=master
    :target: https://travis-ci.org/blowekamp/itkSuperPixel

This ITK module provides an implimentation of the Simple Linear
Iterative Clustering (SLIC) superpixel segmentation algorithm. And it
is planned have more superpixel algoritms.


General
------

This is a module for ITK: The Insight Toolkit for Segmentation and
Registration. It is designed to work with the ITKv4 modular system and
to be placed in ITK/Modules/External.


Getting Started
---------------

This module should be cloned into the ITK reposiory as a subdirectory
in the "Modules/External" directory.

The following is a brief list of instructions to get a external module
into ITK:

cd ITK/Modules/External/
git clone https://github.com/blowekamp/itkSuperPixel.git

Then configure ITK as not make but set "Module_SuperPixel" to "ON" to
enable this module. The external module will need to be manually
updated from the git respository.


License
-------

This software is distributed under the Apache License. Please see
LICENSE for details.


Author
------


Bradley Lowekamp
