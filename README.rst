itkSuperPixel
=================================

.. |CircleCI| image:: https://circleci.com/gh/InsightSoftwareConsortium/itkSuperPixel.svg?style=shield
    :target: https://circleci.com/gh/InsightSoftwareConsortium/itkSuperPixel

.. |TravisCI| image:: https://travis-ci.org/InsightSoftwareConsortium/itkSuperPixel.svg?branch=master
    :target: https://travis-ci.org/InsightSoftwareConsortium/itkSuperPixel

.. |AppVeyor| image:: https://img.shields.io/appveyor/ci/blowekamp/itksuperpixel.svg
    :target: https://ci.appveyor.com/project/blowekamp/itksuperpixel

=========== =========== ===========
   Linux      macOS       Windows
=========== =========== ===========
|CircleCI|  |TravisCI|  |AppVeyor|
=========== =========== ===========

This ITK module provides an implimentation of the Simple Linear
Iterative Clustering (SLIC) superpixel segmentation algorithm.


Getting Started
---------------

This module should be cloned into the ITK reposiory as a subdirectory
in the "Modules/External" directory.

The following is a brief list of instructions to get a external module
into ITK:

cd ITK/Modules/External/
git clone https://github.com/blowekamp/itkSuperPixel.git

Then configure ITK as with cmake and set "Module_SuperPixel" to "ON" to
enable this module. The external module will need to be manually
updated from the git respository.


License
-------

This software is distributed under the Apache License. Please see
LICENSE for details.


Author
------


Bradley Lowekamp
