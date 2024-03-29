# OptogenSim: computational tool for Optogenetics light delivery design

Optogenetics is a powerful tool for in in-vivo brain studies due to its
unparalleled ability to modulate neural activity in living tissues and
animals. Precise temporal and spatial neuronal control can be achieved
through specific cell-type targeting by the use of optogenetic proteins
such as channelrhodopsin and optical stimulation at specific
wavelengths. Given the demonstrated power of this technique, there is
now great need to develop improved light delivery strategies that will
more accurately stimulate the neurons of interest while reducing the
nonspecific effects such as tissue heating or photodamage.

In order to help optimize light delivery strategies, OptogenSIM, a 3D
Monte Carlo simulation platform for optogenetic applications, has been
developed under an extensive multisite collaboration since August 2012.
The development partners in OptogenSIM are: Professor [Steven Jacques's
Lab](http://omlc.ogi.edu/~jacquess/)<sup>1</sup> at Oregon Health and Science
University, the [BIST Lab](https://pantherfile.uwm.edu/pashaie/www/)<sup>2</sup>
at the University of Wisconsin at Milwaukee, and [the LOCI group at
the University of Wisconsin at Madison](https://eliceirilab.org/).
OptogenSIM can be used for simulating light delivery in brain in a wide
variety of optogenetic applications. It aims to provide a 3D simulator
for predicting light distribution in heterogeneous brain tissue with
high resolution, including a voxel-based 3D Monte Carlo tool custom
built for optogentics light delivery, generic optical properties library
for different brain tissues, and a 3D mouse brain atlase segmented with
SPM-based tool<sup>3,4</sup>. This platform allows for delivering light
at any locations of the brain and with commonly used fiber geometries
and light source types, taking into account the light wavelength and
power, optical fiber diameter and numerical aperture(NA), brain tissue
optical properties, light beam type and tissue heterogeneity. Estimated
light density contours can show the region of any specified power
density in the 3D brain space and thus can help optimize the light
delivery settings. A demonstrated validation of mcxyz vs gold standard
MCML can be found
[here](https://loci.wisc.edu/files/loci/software/mcmlVSmcxyz_validation.pdf).

Other rodent brain atlases can also be incorporated into OptogenSIM. We
have tested the simulation on an average 3D mouse brain
atlas<sup>5</sup> and an average 3D rat brain atlas<sup>6</sup> but due
to copyright issues for direct distribution of the atlases with our
tool, these atlases have not been included in the current version.

We would like to thank Dr. Helene Benveniste at the Stony Brook
University for providing the permission to incorporate the 3D mouse
atlas in OptogenSIM. We also would like to thank Vaibhav D. Phad, a LOCI
graduate intern, for discussing and testing the light beam modeling.
This simulator was partially sponsored by the University of Wisconsin
Intercampus award under grant number UDDS B19-2510 and funding from the
Laboratory for Optical and Computational Instrumentation (LOCI).

## Main developers

-   [Yuming Liu](https://loci.wisc.edu/people/yuming-liu) (primary contact)
-   [Steven Jacques](http://omlc.org/~jacquess/)

## Download

- [OptogenSIM standalone for Windows 64](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM_Win64_mcr2014b.zip)
- [OptogenSIM standalone for MAC](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM_Mac_mcr2014b.zip)
- [OptogenSIM quick manual](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM.quickstart.manual10292015.pdf)
- [OptogenSIM source code](https://github.com/uw-loci/optogensim/archive/refs/tags/v1.0.zip)

## Instructions

### Standalone version

1. Make sure the Matlab 2014b MCR is installed properly. The MCR can be
   freely downloaded from the Mathworks, Inc at
  [Windows64 MCR 2014b](https://www.mathworks.com/supportfiles/downloads/R2014b/deployment_files/R2014b/installers/win64/MCR_R2014b_win64_installer.exe)
  or [Mac MCR 2014b](https://www.mathworks.com/supportfiles/downloads/R2014b/deployment_files/R2014b/installers/maci64/MCR_R2014b_maci64_installer.zip).

2. Download the [OptogenSIM standalone for Windows](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM_Win64_mcr2014b.zip)
   or [OptogenSIM standalone for Mac](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM_Mac_mcr2014b.zip).
   Unzip and then run the application `OptogenSIM` (To be noted, before
   opening `OptogenSIM`, if you ever ran the older release, make sure to
   delete/rename/move the old `OSGworking` folder mentioned in step 3 of
   the manual, otherwise the new release could not work properly.):

On Windows: double click on `OptogenSIM`

On Mac: right/control click on `OptogenSIM` &rarr; Show Package Contents
&rarr; Contents &rarr; MacOS &rarr; applauncher (right-click and choose
open)

### Source-code version

Download and unzip [source
code](https://github.com/uw-loci/optogensim/archive/refs/tags/v1.0.zip)
including the Matlab m-files, data, as well as a distribution of
`mcxyz.c`. Open and run the function `OptogenSIM` to launch the GUI.

**Download the [quick
manual](https://github.com/uw-loci/optogensim/releases/download/v1.0/OptogenSIM.quickstart.manual10292015.pdf)
to start.**

### References

<sup>1</sup> [3D Monte Carlo code at Steven Jacques's Lab](http://omlc.ogi.edu/software/mc/mcxyz/)  
<sup>2</sup> [Optogenetics researches at Ramin Pashaie's BIST Lab](https://pantherfile.uwm.edu/pashaie/www/research.htm)  
<sup>3</sup> [Statistic Parametric Mapping (SPM) software](http://www.fil.ion.ucl.ac.uk/spm/)  
<sup>4</sup> [SPMMOUSE software](http://spmmouse.org/)  
<sup>5</sup> [Average mouse brain atlas](http://www.spmmouse.org/)  
<sup>6</sup> [Average Rat Brain Atlas](http://www.idac.tohoku.ac.jp/bir/en/db/rb/)  
<sup>7</sup> [3D Mouse Brain Atlas](http://brainatlas.mbi.ufl.edu/Database/)
