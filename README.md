# xDDx: Numerical Toolbox for Ultrasound Transducer Characterization and Diagnostics

For a quick start, all toolbox files are stored in **xDDx.zip**, available on the [releases page](https://github.com/pavrosni/xDDx/releases) in the **Assets** drop-down list below the release notes.

## Overview
xDDx is an open access transducer characterization toolbox for use in MATLAB or Octave on Windows computers. It consists of two parts. 

The first part, **Holography**, is designed for forward or backward projection of transient acoustic holography measurements to identify transducer surface defects and reveal structural details of the radiated acoustic field. The toolbox includes an automated procedure for correcting axes misalignments to optimize the visualization of transducer surface vibrations.

The second part, **Simulation**, utilizes the same projection algorithm as the first part to simulate fields radiated by user-defined transducer designs in a uniform medium without attenuation. 

The core algorithm is based on the Rayleigh integral implemented in C++ executables. The toolbox was developed in MATLAB and includes both a scripting interface and a graphical user interface (GUI) as a wrapper around the C++ executables. The xDDx algorithms were developed in two versions: for compute unified device architecture (CUDA-compatible graphics processing units – GPUs), and for central processing units (CPUs). The algorithms were specifically optimized for speed and can be used both for post-processing of planar scan data into holograms and for field projection calculations.

A detailed description of the toolbox’s capabilities is given in the xDDx User Manual and the xDDx publication:

> P. B. Rosnitskiy, O. A. Sapozhnikov, V. A. Khokhlova, W. Kreider, S. A. Tsysar, G. P. L. Thomas, K. Contreras, and T. D. Khokhlova, “xDDx: a Numerical Toolbox for Ultrasound Transducer Characterization and Design with Acoustic Holography,” _IEEE Trans. Ultrason., Ferroelectr., Freq. Control_ (Early Access), 2025.

## Installation 
The xDDx toolbox is located in the archive "xDDx.zip", which is available from the link https://github.com/pavrosni/xDDx/releases, "Assets" section

The archive contains a single folder with the same name, "xDDx". Users can unpack the archive into a folder of their choice and refer to the `xDDx\examples` folder to access the tools. 

As the toolbox consists of two parts, the tools for each part are organized into the following folders:
Holography Toolbox (`xDDx\examples\holography_toolbox`)
Simulation Toolbox (`xDDx\examples\simulation_toolbox`) 

The xDDx Toolbox has its tools sorted by categories, with certain subfolders present in the corresponding folder. For example:
`xDDx\examples\holography_toolbox\simple_projection_tools`.

The xDDx tools are accessible through scripts e.g. `xDDx\examples\holography_toolbox\quick_start_spherical.m`, and GUI tools are available to assist with some of the interactive steps. To use a tool, a corresponding script (*.m) should be opened and run in MATLAB/Octave environment. It is good practice to make a local copy of each example before using a tool for a particular case, to keep the original examples unchanged. 

Some users may benefit from reading the comments within a tool script before running it, as the comments summarize the tool's capabilities and might be sufficient to understand its input and output. Note that all examples are available “out of the box” and can be run immediately after unpacking the archive, with no further setup required.

## Documentation and Help
xDDx User Manual is available in the release Assets and the "xDDx.zip" archive. Please feel free to reach out to Pavel Rosnitskiy (pavrosni@gmail.com) with any questions about the manual and the toolbox, as well as to report any bugs.

## Quick Start
Users interested in a quick overview of the software's capabilities may start with:
1. The Quick Start example for the Holography Toolbox `xDDx\examples\holography_toolbox\quick_start_spherical.m`, and the 
2. The Single-Frequency Simulation Tool for the Simulation Toolbox `xDDx\examples\simulation_toolbox\transducer_simulation_sf.m`

These tools provide a general overview of the xDDx interface, its inputs, and outputs.

## License
If you find the xDDx toolbox useful for your work, please consider citing the xDDx paper:

> P. B. Rosnitskiy, O. A. Sapozhnikov, V. A. Khokhlova, W. Kreider, S. A. Tsysar, G. P. L. Thomas, K. Contreras, and T. D. Khokhlova, "xDDx: a Numerical Toolbox for Ultrasound Transducer Characterization and Design with Acoustic Holography," _IEEE Trans. Ultrason., Ferroelectr., Freq. Control_ (Early Access), 2025

and the acoustic holography paper:

> O. A. Sapozhnikov, S. A. Tsysar, V. A. Khokhlova, and W. Kreider, "Acoustic holography as a metrological tool for characterizing medical ultrasound sources and fields," _The Journal of the Acoustical Society of America_, vol. 138, no. 3. Acoustical Society of America (ASA), pp. 1515–1532, Sep. 01, 2015.

General license information is available in `license.txt`
