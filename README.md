pyDAS is developed by [Frederick Massin](https://github.com/FMassin) and Pascal Edme at the [Swiss Seismological Service](http://seismo.ethz.ch) of the ETHZ.

# What is pyDAS?
Currently, pyDAS is a set of Python functions for converting and managing data related to DAS. Here's a brief breakdown:

1. **Dependencies**: The code relies on libraries such as `numpy`, `h5py`, and `obspy`, which are used for numerical operations, handling HDF5 files, and dealing with seismic data formats, respectively.

2. **Default Templates and Schema**: Some default XML schema and binding templates are defined for SeisComP, a seismology software system.

3. **Functions**:
   - **`h5toseed`**: Converts seismic data stored in HDF5 files to the SEED format (a standard for storing seismic data), while also generating associated metadata in XML formats. This function also handles sensitivity scaling, data reading, and organizes data into seismic stations and channels.
   - **`seedtosc`**: Converts metadata from FDSNXML to SCXML format, merges it with existing SeisComP XML files, and manages configurations.
   - **`scexec`**: Executes a SeisComP module command with provided configurations, useful for handling seismic data processing via SeisComP.

4. **Purpose**: The code is intended to automate and streamline the process of managing DAS seismic data, converting formats, and preparing data for analysis and archiving, especially in environments that utilize SeisComP.

# Testing pyDAS
To test pyDAS, you can run [the Jupyter Notebook](pydas_test.ipynb) located in the root directory. This notebook will guide you through an example on how to use pyDAS and its functions.