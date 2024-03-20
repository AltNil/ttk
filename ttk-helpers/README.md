# ttk-helpers

This project contains some helping elements like script when working with ttk

## LUTGenerator
Contains scripts helping to generate lookuptables for ttk. The scripts are:


### verifier.py
Verifier for triangulation comparison. Takes two parameters:

1. Relative path to an pvsm file. This files needs to contain the comparetriangulation filter
2. Option how the output should be formatted. Existing methods are:
    - Long -> Showing full diff (positions and values)
    - Short -> Showing which values differ
    - ExportLong -> Excluding coloring, otherwise same as Long
    - ExportShort -> Excluding coloring, otherwise same as Short

### generateLutByInput.py
Generates the luts by executing specified cases using three parameters

1. csv file coontaining:
    - PathToTestfile;gridsizeX,gridsizeY,gridsizeZ
2. Number of case distinctions per dimension (one number for all)
3. Number of dimension contained


### generateMethods.py
Generates an sceleton for the three methods as well as code for the compareTriangulation. Taking two Parameters:

1. Name of the structure getting
2. Name of the plural of the structure getting


### teststates
Teststates need to be generated new for each computer.
Namings: teststate_{number}.psvm where number has 3 digits and is the max dim coordinate.
Testsets are Cubes (Wavelet 0 till number for each dimension).
