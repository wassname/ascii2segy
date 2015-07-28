usage: ascii2segy_segpy.py [-h] [-o OUTFILE] [-s SAMPLE_RATE] [-m MAX_TIME]
                           [--dtype DTYPE] [-d SEISMIC_DIMENSIONS] [-e]
                           [-v VERBOSITY] [-t]
                           infile

ascii2segy

This script uses segpy-rewrite to convert velf ascii files to segy files. 
Please use 64 bit python for large files as data is loaded into memory.

Run with -h to get more information on parameters.

Outputs
-------
A segy file containing the velocities from the input ascii file.

Inputs
-------
Example input file format:

###############################################################################
#BEGIN HEADER
#X
#Y
#Z (time)
#FLOAT  "Stacking Velocity"
#FLOAT  "XLINE"
#FLOAT  "INLINE"
#END HEADER
   7618446.0      141544.0             0.0          1510.0       10040.0        2060    
   7618446.0      141544.0            40.0          1510.0       10040.0        2060
   7618446.0      141544.0            80.0          1510.0       10040.0        2060

Examples
--------
ascii2segy.py test.csv -o test.sgy
ascii2segy.py -h
ascii2segy.py test.csv -s 40 -m 5000 -c 15 -n 4 -d 2 -v 3 -e 1

Installation
-------
Make sure you have 64bit python 3.3 installation
Install the requirements listed in requirements.txt. Note that this uses 
segpy 2.0.0a2 from github.com/rob-smallshire/segpy which you can install using 
the command below: 
```
pip install git+git://github.com/sixty-north/segpy/tree/91562fddfd6d8424ee4161f4417982243512c150
```
In a command propmt with python64 3.3 install and test:
```
pip install -r requirements.txt
python ascii2segy ./test_data/test_input.tsv
```

License
------
This inherits it's license from segpy

More Info
-------
Created on Mon Aug 04 09:11:11 2014

positional arguments:
  infile                Input ascii velovity file in western (velf) format.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Optional output segy filename, e.g. outfile.sgy.
                        Defaults to infile.sgy. (default: None)
  -s SAMPLE_RATE, --sample-rate SAMPLE_RATE
                        If not provided this is autodetected from input data.
                        The velocities are resampled to this rate in
                        milliseconds. (default: None)
  -m MAX_TIME, --max-time MAX_TIME
                        If not provided this is autodetected from input data.
                        The maximum time of the segy in ms. (default: None)
  --dtype DTYPE         Datatype for segy {8: 'int8', 1: 'ibm', 2: 'int32', 3:
                        'int16', 5: 'float32'}
  -d SEISMIC_DIMENSIONS, --seismic-dimensions SEISMIC_DIMENSIONS
                        Type of seismic, options 2 or 3. (default: 3)
  -e, --use-extended-header
                        Use the segy extended header in the output, this will
                        contain the header from the text file. Defaults to
                        false. (default: False)
  -v VERBOSITY, --verbosity VERBOSITY
                        Print more information, can be 0, 1, or 2. (default:
                        0)
  -t, --test            Run sanity tests. (default: False)
To exit: use 'exit', 'quit', or Ctrl-D.
An exception has occurred, use %tb to see the full traceback.
