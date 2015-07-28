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
```
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
```
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

