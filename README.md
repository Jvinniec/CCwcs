# CCwcs

DESCRIPTION
---------------------
A C++ wrapper library for wcstools. This library introduces a single class call 'WcsHandler' that should function as a c++11 compliant wrapper for the wcstools 'WorldCoor' structure. 

Obtaining the code
---------------------
The code can be checked out from the git repository by submitting the follwoing on the command line:
	
	$ git clone https://github.com/Jvinniec/CCwcs.git

This will create a directory called "CCwcs" in the directory that you execute the above command from.

BUILDING
---------------------
Building should be as simple as:

	 $ ./configure [--with-wcs=/path/to/wcstools] [--prefix=/installation/directory]
	 $ make
	 $ make install

Alternatively, if there is no 'configure' file, you can run the included autogen.sh script.

REPORTING ISSUES
-
You can report issues that you come across in the github repository at https://github.com/Jvinniec/CCwcs/issues

DISCLAIMER
---------------------
This code is developed on top of the wcstools which can be downloaded
at the following url:

http://tdc-www.harvard.edu/software/wcstools/

Note that no rights or guarantees for CCwcs or wcstools are provided by
the CCwcs developer(s). However, please post an issue if you notice
either the code which exists has an error or if you see a feature
missing which you would like to have implemented. 

Or ... I guess you could always clone the code and implement it yourself :D
