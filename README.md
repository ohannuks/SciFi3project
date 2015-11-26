Otto's and Janne's scifi3 project

Otto:
Note; the master branch is no longer supported. This code is meant to be a project for a course at our University, so some of the files are a bit scattered around as we only worked on the project for ~1 month (2 months on SAMRAI). Having checked this 26.11.2015, I think the currently supported branch is the 'doc' branch.

The code solves Poisson equation using the Particle-In-Cell (PIC) approach in a fully parallelized and vectorized mesh framework, employing Uintah framework's novel parallelization techniques. We have ran the code on the Voima supercomputer @FMI up to 70 cores (scales perfectly); we have tested the MPI support but currently PThreads have not been tested.

A few notes: If you're curious about using Uintah these may be helpful example files to demonstrate what one can do in the framework. The code does NOT compile out of the box; please see about installing both UINTAH and VISIT visualization software on your computer before running the code. The hierarchy and overall look has been largely neglected due to time issues I and Janne had, but the coding standard should be clean. Enjoy!

Oh and Ps. we did not remove the licenses from the files, and we claim no credits for Uintah (University of Utah) or VisIt (Lawrence Livermore National Laboratory):
http://uintah.utah.edu/
https://wci.llnl.gov/simulation/computer-codes/visit/

