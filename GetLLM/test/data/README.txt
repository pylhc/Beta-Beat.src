This directory contains the test data for filecheck.py.
Each directory represents one run and will be validated.
A valid directory has following structure:
	-x [dir- name of the directory is indifferent]
	  -input [dir]
	    -model [dir]
	      -twiss.dat
	    -src_files [dir]
	      -'files from analysis'
	  -output [dir]
	    -to_check [dir-output from modified GetLLM.py]
	    -valid [dir-output from origin/valid GetLLM.py]

For new test files create a directory by using this structure.