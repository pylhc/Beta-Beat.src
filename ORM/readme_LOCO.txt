Main script for executing LOCO analysis is LOCO.py. This script will generate the "model" orbit response matrix, generate the Jacobian matrix, and find the model calibrations that minimize the difference between the measured and model orbit response matrices.

The fitting process is iterated some number of times (currently set to 4).

The resulting calibrations for nth iteration are written to results/AllCalib_n.py (results are written for each iteration for troubleshooting/checking convergence).

The example data file for measured orbit response is measuredORM.dat. The columns to be included are:
	NAME = corrector-bpm pair (same device names as used in MADX)
	S = longitudinal position of bpm	
	MX = (measured) x orbit position at bpm
	MY = (measured) y orbit position at bpm
	dMX = (measured) dx/dtheta for bpm/dipole pair
	dMY = (measured) dy/dtheta for bpm/dipole pair
