"""Helper module which represents a SDDS file with file_path, file_name, file_header, bpm_positions and bpm_data.

Note: After creating an instance you will have to call read_in() manually if you are dealing with an input file,
otherwise write_out() for writing data to file. Data is stored as numpy.array().
"""

import os
import numpy

class SDDS_file:
	def __init__(self, sdds_file_path):
		self.file_path = sdds_file_path
		self.file_name = os.path.basename(sdds_file_path)
		self.file_header = []
		
		self.bpm_data = {'X':{}, 'Y':{}}
		self.bpm_positions = {}
	
	def read_in(self):
		"""Reads in data from its file_path."""
		with open(self.file_path) as sdds_file:
			for line in sdds_file:
				if line.startswith('#'): # comment line
					self.file_header.append(line)
				elif line[0] in ['0', '1']: # data line
					line = line.split()
					
					if line[0] == '0': plane = 'X'
					else: plane = 'Y'
					bpm_name = line[1]
					self.bpm_positions[bpm_name] = float(line[2])
					self.bpm_data[plane][bpm_name] = numpy.array([float(d) for d in line[3:]])
				else:
					raise ValueError('Unknown line in SDDS file detected: %s' % line)
	
	def write_out(self):
		"""Writes out data in sdds format to its file_path."""
		with open(self.file_path, 'w') as sdds_file:
			# write header
			sdds_file.write(''.join(self.file_header))
			
			# write data
			for plane in self.bpm_data:
				if plane == 'X': plane_string = '0'
				else: plane_string = '1'
				for bpm in self.bpm_data[plane]:
					sdds_file.write('%s %s      %.5f  %s\n' % (plane_string, bpm, self.bpm_positions[bpm], '  '.join(format(d, '.5f') for d in self.bpm_data[plane][bpm])))
	
	def get_bpm_data(self):
		return self.bpm_data
	
	def get_bpm_position(self):
		return self.bpm_positions
	
	def get_file_header(self):
		return self.file_header
	
	def get_path(self):
		return self.file_path
	
	def get_file_name(self):
		return self.file_name