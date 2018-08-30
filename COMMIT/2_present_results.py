import os
import argparse

import numpy as np

import nibabel as nib

import numpy as np
from dipy.io.streamline import load_trk, save_trk
from dipy.tracking.streamline import Streamlines

def build_argparser():
	DESCRIPTION = "Modelling and filtering with COMMIT"
	p = argparse.ArgumentParser(description=DESCRIPTION)
	p.add_argument('subject', help='ID of subject being processed, required for naming the output')
	p.add_argument('tracks_type', help='Tractography variant used, combination of DT/CSD and seed/select') 
	p.add_argument('model_type', help='Model fitted to data, either SZB or LiFE') 
	p.add_argument('input', help='Location of input folder with results obtained from COMMIT framework')
	p.add_argument('-f', '--force', action='store_true', help='Overwrite existing output files')
	return p
  
def main():
	  
	parser = build_argparser()
	args = parser.parse_args()

	# Load the streamlines from COMMIT
	streamlines_name= args.input + '/dictionary_TRK_fibers.trk'
	print(streamlines_name)

	streams, hdr = load_trk(streamlines_name)
	streamlines = Streamlines(streams)

	# Load weights
	weights_name= args.input + '/Results_StickZeppelinBall/xic.txt'

	# Check lines counts of both text files
	def file_len(fname):
	    with open(fname) as f:
		 for i, l in enumerate(f):
		     pass
	    return i + 1
	print "Number of lines (streamlines) for file with SZB (COMMIT framework) weights is {}".format(file_len(weights_name))

	for j, k in enumerate(streamlines):
		pass
	print j+1
	print "Number of streamlines from COMMIT framework is".format(j)

	# Iterate simultaneously through both files and extract non-zero weights and corresponding streamlines
	from itertools import izip

	filtered_streamlines=[]
	with open(args.input + '/' + args.subject + '_SZB_filtered_weights.txt', "w+") as fil_wgh:
		for line_from_weights, line_from_streamlines in izip(open(weights_name), streamlines):
			if float(line_from_weights)>0:
				filtered_streamlines.append(line_from_streamlines)
				fil_wgh.write(line_from_weights)

	save_trk(args.input + '/' + args.subject + '_' + args.tracks_type + '_' + args.model_type + '_streamlines.trk', filtered_streamlines, np.eye(4))

if __name__ == '__main__':
	main()    

