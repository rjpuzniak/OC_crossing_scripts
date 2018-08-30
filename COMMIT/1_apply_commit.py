import os
import argparse

import numpy as np
import pylab as py

import nibabel as nib

def build_argparser():
    DESCRIPTION = "Modelling and filtering with COMMIT"
    p = argparse.ArgumentParser(description=DESCRIPTION)
    
    p.add_argument('dwi', help='Diffusion-weighted image (.nii)')
    p.add_argument('dwi_scheme', help='Gradient scheme (gradient table) created with fsl2scheme function from Camino (.scheme)')
    p.add_argument('peaks', help='Peaks of ODF extracted from Mrtrixs ODFs with sh2peaks (.nii.gz)')
    p.add_argument('tracks', help='Tractogram converted from Mrtrixs .tck (.trk)')
    p.add_argument('wm', help='White matter mask - i.e. extracted from 5TT image (.nii.gz)')
    p.add_argument('output', help='Location of output folder') # Later can set default value (current folder) and allow changing output location
    p.add_argument("-model", "--model_type", choices=["LiFE_Model", "StickZeppelinBall_Model"], help="Choice of the model that will be fitted to data")
    p.add_argument('-f', '--force', action='store_true', help='Overwrite existing output files')
    return p
  
def main():
  
    parser = build_argparser()
    args = parser.parse_args()

    # check for errors in input 
    try:
        dwi = nib.load(args.dwi)
    except:
        parser.error("Expecting DWI as first image")     
    # and so on

    # Essential part of COMMIT filtering
    
    # Create a dictionary for a tractogram
    from commit import trk2dictionary

    trk2dictionary.run(
	filename_trk = args.tracks, # 'LausanneTwoShell/fibers.trk',
	path_out = args.output, #'LausanneTwoShell/CommitOutput',
	filename_peaks = args.peaks, #'LausanneTwoShell/peaks.nii.gz',
	filename_mask = args.wm, #'LausanneTwoShell/WM.nii.gz',
	fiber_shift = 0.5,
	peaks_use_affine = True
    )
   
    # Precompute the rotation matrices used internally by COMMIT to create the lookup-tables for the response functions
    import commit
    commit.core.setup()  
    
    # Load the data and fit selected model
    if (args.model_type == 'StickZeppelinBall_Model'):
      
      mit = commit.Evaluation( '.', args.output )
      
      mit.CONFIG['doNormalizeSignal'] = False
      mit.CONFIG['doDemean'] = False
      
      mit.load_data( args.dwi, args.dwi_scheme )
      
      mit.set_model( 'StickZeppelinBall' )
      mit.model.set( 1.7E-3, [ 0.7 ], [ 1.7E-3, 3.0E-3 ] )
      mit.generate_kernels( regenerate=True )
      mit.load_kernels()

    
    if (args.model_type == 'LiFE_Model'):
      
      mit = commit.Evaluation( '.', args.output )  
      
      mit.CONFIG['doNormalizeSignal'] = False
      mit.CONFIG['doDemean'] = True
      
      mit.load_data( args.dwi, args.dwi_scheme )
      
      mit.set_model( 'StickZeppelinBall' )
      mit.model.set( 1.7E-3, [], [] )
      mit.generate_kernels( regenerate=True )
      mit.load_kernels()

    
    # Load in memory the sparse data-structure previously created with trk2dicitonary.run():
    mit.load_dictionary('.')

    # Now it's time to build the linear operator A to compute the matrix-vector multiplications for solving the linear system. This operator uses information from the segments loaded in the previous step and the lookup-tables for the response functions; it also needs to know the workload to be assigned to each thread durint the multiplications. To this aim, run the following commands:

    mit.set_threads(8)
    mit.build_operator()

    # fit the model (Stick-Zeppelin-Ball in this case) to the data
    mit.fit( tol_fun = 1e-3, max_iter = 500 )

    # Save results    
    #mit.save_results()
    if (args.model_type == 'StickZeppelinBall_Model'):
	mit.save_results(save_coeff = True)	#make sure it's saved in proper location
    if (args.model_type == 'LiFE_Model'):
	mit.save_results(save_coeff = True)

if __name__ == '__main__':
	main()    

