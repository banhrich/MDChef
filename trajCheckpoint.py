#!/usr/bin/python
################################################################################
#
# General purpose checkpoint file: checks if an 'analysie file' already exists.
# If so, use an associated checkpoint file and obtain the beginning frame (bf).
# Returns beginnining (bf) and end (ef) frame numbers as integers.
# Input: Universe, Folder Name, Base File Name, Name Extension, DataType (npy, dat, txt)
# Output: [return] beginning frame, end frame
#
# Note:
# Load: from trajCheckpoint import checkpoint_frame
# Load: from trajCheckpoint import checkpoint_write
# * Before looping over the trajectory: checkpoint_frame
# * While looping over the trajectory: checkpoint_write
#
# by Richard in 2017 for Python 2.7 MDAnalysis 0.15.0
#
# Changelog (January 26, 2017):
# * Extracted this originally from rmsdvstime.py script and repurposed it.
################################################################################

import os
import sys
import MDAnalysis

# Description: Checks if there is a checkpoint file and returns
# beginning and end frames (bf, ef)
def checkpoint_frame(folder, base_filename, name, datatype):
    # u is the Universe (MDAnalysis)
    # base_filename is the name of the file for this analysis (not -n name argument)

    # Set default beginning and end frames, data file exist bool
    bf, ef, fex = None, None, False

    #Check if folder exists, if not create one
    if not os.path.exists(folder+"/checkpoints"):
        print "Creating directory: %s" % (folder+"/checkpoints")
        os.makedirs(folder+"/checkpoints")

    # Check if file exists; if so, check for checkpoint file
    if os.path.isfile(folder+'/'+base_filename+'.'+name+'.'+datatype):
        fex = True
        # Load checkpoint file and obtain last frame analyzed (+1 for next analysis)
        print "Data file already found. Looking at checkpoint (.cpt) file."
        bf = int(open(folder+'/checkpoints/'+base_filename+'.'+name+'.cpt').readlines()[0]) + 1
    else:
        # If not, a new file will be generated later.
        print "Generating new file. Starting from first frame t=0"
        #f = open(folder+'/'+base_filename+'.'+name+'.'+datatype, 'w')
        #f.close()

    # Print Beginning and End Frames to Console
    if bf != None:
        print "Beg frame: %s" % bf
    else:
        print "Beg frame: 0 (checkpoint frame not found)"
    if ef != None:
        print "End frame: %s" % ef
    else:
        print "End frame: not specified; will run until end of simulation"

    # Test Cases before Proceeding with the Code:
    #if bf == nframes:
    #    print "Terminal frame already reached %s." % nframes
    #    #sys.exit("Stop Signalled: Terminal frame already reached.")

    if bf >= ef and (bf != None and ef != None):
        sys.exit("Stop Signalled: Beginning frame %s is equal or greater to End frame %s" % (bf, ef))
    #print "Checks passed. Running Analysis. Let it rain numbers!"
    print '--------------------------------------'

    # Return Beginning and End Frames, File Exist Bool
    return bf, ef, fex

def checkpoint_write(frame, cpt_name):

    # Write last frame to checkpoint file
    with open(cpt_name, 'w') as o:
        o.write('{0:}'.format(frame))
