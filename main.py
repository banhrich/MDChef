# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""MDAnalysis: One Pass, Checkpoint, Multi-Analysis --- :mod:`mdchef`
=====================================================================

:Author: Customized by Richard Banh
:Year: 2010--2013 (Customized on January 30, 2017)
:Copyright: GNU Public License v3

Changelog
---------------------------------------------------------------------
January 31, 2017 by Richard Banh:
* Created file
---------------------------------------------------------------------

The module contains classes and functions to execute multiple analyses through
one pass of the trajectory while recording checkpoint frames in cpt files
(located in the directory /checkpoints) through each iteration.

The user needs to set a variable equal to 'restaurant' and provide the
appropriate variables. Afterwards, the user adds recipes (analyses) to the menu.
This is followed by serving (running) the analysis!

    >>> MD = restaurant(gro,xtc,folder,name).
    >>> MD.menu.add(name, recipe)
    >>> MD.serve()

The anylysis script is a class with 3 main functions: __init__, loopcode, and
postloopcode. Code written here will be executed before, during, and after
iterating over the trajectory, respectively. View code in recipes.py for
examples.

Custom recipe class structure:
    >>> class general(recipe):
    >>>
    >>>     def __init__(self, u, nframes,
    >>>                  folder, base_filename, name, datatype='pkl',
    >>>                  bin_size=120, debug=False):
    >>>         print "==== %s %s ====" % (base_filename, name)
    >>>
    >>>         # Set Analysis Variables: Core structure between all analyses
    >>>         recipe.__init__(self, folder=folder, base_filename=base_filename, name=name,
    >>>                           datatype=datatype, debug=debug)
    >>>
    >>>         # Code Before Iterating Over the Trajectory
    >>>         ...
    >>>
    >>>     def loopcode(self, frame, time):
    >>>         # Code During Iteration
    >>>         ...
    >>>
    >>>     def postloopcode(self):
    >>>         # Code After Iteration
    >>>         None


Fill Tutorial Below
--------------------

Functions
---------

.. class:: menu
.. class:: recipe
.. class:: restaurant

Helper functions
----------------

The following functions are used by the other functions in this
module. They are probably of more interest to developers than to
normal users.

.. autofunction:: checkpoint_frame
.. autofunction:: checkpoint_write

"""


# Import Modules
import os
import sys
import numpy as np
import pandas as pd
import MDAnalysis

# Custom Modules
from mdchef.align_custom import alignto
from mdchef.trajCheckpoint import checkpoint_frame, checkpoint_write

class menu:
    r"""Returns menu class. Used to add recipes (analyses).

    Parameters
    ----------
    None

    Returns
    -------
    menu class

    Example
    -------
    insert example ...

    Notes
    -----
    insert notes ...

    See Also
    --------
    other scripts to see also ...

    """
    def __init__(self):
        self.d = {}

    def __getitem__(self, i):
        return self.d[i]

    def __iter__(self):
        return self.d.itervalues()

    def add(self, name, recipe):
        '''Description: add recipe (analysis) to the menu'''
        self.d[name] = recipe

    def remove(self, name):
        '''Description: remove recipe (analysis) from the menu. Requires name.'''
        del self.d[name]

    def list_show(self):
        '''Description: print a list of recipes (analyses) on the menu'''

        print "=============================================="
        print "List of Analyses"
        print "=============================================="
        for i, name in enumerate(self.d):
            print "(%s): Analysis Name: %s" % (i+1, name)
        print ""

    def check_cases(self, nframes):
        '''Description: check to see if any recipes(analyses) need to be removed
        due to beginning frame errors/limits. If no analyses are found after
        this check, the script will terminate.'''

        print "=============================================="
        print "Running test cases on each analysis"
        print "=============================================="
        # Check each analysis and delete if bf is over nframes
        del_later = [] # list of names
        for name in self.d:
            #print "%s: (1) bf value: %s  (2) nframes: %s" % (name, self.d[name].bf, nframes)
            if self.d[name].bf >= nframes:
                print "Removed %s from analysis. %s frames. " % (name, self.d[name].bf)
                print " -- Reason: beginning frame exceeded the length of the trajectory (%s frames)" % nframes
                del_later.append(name)
        if len(del_later) > 0:
            for i in range(len(del_later)):
                self.remove(del_later[i])
        print ""

        # Check if there are any menu left, if not stop the script
        if len(self.d) == 0:
            sys.exit("Exit Signalled: There are zero analyses that may run.")

    def min_beg_frame(self):
        '''Description: return the minimum beginning frame (bf) value from all
        analyses.'''
        bf_list = []
        if len(self.d) > 0:
            self.min_bf = None
            for name in self.d:
                bf_list.append(self.d[name].bf)
            self.min_bf = min(bf_list)
            return self.min_bf, bf_list
        else:
            return "Attempted to obtain smallest beginning frame, but there are no analyses."

class recipe:
    def __init__(self, folder, base_filename, name, datatype, on=False, debug=False):
        # Set Variables
        self.folder = folder
        self.base_filename = base_filename
        self.name = name
        self.datatype = datatype
        self.on = on
        self.debug = debug

        # File Directory + Names + File Extensions
        self.data  = self.folder+'/'            +self.base_filename+'.'+self.name+'.'+ datatype
        self.chkpt = self.folder+'/checkpoints/'+self.base_filename+'.'+self.name+'.'+'cpt'

        # Load Cpt File (or set defaults)
        self.cpt_load()

    def cpt_load(self):
        '''Description: check if data file exists, if so check the
        checkpoint(cpt) file associated and return the beginning and end frames
        (bf, ef) and True/False for if the data file exists.'''
        # Check if the file exists, if so check the checkpoint(cpt) file and obtain the beginning frame
        # Beginning frame (bf): int; End frame (ef): int; File Exists (fex): bool
        self.bf, self.ef, self.fex = \
        checkpoint_frame(self.folder, self.base_filename, self.name, self.datatype)

    def cpt_out(self, frame):
        '''Description: write last frame to the checkpoint (cpt) file.'''
        # Write last frame to checkpoint file
        checkpoint_write(frame, self.chkpt)

    def runinloop(self, frame, time, func):
        '''Description: executes code during trajectory for every iteration.
        Writes a checkpoint file for each step.'''
        # Run Function
        func(frame,time)
        # Write last frame to checkpoint file
        self.cpt_out(frame)

    def runafterloop(self, func=None):
        '''Description: executes code after trajectory iteration. Function
        parameter is optional.'''
        if func != None:
            func()
        print "Done %s_%s" % (self.base_filename,self.name)

class restaurant:
    def __init__(self, gro, xtc, folder, name, gro_ref=None):
        # Optional Parameters: Defaults
        if gro_ref is None:
            gro_ref = gro

        # Check and Create Directory
        if not os.path.exists(folder):
            print "Creating directory: %s" % (folder)
            os.makedirs(folder)

        # Folder and Name Variables for Files
        self.folder = folder
        self.name = name

        # Load Gro and Xtc
        self.gro, self.xtc  = gro, xtc
        self.u = MDAnalysis.Universe(gro, xtc) # Mobile
        self.u_ref = MDAnalysis.Universe(gro_ref) # Reference
        self.nframes = len(self.u.trajectory)

        # Default Beginning(b) and End(e) Frame Values:
        self.b, self.e = None, None

        # Set Menu/Analyses
        self.menu = menu()

    def serve(self):
        # Check analyses for permitted bf values
        self.menu.list_show()
        self.menu.check_cases(self.nframes)
        self.menu.list_show()
        self.b, self.bf_list = self.menu.min_beg_frame()

        print "=============================================="
        print "All checks passed. Running all analyses!"
        print "=============================================="
        print "Beg Frame: %s; from bf values of [%s] " % (self.b, ', '.join(sorted(map(str,self.bf_list))))
        print "End Frame: %s" % self.e

        # Loop Over Trajectory
        for ts in self.u.trajectory[self.b:self.e]:
            # Loop over each analysis and perform analysis for this iteration

            # Counter for time frame (ts of 0 = None)
            frame = ts.frame if ts.frame > 0 else None

            for a in self.menu:
                if a.bf == frame: a.on = True # Turn on analysis
                if a.on: a.runinloop(ts.frame, ts.time, a.loopcode)
        print ''

        # Post Processing
        for a in self.menu:
            a.runafterloop(a.postloopcode)

    def get_input_files(self):
        '''Description: prints input files used''''
        print "GRO: " + self.gro; print "XTC: " + self.xtc
        print "GRO_REF: " + self.gro_ref

    def universe(self):
        '''Description: returns MDAnalysis.Universe (for mobile)'''
        return self.u
        
    def universe_ref(self):
        '''Description: returns MDAnalysis.Universe for reference'''
        return self.u_ref
