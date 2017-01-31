class rmsd_vs_time(recipe):
    '''
    Description: Calculate RMSD: mobile vs reference
    '''
    def __init__(self, ref, mob, nframes,
                 folder, base_filename, name, datatype='dat',
                 keep_original=True, selection='protein and name CA',
                 debug=False):
        print "==== %s %s ====" % (base_filename, name)

        # Set Analysis Variables: Core structure between all analyses
        recipe.__init__(self, folder=folder, base_filename=base_filename, name=name,
                          datatype=datatype, debug=debug)

        # Specific analysis variables
        self.keep_original = keep_original
        self.selection = selection # Default: 'protein and name CA'

        # Universe and Atom Groups
        self.ref = ref.select_atoms(self.selection) # Reference: protein backbone CA atoms
        self.mob = mob.select_atoms(self.selection) # Mobile: protein backbone CA atoms

    def loopcode(self, frame, time):
        # Calculate RMSD and Write to File
        if self.debug: print "Before RMSD: ", self.mob.positions[0]
        old, new = alignto(mobile=self.mob, reference=self.ref,
                           select=self.selection, keep_original=self.keep_original)
        if self.debug: print "After RMSD:  ", self.mob.positions[0]
        # Write Frame , Old RMSD, New RMSD to file
        if self.debug: print frame, old, new
        with open(self.data, 'a') as opf:
            opf.write('{0:}\t'.format(frame))
            opf.write('{0:}\t'.format(time))
            opf.write('{0:}\t'.format(old))
            opf.write('{0:}\n'.format(new))

    def postloopcode(self):
        None
