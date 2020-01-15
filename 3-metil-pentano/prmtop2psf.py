import parmed as pmd

# Load an Amber topology file
parm = pmd.load_file('3-metil-pentano.prmtop', xyz='3-metil-pentano.rst7')
parm.write_psf('3-metil-pentano.psf')
