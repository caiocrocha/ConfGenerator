def get_cmd_line():
    parser = argparse.ArgumentParser(description = 'PDB and PSF reader, chain rotator.')
    parser.add_argument('--pdb', action = 'store',     dest = 'pdb', required = True, help = 'PDB file.')
    parser.add_argument('--psf', action = 'store',     dest = 'psf', required = True, help = 'PSF file.')
    parser.add_argument('--new_pdb', action = 'store',     dest = 'new_pdb', required = True, help = 'New PDB file.')
    parser.add_argument('--new_psf', action = 'store',     dest = 'new_psf', required = True, help = 'New PSF file.')
    arguments = parser.parse_args()
    # Assign from cmd line.
    return (arguments.pdb,
            arguments.psf,
            arguments.new_pdb,
            arguments.new_psf)
