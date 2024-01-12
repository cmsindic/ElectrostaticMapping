import os
import ElectricField as ef
import pandas as pd
import requests
from argparse import ArgumentParser
from Residues import *
from fileProcessing.PDB import get_solvated_pdb
from spatial.cartesian import *
from time import time
from scipy.spatial import distance_matrix

parser = ArgumentParser(
    prog='Electrostatic Field Collector',
    description="""Map electrostatic fields in a pdb file
    from amber99sb forcefield without needing gromacs or other
    such software to ascertain charges from forcefield files""")

parser.add_argument('--filename', required=False)
parser.add_argument('--batch', default=False, action='store_true')
parser.add_argument('--interactive', default=False, action='store_true')
parser.add_argument('--quiet', default=False, action='store_true')
parser.add_argument('--exclude', default=False, action='store_true',
    help="Excludes same-molecule atoms from EF calculations")
parser.add_argument('--strict', default=False, action='store_true',
    help="Will not return if any object is not found in the RTP file")
parser.add_argument('--single_results', default=False, action='store_true',
    help="Will generate output for each individual water protons' fields")
parser.add_argument('--cutoff', type=float, required=False)
parser.add_argument('--save', action='store_true', required=False)
parser.add_argument('--solvated', action='store_true',required=False)
args = parser.parse_args()


def handle_output(bulk_solvent, in_enzyme, near_as, pdb_id, charges_coords):
    ''' Function to organize printing or writing of output data.
    bulk_solvent, in_enzyme, near_as are at this point, collections
    of water molecules partitioned by their locations.
    '''
    output_header_1 = "Average electric field strength on water protons:"
    output_header_2 = "Sample Partition  |  Field Strength (a.u.)  |  N Hydrogen"
    output_header_3 = "---------------------------------------------------------"

    partitions = {'Bulk Solvent    ': bulk_solvent,
              'Near Enzyme     ': in_enzyme,
              'Near Active Site': near_as}
              
    
    def impute_individual_fields(partitions):
        ''' To be used for turning 
        bulk_solvent, in_enzyme, near_as,
        which are collections of molecules,
        into collections of field strengths 
        felt by the hydrogen on each. 
        Update partitions accordingly. 
        '''
        for k, v in partitions.items():
            sele = list(v)
            sele_h = [h for water in sele for h in water.hydrogen]
            fields = ef.get_e_field_list(
                targets=sele_h,
                charges_coords=charges_coords,
                cutoff=4
                )
            partitions[k] = fields
            
        return partitions    
            
    
    def get_avg_magnitude(field_data, places=7):
        ''' Get field magnitude across partition,
        rounded to x places. 
        Used if single_results is False.
        '''
        sele = list(field_data)
        sele_h = [h for water in sele for h in water.hydrogen]
        avg_magnitude = ef.get_avg_e_field(
            targets=sele_h,
            charges_coords=charges_coords,
            cutoff=4
            )
            
        return round(avg_magnitude, places)
    
    
    def write_data(key, avg_field_mag):
        ''' Shortcut to write formatted data lines.
        '''
        return str(k + " | " + str(avg_magnitude) + 
               "             | " 
               + str(len(sele)) + '\n')
    
    # If detailed output about each proton is desired
    if args.single_results:
        partitions = impute_individual_fields(partitions)
        
        # If detailed info is to be saved
        if args.save:
            pwd = os.getcwd()
            os.chdir('output')
            outfile = pdb_id + '_verbose.out'
            if os.path.exists(outfile):
                os.chdir(pwd)
                return
            max_len = max([len(v) for v in partitions.values()])
            for k, v in partitions.items():
                for i in range(max_len):
                    try: 
                        v[i] = v[i]
                    except:
                        v.append(1000.000)
                partitions[k] = v
            with open(outfile, 'w') as f:
                f.write('Bulk, Enzyme, Active Site\n')
                i = 0
                B, E, A = partitions.values()
                while i < max_len:
                    b, e, a = B[i], E[i], A[i]
                    f.write(f"{b},{e},{a}\n")
                    i += 1
            os.chdir(pwd)
            
        # If detailed output about each proton is desired
        # but not saved.
        else:
            ...
         
    # If avg field info over partitions is desired and
    # the user wants it to be saved.
    elif args.save:
        pwd = os.getcwd()
        os.chdir('output')
        outfile = pdb_id + '.out'
        if os.path.exists(outfile):
            os.chdir(pwd)
            return  
        with open(outfile, 'w') as f:
            f.write(output_header_1 + "\n")
            f.write(output_header_2 + "\n")
            f.write(output_header_3 + "\n")
            for k, v in partitions.items():
                avg_magnitude = get_avg_magnitude(v)
                f.write(write_data(k, avg_magnitude))
                
        os.chdir(pwd)

    # If avg field info over partitions is desired and
    # the user wants it to be printed, not saved.
    else:
        print(output_header_1)
        print(output_header_2)
        print(output_header_3)
        for k, v in partitions.items():
            avg_magnitude = get_avg_magnitude(v)
            print(write_data(k, avg_magnitude))

    

def single_file_op(pdb=args.filename):
    ''' Get residues from solvated pdb if pdb is solvated using PACKMOL
    or some other method, else get residues from original pdb.
    Primary controller of data collection.
    '''
    t1 = time()
    pdb_id = os.path.split(pdb)[-1][:-4]
    rtp_file = os.path.join(".\\amber99sb-ildn.ff", "aminoacids.rtp")

    pdb = pdb if not args.solvated else get_solvated_pdb(main_pdb)
    model = ResiduesInFile(pdb, pdb_id, rtp_file)
    model.get_residues()

    # (charge, array([x, y, z])) for all atoms in pdb file
    cac = model.get_charges_and_coords(
        model.get_all_molecules(),
        quiet=args.quiet,
        strict=args.strict
    )

    ''' Collect active site information from original file
    given that the residue IDs of the active site are taken from
    the file before solvation.'''
    active_site = ActiveSite(pdb_id)
    active_site.get_residues(
        residues=model.get_all_molecules(),
        interactive=args.interactive
    )
    active_site.get_coords()

    model.convert_solvent_to_water_objects(exclude=args.exclude)

    in_enzyme = list(
        w for w in model.solvent
            if target_near_coords(
                target=w.oxygen.coord,
                coords=model.solute_coords,
                cutoff=4)
    )

    near_as = list(
        w for w in model.solvent
            if target_near_coords(
                target=w.oxygen.coord,
                coords=active_site.coords,
                cutoff=4)
    )

    bulk_solvent = [m for m in model.solvent if not m in in_enzyme]
    in_enzyme = [w for w in in_enzyme if not w in near_as]
    
    handle_output(
        bulk_solvent,
        in_enzyme,
        near_as,
        pdb_id,
        cac
    )

    #print("Run time:", time()-t1)
    
    
if args.batch:
    dirs = [d for d in os.scandir() if os.path.isdir(d.name)]
    model_dirs = [d for d in dirs if 'ases' in d.name]
    input_files = [f for d in model_dirs for f in os.scandir(d)]
    input_files = [os.path.abspath(f.path) for f in input_files]
    input_files = [f for f in input_files if not 'wH' in f]
    for f in input_files:
        try:
            single_file_op(f)
        except (AssertionError, ValueError, ZeroDivisionError,
                requests.exceptions.ConnectionError):
            pass
else:
    single_file_op()



