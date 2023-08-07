############################################################################

############## GET ACTIVE SITE LOCATION FROM ORIGINAL PDB ##################


if args.activesite:
    from ActiveSite import ActiveSite
    active_site = ActiveSite(index_file, pdb)
    active_site.get_residues(residues)
    central_point = active_site.get_central_point()

    sphere_radii = [x/2 for x in range(2, 22, 1)]


    outdir = 'outfiles'
    outfile = os.path.join(outdir, pdbcode + '.csv')
    if args.save:
        with open(outfile, 'w') as f:
            f.write('Plot of E field mag vs distance from AS \n')
            f.write('r(A), E Field Mag (a.u.)\n')
            for r,field in enumerate(spheres.keys()):
                f.write(str(r) + ',' + str(field) + '\n')



############################################################################

############## GET E FIELD AT POINTS AROUND ACTIVE SITE ####################

import ElectricField


'''
from matplotlib import pyplot as plt
plot = plt.figure()
plt.scatter(spheres.keys(),avg_mags)
plt.show()
'''
