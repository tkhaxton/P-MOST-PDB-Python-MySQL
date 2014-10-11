#!/usr/bin/python

import MySQLdb as mdb
import sys
import os
from numpy import matrix
from numpy import linalg

amino_acid_directory = '../amino_acid_pdb_files'
site_atoms_file='../parameter_files/site_atoms.txt'
site_director_atoms_file='../parameter_files/site_director_atoms.txt'
cterminus_atoms_file='../parameter_files/Cterminus_atoms.txt'
elements_file='../parameter_files/elements.txt'
corrections_atoms_file='../parameter_files/corrections_atoms.txt'

#pdb_directory='../../PDB_files/Examples'
pdb_directory=sys.argv[1]

def create_site_director_table(cur, filename):
    #   Create table of atom names defining each coarse-grained site
    print "Creating site director table..."
    cur.execute("CREATE TABLE site_director_atoms(resname VARCHAR(3), siteid INTEGER, sitetype INTEGER, atom1 VARCHAR(4), atom2 VARCHAR(4), atom3 VARCHAR(4))")
    with open(filename) as f:
        for x in f:
            resname = x[0:3]
            nsites = int(x[3])
            for i in range(nsites):
                sitetype = int(x[4 + i * 13])
                atom1 = x[5 + i * 13:9 + i *13]
                atom2 = x[9 + i * 13:13 + i *13]
                atom3 = x[13 + i * 13:17 + i *13]
                cur.execute("INSERT INTO site_director_atoms (resname, siteid, sitetype, atom1, atom2, atom3) VALUES(%s, %s, %s, %s, %s, %s)", (resname, i, sitetype, atom1, atom2, atom3))

def create_site_table(cur, filename):
    #   Create table of atoms defined by residue, associated coarse-grained site, and atom name
    print "Creating site table..."
    cur.execute("CREATE TABLE site_atoms(resname VARCHAR(3), atomname VARCHAR(4), siteid INTEGER)")
    with open(filename) as f:
        for x in f:
            resname = x[0:3]
            natoms = (len(x) - 3) / 5
            for i in range(natoms):
                atomname = x[3 + 5 * i:7 + 5 * i]
                siteid = int(x[7 + 5 * i])
                cur.execute("INSERT INTO site_atoms (resname, atomname, siteid) VALUES(%s, %s, %s)", (resname, atomname, siteid))

def create_cterminus_atoms_table(cur, filename):
    #   Create table of atom names that may be associated with a C-terminal site
    print "Creating C terminus table..."
    cur.execute("CREATE TABLE cterminus_atoms(atom1 VARCHAR(4), atom2 VARCHAR(4), atom3 VARCHAR(4))")
    with open(filename) as f:
        for i, x in enumerate(f):
            if i==0:
                atom1 = x[0:4]
            if i==1:
                atom2 = x[0:4]
            if i==2:
                atom3 = x[0:4]
    cur.execute("INSERT INTO cterminus_atoms (atom1, atom2, atom3) VALUES(%s, %s, %s)", (atom1, atom2, atom3))

def create_amino_acid_atom_table(cur, directory):
    #   Create table of atoms associated with their amino acid residue
    print "Creating amino acid table..."
    cur.execute("CREATE TABLE amino_acid_atom(resname VARCHAR(3), atomname VARCHAR(4), serial INT)")
    cur.execute("CREATE TABLE connect(resname VARCHAR(3), selfserial INT, serial1 INT, serial2 INT, serial3 INT, serial4 INT)")
    for root, dirs, filenames in os.walk(directory):
        for file in filenames:
            with open(directory + '/' + file,'r') as f:
                for x in f:
                    if(len(x) >= 6):
                        record = x[0:6]
                        if record == 'ATOM  ':
                            serial=int(x[6:11])
                            resname=x[17:20]
                            atomname=x[12:16]
                            cur.execute("INSERT INTO amino_acid_atom (resname, atomname, serial) VALUES(%s, %s, %s)", (resname, atomname, serial))
                        elif record == 'CONECT':
                            linelength = len(x)
                            if (linelength - 12) % 5 != 0:
                                print 'CONECT line length %i!' % (linelength)
                                exit(1)
                            selfserial = int(x[6:11])
                            serial1 = int(x[11:16])
                            if linelength > 17:
                                serial2 = int(x[16:21])
                                if linelength > 22:
                                    serial3 = int(x[21:26])
                                    if linelength > 27:
                                        serial4 = int(x[26:31])
                                    else:
                                        serial4 = None
                                else:
                                    serial3, serial4 = None, None
                            else:
                                serial2, serial3, serial4 = None, None, None
                            cur.execute("INSERT INTO connect (resname, selfserial, serial1, serial2, serial3, serial4) VALUES(%s, %s, %s, %s, %s, %s)", (resname, selfserial, serial1, serial2, serial3, serial4))
                #   Extra H1 and H3 atoms for protonated amino terminus
                cur.execute("INSERT INTO amino_acid_atom (resname, atomname, serial) VALUES(%s, %s, %s)", (resname, ' H1 ', serial + 1))
                cur.execute("INSERT INTO amino_acid_atom (resname, atomname, serial) VALUES(%s, %s, %s)", (resname, ' H3 ', serial + 2))
                #   Connect to H1 and H3 to N
                cur.execute("INSERT INTO connect (resname, selfserial, serial1) VALUES(%s, (SELECT serial FROM amino_acid_atom WHERE resname = %s AND atomname = %s), %s)", (resname, resname, ' N  ', serial + 1))
                cur.execute("INSERT INTO connect (resname, selfserial, serial1) VALUES(%s, (SELECT serial FROM amino_acid_atom WHERE resname = %s AND atomname = %s), %s)", (resname, resname, ' N  ', serial + 2))
                cur.execute("INSERT INTO connect (resname, selfserial, serial1) VALUES(%s, %s, (SELECT serial FROM amino_acid_atom WHERE resname = %s AND atomname = %s))", (resname, serial + 1, resname, ' N  '))
                cur.execute("INSERT INTO connect (resname, selfserial, serial1) VALUES(%s, %s, (SELECT serial FROM amino_acid_atom WHERE resname = %s AND atomname = %s))", (resname, serial + 2, resname, ' N  '))
    cur.execute("CREATE TABLE amino_acid(resname VARCHAR(3)) SELECT DISTINCT resname FROM amino_acid_atom")

def create_pdb_table(cur, directory):
    #   Create table of atomic coordinates input from Protein Data Bank (PDB) files
    print "Creating PDB table..."
    cur.execute("CREATE TABLE pdb_onefile(entry VARCHAR(11), chainid VARCHAR (1), resseq INT, icode VARCHAR(1), rescounter INT, resname VARCHAR(3), atomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, siteid INT, PRIMARY KEY(entry, chainid, rescounter, atomname))")
    cur.execute("CREATE TABLE pdb(entry VARCHAR(11), chainid VARCHAR (1), resseq INT, icode VARCHAR(1), rescounter INT, resname VARCHAR(3), atomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, siteid INT, PRIMARY KEY(entry, chainid, rescounter, atomname))")
    cur.execute("CREATE TABLE pdb_file_error_table(entry VARCHAR(11), reason VARCHAR(40))")
    for root, dirs, filenames in os.walk(directory):
        for file in filenames:
            if file[-3:] == 'pdb' or file[-3:] == 'ent' or file[-3:] == 'brk':
                print "Reading file", file, "..."
                breakflag=0
                add_entry_to_pdb_table(directory, file)
                
def add_entry_to_pdb_table(directory, file):
    #   Add coordinates from a single PDB file
    cur.execute("TRUNCATE TABLE pdb_onefile")
    founddbref = 0
    foundatom = 0
    with open(directory + '/' + file,'r') as f:
        for x in f:
            if(len(x) >= 6):
                record = x[0:6]
                if record == 'DBREF ' or record == 'DBREF1':
                    founddbref=1
                if record == 'MODEL ':
                    cur.execute("INSERT INTO pdb_file_error_table (entry, reason) VALUES(%s, %s)", (file, 'Ensemble of models'))
                    return
                elif record == 'ATOM  ':
                    altloc = x[16]
                    if altloc == ' ':
                        #   Only include atoms with no alternate location, because atoms in PDB are often not grouped correctly by altloc
                        atomname = x[12:16]
                        chainid = x[21]
                        resseq = int(x[22:26])
                        icode = x[26]
                        if foundatom == 0:
                            if founddbref == 0:
                                cur.execute("INSERT INTO pdb_file_error_table (entry, reason) VALUES(%s, %s)", (file, 'No DBREF entry found'))
                                return
                            foundatom = 1
                            rescounter = 0
                        else:
                            if chainid==lastchainid:
                                if resseq < lastresseq:
                                    cur.execute("INSERT INTO pdb_file_error_table (entry, reason) VALUES(%s, %s)", (file, 'Residue indices out of order'))
                                    return
                                elif resseq > lastresseq:
                                    #   New residue on same chain
                                    rescounter += (resseq - lastresseq)
                                elif icode != lasticode:
                                        rescounter += 1
                            else:
                                #   New residue on new chain
                                rescounter += 1
                                #   Nterminusflag=1
                        resname = x[17:20]
                        cur.execute("SELECT %s IN (SELECT resname FROM amino_acid)", (resname))
                        testval=cur.fetchone()[0]
                        if testval!=1:
                            cur.execute("INSERT INTO pdb_file_error_table (entry, reason) VALUES(%s, %s)", (file, 'Residue not among the known amino acids'))
                            return
                        cur.execute("SELECT %s IN (SELECT atomname FROM amino_acid_atom WHERE resname = %s)", (atomname, resname))
                        testval=cur.fetchone()[0]
                        if testval!=1:
                            cur.execute("INSERT INTO pdb_file_error_table (entry, reason) VALUES(%s, %s)", (file, 'Unknown atom'))
                            return
                        rx = float(x[30:38])
                        ry = float(x[38:46])
                        rz = float(x[46:54])
                        cur.execute("INSERT INTO pdb_onefile (entry, chainid, resseq, icode, rescounter, resname, atomname, x, y, z) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", (file, chainid, resseq, icode, rescounter, resname, atomname, rx, ry, rz))
                        lastchainid = chainid
                        lastresseq = resseq
                        lasticode = icode
    cur.execute("INSERT INTO pdb SELECT * FROM pdb_onefile")

def map_to_cgmodel(cur):
    #   Map atomic coordinates to the coarse-grained representation
    print "Mapping to coarse-grained model..."
    cur.execute("CREATE TABLE cgatom1(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, resname VARCHAR(3), siteid INT, x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, siteid))")
    cur.execute("CREATE TABLE cgatom2(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, siteid INT, x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, siteid))")
    cur.execute("CREATE TABLE cgatom3(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, siteid INT, x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, siteid))")
    cur.execute("CREATE TABLE cgtemp(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, resname VARCHAR(3), siteid INT, x DOUBLE, y DOUBLE, z DOUBLE, dif1x DOUBLE, dif1y DOUBLE, dif1z DOUBLE, dif1norm DOUBLE, dif2x DOUBLE, dif2y DOUBLE, dif2z DOUBLE, projection DOUBLE, dif2norm DOUBLE, PRIMARY KEY(entry, chainid, rescounter, siteid))")
    cur.execute("CREATE TABLE cg(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, resname VARCHAR(3), siteid INT, x DOUBLE, y DOUBLE, z DOUBLE, exx DOUBLE, exy DOUBLE, exz DOUBLE, eyx DOUBLE, eyy DOUBLE, eyz DOUBLE, ezx DOUBLE, ezy DOUBLE, ezz DOUBLE, PRIMARY KEY(entry, chainid, rescounter, siteid))")

    #   Identify atoms defining the coarse-grained sites
    cur.execute("INSERT INTO cgatom1 SELECT pdb.entry, pdb.chainid, pdb.rescounter, pdb.resname, site_director_atoms.siteid, pdb.x, pdb.y, pdb.z FROM pdb JOIN site_director_atoms ON (pdb.resname = site_director_atoms.resname and pdb.atomname = site_director_atoms.atom1)")
    cur.execute("INSERT INTO cgatom2 SELECT pdb.entry, pdb.chainid, pdb.rescounter, site_director_atoms.siteid, pdb.x, pdb.y, pdb.z FROM pdb JOIN site_director_atoms ON (pdb.resname = site_director_atoms.resname and pdb.atomname = site_director_atoms.atom2)")
    cur.execute("INSERT INTO cgatom3 SELECT pdb.entry, pdb.chainid, pdb.rescounter, site_director_atoms.siteid, pdb.x, pdb.y, pdb.z FROM pdb JOIN site_director_atoms ON (pdb.resname = site_director_atoms.resname and pdb.atomname = site_director_atoms.atom3)")

    #   Add atoms for C-terminal sites, labeling them with siteid = -1
    cur.execute("CREATE TABLE Cterminal_residues (entry VARCHAR(11), chainid VARCHAR (1), rescounter INT)")
    cur.execute("INSERT INTO Cterminal_residues (entry, chainid, rescounter) SELECT entry, chainid, max(rescounter) FROM pdb GROUP BY entry, chainid;")
    cur.execute("INSERT INTO cgatom1 (entry, chainid, rescounter, resname, siteid, x, y, z) SELECT pdb.entry, pdb.chainid, pdb.rescounter, pdb.resname, %s, pdb.x, pdb.y, pdb.z FROM pdb JOIN Cterminal_residues JOIN Cterminus_atoms ON (pdb.entry = Cterminal_residues.entry and pdb.chainid = Cterminal_residues.chainid and pdb.rescounter = Cterminal_residues.rescounter and pdb.atomname = Cterminus_atoms.atom1)", (-1))
    cur.execute("INSERT INTO cgatom2 (entry, chainid, rescounter, siteid, x, y, z) SELECT pdb.entry, pdb.chainid, pdb.rescounter, %s, pdb.x, pdb.y, pdb.z FROM pdb JOIN Cterminal_residues JOIN Cterminus_atoms ON (pdb.entry = Cterminal_residues.entry and pdb.chainid = Cterminal_residues.chainid and pdb.rescounter = Cterminal_residues.rescounter and pdb.atomname = Cterminus_atoms.atom2)", (-1))
    cur.execute("INSERT INTO cgatom3 (entry, chainid, rescounter, siteid, x, y, z) SELECT pdb.entry, pdb.chainid, pdb.rescounter, %s, pdb.x, pdb.y, pdb.z FROM pdb JOIN Cterminal_residues JOIN Cterminus_atoms ON (pdb.entry = Cterminal_residues.entry and pdb.chainid = Cterminal_residues.chainid and pdb.rescounter = Cterminal_residues.rescounter and pdb.atomname = Cterminus_atoms.atom3)", (-1))

    #   Define orthonormal triad for each site
    cur.execute("INSERT INTO cgtemp (entry, chainid, rescounter, resname, siteid, x, y, z, dif1x, dif1y, dif1z) SELECT cgatom1.entry, cgatom1.chainid, cgatom1.rescounter, cgatom1.resname, cgatom1.siteid, cgatom1.x, cgatom1.y, cgatom1.z, cgatom2.x-cgatom1.x, cgatom2.y-cgatom1.y, cgatom2.z-cgatom1.z FROM cgatom1 JOIN cgatom2 ON (cgatom1.entry = cgatom2.entry and cgatom1.chainid = cgatom2.chainid and cgatom1.rescounter = cgatom2.rescounter and cgatom1.siteid = cgatom2.siteid)")
    cur.execute("UPDATE cgtemp SET dif1norm = SQRT(POW(dif1x, 2) + POW(dif1y, 2) + POW(dif1z, 2))")
    cur.execute("DELETE FROM cgtemp WHERE dif1norm = 0")
    cur.execute("UPDATE cgtemp JOIN cgatom3 ON (cgtemp.entry = cgatom3.entry and cgtemp.chainid = cgatom3.chainid and cgtemp.rescounter = cgatom3.rescounter and cgtemp.siteid = cgatom3.siteid) SET cgtemp.dif1x = cgtemp.dif1x / cgtemp.dif1norm, cgtemp.dif1y = cgtemp.dif1y / cgtemp.dif1norm, cgtemp.dif1z = cgtemp.dif1z / cgtemp.dif1norm, cgtemp.dif2x = cgatom3.x - cgtemp.x, cgtemp.dif2y = cgatom3.y - cgtemp.y, cgtemp.dif2z = cgatom3.z - cgtemp.z")
    cur.execute("UPDATE cgtemp SET projection = dif1x * dif2x + dif1y * dif2y + dif1z * dif2z")
    cur.execute("UPDATE cgtemp SET dif2x = dif2x - dif1x * projection, dif2y = dif2y - dif1y * projection, dif2z = dif2z - dif1z * projection")
    cur.execute("UPDATE cgtemp SET dif2norm = SQRT(POW(dif2x, 2) + POW(dif2y, 2) + POW(dif2z, 2))")
    cur.execute("DELETE FROM cgtemp WHERE dif2x IS NULL OR dif2norm = 0")
    cur.execute("UPDATE cgtemp SET dif2x = dif2x / dif2norm, dif2y = dif2y / dif2norm, dif2z = dif2z / dif2norm")
    cur.execute("INSERT INTO cg (entry, chainid, rescounter, resname, siteid, x, y, z, exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz) SELECT entry, chainid, rescounter, resname, siteid, x, y, z, dif1x, dif1y, dif1z, dif2x, dif2y, dif2z, dif1y * dif2z - dif1z * dif2y, dif1z * dif2x - dif1x * dif2z, dif1x * dif2y - dif1y * dif2x FROM cgtemp")

def record_moments_relative_to_oriented_cg_sites(cur):
    #   Add a siteid column to pdb to handle different sites for C-terminal and non-C-terminal residues
    cur.execute("UPDATE pdb JOIN site_atoms ON (pdb.resname = site_atoms.resname and pdb.atomname = site_atoms.atomname) SET pdb.siteid = site_atoms.siteid")
    cur.execute("UPDATE pdb JOIN Cterminus_atoms JOIN Cterminal_residues ON (pdb.entry = Cterminal_residues.entry and pdb.chainid = Cterminal_residues.chainid and pdb.rescounter = Cterminal_residues.rescounter and (pdb.atomname = Cterminus_atoms.atom1 or pdb.atomname = Cterminus_atoms.atom2 or pdb.atomname = Cterminus_atoms.atom3)) SET pdb.siteid = %s", (-1))

    #   Create table of relative positions
    print "Calculating positions relative to oriented sites..."
    cur.execute("CREATE TABLE relative(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, resname VARCHAR(3), siteid INT, atomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, orientedx DOUBLE, orientedy DOUBLE, orientedz DOUBLE, PRIMARY KEY(entry, chainid, rescounter, atomname))")
    cur.execute("INSERT INTO relative (entry, chainid, rescounter, resname, siteid, atomname, x, y, z) SELECT pdb.entry, pdb.chainid, pdb.rescounter, pdb.resname, cg.siteid, pdb.atomname, pdb.x - cg.x, pdb.y - cg.y, pdb.z - cg.z FROM pdb JOIN cg ON (pdb.entry = cg.entry and pdb.chainid = cg.chainid and pdb.rescounter = cg.rescounter and pdb.siteid = cg.siteid)")
    cur.execute("UPDATE relative JOIN cg ON (relative.entry = cg.entry and relative.chainid = cg.chainid and relative.rescounter = cg.rescounter and relative.siteid = cg.siteid) SET relative.orientedx = relative.x * cg.exx + relative.y * cg.exy + relative.z * cg.exz, relative.orientedy = relative.x * cg.eyx + relative.y * cg.eyy + relative.z * cg.eyz, relative.orientedz = relative.x * cg.ezx + relative.y * cg.ezy + relative.z * cg.ezz")

    #   Create table of average values and standard deviations
    print "Calculating moments of relative positions..."
    cur.execute("CREATE TABLE moments(resname VARCHAR(3), siteid INT, atomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, dx DOUBLE, dy DOUBLE, dz DOUBLE, rmsd DOUBLE, count INT, PRIMARY KEY(resname, siteid, atomname))")
    cur.execute("INSERT INTO moments (resname, siteid, atomname, x, y, z, dx, dy, dz, count) SELECT resname, siteid, atomname, avg(orientedx), avg(orientedy), avg(orientedz), stddev(orientedx), stddev(orientedy), stddev(orientedz), count(orientedx) FROM relative GROUP BY resname, siteid, atomname")
    cur.execute("UPDATE moments SET rmsd = SQRT(POW(dx, 2) + POW(dy, 2) + POW(dz, 2))")

def calculate_covariances(cur, filename):
    #   Read corrections file
    cur.execute("CREATE TABLE corrections(resname VARCHAR(3), atom1 VARCHAR(4), atom2 VARCHAR(4))")
    with open(filename) as f:
        for x in f:
            cur.execute("INSERT INTO corrections (resname, atom1, atom2) VALUES(%s, %s, %s)", (x[0:3], x[3:7], x[7:11]))

    #   Create lists of atoms according to whether they predict or are predicted by linear correction
    cur.execute("CREATE TABLE correcting_atom_list(resname VARCHAR(3), atomname VARCHAR(4))")
    cur.execute("CREATE TABLE corrected_atom_list(resname VARCHAR(3), atomname VARCHAR(4))")
    cur.execute("CREATE TABLE direct_list(resname VARCHAR(3), atomname VARCHAR(4))")
    cur.execute("INSERT INTO correcting_atom_list (resname, atomname) SELECT DISTINCT resname, atom2 FROM corrections")
    cur.execute("INSERT INTO corrected_atom_list (resname, atomname) SELECT DISTINCT resname, atom1 FROM corrections")
    cur.execute("INSERT INTO direct_list (resname, atomname) SELECT resname, atomname FROM amino_acid_atom WHERE (resname, atomname) NOT IN (SELECT resname, atomname FROM corrected_atom_list)")

    print "Calculating lab-frame predicted positions for atoms used in linear correction..."
    cur.execute("CREATE TABLE correcting_labframe(entry VARCHAR(11), chainid VARCHAR(1), rescounter INT, resname VARCHAR(3), correctingatomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, correctingatomname))")
    cur.execute("INSERT INTO correcting_labframe (entry, chainid, rescounter, resname, correctingatomname, x, y, z) SELECT cg.entry, cg.chainid, cg.rescounter, moments.resname, moments.atomname, cg.x + moments.x  * cg.exx + moments.y  * cg.eyx + moments.z  * cg.ezx, cg.y + moments.y  * cg.exy + moments.y  * cg.eyy + moments.z  * cg.ezy, cg.z + moments.x  * cg.exz + moments.y  * cg.eyz + moments.z  * cg.ezz FROM moments JOIN cg JOIN correcting_atom_list ON (cg.resname = moments.resname and cg.siteid = moments.siteid and cg.siteid >= 0 and cg.resname = correcting_atom_list.resname and moments.atomname = correcting_atom_list.atomname)")
                    
    print "Calculating predicted positions in frame of corrected atom..."
    cur.execute("CREATE TABLE correcting_correctedframe(entry VARCHAR(11), chainid VARCHAR(1), rescounter INT, resname VARCHAR(3), siteid INT, correctedatomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, correctedatomname))")
    cur.execute("INSERT INTO correcting_correctedframe(entry, chainid, rescounter, resname, siteid, correctedatomname, x, y, z) SELECT cg.entry, cg.chainid, cg.rescounter, cg.resname, cg.siteid, corrections.atom1, (correcting_labframe.x - cg.x) * cg.exx + (correcting_labframe.y - cg.y) * cg.exy + (correcting_labframe.z - cg.z) * cg.exz, (correcting_labframe.x - cg.x) * cg.eyx + (correcting_labframe.y - cg.y) * cg.eyy + (correcting_labframe.z - cg.z) * cg.eyz, (correcting_labframe.x - cg.x) * cg.ezx + (correcting_labframe.y - cg.y) * cg.ezy + (correcting_labframe.z - cg.z) * cg.ezz FROM correcting_labframe JOIN cg JOIN corrections JOIN site_atoms ON (corrections.resname = cg.resname and corrections.atom2 = correcting_labframe.correctingatomname and cg.entry = correcting_labframe.entry and cg.chainid = correcting_labframe.chainid and cg.rescounter = correcting_labframe.rescounter and cg.siteid = site_atoms.siteid and cg.resname = site_atoms.resname and corrections.atom1 = site_atoms.atomname)")

    #   Update relative table to be relative to average
    cur.execute("UPDATE relative JOIN moments ON (relative.resname = moments.resname AND relative.atomname = moments.atomname) SET relative.orientedx = relative.orientedx - moments.x, relative.orientedy = relative.orientedy - moments.y, relative.orientedz = relative.orientedz - moments.z")
                
    print "Calculating covariances..."
    cur.execute("CREATE TABLE cov(resname VARCHAR(3), atomname VARCHAR(4), pos_difference_av_x DOUBLE, pos_difference_av_y DOUBLE, pos_difference_av_z DOUBLE, predicted_pos_av_x DOUBLE, predicted_pos_av_y DOUBLE, predicted_pos_av_z DOUBLE, cov_cross_xx DOUBLE, cov_cross_xy DOUBLE, cov_cross_xz DOUBLE, cov_cross_yx DOUBLE, cov_cross_yy DOUBLE, cov_cross_yz DOUBLE, cov_cross_zx DOUBLE, cov_cross_zy DOUBLE, cov_cross_zz DOUBLE, cov_self_xx DOUBLE, cov_self_xy DOUBLE, cov_self_xz DOUBLE, cov_self_yx DOUBLE, cov_self_yy DOUBLE, cov_self_yz DOUBLE, cov_self_zx DOUBLE, cov_self_zy DOUBLE, cov_self_zz DOUBLE, interceptx DOUBLE, intercepty DOUBLE, interceptz DOUBLE, slopexx DOUBLE, slopexy DOUBLE, slopexz DOUBLE, slopeyx DOUBLE, slopeyy DOUBLE, slopeyz DOUBLE, slopezx DOUBLE, slopezy DOUBLE, slopezz DOUBLE, PRIMARY KEY(resname, atomname))")
    cur.execute("INSERT INTO cov (resname, atomname, pos_difference_av_x, pos_difference_av_y, pos_difference_av_z, predicted_pos_av_x, predicted_pos_av_y, predicted_pos_av_z, cov_cross_xx, cov_cross_xy, cov_cross_xz, cov_cross_yx, cov_cross_yy, cov_cross_yz, cov_cross_zx, cov_cross_zy, cov_cross_zz, cov_self_xx, cov_self_xy, cov_self_xz, cov_self_yx, cov_self_yy, cov_self_yz, cov_self_zx, cov_self_zy, cov_self_zz) SELECT relative.resname, relative.atomname, avg(relative.orientedx), avg(relative.orientedy), avg(relative.orientedz), avg(correcting_correctedframe.x), avg(correcting_correctedframe.y), avg(correcting_correctedframe.z), avg(relative.orientedx * correcting_correctedframe.x), avg(relative.orientedx * correcting_correctedframe.y), avg(relative.orientedx * correcting_correctedframe.z), avg(relative.orientedy * correcting_correctedframe.x), avg(relative.orientedy * correcting_correctedframe.y), avg(relative.orientedy * correcting_correctedframe.z), avg(relative.orientedz * correcting_correctedframe.x), avg(relative.orientedz * correcting_correctedframe.y), avg(relative.orientedz * correcting_correctedframe.z), avg(correcting_correctedframe.x * correcting_correctedframe.x), avg(correcting_correctedframe.x * correcting_correctedframe.y), avg(correcting_correctedframe.x * correcting_correctedframe.z), avg(correcting_correctedframe.y * correcting_correctedframe.x), avg(correcting_correctedframe.y * correcting_correctedframe.y), avg(correcting_correctedframe.y * correcting_correctedframe.z), avg(correcting_correctedframe.z * correcting_correctedframe.x), avg(correcting_correctedframe.z * correcting_correctedframe.y), avg(correcting_correctedframe.z * correcting_correctedframe.z) FROM relative JOIN correcting_correctedframe ON (relative.entry = correcting_correctedframe.entry AND relative.chainid = correcting_correctedframe.chainid AND relative.rescounter = correcting_correctedframe.rescounter AND relative.atomname = correcting_correctedframe.correctedatomname) GROUP BY relative.resname, relative.atomname")
    cur.execute("UPDATE cov SET cov_cross_xx = cov_cross_xx - pos_difference_av_x * predicted_pos_av_x, cov_cross_xy = cov_cross_xy - pos_difference_av_x * predicted_pos_av_y, cov_cross_xz = cov_cross_xz - pos_difference_av_x * predicted_pos_av_z, cov_cross_yx = cov_cross_yx - pos_difference_av_y * predicted_pos_av_x, cov_cross_yy = cov_cross_yy - pos_difference_av_y * predicted_pos_av_y, cov_cross_yz = cov_cross_yz - pos_difference_av_y * predicted_pos_av_z, cov_cross_zx = cov_cross_zx - pos_difference_av_z * predicted_pos_av_x, cov_cross_zy = cov_cross_zy - pos_difference_av_z * predicted_pos_av_y, cov_cross_zz = cov_cross_zz - pos_difference_av_z * predicted_pos_av_z, cov_self_xx = cov_self_xx - predicted_pos_av_x * predicted_pos_av_x, cov_self_xy = cov_self_xy - predicted_pos_av_x * predicted_pos_av_y, cov_self_xz = cov_self_xz - predicted_pos_av_x * predicted_pos_av_z, cov_self_yx = cov_self_yx - predicted_pos_av_y * predicted_pos_av_x, cov_self_yy = cov_self_yy - predicted_pos_av_y * predicted_pos_av_y, cov_self_yz = cov_self_yz - predicted_pos_av_y * predicted_pos_av_z, cov_self_zx = cov_self_zx - predicted_pos_av_z * predicted_pos_av_x, cov_self_zy = cov_self_zy - predicted_pos_av_z * predicted_pos_av_y, cov_self_zz = cov_self_zz - predicted_pos_av_z * predicted_pos_av_z, interceptx = pos_difference_av_x, intercepty = pos_difference_av_y, interceptz = pos_difference_av_z")

    #   Export matrices to python, calculate inverse matrix in python, then import back into MySQL
    cur.execute("SELECT resname, atomname, cov_self_xx, cov_self_xy, cov_self_xz, cov_self_yx, cov_self_yy, cov_self_yz, cov_self_zx, cov_self_zy, cov_self_zz, cov_cross_xx, cov_cross_xy, cov_cross_xz, cov_cross_yx, cov_cross_yy, cov_cross_yz, cov_cross_zx, cov_cross_zy, cov_cross_zz, predicted_pos_av_x, predicted_pos_av_y, predicted_pos_av_z FROM cov")
    rows = cur.fetchall()
    for row in rows:
        resname = row[0]
        atomname = row[1]
        self = matrix([[row[2],row[3],row[4]],[row[5],row[6],row[7]],[row[8],row[9],row[10]]])
        cross = matrix([[row[11],row[12],row[13]],[row[14],row[15],row[16]],[row[17],row[18],row[19]]])
        predicted_pos = matrix([[row[20]], [row[21]], [row[22]]])
        slope = cross * (self.I).T
        interceptsubtract = slope * predicted_pos
        cur.execute("UPDATE cov SET slopexx = %s, slopexy = %s, slopexz = %s, slopeyx = %s, slopeyy = %s, slopeyz = %s, slopezx = %s, slopezy = %s, slopezz = %s, interceptx = interceptx - %s, intercepty = intercepty - %s, interceptz = interceptz - %s WHERE resname = %s AND atomname = %s", (slope[0, 0], slope[0, 1], slope[0, 2], slope[1, 0], slope[1, 1], slope[1, 2], slope[2, 0], slope[2, 1], slope[2, 2], interceptsubtract[0, 0], interceptsubtract[1, 0], interceptsubtract[2, 0], resname, atomname))

    print "Predicting positions using corrected backmapping..."
    #   Join with pdb so I only include atoms that exist in original configuration
    cur.execute("CREATE TABLE relative_corrected(entry VARCHAR(11), chainid VARCHAR (1), rescounter INT, resname VARCHAR(3), siteid INT, atomname VARCHAR(4), x DOUBLE, y DOUBLE, z DOUBLE, PRIMARY KEY(entry, chainid, rescounter, atomname))")
    cur.execute("INSERT INTO relative_corrected (entry, chainid, rescounter, resname, siteid, atomname, x, y, z) SELECT correcting_correctedframe.entry, correcting_correctedframe.chainid, correcting_correctedframe.rescounter, correcting_correctedframe.resname, correcting_correctedframe.siteid, correcting_correctedframe.correctedatomname, relative.orientedx - (cov.interceptx + cov.slopexx * correcting_correctedframe.x + cov.slopexy * correcting_correctedframe.y + cov.slopexz * correcting_correctedframe.z), relative.orientedy - (cov.intercepty + cov.slopeyx * correcting_correctedframe.x + cov.slopeyy * correcting_correctedframe.y + cov.slopeyz * correcting_correctedframe.z), relative.orientedz - (cov.interceptz + cov.slopezx * correcting_correctedframe.x + cov.slopezy * correcting_correctedframe.y + cov.slopezz * correcting_correctedframe.z) FROM correcting_correctedframe JOIN relative JOIN cov ON (relative.entry = correcting_correctedframe.entry and relative.chainid = correcting_correctedframe.chainid and relative.rescounter = correcting_correctedframe.rescounter and relative.siteid >= 0 and relative.atomname = correcting_correctedframe.correctedatomname and cov.resname = relative.resname and cov.atomname = relative.atomname)")

def calculate_rmsd(cur):
    #   Calculate root-mean-square displacement by residue and atom between original and backmapped atomic configurations
    print "Calculating rmsd..."
    cur.execute("CREATE TABLE rmsd(resname VARCHAR(3), siteid INT, atomname VARCHAR(4), dx DOUBLE, dy DOUBLE, dz DOUBLE, rmsd DOUBLE, count INT, PRIMARY KEY(resname, siteid, atomname))")
    cur.execute("INSERT INTO rmsd (resname, siteid, atomname, dx, dy, dz, rmsd, count) SELECT moments.resname, moments.siteid, moments.atomname, moments.dx, moments.dy, moments.dz, moments.rmsd, moments.count FROM moments JOIN direct_list ON (moments.resname = direct_list.resname and moments.atomname = direct_list.atomname)")
    
    cur.execute("INSERT INTO rmsd (resname, siteid, atomname, dx, dy, dz, rmsd, count) SELECT resname, siteid, atomname, sqrt(avg(x*x)), sqrt(avg(y*y)), sqrt(avg(z*z)), sqrt(avg(x*x) + avg(y*y) + avg(z*z)), count(x) FROM relative_corrected GROUP BY resname, siteid, atomname")

con = mdb.connect('localhost', 'root', )
with con:
    cur = con.cursor()
    try:
        cur.execute("DROP DATABASE proteindb")
    except:
        pass
    cur.execute("CREATE DATABASE proteindb")
    cur.execute("USE proteindb")
    create_amino_acid_atom_table(cur, amino_acid_directory)
    create_site_director_table(cur, site_director_atoms_file)
    create_site_table(cur, site_atoms_file)
    create_cterminus_atoms_table(cur, cterminus_atoms_file)
    create_pdb_table(cur, pdb_directory)
    map_to_cgmodel(cur)
    record_moments_relative_to_oriented_cg_sites(cur)
    calculate_covariances(cur, corrections_atoms_file)
    calculate_rmsd(cur)

