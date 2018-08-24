from __future__ import division, print_function
import sys
inputpdb = sys.argv[1]
# Usage:
# python Center.py 1xyz.pdb

#===============================================================================================================================================================
# From mail:
# 1. determine the "center of gravity", that is, the average x,y,z tuple of all coordinates -> (Center of geometry)
# 2. determine the average distance of all coordinates from this center.

# Lets read the pdb and hold X,Y,Z coordinates in holder (dicitonary)
# HOLDER (Key -> name of atom (for future reference), values -> xyz
#===============================================================================================================================================================
HOLDER = {}
with open(inputpdb) as pdb:
    for lines in pdb:
        if lines.startswith("ATOM") : # or lines.startswith("HETATM") ## in case you want to include the non-protein/modified residues also.
            
            # line in a protein:
            # ATOM      1  N   ASN A 516      41.511  25.152  36.876  1.00 22.29      A    N  
            # 0         1  2   3   4  5       6       7       8       9    10         11   12 # index
            lines_split = lines.split() #splitted by space
            
            # Holder_Key = Atom type (name of atom, index->1,2,3,4,5 or [1:6]) values-> X,Y,Z (index->6,7,8 or [6:9])
            Name_of_atom = "{}_{}_{}_{}_{}".format(*lines_split[1:6]) # fancy way to write "1_N_ASN_A_516"  
            HOLDER[Name_of_atom] = map(float, lines_split[6:9]) #Change all xyz from string to float 

#===============================================================================================================================================================
# Calculate Center of Geometry of protein:

# 1. Get all xyz values from HOLDER(atom name is irrelevant)
ALL_X, ALL_Y, ALL_Z = zip(*HOLDER.values())

# 2. Total number of atoms :)
SIZE_X, SIZE_Y, SIZE_Z = len(ALL_X), len(ALL_Y), len(ALL_Z)

# 3. Coordinates of center of geometry = (Mean_X, Mean_Y, Mean_Z)
Center_of_geometry = [sum(ALL_X) / SIZE_X, 
                      sum(ALL_Y) / SIZE_Y, 
                      sum(ALL_Z) / SIZE_Z ]

#Output (for reference):
print("Center_of_geometry for {}: {}".format(inputpdb, Center_of_geometry))
#===============================================================================================================================================================

#Code for measuring Euclidean distance of two 3D coordinates [ https://en.wikipedia.org/wiki/Euclidean_distance#Three_dimensions ]
Distance = lambda ATOM_A, ATOM_B : sum([(ATOM_A[i] - ATOM_B[i])**2 for i in range(3)]) ** 0.5

#Start new HOLDER for distance
Distance_HOLDER = {}

#To get maximum distance or diameter 
MAXIMUM_DISTANCE = 0.0
MAXIMUM_ATOMSET  = ""

#Print distances in a log file, for reference
with open("Distance_from_COG.txt", "w") as OUTFILE:
    OUTFILE.write("ATOM_NAME\tDistance(in_Angstrom)\n")
    for each_atom in HOLDER:
        Distance_from_COG = Distance(Center_of_geometry, HOLDER[each_atom])
        
        Distance_HOLDER[each_atom] = Distance_from_COG
        
        OUTFILE.write("{}\t{}\n".format(each_atom, Distance_from_COG))
        
        if Distance_from_COG > MAXIMUM_DISTANCE:
            MAXIMUM_DISTANCE = Distance_from_COG
            MAXIMUM_ATOMSET  = each_atom

print("MAXIMUM_DISTANCE: {} for ATOM {} from Center_of_geometry".format(MAXIMUM_DISTANCE, MAXIMUM_ATOMSET))




