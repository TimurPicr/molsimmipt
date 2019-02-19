# create 2d lattice, low density so diffuses

dimension	2

lattice		sq2 0.5 origin 0.25 0.25 0.0
region		box block 0 31 0 7 -0.5 0.5
create_box	2 box
create_atoms	1 box

# LAMMPS letters via regions, convert to type 2 atoms

region	        L1 block 2 3 1 6 -0.5 0.5
region	        L2 block 2 5 1 2 -0.5 0.5
region		L union 2 L1 L2

region	        A1 block 6 7 1 6 -0.5 0.5
region	        A2 block 8 9 1 6 -0.5 0.5
region	        A3 block 6 9 3 4 -0.5 0.5
region	        A4 block 6 9 5 6 -0.5 0.5
region		A union 4 A1 A2 A3 A4

region	        1M1 block 10 11 1 6 -0.5 0.5
region	        1M2 block 14 15 1 6 -0.5 0.5
region	        1M3 block 10 15 5 6 -0.5 0.5
region	        1M4 block 12 13 3 6 -0.5 0.5
region		1M union 4 1M1 1M2 1M3 1M4

region	        2M1 block 16 17 1 6 -0.5 0.5
region	        2M2 block 20 21 1 6 -0.5 0.5
region	        2M3 block 16 21 5 6 -0.5 0.5
region	        2M4 block 18 19 3 6 -0.5 0.5
region		2M union 4 2M1 2M2 2M3 2M4

region	        P1 block 22 23 1 6 -0.5 0.5
region	        P2 block 24 25 3 6 -0.5 0.5
region	        P3 block 22 25 3 4 -0.5 0.5
region	        P4 block 22 25 5 6 -0.5 0.5
region		P union 4 P1 P2 P3 P4

region	        S1 block 26 29 5 6 -0.5 0.5
region	        S2 block 26 27 3 6 -0.5 0.5
region	        S3 block 26 29 3 4 -0.5 0.5
region	        S4 block 28 29 1 4 -0.5 0.5
region	        S5 block 26 29 1 2 -0.5 0.5
region		S union 5 S1 S2 S3 S4 S5

region          LAMMPS union 6 L A 1M 2M P S
set		region LAMMPS type 2
