#!/usr/bin/python

import sys

cli_args = sys.argv

if len(cli_args) != 2:
    print("Usage: ", cli_args[0], " <log to process>")
    sys.exit()

with open(cli_args[1]) as inp_file:
    nextline = ""
    while "MaxReplicaForce" not in nextline:
        nextline = inp_file.readline()

    while True:
        nextline = inp_file.readline()
        if "Climbing" in nextline:
            break
        log_data = nextline.split()
        tstep = log_data[0]
        with open("neb_path."+tstep, "w") as neb_file:
            for i in range(9, len(log_data), 2):
                neb_file.write(log_data[i] + " " + log_data[i+1] + "\n")

    while "MaxReplicaForce" not in nextline:
        nextline = inp_file.readline()

    while True:
        nextline = inp_file.readline()
        if not nextline:
            break
        log_data = nextline.split()
        tstep = log_data[0]
        with open("neb_path."+tstep, "w") as neb_file:
            for i in range(9, len(log_data), 2):
                neb_file.write(log_data[i] + " " + log_data[i+1] + "\n")
    

sys.exit()
