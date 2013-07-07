"""
This program is designed to parse SNP files from OpenSNP.  

Parameters are stored in an text file whose location is passed in as a parameter to the 
program.  The file is tab-separated and contains identifier/value pairs.  See parsefiles.txt 
for an example.
"""

import sys
import fnmatch
import os
import parser
from datetime import datetime

# import re

true_values = ["TRUE", "T", "1", "YES", "Y"]
start_time = datetime.now()
DIR = "DIR"
FILES = "FILES"
RSID  = "RSID"
CHROMSTART = "CHROMSTART"
CHROMEND = "CHROMEND"
POSSTART = "POSSTART"
POSEND = "POSEND"
SHOWFILEPROGRESS = "SHOWFILEPROGRESS"
SHOWPROGRESSLINES = "SHOWPROGRESS#LINES"
SHOWSELECTEDFILES = "SHOWSELECTEDFILES"

def bool_param(param):
    if (param.upper() in true_values):
        value = True
    else:
        value = False
    return value
    
    
# Process the SNPs according to the params passed in.
def cleanup(params):
    if (DIR not in params):
        params[DIR] = "."
    if (FILES not in params):
        params[FILES] = {{0,"Default"}, ["*"]}
    if (RSID in params):
        params[RSID] = params[RSID].upper()
    else:
        params[RSID] = "*"        
    if (CHROMSTART in params):
        params[CHROMSTART] = int(params[CHROMSTART])
    else:
        params[CHROMSTART] = 0
    if (CHROMEND in params):
        params[CHROMEND] = int(params[CHROMEND])
    else:
        params[CHROMEND] = 23
    if (POSSTART in params):
        params[POSSTART] = int(params[POSSTART])
    else:
        params[POSSTART] = 0
    if (POSEND in params):
        params[POSEND] = int(params[POSEND])
    else:
        params[POSEND] = sys.maxint
    if (SHOWFILEPROGRESS in params):
        params[SHOWFILEPROGRESS] = bool_param(params[SHOWFILEPROGRESS])
    else:
        params[SHOWFILEPROGRESS] = False
    if (SHOWPROGRESSLINES in params):
        params[SHOWPROGRESSLINES] = int(params[SHOWPROGRESSLINES])
    else:
        params[SHOWPROGRESSLINES] = 0
    if (SHOWSELECTEDFILES in params):
        params[SHOWSELECTEDFILES] = bool_param(params[SHOWSELECTEDFILES])
    else:
        params[SHOWSELECTEDFILES] = false
    return

# Remove leading and trailing double quotes
def strip_quotes(val):
    if val[0:1] == "\"":
        val = val[1:]
    if val[len(val)-1:] == "\"":
        val = val[:len(val)-1]
    return val
# Increment a counter in a dictionary
def increment(counter, value):
    if (value in counter):
        counter[value] += 1
    else:
        counter[value] = 1
    return

# If the file matches one of the file patterns in the parm, return the label to use
def file_name_label(filename, parm):
    label = ""
    for entry in sorted(parm.items()):
        for file_match in entry[1]:
            if fnmatch.fnmatch(filename, file_match):
                label = entry[0][1]
                break
        if label != "":
            break
    return label

def bypass(filename):
    return "-exome-" in filename

# Process the SNPs according to the params passed in.
def process(params):
    cleanup(params)
    rsids = {}
    files = 0
    print "Processing files"
#    rsid_re = re.compile(params["RSID"])
    chrom_start = params[CHROMSTART]
    chrom_end = params[CHROMEND]
    pos_start = params[POSSTART]
    pos_end = params[POSEND]
    d = params[DIR]
    directory = os.listdir(d)
    dir_count = 0
    fileTypeCounts={}
    skipped_files = []
    
    # Count number of files to process for progress reporting
    for filename in directory:
        if (bypass(filename)):
            continue
        if (file_name_label(filename, params[FILES]) != ""):
            if params[SHOWSELECTEDFILES]: 
                print "Selected:", filename
            dir_count += 1
            if (fnmatch.fnmatch(filename, "*23andme.txt")):
                increment(fileTypeCounts, "23andme")
            elif (fnmatch.fnmatch(filename, "*illumina.txt")):
                    increment(fileTypeCounts, "illumina")
            elif (fnmatch.fnmatch(filename, "*IYG.txt")):
                    increment(fileTypeCounts, "iyg")
            elif (fnmatch.fnmatch(filename, "*decodeme.txt")):
                    increment(fileTypeCounts, "decodeme")
                    
    # Summarize counts by file type
    print
    for count in fileTypeCounts.iteritems():
        print count[0], ":", count[1]
    print
            
    # Process each selected file.  Returns a dictionary structured as {rsid:{(chromosome, position):{genotype:count_of_occurences}}}
    for filename in directory:
        label = file_name_label(filename, params[FILES])
        if (label != ""):
            if (bypass(filename)):
                skipped_files.append(filename)
                continue
            with open(os.path.join(d, filename)) as f:
                files += 1
                if params[SHOWFILEPROGRESS] and params[SHOWPROGRESSLINES] <= 0:
                    elapsed = divmod((datetime.now() - start_time).total_seconds(), 60)
                    print "files: %d/%d  elapsed: %d minutes, %d seconds\r" % (files, dir_count, elapsed[0], elapsed[1])
                sys.stdout.flush()
                lines_processed = 0
                lines_read = 0
                line = f.readline()
                exit_file = False
                while (line != "" and not exit_file):
                    lines_read += 1
                    skip = False;
                    if (line[0:1] == "#"):
                        skip = True
                    else:
                        # 23andme structure (tab separated):
                        # rsid      chromosome    pos    genotype
                        # rs4477212    1         72017    AA
                        if (fnmatch.fnmatch(filename, "*23andme.txt")):
                            data = line.strip().split("\t")
                            if (len(data) != 4):
                                skip = True
                                
                        # illumina can have two structures:
                        elif (fnmatch.fnmatch(filename, "*illumina.txt")):
                            if (line[0:1] == "\""):
                                # RSID,CHROMOSOME,POSITION,RESULT
                                # "rs3094315","1","742429","AA"
                                data = map(strip_quotes, line.strip().split(",")) 
                                if (len(data) != 4):
                                    skip = True
                            else:
                                #rsid        chromosome   position    allele1    allele2
                                #rs4477212        1        82154    T    T
                                data = line.strip().split("\t")
                                if (len(data) == 5):
                                    data = data[0], data[1], data[2], data[3] + data[4]
                                else:
                                    skip = True
                                    
                        # IYG structure:
                        # RSID,RESULT
                        # rs2131925    TT
                        elif (fnmatch.fnmatch(filename, "*IYG.txt")):
                            data = strip_quotes(line.strip()).split("\t")
                            if (len(data) == 2):
                                data = (data[0], "0", "0", data[1])
                            else:
                                skip = True
                                
                        # decodeme structure:
                        # Name,Variation,Chromosome,Position,Strand,YourCode
                        # rs4345758,C/T,1,28663,+,TT
                        elif (fnmatch.fnmatch(filename, "*decodeme.txt")):
                            data = line.strip().split(",")
                            if (len(data) == 6):
                                data = (data[0], data[2], data[3], data[5])
                            else:
                                skip = True
                        else: 
                            skip = True
                            skipped_files.append(filename)
                            exit_file = True
                            
                        if (not skip):
                            try:
                                chrom = int(data[1])
                            except ValueError:
                                skip = True
                            try:
                                pos = int(data[2])
                            except ValueError:
                                skip = True
                        if len(data[0]) > 20:
                            skipped_files.append(filename)
                            exit_file = True
                    if not exit_file and not skip:
                        rsid = data[0].upper()
                        genotype = data[3]
                        
                        # Optimization.  Relies on files being in chromosome # sequence
                        if chrom_end > 0 and chrom > chrom_end:
                            exit_file = True
                            continue
                        if (fnmatch.fnmatch(rsid, params[RSID]) and chrom >= chrom_start 
                            and pos >= pos_start and pos <= pos_end):
                            lines_processed += 1
                            chrom_position = (chrom,pos)
                            if (rsid in rsids):
                                labels = rsids[rsid]
                                if (label in labels):
                                    chrom_positions = labels[label]
                                    if chrom_position in chrom_positions:
                                        genotypes = chrom_positions[chrom_position]
                                        if (genotype in genotypes):
                                            genotypes[genotype] = genotypes[genotype] + 1
                                        else:
                                            genotypes[genotype] = 1
                                        chrom_positions[chrom_position] = genotypes
                                    else:
                                        chrom_positions[chrom_position] = {genotype:1}
                                    labels[label] = chrom_positions
                                else:
                                    labels[label] = {(chrom_position):{genotype:1}}
                                rsids[rsid] = labels
                            else:
                                rsids[rsid] = {label:{(chrom_position):{genotype:1}}}
                    if params[SHOWPROGRESSLINES] > 0 and lines_read % params[SHOWPROGRESSLINES] == 0:
                        elapsed = divmod((datetime.now() - start_time).total_seconds(), 60)
                        print "file: {:,d}/{:,d}  lines read: {:,d}  processed {:,d}  elapsed: {:,d} minutes, {:,d} seconds\r".format(files, dir_count, lines_read, lines_processed, int(elapsed[0]), int(elapsed[1])) 
                        sys.stdout.flush()
                    line = f.readline()
                if (params[SHOWSELECTEDFILES] and lines_processed == 0):
                    print "    No lines to process in ", filename
    if(len(skipped_files) > 0):
        print 
        print "Skipped Files"
        for file in skipped_files:
            print file
    return rsids

"""
# This doesn't work but I'm keeping it for now as an example of how to manipulate the keys in sorting a dictionary for output
def rsid_key_seq( item ):
    # RS numbers can begin with "RS", "VG" or "I".  We want to sequence all the "I"s first, then the "R"s, then the "V"s.
    if item[0:1] == "R":
        return int(item[2:]) + 1000000000
    elif item[0:1] == "V":
        return int(item[2:]) + 5000000000
    else:
        return int(item[1:])
"""        
if __name__=="__main__":
    filename = sys.argv[1]
    params = {}
    f=file(filename,"r")
    for line in f:
            if line[0:1] != "#":
                values = line.strip().split("\t")
                if (len(values) > 1 and values[0].strip()):
                    name = values[0].strip()
                    val = values[1].strip()
                    if (name[0:len(FILES)] == FILES):
                        # If this is a files param, assemble all the file values as a 
                        # dictionary in the form {(priority,label): [fileName, fileName, fileName]}
                        # so if multiple lines have the same (priority,label) combination, they
                        # will be consolidated
                        priority_label = (sys.maxint,"No Label") # default
                        priority_label_vals = name.split(":")
                        priority_label_len = len(priority_label_vals)
                        if (priority_label_len > 2):
                            priority_label = (int(priority_label_vals[2]), priority_label_vals[1])
                        priority_labels = {} # default
                        if FILES in params:
                            priority_labels = params[FILES]  # Get previously parsed files
                        files_list = [] # default
                        if priority_label in priority_labels:
                            files_list = priority_labels[priority_label]                            
                        files_list.extend([x.strip() for x in val.split(",")])
                        priority_labels[priority_label] = files_list
                        name = FILES
                        val = priority_labels
                    params[name.upper()]=val
    f.close()
    results = process(params)

print "\n"

if (len(results) > 0):
#    for entry in sorted(results.items(), key=lambda t: rsid_key_seq(t[0])):
    for entry in sorted(results.items()):
        print entry[0], ":", entry[1], "           " # Add whitespace to ensure there are no leftover chars when we overprint the prior line.
else:
    print "Nothing matched selections"

elapsed = divmod((datetime.now() - start_time).total_seconds(), 60)
print
print "Elapsed: %d minutes, %d seconds" % (elapsed[0], elapsed[1])


            