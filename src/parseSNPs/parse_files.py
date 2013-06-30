"""
This program is designed to parse SNP files from OpenSNP.  

Parameters are stored in an text file whose location is passed in as a parameter to the 
program.  The file is tab-separated and contains identifier/value pairs.  Allowed 
identifiers include:
DIR - The directory the SNP files are found in.  If "." or unspecified, the current directory
    is used.
FILE - The name of the SNP file to process.  If "*" or unspecified, all files are processed.
RSID - The RSID to search for int The file.  If unspecified or "*", all RSIDs are processed.
CHROMSTART - The starting chromosome to process.  If unspecified, all up to the CHROMEND# 
    are processed
CHROMEND - The last chromosome to process.  If unspecified, all starting with CHROMSTART# 
    are processed
POSSTART - The starting position on the chromosome to process.  If unspecified, all positions 
    up to the END# are processed
POSEND - The last position on the chromosome to process.  If unspecified, all positions starting 
    with START# are processed
SHOWFILEPROGRESS - If true, show the progress of files as we process each
SHOWLINEPROGRESS - If true, show the progress of files and the lines in them as we process each
"""

import sys
import fnmatch
import os
from datetime import datetime

# import re

true_values = ["TRUE", "T", "1", "YES", "Y"]
start_time = datetime.now()

def bool_param(param):
    if (param.upper() in true_values):
        value = True
    else:
        value = False
    return value
    
    
# Process the SNPs according to the params passed in.
def cleanup(params):
    if ("DIR" not in params):
        params["DIR"] = "."
    if ("FILE" not in params):
        params["FILE"] = "*"
    if ("RSID" in params):
        params["RSID"] = params["RSID"].upper()
    else:
        params["RSID"] = "*"        
    if ("CHROMSTART" in params):
        params["CHROMSTART"] = int(params["CHROMSTART"])
    else:
        params["CHROMSTART"] = 0
    if ("CHROMEND" in params):
        params["CHROMEND"] = int(params["CHROMEND"])
    else:
        params["CHROMEND"] = 23
    if ("POSSTART" in params):
        params["POSSTART"] = int(params["POSSTART"])
    else:
        params["POSSTART"] = 0
    if ("POSEND" in params):
        params["POSEND"] = int(params["POSEND"])
    else:
        params["POSEND"] = sys.maxint
    if ("SHOWFILEPROGRESS" in params):
        params["SHOWFILEPROGRESS"] = bool_param(params["SHOWFILEPROGRESS"])
    else:
        params["SHOWFILEPROGRESS"] = False
    if ("SHOWLINEPROGRESS" in params):
        params["SHOWLINEPROGRESS"] = bool_param(params["SHOWLINEPROGRESS"])
    else:
        params["SHOWLINEPROGRESS"] = False
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

# Process the SNPs according to the params passed in.
def process(params):
    cleanup(params)
    rsids = {}
    files = 0
    print "Processing files"
#    rsid_re = re.compile(params["RSID"])
    chrom_start = params["CHROMSTART"]
    chrom_end = params["CHROMEND"]
    pos_start = params["POSSTART"]
    pos_end = params["POSEND"]
    d = params["DIR"]
    directory = os.listdir(d)
    dir_count = 0
    fileTypeCounts={}
    skipped_files = []
    
    # Count number of files to process for progress reporting
    for filename in directory:
        if (fnmatch.fnmatch(filename, params["FILE"])):
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
        if (fnmatch.fnmatch(filename, params["FILE"])):
            with open(os.path.join(d, filename)) as f:
                files += 1
                if params["SHOWFILEPROGRESS"] and not params["SHOWLINEPROGRESS"]:
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
                    elif not skip:
                        rsid = data[0].upper()
                        genotype = data[3]
                        
                        # Optimization.  Relies on files being in chromosome # sequence
                        if (rsid == "RS10403190"):
                            i = 1
                        if chrom_end > 0 and chrom > chrom_end:
                            exit_file = True
                            continue
                        if (fnmatch.fnmatch(rsid, params["RSID"]) and chrom >= chrom_start 
                            and pos >= pos_start and pos <= pos_end):
                            lines_processed += 1
                            chrom_position = (chrom,pos)
                            if (rsid in rsids):
                                chrom_positions = rsids[rsid]
                                if chrom_position in chrom_positions:
                                    genotypes = chrom_positions[chrom_position]
                                    if (genotype in genotypes):
                                        genotypes[genotype] = genotypes[genotype] + 1
                                    else:
                                        genotypes[genotype] = 1
                                    chrom_positions[chrom_position] = genotypes
                                else:
                                    chrom_positions[chrom_position] = {genotype:1}
                                rsids[rsid] = chrom_positions
                            else:
                                rsids[rsid] = {(chrom_position):{genotype:1}}
                    if params["SHOWLINEPROGRESS"]:
                        elapsed = divmod((datetime.now() - start_time).total_seconds(), 60)
                        print "file: %d/%d  lines read: %d, processed %d  elapsed: %d minutes, %d seconds\r" % (files, dir_count, lines_read, lines_processed, elapsed[0], elapsed[1]) 
                        sys.stdout.flush()
                    line = f.readline();
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
    with open(filename) as f:
        line = f.readline()
        while (line != ""):
            if line[0:1] != "#":
                values = line.strip().split("\t")
                if (len(values) > 1 and values[0].strip()):
                    params[values[0].strip().upper()]=values[1].strip()
            line = f.readline()

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


            