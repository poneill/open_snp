"""
This program is designed to parse SNP files from OpenSNP.  The file contents are tab-separated and look 
generally like this (minus the column headings):

   (1)       (2)     (3)    (4)
rs3131972     1    742584    AG
rs12184325    1    743968    CC
rs3131969     1    744045    AG
rs12562034    1    758311    GG
rs12124819    1    766409    AA
rs2518996     1    782397    GG
rs12132517    1    788664    GG
rs11240777    1    788822    GG

The columns are:
(1) RSID - The SNP identifier (an rsid or an internal id)
(2) Chromosome #
(3) Position on the chromosome on the reference human genome
(3) The genotype call oriented with respect to the plus strand on the human reference 
sequence using reference human assembly build 36.

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
import operator

# import re

true_values = ["TRUE", "T", "1", "YES", "Y"]
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
    # Count number of files to process for progress reporting
    for filename in directory:
        if (fnmatch.fnmatch(filename, params["FILE"])):
            dir_count += 1
            
    # Process each selected file.  Returns a dictionary structured as {rsid:{(chromosome, position):{genotype:count_of_occurences}}}
    for filename in directory:
        if (fnmatch.fnmatch(filename, params["FILE"])):
            with open(os.path.join(d, filename)) as f:
                files += 1
                if params["SHOWFILEPROGRESS"] and not params["SHOWLINEPROGRESS"]:
                    print "file: " + str(files) + "/" + str(dir_count) + "\r",
                sys.stdout.flush()
                lines = 0
                line = f.readline()
                while (line != ""):
                    skip = False;
                    if (line[0:1] == "#"):
                        skip = True
                    else:
                        data = line.strip().split("\t")
                        if (len(data) != 4):
                            skip = True
                        else:
                            try:
                                chrom = int(data[1])
                            except ValueError:
                                skip = True
                            try:
                                pos = int(data[2])
                            except ValueError:
                                skip = True
                    if (not skip):
                        rsid = data[0].upper()
                        genotype = data[3]
    
#                        if (rsid_re.match(rsid)                          
                        if (fnmatch.fnmatch(rsid, params["RSID"]) 
                            and chrom >= chrom_start and chrom <= chrom_end 
                            and pos >= pos_start and pos <= pos_end):
                            lines += 1
                            chrom_position = (chrom,pos)
                            if (rsid in rsids):
                                chrom_positions = rsids[rsid]
                                if [chrom_position in chrom_positions]:
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
                                print "file: " + str(files) + "/" + str(dir_count) \
                                    + "  line: " + str(lines) + "\r", 
                                sys.stdout.flush()
                    line = f.readline();
    return rsids

def rsid_key_seq( item ):
    f=operator.itemgetter(0)
    return int(f(item)[2:])
        
if __name__=="__main__":
    filename = sys.argv[1]
    params = {}
    with open(filename) as f:
        line = f.readline()
        while (line != ""):
            values = line.strip().split("\t")
            if (len(values) > 1 and values[0].strip()):
                params[values[0].strip().upper()]=values[1].strip()
            line = f.readline()

    results = process(params)

print "\n"

for entry in sorted(results.items(), key=lambda t: int(t[0][2:])):
    print entry
            