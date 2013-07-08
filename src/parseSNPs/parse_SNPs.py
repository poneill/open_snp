"""
This program is designed to parse SNP files from OpenSNP.  

Parameters are stored in an text file whose location is passed in as a parameter to the 
program.  The file is tab-separated and contains identifier/value pairs.  See parsefiles.txt 
for an example.
"""

import sys
import fnmatch
import os
from snp_utils import get_elapsed, string_to_bool, inc_dict_counter
from datetime import datetime

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
        params[SHOWFILEPROGRESS] = string_to_bool(params[SHOWFILEPROGRESS])
    else:
        params[SHOWFILEPROGRESS] = False
    if (SHOWPROGRESSLINES in params):
        params[SHOWPROGRESSLINES] = int(params[SHOWPROGRESSLINES])
    else:
        params[SHOWPROGRESSLINES] = 0
    if (SHOWSELECTEDFILES in params):
        params[SHOWSELECTEDFILES] = string_to_bool(params[SHOWSELECTEDFILES])
    else:
        params[SHOWSELECTEDFILES] = false
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

# Main processing method.  The one parameter, "parms", is a dictionary.  All selections are optional.  
# 
# Valid keys are as follows:
# DIR - The directory containing the SNP files.  E.g. 'C:\\OpenSNP'
# FILES - Used to select which files should be processed and allow these to be grouped under
#        A common label.  it will likely be common to create multiple groups that we want to 
#        compare with each group matching a phenotype.  The value for this parameter is a 
#        Dictionary.  The key to the dictionary is a tuple containing the group priority and
#        group name.  The value is a list of expressions that match the files to be included
#        in the group, with wildcards allowed.  So the parameter might look like this:
#          {(1, 'Tongue rollers'): ['user10_*.txt', 'user11_*.txt', 'user13_*.txt', 'user14_*.txt'],
#           (2, 'Non-tongue rollers'): ['user20_*.txt', 'user22_*.txt', 'user26_*.txt']} 
#        This describes two groups.  The first, with four file selections identifies users with
#        the "tongue roller" phenotype.  The second identifies users with the "non-tongue-roller"
#        phenoptype.  The reason each group has an associated priority is that this allows for a 
#        us to control which matches get checked first.  In the above case, we could have a third
#        group that matched "*.txt" files for all users who did not specify whether or not they
#        are tongue rollers (i.e. they did not match the first or second group).
# RSID - Allows the selections to be limited to one or more specific or all RSIDs.  
#        The value can be specified with wild cards. E.g. RSID10403190 or RSID104*
# CHROMSTART - Ignore all chromosomes numbers below this value
# CHROMEND - Ignore all chromosome numbers greater than this value
# POSSTART - Ignore all positions in chromosomes less than this value
# POSEND - Ignore all positions in chromosomes greater than this value
# SHOWFILEPROGRESS - If true and SHOWPROGRESSLINES is false, show progress information
#        as each new file is processed*
# SHOWPROGRESSLINES - An integer (e.g. 100000).  If > 0, write out a progress line when every 
#        nth record is processed, where "n" matches this integer.
# SHOWSELECTEDFILES - If true, show the names of files that are processed and also any processed
#        files that had no lines matching the limits above.*
# * True/False values try to match any reasonable value.  True can be represented by "TRUE", "T", 
#        "1", "YES" or "Y" in any case.
#
# The returned value is a structure that summarizes what was read. Its layout is:
#    {rsids:{(chromosome number, position):{label:{genotype:count}}}}
# That is for each RSID, it lists
#    ... for each chromosome number and position combination,
#        ... for each file group label
#            ... for each genotype
#                ...it lists the count
def parse_snps(params):
    cleanup(params)
    rsids = {}
    files = 0
    print "Processing files"
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
                inc_dict_counter(fileTypeCounts, "23andme")
            elif (fnmatch.fnmatch(filename, "*illumina.txt")):
                    inc_dict_counter(fileTypeCounts, "illumina")
            elif (fnmatch.fnmatch(filename, "*IYG.txt")):
                    inc_dict_counter(fileTypeCounts, "iyg")
            elif (fnmatch.fnmatch(filename, "*decodeme.txt")):
                    inc_dict_counter(fileTypeCounts, "decodeme")
                    
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
                    elapsed = snp_utils.get_elapsed();
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
                                chrom_positions = rsids[rsid]
                                if chrom_position in chrom_positions:
                                    labels = chrom_positions[chrom_position]
                                    if (label in labels):
                                        genotypes = labels[label]
                                        if (genotype in genotypes):
                                            genotypes[genotype] = genotypes[genotype] + 1
                                        else:
                                            genotypes[genotype] = 1
                                        labels[label] = genotypes
                                    else:
                                        labels[label] = {genotype:1}
                                    chrom_positions[chrom_position] = labels
                                else:
                                    chrom_positions[chrom_position] = {(label):{genotype:1}}
                                rsids[rsid] = chrom_positions
                            else:
                                rsids[rsid] = {(chrom_position):{label:{genotype:1}}}
                    if params[SHOWPROGRESSLINES] > 0 and lines_read % params[SHOWPROGRESSLINES] == 0:
                        elapsed = get_elapsed();
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
            