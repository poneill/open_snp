"""
This program is designed to parse SNP files from OpenSNP.  

Parameters are stored in an text file whose location is passed in as a parameter to the 
program.  The file is tab-separated and contains identifier/value pairs.  See parsefiles.txt 
for an example.
"""
import sys
import fnmatch
import os
from snp_classes import *
from snp_utils import get_elapsed, string_to_bool, inc_dict_counter
from datetime import datetime

def bypass(filename):
    return "-exome-" in filename

# Main processing method.  The one parameter, "parms" is an instance of the Params class.
# The returned value is a ResultsSet instance
def parse_snps(params):
    files = 0
    print "Processing files"
    directory = os.listdir(params.getDir())
    dir_count = 0
    fileTypeCounts={}
    skipped_files = []
    results_set = ResultsSet()   
    # Count number of files to process for progress reporting
    for filename in directory:
        if (bypass(filename)):
            continue
        if (params.getFileGroupLabel(filename) != None):
            if params.getShowSelectedFiles(): 
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
        label = params.getFileGroupLabel(filename)
        if (label != None):
            if (bypass(filename)):
                skipped_files.append(filename)
                continue
            with open(os.path.join(params.getDir(), filename)) as f:
                files += 1
                if params.getShowFileProgress() and params.getShowLinesProgressInterval() <= 0:
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
                        if chrom > params.getChromEnd:
                            exit_file = True
                            continue
                        if (fnmatch.fnmatch(rsid, params.getRSID()) and chrom >= params.getChromStart()
                            and pos >= params.getPosStart() and pos <= params.getPosEnd()):
                            lines_processed += 1
                            
                            # Add to the results
                            result = results_set.getResultAlways(rsid, chrom, pos)
                            result.addOne(label, genotype)
                            
                    if params.getShowLinesProgressInterval() > 0 and lines_read % params.getShowLinesProgressInterval() == 0:
                        elapsed = get_elapsed();
                        print "file: {:,d}/{:,d}  lines read: {:,d}  processed {:,d}  elapsed: {:,d} minutes, {:,d} seconds\r".format(files, dir_count, lines_read, lines_processed, int(elapsed[0]), int(elapsed[1])) 
                        sys.stdout.flush()
                    line = f.readline()
                if (params.getShowSelectedFiles() and lines_processed == 0):
                    print "    No lines to process in ", filename
    if(len(skipped_files) > 0):
        print 
        print "Skipped Files"
        for file in skipped_files:
            print file
    return results_set
            