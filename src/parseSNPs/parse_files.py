"""
This program is designed to parse SNP files from OpenSNP.  It reads parameters to control
the parsing from a text file whose location is passed in as a parameter to the 
program.  The results are output to the console.  The heavy lifting is done by parser.py. 

Parameters are stored in a text file.  See parsefiles.txt 
for an example.
"""

import sys
from snp_utils import *
from parse_SNPs import DIR, FILES, RSID, CHROMSTART, CHROMEND, POSSTART, POSEND, SHOWFILEPROGRESS, SHOWPROGRESSLINES, SHOWSELECTEDFILES, parse_snps
from datetime import datetime

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
    results = parse_snps(params)

print "\n"

if (len(results) > 0):
#    for entry in sorted(results.items(), key=lambda t: rsid_key_seq(t[0])):
    for entry in sorted(results.items()):
        print entry[0], ":", entry[1], "           " # Add whitespace to ensure there are no leftover chars when we overprint the prior line.
else:
    print "Nothing matched selections"

elapsed = elapsed = get_elapsed();
print
print "Elapsed: %d minutes, %d seconds" % (elapsed[0], elapsed[1])


            