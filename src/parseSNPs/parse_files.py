"""
This program is designed to parse SNP files from OpenSNP.  It reads parameters to control
the parsing from a text file whose location is passed in as a parameter to the 
program.  The results are output to the console.  The heavy lifting is done by parser.py. 

Parameters are stored in a text file.  See parsefiles.txt 
for an example.
"""
from snp_classes import *
from snp_utils import *
from parse_SNPs import parse_snps
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
    params = Params()
    filename = sys.argv[1]
    f=file(filename,"r")
    for line in f:
        if line[0:1] != "#": # Skip comment lines
            values = line.strip().split("\t")
            if (len(values) > 1 and values[0].strip()): # Must have a non-blank name and a value
                name = values[0].strip()
                val = values[1].strip()
                if (name.startswith("FILES")):
                    # The "FILES" line int the parameter file will look like this:
                    #   FILES:Group 1:1    user10_*.txt, user11_*.txt, user13_*.txt, user14_*.txt
                    # The line is in two sections, separated by a tab.
                    # The first section containers the group label and priority, separated by a colon
                    # The second section lists file selectors for the group, separated by commas.
                    file_group = FileGroup()
                    label_and_priority = name.split(":")
                    if len(label_and_priority) > 2:
                        file_group.set_label(label_and_priority[1])         
                        file_group.set_priority(int(label_and_priority[2]))                  
                    for file_selector in val.split(","):
                        file_group.add_file_selector(file_selector)
                    params.add_file_group(file_group)
                else:
                    # Most lines in the file are simple name/value pairs separated by a tab.  E.g.:
                    # RSID    RS10403190
                    name = name.upper()
                    if(name == "DIR"):
                        params.set_directory_location(val)
                    elif( name == "RSID"):
                        params.set_rsid(val)
                    elif( name == "CHROMOSOMES"):
                        if ( len( val ) > 0 ):
                            for chromosome in map(strip, val.split(",")):
                                params.add_chromosome( chromosome )
                    elif( name == "POSSTART"):
                        params.set_position_start(int(val))
                    elif( name == "POSEND"):
                        params.set_position_end(int(val))
                    elif( name == "SHOWFILEPROGRESS"):
                        # True can be represented by "TRUE", "T", "1", "YES" or "Y" in any case.
                        params.set_show_file_progress(string_to_bool(val))
                    elif( name == "SHOWPROGRESS#LINES"):
                        params.set_show_lines_progress_interval(int(val))
                    elif( name == "SHOWSELECTEDFILES"):
                        # True can be represented by "TRUE", "T", "1", "YES" or "Y" in any case.    
                        params.set_show_selected_files(string_to_bool(val))
    f.close()
    results_set = parse_snps(params)

    print "\n"
    print str(params)
    print "\n" 
    if (len(results_set) > 0):
    #    for entry in sorted(results.items(), key=lambda t: rsid_key_seq(t[0])):
        print str( results_set )
    else:
        print "Nothing matched selections"
    
    elapsed = elapsed = get_elapsed();
    print
    print "Elapsed: %d minutes, %d seconds" % (elapsed[0], elapsed[1])
