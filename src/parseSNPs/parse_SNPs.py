"""
This program is designed to parse SNP files from OpenSNP.  

Parameters that control what data is selected are passed to parse_snps in a Params instance.  
(See snp_classes)
"""
import sys
import fnmatch
import os
from snp_classes import *
from snp_utils import get_elapsed, string_to_bool, increment_dictionary_counter
from datetime import datetime

def bypass(filename):
    return "-exome-" in filename

# Main processing method.  The one parameter, "parms" is an instance of the Params class.
# The returned value is a ResultsSet instance
def parse_snps(params):
    print "Processing files"
    files = 0
    directory_files = os.listdir(params.get_directory_location())
    dir_count = 0
    fileTypeCounts = {}
    skipped_files = []
    
    # Count number of files to process for progress reporting
    for filename in directory_files:
        if (bypass(filename)):
            continue
        if (params.get_file_group_label(filename) != None):
            if params.get_show_selected_files(): 
                print "Selected:", filename
            dir_count += 1
            processor = AbstractSNPProcessor.get_processor(filename);
            if (processor):
                counter_name = processor.get_file_type_label()
                increment_dictionary_counter( fileTypeCounts, counter_name )
                    
    # Summarize counts by file type
    print
    for count in fileTypeCounts.iteritems():
        print count[0], ":", count[1]
    print
            
    # Process each selected file.  Returns a ResultsSet instance
    results_set = ResultsSet()
    for filename in directory_files:
        label = params.get_file_group_label(filename)
        if (label != None):
            if (bypass(filename)):
                skipped_files.append(filename)
                continue
            with open(os.path.join(params.get_directory_location(), filename)) as f:
                files += 1
                if params.get_show_file_progress() and params.get_show_lines_progress_interval() <= 0:
                    elapsed = snp_utils.get_elapsed();
                    print "files: %d/%d  elapsed: %d minutes, %d seconds          \r" % (files, dir_count, elapsed[0], elapsed[1]),
                    sys.stdout.flush()
                lines_processed = lines_read = 0
                line = f.readline()
                exit_file = False
                while (line != "" and not exit_file):
                    lines_read += 1
                    snp_values = None
                    if not line.startswith( "#" ):
                        snp_values = AbstractSNPProcessor.get_processor(filename).parse_line(line)
                        if ( snp_values == None ) or ( len(snp_values.get_rsid()) > 20 ):
                            skipped_files.append(filename)
                            exit_file = True
                            
                        if not exit_file:
                            
                            if ( params.process( snp_values )):
                                lines_processed += 1
                                
                                # Add to the results
                                result = results_set.get_result_always(snp_values)
                                result.add_one(label, snp_values.get_genotype())
                            
                    if params.get_show_lines_progress_interval() > 0 and lines_read % params.get_show_lines_progress_interval() == 0:
                        elapsed = get_elapsed();
                        print "file: {:,d}/{:,d}  lines read: {:,d}  processed {:,d}  elapsed: {:,d} minutes, {:,d} seconds          \r".format(files, dir_count, lines_read, lines_processed, int(elapsed[0]), int(elapsed[1])), 
                        sys.stdout.flush()
                    line = f.readline()
                if (params.get_show_selected_files() and lines_processed == 0):
                    print "    No lines to process in ", filename
    if(len(skipped_files) > 0):
        print 
        print "Skipped Files"
        for file in skipped_files:
            print file
    return results_set
            
