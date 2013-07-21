"""
These classes are used by parse_SNPs for parsing SNP files from OpenSNP and structuring
the results.

Author: David Gray
"""

import sys
import fnmatch
from snp_utils import *

####################################################################################
#
# A container for holding all the results
#
####################################################################################
class ResultsSet:
    # Constructor
    def __init__(self):
        self.results = {}

    # Get the length of the results    
    def __len__(self):
        return(len(self.results))
    
    # Convert contents to string
    def __str__(self):
        string_out = "( Results:"
        for result in sorted(self.get_results_iterator()):
            string_out += str(result)
        string_out += " )\n"
        return string_out
    
    # Get a result if found or create a result for the SnpValues instance passed in and return the result
    def get_or_create_result(self, rsid, chromosome, position):
        key = Result.get_key_static(rsid, chromosome, position)
        if (key in self.results):
            result = self.results[key]
        else:
            result = Result(rsid, chromosome, position)
            self.results[key] = result
        return result
    
    # return an iterator of the results
    def get_results_iterator(self):
        return self.results.itervalues()
        
####################################################################################
#
# Container for a result pertaining to an rsid, chromosome number and position
#
####################################################################################
class Result:
    
    # Constructor
    def __init__(self, rsid, chromosome, position):
        self.rsid = rsid
        self.chromosome = chromosome
        self.position = position
        self.groups = {}
    
    # Convert the contents to a String
    def __str__(self):
        string_out = "\n   "
        string_out += "( Result: rsid " + str(self.rsid) 
        string_out += ", chromosome " + str(self.chromosome) 
        string_out += ", position " + str(self.position)
        for entry in sorted(self.groups.items()):
            string_out += "\n      " + str(entry[1])   
        string_out += "\n   " + ")\n"
        return string_out
    
    # Statically get a key so it can be accessed without a Result instance
    @staticmethod
    def get_key_static(rsid, chromosome, position):
        return (rsid, chromosome, position)
    
    # Get a key that can be used when placing these in a dictionary
    def get_key(self):
        return Result.get_key_static(self.rsid, self.chromosome, self.position)
    
    # Get the RSID for these results
    def get_rsid (self):
        return self.rsid
    
    # Get the chromosome number for these results
    def get_chromosome (self):
        return self.chromosome
    
    # Get the position for these results
    def get_position (self):
        return self.position
    
    # Add one genotype for a given group label
    def add_one (self, label, gtype):
        if label not in self.groups:
            self.groups[ label ] = Group(label)
        self.groups[ label ].add_genotype(gtype)
    
    # Get the results for a group given the label
    def get_group (self, label):
        if label in self.groups:
            return self.groups[ label ]
        else:
            return Group(label)
    
    # Get all the groups as a dictionary of {label : group}
    def get_groups(self):
        return self.groups

####################################################################################
#
# Class to manage results for a file group.  
#
####################################################################################
class Group:
    # Constructor
    def __init__(self, label):
        self.label = label
        self.gtypes = {}
    
    # Convert the contents to a string
    def __str__(self):
        string_out = "( Group: label '" + self.label + "': "
        group_info = ""
        for entry in sorted(self.gtypes.items()):
            if len(group_info) > 0:
                group_info += ", "
            group_info += entry[0] + "=" + str(entry[1])   
        string_out += group_info + " )"
        return string_out
    
    # Get the group label    
    def get_label (self):
        return self.label
    
    # Add a single genotype instance for the group    
    def add_genotype (self, gtype):
        if gtype in self.gtypes:
            self.gtypes[gtype] += 1
        else:
            self.gtypes[gtype] = 1
    
    # Get the count of a single genotype in the group
    def get_count (self, gtype):
        if gtype in self.gtypes:
            return self.gtypes[gtype]
        else:
            return 0
            
    # Get all the genotype counts in the group as a dictionary in the format {genotype : count}
    def get_counts(self):
        return self.gtypes
    
####################################################################################
#
# Class to manage the execution parameters for parse_SNPs
#
####################################################################################
class Params:
    # Constructor
    def __init__(self):
        self.dir = "."
        self.rsid = "*"
        self.chromosomes = None
        self.pos_start = 0
        self.pos_end = sys.maxint
        self.show_lines_progress_interval = 0
        self.show_file_progress = False
        self.show_selected_files = False
        self.file_groups = []
        # {{0,"Default"}, ["*"]}
    
    # Convert the contents to a string
    def __str__(self):
        string_out = "( Params: \n   rsid " + self.rsid
        string_out += "\n   chromosomes " + str(self.chromosomes)
        string_out += "\n   pos_start " + str(self.pos_start) 
        string_out += "\n   pos_end " + str(self.pos_end) 
        string_out += "\n   show_file_progress " + str(self.show_file_progress) 
        string_out += "\n   show_selected_files " + str(self.show_selected_files) 
        string_out += "\n   show_lines_progress_interval " + str(self.show_lines_progress_interval) 
        string_out += "\n   dir " + self.dir
        string_out += "\n   file groups:"
        for file_group in self.file_groups:
            string_out += str(file_group)  
        string_out += "\n)"
        return string_out
                
    # Add a chromosome to include.  Can be a chromosome number or letter such as X or Y.  
    # If no chromosomes are added, all are processed.    
    def add_chromosome (self, chromosome):
        if self.chromosomes == None:
            self.chromosomes = [ chromosome ]
        else:
            self.chromosomes.append(chromosome)
    
    # Add a file group
    def add_file_group (self, file_group):
        i = 0
        inserted = False
        while i < len (self.file_groups):
            # Keep groups in priority sequence
            if file_group.get_priority_seq() < self.file_groups[i].get_priority_seq():
                self.file_groups.insert(i, file_group)
                inserted = True
                break
            i += 1
        if not inserted:
            self.file_groups.append(file_group)
                
    # Get the chromosomes to include.  If None, all are included
    def get_chromosomes (self):
        return self.chromosomes
    
    # Get the directory containing the SNP files.  E.g. 'C:\\OpenSNP'
    def get_directory_location (self):
        return self.dir
    
    # If the filename passed in matches a file group, return the group label    
    def get_file_group_label (self, file_name):
        label = "Default"
        if len(self.file_groups) > 0:
            label = None
            for file_group in self.file_groups:
                if file_group.matches(file_name):
                    label = file_group.get_label()
                    break
        return label
    
    # Get the ending chromosome position to include. Ignore all positions in chromosomes greater than this value.
    def get_position_end (self):
        return self.pos_end
    
    # Get the starting chromosome position to include. Ignore all chromosome numbers below this value.
    def get_position_start (self):
        return self.pos_start
    
    # Get the pattern of RSIDs to process.  Allows the selections to be limited 
    # to one or more specific or all RSIDs. The value can be specified with 
    # wild cards. E.g. RSID10403190 or RSID104*
    def get_rsid (self):
        return self.rsid
    
    # If True and get_show_lines_progress_interval is zero, show progress information as each new file is processed
    def get_show_file_progress (self):
        return self.show_file_progress
    
    # Get the line progress interval.  For example, if this value is 100, update the line progress every 100th line.
    def get_show_lines_progress_interval (self):
        return self.show_lines_progress_interval
    
    # Get an indicator whether selected files should be listed to the console.
    # If True, also list files with no data to process that passes the selections in this file    
    def get_show_selected_files (self):
        return self.show_selected_files
    
    # Determine whether a SnpValues instance should be processed
    def process (self, snp_values):
        chromosome_ok = False
        if(self.chromosomes):
            for chromosome in self.chromosomes:
                chromosome_ok = fnmatch.fnmatch(snp_values.get_chromosome(), chromosome)
                if chromosome_ok:
                    break;
        else:
            chromosome_ok = True
        return fnmatch.fnmatch(snp_values.get_rsid(), self.rsid) \
            and chromosome_ok \
            and snp_values.get_position() >= self.pos_start \
            and snp_values.get_position() <= self.pos_end
    
    # Set the directory containing the SNP files.  E.g. 'C:\\OpenSNP'
    def set_directory_location (self, dir):
        self.dir = dir.strip()
    
    # Set the ending chromosome position to include. Ignore all positions in chromosomes greater than this value.
    def set_position_end (self, pos_end):
        self.pos_end = pos_end
    
    # Set the starting chromosome position to include.  Ignore all positions in chromosomes less than this value.
    def set_position_start (self, pos_start):
        self.pos_start = pos_start
    
    # Set pattern of RSIDs to process.  Allows the selections to be limited 
    # to one or more specific or all RSIDs. The value can be specified with 
    # wild cards. E.g. RSID10403190 or RSID104*
    def set_rsid (self, rsid):
        self.rsid = rsid.strip().upper()
    
    # If True and get_show_lines_progress_interval is zero, show progress information as each new file is processed
    def set_show_file_progress (self, show_file_progress):
        self.show_file_progress = show_file_progress
    
    # Set the line progress interval.  For example, if this value is 100, update the line progress every 100th line.
    def set_show_lines_progress_interval (self, show_lines_progress_interval):
        self.show_lines_progress_interval = show_lines_progress_interval
    
    # Set an indicator whether selected files should be listed to the console.
    # If True, also list files with no data to process that passes the selections in this file    
    def set_show_selected_files (self, show_selected_files):
        self.show_selected_files = show_selected_files

####################################################################################
#
# Class to manage file groups passed as parameters to parse_SNPs
#
####################################################################################
class FileGroup:
    # Constructor
    def __init__(self, label, priority_seq):
        self.label = label
        self.priority_seq = priority_seq
        self.files = []
    
    # Convert the contents to a string
    def __str__(self):
        string_out = "\n      ( FileGroup: label " + self.label
        string_out += ", priority_seq " + str(self.priority_seq) 
        string_out += ", file selectors:"
        for entry in sorted(self.files):
            string_out += "\n         " + entry   
        string_out += "\n      )"
        return string_out
    
    # Add a file selector to the group
    def add_file_selector (self, file_selector):
        self.files.append(file_selector.strip())
    
    # Get the file group label
    def get_label (self):
        return self.label
    
    # Get the file group priority sequence    
    def get_priority_seq (self):
        return self.priority_seq
    
    # Return True if the file name matches a file selector in this group
    def matches(self, file_name):
        matches = False;
        for file_selector in self.files:
            if (fnmatch.fnmatch(file_name, file_selector)):
                matches = True
                break
        return matches

####################################################################################
#
# Immutable class to manage file groups passed as parameters to parse_SNPs
#
####################################################################################
class SnpValues:
    
    # Constructor
    def __init__(self, rsid, chromosome, position, genotype):
        self.rsid = rsid.upper()
        self.chromosome = chromosome
        self.position = position
        self.genotype = genotype.upper()
    
    # Convert the contents to a string
    def __str__(self):
        string_out = "( SnpValues: RSID " + self.rsid
        string_out += ", chromosome " + str(self.chromosome) 
        string_out += ", position " + str(self.position) 
        string_out += ", genotype " + str(self.genotype) 
        string_out += " )"
        return string_out
    
    # Get the RSID    
    def get_rsid (self):
        return self.rsid
    
    # Get the chromosome number    
    def get_chromosome (self):
        return self.chromosome 
    
    # Get the position on the chromosome    
    def get_position (self):
        return self.position
    
    # Get the genotype
    def get_genotype (self):
        return self.genotype

####################################################################################
#
# Abstract class defining methods to process SNP files
#
####################################################################################
class AbstractSNPProcessor(object):
    
    # Return true if this processor is appropriate for the file passed in
    def handles_file(self, filename):
        raise NotImplementedError("Should have implemented parse_line")
    
    # Get a label to use when summarizing counts for this file type
    def get_file_type_label(self):
        raise NotImplementedError("Should have implemented parse_line")
    
    # Parse a SNP file line of data for this file type
    def parse_line(self, line):
        raise NotImplementedError("Should have implemented parse_line")
    
    # Get a processor given the file name
    @staticmethod
    def get_processor(filename):
        processor_out = None
        for processor in processors:
            if (processor.handles_file(filename)):
                processor_out = processor
                break
        return processor_out

####################################################################################
#
# Concrete class to process 23andme files
#
####################################################################################
class TwentyThreeAndMeSNPProcessor(AbstractSNPProcessor):
    
    # Constructor
    def __init__(self):
        pass
    
    # Return true if this processor is appropriate for the file passed in
    def handles_file(self, filename):
        return fnmatch.fnmatch(filename, "*23andme.txt")
    
    # Get a label to use when summarizing counts for this file type
    def get_file_type_label(self):
        return "23andme"
    
    # Parse a SNP file line of data for this file type
    # 23andme structure (tab separated):
    # rsid      chromosome    pos    genotype
    # rs4477212    1         72017    AA
    def parse_line(self, line):
        snp_values = None
        data = map(strip, line.split("\t"))
        if (len(data) == 4):
            try:
                snp_values = SnpValues(data[0], data[1], int(data[2]), data[3])
            except ValueError:
                pass  # Nothing to do.  Returning None handles the issue
        return snp_values

####################################################################################
#
# Concrete class to process illumina files
#
####################################################################################
class IlluminaSNPProcessor(AbstractSNPProcessor):
    
    # Constructor
    def __init__(self):
        pass
    
    # Return true if this processor is appropriate for the file passed in
    def handles_file(self, filename):
        return fnmatch.fnmatch(filename, "*illumina.txt")
    
    # Get a label to use when summarizing counts for this file type
    def get_file_type_label(self):
        return "illumina"
    
    # Parse a SNP file line of data for this file type
    def parse_line(self, line):
        snp_values = None
        try:
            # illumina can have two structures:
            if (line[0:1] == "\""):
                # RSID,CHROMOSOME,POSITION,RESULT
                # "rs3094315","1","742429","AA"
                data = map(strip_quotes, line.strip().split(",")) 
                if (len(data) == 4):
                    snp_values = SnpValues(data[0], data[1], int(data[2]), data[3])
            else:
                # rsid        chromosome   position    allele1    allele2
                # rs4477212        1        82154    T    T
                data = line.strip().split("\t")
                if (len(data) == 5):
                    # Example user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt
                    snp_values = SnpValues(data[0], data[1], int(data[2]), data[3] + data[4])
                else:
                    # rsid chromosome position allele1 allele2
                    # rs11240777 1 788822 AA
                    data = line.strip().split(" ")
                    if (len(data) == 4):
                        # Example user981_file487_yearofbirth_1966_sex_unknown.ftdna-illumina.txt
                        snp_values = SnpValues(data[0], data[1], int(data[2]), data[3])
        except ValueError:
            pass  # Nothing to do.  Returning None handles the issue
        return snp_values

####################################################################################
#
# Concrete class to process IYG files
#
####################################################################################
class IYGSNPProcessor(AbstractSNPProcessor):
    
    # Constructor
    def __init__(self):
        pass
    
    # Return true if this processor is appropriate for the file passed in
    def handles_file(self, filename):
        return fnmatch.fnmatch(filename, "*IYG.txt")
    
    # Get a label to use when summarizing counts for this file type
    def get_file_type_label(self):
        return "iyg"
    
    # Parse a SNP file line of data for this file type
    def parse_line(self, line):
        # IYG structure:
        # RSID,RESULT
        # rs2131925    TT
        snp_values = None
        data = strip_quotes(line.strip()).split("\t")
        if (len(data) == 2):
            try:
                snp_values = SnpValues(data[0], "", 0, data[1])
            except ValueError:
                pass  # Nothing to do.  Returning None handles the issue
        return snp_values

####################################################################################
#
# Concrete class to process decodeme files
#
####################################################################################
class DecodeMeSNPProcessor(AbstractSNPProcessor):
    
    # Constructor
    def __init__(self):
        pass
    
    # Return true if this processor is appropriate for the file passed in
    def handles_file(self, filename):
        return fnmatch.fnmatch(filename, "*decodeme.txt")
    
    # Get a label to use when summarizing counts for this file type
    def get_file_type_label(self):
        return "decodeme"
    
    # Parse a SNP file line of data for this file type
    def parse_line(self, line):
        # decodeme structure:
        # Name,Variation,Chromosome,Position,Strand,YourCode
        # rs4345758,C/T,1,28663,+,TT        
        snp_values = None
        data = line.strip().split(",")
        if (len(data) == 6):
            try:
                snp_values = SnpValues(data[0], data[2], int(data[3]), data[5])
            except ValueError:
                pass  # Nothing to do.  Returning None handles the issue
        return snp_values
    

####################################################################################
#
# Globals    
#
####################################################################################
processors = [ TwentyThreeAndMeSNPProcessor(), IlluminaSNPProcessor(), IYGSNPProcessor(), DecodeMeSNPProcessor() ]

def strip(string):
    if string != None:
        string = string.strip()
    return string
