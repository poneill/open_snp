import sys
import fnmatch

####################################################################################
#
# A container for holding all the results
#
####################################################################################
class ResultsSet:
	# Constructor
	def __init__( self ):
		self.results = {}

	# Get the length of the results	
	def __len__( self ):
		return(len(self.results))
	
	# Convert contents to string
	def __str__( self ):
		string_out = "( Results:"
		for result in sorted( self.getResultsIter() ):
			string_out += str( result )
		string_out += " )\n"
		return string_out
	
	# Get a result if it exists.  Otherwise return None
	def getResultOrNone( self, rsid, chromosome, position):
		key = Result.getKeyStatic(rsid, chromosome, position)
		result = None
		if (key in self.results):
			result = results[key]
		return result	
	
	# Get a result if found or create a result for the parms passed in and return it
	def getResultAlways( self, rsid, chromosome, position):
		key = Result.getKeyStatic(rsid, chromosome, position)
		if (key in self.results):
			result = self.results[key]
		else:
			result = Result(rsid, chromosome, position)
			self.results[key] = result
		return result
	
	# return an iterator of the results
	def getResultsIter( self ):
		return self.results.itervalues()
		
####################################################################################
#
# Container for a result pertaining to an rsid, chromosome number and position
#
####################################################################################
class Result:
	
	# Constructor
	def __init__( self, rsid, chromosome, position):
		self.rsid = rsid
		self.chromosome = chromosome
		self.position = position
		self.groups = {}
	
	# Convert the contents to a String
	def __str__( self ):
		string_out = "\n   " + "( Result: rsid " + str(self.rsid) + ", chromosome " + str(self.chromosome) + ", position " + str(self.position) + ": "
		for entry in sorted(self.groups.items()):
			string_out += "\n      " + str(entry[1])   
		string_out += "\n   " + ")\n"
		return string_out
	
	# Statically get a key so it can be accessed without a Result instance
	@staticmethod
	def getKeyStatic( rsid, chromosome, position):
		return (rsid, chromosome, position)
	
	# Get a key that can be used when placing these in a dictionary
	def getKey( self ):
		return getKey(self.rsid, self.chromosome, self.position)
	
	# Get the RSID for these results
	def getRsid ( self ):
		return rsid
	
	# Get the chromosome number for these results
	def getChromosome ( self ):
		return chromosome
	
	# Get the postition for these results
	def getPosition ( self ):
		return position
	
	# Add one genotype for a given group label
	def addOne ( self, label, gtype ):
		if label not in self.groups:
			self.groups[ label ] = Group( label )
		self.groups[ label ].addGtype( gtype )
	
	# Get the results for a group given the label
	def getGroup ( self, label ):
		if label in self.gtypes:
			return self.gtypes[ label ]
		else:
			return Group( label )
	
	# Get all the groups as a dictionary of {label : group}
	def getGroups( self ):
		return self.groups

####################################################################################
#
# Class to manage results for a file group.  
#
####################################################################################
class Group:
	# Constructor
	def __init__( self, label):
		self.label = label
		self.gtypes = {}
	
	# Convert the contents to a string
	def __str__( self ):
		string_out = "( Group: label '" + self.label + "': "
		first = True
		for entry in sorted(self.gtypes.items()):
			if not first:
				string_out += ", "
			first = False
        	string_out += entry[0] + "=" + str(entry[1])   
		string_out += " )"
		return string_out
	
	# Get the group label	
	def getLabel ( self ):
		return label
	
	# Add a single genotype instance for the group	
	def addGtype ( self, gtype ):
		if gtype in self.gtypes:
			self.gtypes[gtype] += 1
		else:
			self.gtypes[gtype] = 1
	
	# Get the count of a single genotype in the group
	def getCount ( self, gtype ):
		if gtype in self.gtypes:
			return self.gtypes[gtype]
		else:
			return 0
			
	# Get all the genotype counts in the group as a dictionary in the format {genotype : count}
	def getCounts( self ):
		return self.gtypes
	
####################################################################################
#
# Class to manage the execution parameters for parse_SNPs
#
####################################################################################
class Params:
	# Constructor
	def __init__( self ):
		self.dir = "."
		self.rsid = "*"
		self.chrom_start = 0
		self.chrom_end = 23
		self.pos_start = 0
		self.pos_end = sys.maxint
		self.show_progress_lines = 0
		self.show_file_progress = False
		self.show_selected_files = False
		self.file_groups = []
		# {{0,"Default"}, ["*"]}
	
	# Convert the contents to a string
	def __str__( self ):
		string_out = "( Params: \n   rsid " + self.rsid
		string_out += "\n   chrom_start " + str( self.chrom_start ) 
		string_out += "\n   chrom_end " + str( self.chrom_end ) 
		string_out += "\n   pos_start " + str( self.pos_start ) 
		string_out += "\n   pos_end " + str( self.pos_end ) 
		string_out += "\n   show_file_progress " + str( self.show_file_progress ) 
		string_out += "\n   show_selected_files " + str( self.show_selected_files ) 
		string_out += "\n   show_lines_progress_interval " + str( self.show_lines_progress_interval ) 
		string_out += "\n   dir " + self.dir
		string_out += "\n   file groups:"
		for file_group in self.file_groups:
			string_out += str( file_group )  
		string_out += "\n)"
		return string_out
	
	# Add a file group
	def addFileGroup ( self, file_group ):
		i = 0
		inserted = False
		while i < len ( self.file_groups ):
			# Keep groups in priority sequence
			if file_group.getPriority() < self.file_groups[i].getPriority():
				self.file_groups.insert(i, file_group)
				inserted = True
				break
			i += 1
		if not inserted:
			self.file_groups.append(file_group)
				
	# Get the ending chromosome number to include.  Ignore all chromosome numbers greater than this value.
	def getChromEnd ( self ):
		return self.chrom_end
	
	# Get the starting chromosome number to include.  Ignore all chromosome numbers below this value.
	def getChromStart ( self ):
		return self.chrom_start
	
	# Get the directory containing the SNP files.  E.g. 'C:\\OpenSNP'
	def getDir ( self ):
		return self.dir
	
	# If the filename passed in matches a file group, return the group label	
	def getFileGroupLabel ( self, file_name ):
		label = "Default"
		if len(self.file_groups) > 0:
			label = None
			for file_group in self.file_groups:
				if file_group.matches( file_name ):
					label = file_group.getLabel()
					break
		return label
	
	# Get the ending chromosome position to include. Ignore all positions in chromosomes greater than this value.
	def getPosEnd ( self ):
		return self.pos_end
	
	# Get the starting chromosome position to include. Ignore all chromosome numbers below this value.
	def getPosStart ( self ):
		return self.pos_start
	
	# Get the pattern of RSIDs to process.  Allows the selections to be limited 
	# to one or more specific or all RSIDs. The value can be specified with 
	# wild cards. E.g. RSID10403190 or RSID104*
	def getRSID ( self ):
		return self.rsid
	
	# If True and getShowLinesProgressInterval is zero, show progress information as each new file is processed
	def getShowFileProgress ( self ):
		return self.show_file_progress
	
	# Get the line progress interval.  For example, if this value is 100, update the line progress every 100th line.
	def getShowLinesProgressInterval ( self ):
		return self.show_lines_progress_interval
	
	# Get an indicator whether selected files should be listed to the console.
	# If True, also list files with no data to process that passes the selections in this file	
	def getShowSelectedFiles ( self ):
		return self.show_selected_files
				
	# Set the ending chromosome number to include.Ignore all chromosome numbers greater than this value.	
	def setChromEnd ( self, chrom_end ):
		self.chrom_end = chrom_end
	
	# Set the starting chromosome number to include.  Ignore all chromosome numbers below this value.
	def setChromStart ( self, chrom_start ):
		self.chrom_start = chrom_start
	
	# Set the directory containing the SNP files.  E.g. 'C:\\OpenSNP'
	def setDir ( self, dir ):
		self.dir = dir.strip().upper()
	
	# Set the ending chromosome position to include. Ignore all positions in chromosomes greater than this value.
	def setPosEnd ( self, pos_end ):
		self.pos_end = pos_end
	
	# Set the starting chromosome position to include.  Ignore all positions in chromosomes less than this value.
	def setPosStart ( self, pos_start ):
		self.pos_start = pos_start
	
	# Set pattern of RSIDs to process.  Allows the selections to be limited 
	# to one or more specific or all RSIDs. The value can be specified with 
	# wild cards. E.g. RSID10403190 or RSID104*
	def setRSID ( self, rsid ):
		self.rsid = rsid.strip().upper()
	
	# If True and getShowLinesProgressInterval is zero, show progress information as each new file is processed
	def setShowFileProgress ( self, show_file_progress ):
		self.show_file_progress = show_file_progress
	
	# Set the line progress interval.  For example, if this value is 100, update the line progress every 100th line.
	def setShowLinesProgressInterval ( self, show_lines_progress_interval ):
		self.show_lines_progress_interval = show_lines_progress_interval
	
	# Set an indicator whether selected files should be listed to the console.
	# If True, also list files with no data to process that passes the selections in this file	
	def setShowSelectedFiles ( self, show_selected_files ):
		self.show_selected_files = show_selected_files

####################################################################################
#
# Class to manage file groups passed as parameters to parse_SNPs
#
####################################################################################
class FileGroup:
	# Constructor
	def __init__( self ):
		self.label = "label"
		self.priority = 1
		self.files = []
	
	# Convert the contents to a string
	def __str__( self ):
		string_out = "\n      ( FileGroup: label " + self.label
		string_out += ", priority " + str(self.priority) 
		string_out += ", file selectors:"
		for entry in sorted(self.files):
			string_out += "\n         " + entry   
		string_out += "\n      )"
		return string_out
	
	# Add a file selector to the group
	def addFileSelector ( self, file_selector ):
		self.files.append( file_selector.strip() )
	
	# Get the file group label
	def getLabel ( self ):
		return self.label
	
	# Get the file group priority	
	def getPriority ( self ):
		return self.priority
	
	# Return True if the file name matches a file selector in this group
	def matches( self, file_name ):
		matches = False;
		for file_selector in self.files:
			if (fnmatch.fnmatch(file_name, file_selector)):
				matches = True
				break
		return matches
	
	# Set the file group label
	def setLabel ( self, label ):
		self.label = label.strip()
	
	# Set the file group priority	
	def setPriority ( self, priority ):
		self.priority = priority
