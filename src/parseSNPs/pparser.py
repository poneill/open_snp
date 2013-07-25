import sys
import json
import re #regular expressions
import csv
from string import maketrans

def Rev(s): 
    print(s[::-1].translate(maketrans('ACGT', 'TGCA')))

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)


def cleanString(s): return "".join(i for i in s if ord(i)<128)


if __name__ == "__main__":
    print sys.argv

#filename = sys.argv[1]


#json data that will be written to file
output = []
i = 0

#read CSV column by column, then row by row and create a JSON file from it, containing only phenotypes for which
#a value was provided

with open('phenotypes.csv', 'rb') as csvfile:
        csvreader = csv.DictReader(csvfile,dialect='excel', delimiter=';', quotechar='"')
        #phenotype is a column
        for phenotype in csvreader.fieldnames:
            
            phenotypeDictionary = []
            
            for row in csvreader:

                if(len(row[phenotype].replace("-",""))>0):
                    dictionary = {'user_id': row["user_id"], 'value': row[phenotype]}
                    phenotypeDictionary.append(dictionary)

            
            print "processing "+ str(i)+"/" +str(len(csvreader.fieldnames)) +" phenotype: "+ phenotype
            output.append({"phenotype":phenotype,"phenotype_id":i,"data": phenotypeDictionary})
            i = i+1
            csvfile.seek(0) #reset the file, otherwise csvdata iterates over rows only once
                

        


#create json without special formatting
with open('phenotypes_plain.json', 'w') as outfile:
  json.dump(output, outfile)
  outfile.close()

#create json with pretty print for readability
with open('phenotypes_pretty.json', 'w') as outfile:
    outfile.write(json.dumps(output, sort_keys=True, indent=4, separators=(',', ': ')))
    outfile.close()

print "Created henotypes_plain.json and phenotypes_pretty.json"


#now create a json dictionary for each user

output = []

with open('phenotypes.csv', 'rb') as csvfile:
            csvreader = csv.DictReader(csvfile,dialect='excel', delimiter=';', quotechar='"')
            #phenotype is a column
            
            fieldNamesList = list(csvreader.fieldnames)
            for row in csvreader:
                userDictionary = []
                i = 0

                for phenotype in fieldNamesList:

                    if len(row[phenotype].replace("-",""))>0:
                        dictionary = {phenotype: row[phenotype],"phenotype_id":i}
                        userDictionary.append(dictionary)
                    i = i + 1
                    
                print "processing user_id: "+ row["user_id"]
                output.append({"user_id":row["user_id"],"data": userDictionary})
 #           csvfile.seek(0) #reset the file, otherwise csvdata iterates over rows only once
                

 #create json without special formatting
with open('users_plain.json', 'w') as outfile:
  json.dump(output, outfile)
  outfile.close()

#create json with pretty print for readability
with open('users_pretty.json', 'w') as outfile:
    outfile.write(json.dumps(output, sort_keys=True, indent=4, separators=(',', ': ')))
    outfile.close()

print "Created users_plain.json and users_pretty.json"


#writes pretty printed json string to file
#data = json.dumps({'4': 5, '6': 7}, sort_keys=True,indent=4, separators=(',', ': '))
#with open('data.txt', 'w') as outfile:
#  json.dump(data, outfile)



 
