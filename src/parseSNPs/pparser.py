import sys
import json
import re #regular expressions


from string import maketrans

def Rev(s): 
    print(s[::-1].translate(maketrans('ACGT', 'TGCA')))

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)


def cleanString(s): return "".join(i for i in s if ord(i)<128)


if __name__ == "__main__":
    print sys.argv

#filename = sys.argv[1]
#get the lines out of the file

file=open("phenotypes.csv")
text=file.read()
file.close()
lines = text.split("\n")


header = lines[0].split(";")

separatedLine = {}

output = []


for i in range(len(header)):
    column = header[i]
    #replace all non whitespace \W+ with underscores
    #phenotype = (re.sub(r'[\W+]', '_', removeNonAscii(column.strip()))).lower()
    phenotype = column

    print "processing "+ str(i)+"/" +str(len(header)) +" phenotype: "+ phenotype

#create a new phenotype dictionary for each phenotype
    phenotypeDictionary = []

    for lineIndex in range(len(lines)):
        separatedLine = lines[lineIndex].split(";")
        if(i<len(separatedLine) and len(separatedLine[i].replace("-",""))>0):
#populate data if it is not blank
          dictionary = {'user_id': separatedLine[0], 'value': separatedLine[i]}
          phenotypeDictionary.append(dictionary)

    output.append({"phenotype":phenotype,"phenotype_id":i,"data": phenotypeDictionary})
          #print phenotype+": user: "+ separatedLine[0]+" value: "+separatedLine[i]
        


#create json without special formatting
with open('phenotypes_plain.json', 'w') as outfile:
  json.dump(output, outfile)
  outfile.close()

#create json with pretty print for readability
with open('phenotypes_pretty.json', 'w') as outfile:
    outfile.write(json.dumps(output, sort_keys=True, indent=4, separators=(',', ': ')))
    outfile.close()

print "Created henotypes_plain.json and phenotypes_pretty.json"


#writes pretty printed json string to file
#data = json.dumps({'4': 5, '6': 7}, sort_keys=True,indent=4, separators=(',', ': '))
#with open('data.txt', 'w') as outfile:
#  json.dump(data, outfile)



 
