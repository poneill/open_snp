import sys
import json
import re #regular expressions

from string import maketrans
from pprint import pprint



def Rev(s): 
    print(s[::-1].translate(maketrans('ACGT', 'TGCA')))

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)


def cleanString(s): return "".join(i for i in s if ord(i)<128)


if __name__ == "__main__":
    print sys.argv

#filename = sys.argv[1]
#get the lines out of the file


json_data = open('phenotypes_plain.json')
data = json.load(json_data)

json_data.close()

desiredPhenotypeID = 138 #138 is tongue roller


group1 ={"tongue roller","roller","yes"}
group2 = {"no","non-roller"}

#iterating over json file
for dictionary in data:
    
    if dictionary["phenotype_id"] == desiredPhenotypeID:
        print dictionary["phenotype"]
        
        for userDictionary in dictionary["data"]:

            for value in group1:
                if userDictionary["value"].lower() == value.lower():
                    print userDictionary["user_id"] + " has "+ value
                    #here we would assign this id to a group we are interested in
                    #group.add_file_selector('user' + userDictionary["user_id"] + '_*.txt')

            for value in group2:
                if userDictionary["value"].lower() == value.lower():
                    print userDictionary["user_id"] + " has "+ value
                    #group.add_file_selector('user' + userDictionary["user_id"] + '_*.txt')
    


 
