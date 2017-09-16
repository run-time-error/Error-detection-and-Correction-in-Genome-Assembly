#0 ->reference position
#1 ->reference position e ki paise
#2 ->query position e ki paise
#3 ->query position
#14 -> scaffold number

deleteCount = 0
insertCount = 0
singleCharacterInsertion = 0
singleCharacterDeletion = 0
SNPCount = 0

insertedString = ""
deletedString = ""

def checkDelete ( words2 , words3 ):
    global deletedString,deleteCount
    if words2[0] != words3[0] and words2[3] != words3[3]:
    #deletion ends here
        deletedString+=words2[1]
        fout.write(deletedString + " " + str(len(deletedString)) + "\n" + 
        "deletion ends reference " + words2[0] + " and query " + words2[3] + " " + words2[14] + "\n\n")
        deleteCount+=1
        deletedString=""
    elif words2[0] != words3[0] and words2[3] == words3[3]:
        #middle of deletion 
        deletedString+=words2[1]
    elif words2[0] == words3[0] and words2[3] != words3[3]:
        #deletion end 
        deletedString += words2[1]
        fout.write(deletedString + " " + str(len(deletedString)) + "\n" + 
        "deletion ends reference " + words2[0] + " and query " + words2[3] + " " + words2[14] + "\n\n")
        deleteCount+=1
        deletedString=""     
        return 
        
def checkInsert ( words2, words3 ):
    global insertedString,insertCount
    if words2[0] != words3[0] and words2[3] != words3[3]:
        #end of insert
        insertedString+=words2[2]
        fout.write(insertedString + " " + str(len(insertedString)) + "\n" + 
        "insertion ends reference " + words2[0] + " query " + words2[3] +  " " + words2[14] + "\n\n")
        insertCount+=1
        insertedString=""
    
    elif words2[0] != words3[0] and words2[3] == words3[3]:
        #insert end
        insertedString+=words2[2]
        fout.write(insertedString + " " + str(len(insertedString)) + "\n" + 
        "insertion ends reference " + words2[0] + " query " + words2[3] +  " " + words2[14] + "\n\n")
        insertCount+=1
        insertedString=""
    elif words2[0] == words3[0] and words2[3] != words3[3]:
        #insert middle
        insertedString+=words2[2]
    return

def checkSNPs (words2, words3 ):
    global insertedString,SNPCount
    global deletedString,singleCharacterDeletion,singleCharacterInsertion
    if(words2[0] != words3[0] and words2[3] != words3[3]):
        if(words2[1].isalpha() and words2[2].isalpha()): 
            fout.write("SNP at reference " + str(words2[0]) + " and query " + str(words2[3]) + " {" + 
            words2[1] + " | " + words2[2] + "}" +  " " + words2[14] + "\n\n" )
            SNPCount+=1
        elif(words2[1].isalpha()==False and words2[2].isalpha()):
            fout.write("Single character insertion at reference " + 
            str(words2[0]) + " and query " + str(words2[3]) + " inserted char is " + 
            str(words2[2]) + " {" + words2[1] + " | " + words2[2] + "}" +  " " + words2[14] + "\n\n" )
            singleCharacterInsertion+=1    
        elif(words2[1].isalpha() and words2[2].isalpha()==False):
            fout.write("Single character deletion at reference " + 
            str(words2[0]) + " and query " + str(words2[3]) + " deleted char is " + 
            str(words2[1]) + " {" + words2[1] + " | " + words2[2] + "}" +  " " + words2[14] + "\n\n" )
            singleCharacterDeletion+=1
    elif words2[0] != words3[0] and words2[3] == words3[3]:
        #next deleted string start
        deletedString += words2[1]
        fout.write("deletion begin at reference " + words2[0] + " and query " + words2[3] +  " " + words2[14] + "\n")
    elif words2[0] == words3[0] and words2[3] != words3[3]:
        #next inserted string start
        insertedString += words2[2]
        fout.write("insertion begin at reference " + words2[0] + " and query " + words2[3] +  " " + words2[14] + "\n")    
    return
    
def moderatorFunction (words1, words2, words3):
    if words1[0] != words2[0] and words1[3] == words2[3]:
        checkDelete(words2, words3)            
    elif words1[0] == words2[0] and words1[3] != words2[3]:
        checkInsert(words2, words3)
    elif words1[0] != words2[0] and words1[3] != words2[3]:
        checkSNPs(words2, words3)   
    return 
    
    
     
with open ('ref_qry.snps') as fin, open('out.txt', 'w') as fout:
    for _ in xrange(5):
        next(fin)
        
    line1 = fin.next()
    line2 = fin.next()
    
    w1 = line1.split()
    w2 = line2.split()
    
    if w1[0] != w2[0] and w1[3] != w2[3]:
        checkSNPs(w1,w2)      
    elif w1[0] != w2[0] and w1[3] == w2[3]:
        checkDelete(w1,w2)
    elif w1[0] == w2[0] and w1[3] != w2[3]:
        checkInsert(w1, w2)
        
    for line in fin:
        line3 = line
        words1 = line1.split()
        words2 = line2.split()
        words3 = line3.split()
    
        moderatorFunction(words1, words2, words3)    

        #fout.write(str(words1[0]) + "\t" + str(words1[1]) + "\t" + str(words1[2]) + "\t" + str(words1[3]) + "\t" + str(words1[14]) + "\t")
        #fout.write(str(words2[0]) + "\t" + str(words2[1]) + "\t" + str(words2[2]) + "\t" + str(words2[3]) + "\t" + str(words2[14]) + "\n")
        line1 = line2
        line2 = line3
    
    w2 = line1.split()
    w1 = line2.split()
    if w1[0] != w2[0] and w1[3] == w2[3]:
        checkDelete(w1, w2)            
    elif w1[0] == w2[0] and w1[3] != w2[3]:
        checkInsert(w1, w2)
    elif w1[0] != w2[0] and w1[3] != w2[3]:
        checkSNPs(w1, w2) 
    fout.write("Deleted : " + str(deleteCount) + "\n" + "Inserted: " + str(insertCount) + "\n" + "singleCharacterDeletion: " + str(singleCharacterDeletion) + "\n" + "singleCharacterInsertion: " + str(singleCharacterInsertion) + "\n" + "SNP: " + str(SNPCount) + "\n")

