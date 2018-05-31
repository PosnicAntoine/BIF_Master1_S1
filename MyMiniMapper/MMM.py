import os, sys, fasta, time, tools_karkkainen_sanders as tks

"""""""""
Functions to display elapsed time and reset counter
"""""""""
last_time = 0

def startChrono():
	global last_time
	last_time = time.time()

def endChrono():
	return (time.time() - last_time) #Retour en secondes

"""""""""""""""
biopalind
This fuction returns the reverse complement of a strand
(biological Palyndrome)
"""""""""""""""

def biologicalPalyndrome(s):
	#dictionnary for complementary bases
	dictinv = {"a":"t","c":"g","g":"c","t":"a"}
	output=""
	for i in range(len(s)):
		#We take all the bases in the reverse order and take the complementary bases
		output+=dictinv[s[len(s)-i-1]]
	return output


"""""""""""""""
usefull functions
"""""""""""""""
"""""""""
returns the Burros-Wheeler-Transform of a string, with an added dollar to find the end.
"""""""""
def getBWT(s, refSA):
	bwt= ""
	for i in refSA:
		if i == 0:
			#If i=0, it means we are at the end ()
			bwt = bwt + "$"
		else:
			#else, we just take the character at pos i-1
			bwt = bwt + s[i-1]
	return bwt

#dictionnary to find corresponding places in F
posdict = {"$":0,"a":1,"c":2,"g":3,"t":4}
"""""""""
returns f, which is the list of places of first appearance of letters (when in alphabetical order). (list of 5, based on posdict)
"""""""""
def getF(bwt):
	#First, we need to count the number of occurence of each letter
	count = [0]*5
	for c in bwt:
		count[posdict[c]] = count[posdict[c]]+1
	#Then, letters in an alphabetical order, the place is the place of the last previous letter.
	f = [0]*5
	curr=0
	for i in range(len(count)):
		f[i] = curr
		curr = curr + count[i]
	return f

"""""""""
returns the list of ranks of all letters.
No subsampling is used here, so it can get pretty heavy for long char sequences
"""""""""
def getRank(bwt):
	ranks = []
	count = [0]*5
	#Same letters are in their order of appearance for the final string
	for c in bwt:
		#So we can remember what their last position was and add it to the result.
		ranks.append(count[posdict[c]])
		#We add one to the letter we just met to have the correct rank for the next one
		count[posdict[c]]+=1
	return ranks


#Gather arguments from the user
if(len(sys.argv)<4):
	#If in the incorrect form, return an error message
	print("Arguments must be of the form : referencefile, readsfile, k, dmax.")
	exit(0)
#referencefile and readsfile must be file names, k and dmax integers.
referencefile, readsfile, kmerLength, dmax = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

#Initialization of the reference file
#readFasta only take the sequence of bases, and the $ is for the BWT, to mark the end of the string.
reference = (fasta.readFasta(referencefile)).lower()+"$"

#Initialization of reads and readsInv, its reverse complementary
reads, readsBioPalind = [],[]
for line in open(readsfile, "r"):
	if line[0] != ">":#lines with > do not contain sequences, but merely comments about the sequences.

		reads.append(line[:-1].lower())#-1 to remove \n. To lower case for practical reasons when calling posdict.
		readsBioPalind.append(biologicalPalyndrome(line[:-1].lower()))#We also stock the biological palyndromes

#We create SA, BWT, Rank and F from reference

print("generating SA")
startChrono()
refSA = tks.simple_kark_sort(reference)
print(" done in "+str(endChrono())+" s")
print("generating BWT")
startChrono()
refBWT = getBWT(reference, refSA)
print(" done in "+str(endChrono())+" s")
print("generating ranks")
startChrono()
refRank = getRank(refBWT)
print(" done in "+str(endChrono())+" s")
print("generating F")
startChrono()
refF = getF(refBWT)
print(" done in "+str(endChrono())+" s")

"""""""""""""""
usefull functions
"""""""""""""""

"""""""""
Recursive search in BWT.
Deprecated, unoptimized.

def searchBWTrec(stri, i):
	if len(stri)==0:
		return refSA[i]
	else :
		currChar = stri[len(stri)-1]
		if refBWT[i] == currChar:
			#We continue searching for the word without its first letter
			return searchBWTrec(stri[:-1], refF[posdict[currChar]] + refRank[i])
		else :
			return



"""""""""

"""""""""
def searchBWT(stri):
	if len(stri)==0:
		return []
	else :
		res = []
		initChar = stri[len(stri)-1]#On prend le dernier caractère
		begin = refF[posdict[initChar]]
		if posdict[initChar] == len(posdict)-1:
			end = len(refBWT)-1
		else:
			end = refF[posdict[initChar]+1]-1
		for i in range(begin, end+1):
			currFound = searchBWTrec(stri[:-1], i)
			if currFound != None:
				res.append(currFound)
		return res

"""""""""

"""""""""
Research function for the BWT
"""""""""
def searchBWT(stri):
	i = len(stri)-1 #On part de la dernière lettre
	currChar = stri[i]#On prend le ième caractère
	begin = refF[posdict[currChar]]#On prend la première occurence de ce caractère
	end = begin #Initialization
	if posdict[currChar] == len(posdict)-1:#Dans le cas où le caractère est un t, le dernier emplacement valide est le dernier de la BWT entière
		end = len(refBWT)-1
	else:
		end = refF[posdict[currChar]+1]-1#On prend l'emplacement juste avant la première occurence du caractère suivant
	i-=1
	while(i >= 0):
		currChar = stri[i]#On prend le ième caractère du read
		#On cherche le premier et dernier caractère de bwt correspondant, pour avoir leurs indices
		while(refBWT[begin]!=currChar):#Tant que dans BWT on ne trouve pas le charactère, on va chercher plus loin
			begin+=1
			if begin>end:
		#		print("AAAAA")
				return [] #Si le début dépasse la fin, c'est fichu, on ne trouvera pas le mot
		while(refBWT[end]!=currChar):#Même chose avec le dernier charactère, on recule jusqu'à trouver sa dernière occurence
			end-=1
			if begin>end:
		#		print("BBBBB")
				return [] #Si le début dépasse la fin, c'est fichu, on ne trouvera pas le mot
		#refF de posdict de currChar donne la première occurence de ce char. refRank donne le rang au sein de ce type de
		begin = refF[posdict[currChar]] + refRank[begin] #refRank récupère le rang du char, on va le trouver dans sa en l'ajoutant au rang de la première occurence.
		end = refF[posdict[currChar]] + refRank[end]
		i-=1
	res=[]
	#Pour le cas où on en trouve plusieurs occurences. En général, il y a juste begin
	for j in range(begin, end+1):
		res.append(refSA[j])
	return res


"""""""""""""""
 delete the seeds that are already in res from listOfFound.
For example : the seed [2,6] is equivalent to [2+n, 6+n], with n<kmersize
"""""""""""""""
def interestingSeeds(kmerI, listOfFound, res):
	interSeeds=[]
	if res==[]:
		return listOfFound
	for currFound in listOfFound:
		#By default, we consider each seed interesting
		interestingFound=True
		for currRes in res:
			currResIndiceKmer=currRes[0]
			for currResIndiceRef in currRes[1]:
				# Except if it's redundant with another already existing seed
				if currFound-kmerI == currResIndiceRef-currResIndiceKmer:
					interestingFound=False
					break
			if not interestingFound:break
		#Et on n'ajoute la seed que si interestinFound est à True
		if interestingFound:
			interSeeds.append(currFound)
	return interSeeds


"""""""""
Looking for every kmer in list readKmers
"""""""""
def searchKmers(readKmers):
	res = []
	for kmer in readKmers:
		currFound = searchBWT(kmer[0])
		currFound = interestingSeeds(kmer[1], currFound, res)
		if currFound != []:
			res.append([kmer[1],currFound])

	return res

"""""""""
Get all possible kmers of size kmerLength
"""""""""
def getKmers(read):
	kmers=[]
	for i in range(len(read)-kmerLength+1):
		kmers.append([read[i:i+kmerLength], i])
	return kmers


"""""""""
Search all possible seeds in ref from multiple reads.
"""""""""
def searchMultipleSeeds(reads):
	kmers=[]
	print("  creating kmers from reads")
	for read in reads :
		kmers.append(getKmers(read))
	res=[]
	print("  searching kmers")
	for i in range(len(reads)):
		res.append(searchKmers
(kmers[i]))
	return res

"""""""""""
From a seed (referencing begining in the read and the reference), extend to the left and the right while score < dmax
"""""""""""
def extendSeed(seed, read):
	score = 0
	if seed[0] > seed[1] or seed[0] + len(read) > len(reference) -seed[1]:
		return

	if seed[0] > 0:
		#extendLeft()
		for i in range(seed[0]):
			#If char is different, increase the number of found mismatches.
			score += 1 if reference[seed[1]-i] != read[seed[0]-1] else 0
			#Return None if score too high.
			if score > dmax:
				return

	if seed[0]+kmerLength < len(reference) and score <= dmax:
		#extendRight()
		for i in range(len(read) - (seed[0]+kmerLength)):
			score += 1 if reference[seed[1]+kmerLength+i] != read[seed[0]+kmerLength+i] else 0
			#Return None if score too high.
			if score > dmax:
				return

	#Return None if score too high.
	if score > dmax:
		return
	return [seed[1] - seed[0], seed[0], score]

"""""""""
Extend multiple seeds
"""""""""
def extendSeeds(seeds, read):
	res = []
	for seed in seeds:
		seedIndiceRead = seed[0]
		for seedIndiceRef in seed[1]:
			currFound = extendSeed([seedIndiceRead, seedIndiceRef], read)
			if currFound != None and not currFound in res:
				res.append([currFound[0], currFound[2]])#Add place found and nb of mismatches
	return res

"""""""""
Extend multiple seeds from multiple reads
"""""""""
def extendMultipleSeeds(multipleSeeds, reads):
	res = []
	for i in range(len(reads)):
		res.append(extendSeeds(multipleSeeds[i], reads[i]))
	return res
"""""""""
Print 2 sequences from read and reference with | when match and : when mismatch
"""""""""
def printmatch(read, pos):
	res ="  "+read.upper()+"\n"
	res+="  "
	for i in range(len(read)):
		res+=("|" if read[i]==reference[i+pos] else ":")
	res+="\n  "+reference[pos:pos+len(read)].upper()+"\n"
	return res

"""""""""
Pretty-print the results from a single read
"""""""""
def prettyPrintSingle(extended, extendedRev, i):
	res = ">read"+str(i+1)+":\n"
	align_counter = 0
	#Dans le cas où on ne trouve rien, on ne retourne rien
	if extended==[] and extendedRev==[]:
		return ""

	for match in extended:
		res+="  >>alignment "+str(align_counter)+"\n"
		res+="  #pos="+str(match[0])+"\n"
		res+="  #strand=+1\n"
		res+="  #d="+str(match[1])+"\n"
		res+=printmatch(reads[i], match[0])
		align_counter+=1

	for match in extendedRev:
		res+="  >>alignment "+str(align_counter)+"\n"
		res+="  #pos="+str(match[0])+"\n"
		res+="  #strand=-1\n"
		res+="  #d="+str(match[1])+"\n"
		res+=printmatch(readsBioPalind[i], match[0])
		align_counter+=1
	return res

"""""""""
Pretty-print the results from a multiple reads
"""""""""
def prettyPrint(extendedSeeds, extendedRevs):
	res = ""
	for i in range(len(extendedSeeds)):
		res += prettyPrintSingle(extendedSeeds[i], extendedRevs[i],i)
	return res



"""""""""""""""
 MAIN FUNCTION
"""""""""""""""

#Get all seeds for all reads, including reverse complementary

begin_all_search = time.time()#To know total time to find everything

print("searching Seeds")
startChrono()
multipleSeeds = searchMultipleSeeds(reads)
print(" done in "+str(endChrono())+" s")
print("searching Seeds Reverse Complement")
startChrono()
multipleSeedsBioPalind = searchMultipleSeeds(readsBioPalind)
print(" done in "+str(endChrono())+" s")
#Extend all seeds
print("extending Seeds")
startChrono()
extendedSeeds = extendMultipleSeeds(multipleSeeds, reads)
print(" done in "+str(endChrono())+" s")
print("extending Seeds Reverse Complement")
startChrono()
extendedSeedsBioPalind = extendMultipleSeeds(multipleSeedsBioPalind, readsBioPalind)
print(" done in "+str(endChrono())+" s")

print("Total time : "+str(time.time()-begin_all_search))

#Count matches
nb_matches=0
for seeds in extendedSeeds:
	for seed in seeds:
		if seed != []:
			nb_matches+=1

for seeds in extendedSeedsBioPalind:
	for seed in seeds:
		if seed != []:
			nb_matches+=1

print("Nb of matches found : "+str(nb_matches))

#Generate text_output from successful extended seeds
text_output = prettyPrint(extendedSeeds, extendedSeedsBioPalind)
#Write in a file and display
text_file = open("output.txt", "w")
text_file.write(text_output)
text_file.close()
#print(text_output)
