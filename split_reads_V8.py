#Algorithm for identification of reads overlaping translocation breakpoint junctions - Split-reads. This tool was specificly design for long insert genome sequencing libraries (liGS), with small reads (<50bp).

#it starts by searching unpaired reads on the regions delimited by the read clusters
#if their is no reads unpaired, the script shuts down

#on the side, the fragment between the clusters, i.e. the delimited regions where the breakpoint is retrived from the NCBI database

#this fragments are bwa indexed
#then, each read is fragmented, in each side, as the example bellow:

#ATGTGTAAACACGTGCACCCACGTG
#			|
#ATGTG				 ACGTG

#This fragments are then mapped against the indexed fragments, and their size
#incresed until their is no more sequence or no more homology.
#Outputs the sam line of the final sequence

#Dependences:
#os
#sys
#Collections
#Biopython, Entrez, SeqIO			

#Run:
#python split_reads_V8.py regionA regionB samfile hg19/hg38 type > results
#########################################################################
import sys
from sys import argv
import os
from Bio import Entrez, SeqIO
import collections

#chromossomes ncbi acession, hg19
chrs={"1":"NC_000001.10", "2":"NC_000002.11", "3":"NC_000003.11", "4":"NC_000004.11", "5":"NC_000005.9", "6":"NC_000006.11", "7":"NC_000007.13", "8":"NC_000008.10", "9": "NC_000009.11", "10":"NC_000010.10", "11":"NC_000011.9", "12": "NC_000012.11","13":"NC_000013.10", "14":"NC_000014.8", "15":"NC_000015.9", "16":"NC_000016.9", "17":"NC_000017.10", "18":"NC_000018.9", "19":"NC_000019.9", "20":"NC_000020.10", "21":"NC_000021.8", "22":"NC_000022.10", "X":"NC_000023.10", "Y":"NC_000024.9"}
#chromossomes ncbi acession, hg38
chrshg38={"1":"NC_000001.11", "2":"NC_000002.12", "3":"NC_000003.12", "4":"NC_000004.12", "5":"NC_000005.10", "6":"NC_000006.12", "7":"NC_000007.14", "8":"NC_000008.11", "9": "NC_000009.12", "10":"NC_000010.11", "11":"NC_000011.10", "12": "NC_000012.12","13":"NC_000013.11", "14":"NC_000014.9", "15":"NC_000015.10", "16":"NC_000016.10", "17":"NC_000017.11", "18":"NC_000018.10", "19":"NC_000019.10", "20":"NC_000020.11", "21":"NC_000021.9", "22":"NC_000022.11", "X":"NC_000023.11", "Y":"NC_000024.10"}

def read_is_unmaped(value):
	"""Read the sam flag read unmapped
	return true if mate is unmapped, return
	false if mate is not unmaped"""
	v=int(value)
	if v & 4:
		"""Its mate is unmaped"""
		return True
	else:
		return False

def read_unmaped(samfile, chrr, start1, end1, dic_unm):
	"""read the cluster zone and check for unpared
	mapped reads. Return the reads names"""
	sam=open(samfile)
	for i in sam:
		line=i.split("\t")
		if line[2]==chrr:#checks the cromossome of the cluster
			if (int(line[3])>=start1 and int(line[3])<=end1): #checks the position fo the read on the cluster
				if read_is_unmaped(line[1])==True:
					l=len(line[9])
					dic_unm["@"+line[0]+"_"+str(l)]=[line[9].strip(), line[10].strip()]#add to the dictionary [readname_readsize]=[read sequence, read quality]
	sam.close()
	return dic_unm


def make_file_and_index(chh, limitA, limitB, outfile, strandd, version):
	"""uses the output limits of the read clusters function, and retrive from
	the ncbi database the portion of chromossome which interest us. 
	It as the possibility of retriving the reverse complement of the sequence.
	Then, writes to a file as:
	>chr_name:start-stop:1/2
	where 1 is forward strand and 2 is reverse strand"""
	if version=="hg19":
		aa=chrs
	if version =="hg38":
		aa=chrshg38
	out=open(outfile, "a")
	Entrez.email="joana.fino@insa.min-saude.pt"
	handle=Entrez.efetch(db="nucleotide", id=aa[chh], rettype="fasta", strand=strandd, seq_start=limitA, seq_stop=limitB)
	record=SeqIO.read(handle, "fasta")
	handle.close()
	out.write(">"+chh+":"+str(limitA)+"-"+str(limitB)+"\n")#":"+str(strandd)+
	out.write(str(record.seq)+"\n")
	out.close()


def split_read(list_fq, chunk, outfile, names):
	"""Split the reads from fastq file, into fragments of
	chunk size, and acording to forward or reverse.
	Write the fragments in a file with the sufix start and end
	acording with the fragment of the read in the begining or the end"""
	t=""
	d=""
	temp=0
	out=open(outfile, "w")
	for key, value in list_fq.items():
		if key in names or len(names)==0:
			out.write(key.split("\t")[0]+".start\n")
			out.write(value[0][:chunk]+"\n")
			out.write("+\n")
			out.write(value[1][:chunk]+"\n")
			
			out.write(key.split("\t")[0]+".end\n")
			out.write(value[0][(chunk+1)*-1:]+"\n")
			out.write("+\n")
			out.write(value[1][(chunk+1)*-1:]+"\n")	
	out.close()

def read_sam_final(samfile):
	"""Read the sam file outputed by the alignment, and checks if their 
	is mapped reads or not. return a dictionary with the hits from the
	alignment"""
	sam=open(samfile)
	dic={}
	to_continue=set()
	for i in sam:
		if i.startswith("@")==False:
			line=i.split("\t")
			if line[1]!="4":
				tm=line[0].split(".")[0]
				sz=int(tm.split("_")[1])
				if sz!=len(line[9]):
					to_continue.add("@"+line[0].split(".")[0])
				if line[0] not in dic:
					dic[line[0]]=[i]
				else:
					dic[line[0]].append(i)
	sam.close()
	return dic, to_continue



def two_dics_2_1(dicfinal, intermediate):
	"""update the dicfinal dictionary with the information 
	from the intermediate dictionary."""
	for key, value in intermediate.items():
		if key not in dicfinal:
			dicfinal[key]=value
		else:
			gg=dicfinal[key]+value
			dicfinal[key]=gg
	return dicfinal

def search_pairs(dickeys):
	newd=[]
	for el in dickeys:
		newd.append(el.split(".")[0])
	newr=[item for item, count in collections.Counter(newd).items() if count > 1]
	oth=[item for item, count in collections.Counter(newd).items() if count == 1]
	return newr, oth


def test_others(list_start, list_end, sz):#used in new_parse_output
	aa=list_start[-1].split("\t")
	l=-2
	while l*-1<=len(list_end):
		bb=list_end[l].split("\t")
		if aa[2]!=bb[2] and len(aa[9])+len(bb[9])>=sz:
			return "end", bb
		l-=1
	aa=list_end[-1].split("\t")
	l=-2
	while l*-1<=len(list_start):
		bb=list_start[l].split("\t")
		if aa[2]!=bb[2] and len(aa[9])+len(bb[9])>=sz:
			return "start", bb
		l-=1
	return "none", []
			

def new_parse_output(dic_final, t):
	readn, oth=search_pairs(dic_final.keys())
	posp={}
	for readname in readn:
		sz=int(readname.split("_")[1])#readsize
		st=dic_final[readname+".start"][-1].split("\t")#start result
		end=dic_final[readname+".end"][-1].split("\t")#end result
		if len(st[9])+len(end[9])>=0.8*sz:
			if (st[2]!=end[2]) or (st[2]==end[2] and t=="del"):
				posp[readname]=[[readname+".start", st[1], len(st[9]), st[2], st[3], st[9]], [readname+".end", end[1], len(end[9]), end[2], end[3], end[9]]]
			else:
				t, l=test_others(dic_final[readname+".start"], dic_final[readname+".end"], sz)#test if other hit gives a diferent chromosome.
				if t=="start":
					posp[readname]=[[readname+".start", l[1], len(l[9]), l[2], l[3], l[9]], [readname+".end", end[1], len(end[9]), end[2], end[3], end[9]]]
				elif t=="end":
					posp[readname]=[[readname+".start", st[1], len(st[9]), st[2], st[3], st[9]], [readname+".end", l[1], len(l[9]), l[2], l[3], l[9]]]
	print oth
	for r in oth:
		if r+".start" in dic_final:
			if len(dic_final[r+".start"][-1].split("\t")[9])>=0.68*sz:
				st=dic_final[r+".start"][-1].split("\t")
				posp[readname]=[[r+".start", st[1], len(st[9]), st[2], st[3], st[9]], [readname+".end", "insuficient alignment for the end chunk."]]
		if r+".end" in dic_final:
			if len(dic_final[r+".end"][-1].split("\t")[9])>=0.68*sz:
				end=dic_final[r+".end"][-1].split("\t")
				posp[readname]=[[r+".start", "insuficient alignment for the start chunk."], [readname+".end", end[1], len(end[9]), end[2], end[3], end[9]]]

	return posp


def read_orientation(flag):
	if int(flag)&16:
		return "reverse"
	else:
		return "forward"
		
def search_pair(samfile, readname):
	f=open(samfile)
	for i in f:
		if i.startswith(readname):
			if int(i.split("\t")[1])&8:
				line=i.split("\t")
				loc=int(line[3])+int(line[5][:-1])
				a=line[2]+":"+line[3]+"-"+str(loc)+"\t"+read_orientation(line[1])
				break
	f.close()
	return a

def calc_pos(l):
	ref=l[3].split(":")
	s=ref[1].split("-")[0]
	start=int(s)+int(l[4])
	end=start+int(l[2])
	return ref[0]+":"+str(start)+"-"+str(end)

def write_output(samfile, posp):
	print("Original read name\tChunk Read Name\tOrientation\tSize\tPosition\tSequence\t Mate Position\t Mate orientation")
	for key, value in posp.items():
		r=key.split("_")[0]
		pair_pos=search_pair(samfile, r)
		if len(value[0])>2:
			startpos=calc_pos(value[0])
			print(r+"\t"+value[0][0]+"\t"+read_orientation(value[0][1])+"\t"+str(value[0][2])+"\t"+startpos+"\t"+value[0][5]+"\t"+pair_pos)
		if len(value[0])==2:
			print(r+"\t"+value[0][0]+"\t"+value[0][1])
		if len(value[1])>2:
			endpos=calc_pos(value[1])
			print(r+"\t"+value[1][0]+"\t"+read_orientation(value[1][1])+"\t"+str(value[1][2])+"\t"+endpos+"\t"+value[1][5]+"\t"+pair_pos)
		if len(value[1])==2:
			print(r+"\t"+value[1][0]+"\t"+value[1][1])
		print("\n")

def organize(unmappeds, ttt, samfile):
	t=5#start size of reads
	os.system("bwa index -a bwtsw junctions.fasta")#index with bwa
	final_pairs={}
	print("all prepared")
	split_read(unmappeds,t, "chunk.fastq", []) #split the read to chunks of 5bp
	os.system("bwa aln -n 0 -o 0 -e 0 -O 10 -M 1000 -t 3 junctions.fasta chunk.fastq > maped.sai") #align by bwa align the parameters are seted for not alowing any mismatches or gaps
	os.system("bwa samse  -n 20 junctions.fasta maped.sai chunk.fastq > map.sam")#make the final part of the alignment
	b, to_continue=read_sam_final("map.sam")#read the samfile, and retrive a dictionary with the aligned reads
	#os.system("rm chunk.fastq")#remove the chunk file, to not mess with the next steps
	while(len(to_continue)!=0):# iterate by the rest of the read
		t+=1
		final_pairs=two_dics_2_1(final_pairs, b)#pass the entries to list
		split_read(unmappeds,t, "chunk.fastq", to_continue)#split with one more bp
		os.system("bwa aln -n 0 -o 0 -e 0 -O 10 -M 1000 -t 3 junctions.fasta chunk.fastq > maped.sai")#align
		os.system("bwa samse -n 20 junctions.fasta maped.sai chunk.fastq > map.sam")
		b, to_continue=read_sam_final("map.sam")#read sam
	#	os.system("rm chunk.fastq")#remove read chunk
	write_output(samfile, new_parse_output(final_pairs, ttt))


########################MAIN RUN############################################################################

def trans(regA,regB, samfile, version, ttt): #TRANS inv
	rrA=regA.split(":")
	rrrA=rrA[1].split("-")
	limitA=[rrA[0],int(rrrA[0]), int(rrrA[1]), int(rrrA[2]),int(rrrA[3])]
	rrB=regB.split(":")
	rrrB=rrB[1].split("-")
	limitB=[rrB[0],int(rrrB[0]), int(rrrB[1]), int(rrrB[2]),int(rrrB[3])]
	print("limits done")
	unmappedsA= read_unmaped(samfile, limitA[0], limitA[1], limitA[2],{})
	unmappedsAA=read_unmaped(samfile, limitA[0], limitA[3], limitA[4],unmappedsA)
	unmappedsB= read_unmaped(samfile, limitB[0], limitB[1], limitB[2],unmappedsAA)
	unmappeds=read_unmaped(samfile, limitB[0], limitB[3], limitB[4],unmappedsB)
	print("unmapeds done")
	if len(unmappeds)==0:#if their is no unmapped reads, kill the program with a message
		print("No unmaped reads! Shuting down...")
		sys.exit()
	make_file_and_index(limitB[0],limitB[2]-25, limitB[3]+25, "junctions.fasta", 1, version)
	make_file_and_index(limitA[0],limitA[2]-25, limitA[3]+25, "junctions.fasta", 1, version)
	print("file index done")
	organize(unmappeds, ttt, samfile)

def orgdel(regA,regB, samfile, version, ttt):#dell dup
	rrA=regA.split(":")
	rrrA=rrA[1].split("-")
	limitA=[rrA[0],int(rrrA[0]), int(rrrA[1])]
	rrB=regB.split(":")
	rrrB=rrB[1].split("-")
	limitB=[rrB[0],int(rrrB[0]), int(rrrB[1])]
	print("limits done")
	unmappedsA= read_unmaped(samfile, limitA[0], limitA[1], limitA[2],{})
	unmappeds= read_unmaped(samfile, limitB[0], limitB[1], limitB[2],unmappedsA)
	print("unmapeds done")
	if len(unmappeds)==0:#if their is no unmapped reads, kill the program with a message
		print("No unmaped reads! Shuting down...")
		sys.exit()
	if ttt=="del":
		make_file_and_index(limitA[0],limitA[2]-25, limitB[1]+25, "junctions.fasta", 1, version)
	else:
		make_file_and_index(limitA[0],limitA[1]-5000, limitA[1], "junctions.fasta", 1, version)
		make_file_and_index(limitB[0],limitB[2], limitB[2]+5000, "junctions.fasta", 1, version)
	print("file index done")
	organize(unmappeds, ttt, samfile)

def orgins(regA,regB, samfile, version, ttt):#ins
	rrA=regA.split(":")
	rrrA=rrA[1].split("-")
	limitA=[rrA[0],int(rrrA[0]), int(rrrA[1]), int(rrrA[2]),int(rrrA[3])]
	rrB=regB.split(":")
	rrrB=rrB[1].split("-")
	limitB=[rrB[0],int(rrrB[0]), int(rrrB[1]), int(rrrB[2]),int(rrrB[3])]
	print("limits done")
	unmappedsA= read_unmaped(samfile, limitA[0], limitA[1], limitA[2],{})
	unmappedsAA=read_unmaped(samfile, limitA[0], limitA[3], limitA[4],unmappedsA)
	unmappedsB= read_unmaped(samfile, limitB[0], limitB[1], limitB[2],unmappedsAA)
	unmappeds=read_unmaped(samfile, limitB[0], limitB[3], limitB[4],unmappedsB)
	print("unmapeds done")
	if len(unmappeds)==0:#if their is no unmapped reads, kill the program with a message
		print("No unmaped reads! Shuting down...")
		sys.exit()
	make_file_and_index(limitA[0],limitA[1]-5000, limitA[1], "junctions.fasta", 1, version)#limite da zona de excisao a 5'
	make_file_and_index(limitA[0],limitA[4], limitA[4]+5000, "junctions.fasta", 1, version)#limite da zona de excisao a 3'
	make_file_and_index(limitB[0],limitB[2]-25, limitB[3]+25, "junctions.fasta", 1, version)#limite da regiao da insercao
	print("file index done")
	organize(unmappeds, ttt, samfile)


#########################################################################################################################  

#python split_reads_V8.py regionA regionB samfile hg19/hg38 type > results
#where regionA and regionB= chr:start_cluster - start breakpoint region - end breapoint region - end cluster
if argv[5]=="trans" or argv[5]=="inv":
	trans(argv[1], argv[2], argv[3], argv[4], argv[5])
elif argv[5]=="del" or argv[5]=="dup":
	orgdel(argv[1], argv[2], argv[3], argv[4], argv[5])
elif argv[5]=="ins":
	orgins(argv[1], argv[2], argv[3], argv[4], argv[5])

