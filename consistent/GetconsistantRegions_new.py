#!/homes/dthybert/software/Python-2.7.5/python
import pysam
import scipy.stats
import sys
import argparse



class GenomeSegment:

	def __init__(self,size, chr,start):
		self.chr=chr
		self.start=start
		self.windowsSmoothing=5
		self.lstPos=[0]*size
		self.lstNbRead=[0]*size
		self.lstFraction=[0.0]*size
		self.lstNormFraction=[0.0]*size
		self.lstOtherInformation=[[]]*size
		self.smoothedNbReads=[0.0]*size
		self.smoothedFraction=[0.0]*size
	
	def addPosition(self, position, index):
		tabLine=position.split("\t")
		self.lstPos[index]=int(tabLine[1])
		self.lstNbRead[index]=int(tabLine[2])
		self.lstFraction[index]=float(tabLine[3])
		self.lstNormFraction[index]=float(tabLine[4])
		self.lstOtherInformation[index]=tabLine[5:]
		
	def _average(self,lst):
		sum=0
		for v in lst:
			sum=v+sum
		return float(sum)/float(len(lst))
	
	def smooth(self,size):
		i=0
		size=len(self.lstPos)
		while i < size:
			smoothNBRead=0.0
			smmothFraction=0.0
			if i < 5:
				smoothNBRead=self._average(self.lstNbRead[:i+self.windowsSmoothing])
				smmothFraction=self._average(self.lstFraction[:i+self.windowsSmoothing])
			elif i > size-5:
			 	smoothNBRead=self._average(self.lstNbRead[i-self.windowsSmoothing:])
			 	smmothFraction=self._average(self.lstFraction[i-self.windowsSmoothing:])
			else:
				smoothNBRead=self._average(self.lstNbRead[i-self.windowsSmoothing:i+self.windowsSmoothing])
				smmothFraction=self._average(self.lstFraction[i-self.windowsSmoothing:i+self.windowsSmoothing])
			
			self.smoothedNbReads[i]=smoothNBRead
			self.smoothedFraction[i]=smmothFraction
			i=i+1
	
	
	def IdentifyGoodRegion(self, nbReadMini, FreqThreshold):
		lstRegions=[]
		start=self.start
		end=self.start
		i=0
		while i < len(self.smoothedNbReads):
			if self.smoothedNbReads[i] < nbReadMini and self.smoothedFraction[i] <FreqThreshold:
				if start!=end:
					lstRegions.append([self.chr, start,end])
				start=self.start+i
				end=self.start+i
			else:
				end=end+1
			i=i+1
		
		return lstRegions


		
def Z_score(val, mean,std):
	return (float(val)-float(mean))/float(std)

def loadStatistics(strconfigFile):
	statByFile={}
	objFile=open(strconfigFile)
	for line in objFile:
		if line[0]=="#":
			continue
		tabLine=line.split()
		file=tabLine[0]
		mean=float(tabLine[1])
		std=float(tabLine[2])
		statByFile[file]=[mean,std]
	return statByFile


def getString(dico, file,pos):
	#print pos
	lsttag=dico[file][pos]
	stringTag="-"
	for tag in lsttag:
		if stringTag=="-":
			stringTag=str(tag)
		else:
			stringTag=stringTag+","+str(tag)
	return 	stringTag



def getLineToPrint(dico,index,pos,chr):
	
	nbTotalOK=0
	nbTotal=0
	fractionOk=0.0
	correctoedFractionOk=0.0
	lstTotal=[]
	lstFraction=[]
	i=0
	for sample in dico.keys():
		lstTag=dico[sample][index]
		nbTagOK=0	
		nbTagMQbad=0
		nbOrphanTag=0
		for tag in lstTag:
			if tag==1:
				nbTagOK=nbTagOK+1				
			if tag==4:
				nbTagMQbad=nbTagMQbad+1
			if tag==0:
				nbOrphanTag=nbOrphanTag+1

		lstTotal.append(nbTagOK)		
		sizeSample=len(lstTag)-nbTagMQbad-nbOrphanTag
		#print sizeSample,len(lstTag)
		if sizeSample==0:
			fraction=0
		else:
			fraction=float(nbTagOK)/float(sizeSample)
		
		lstFraction.append(fraction)
		nbTotal=nbTotal+sizeSample
		nbTotalOK=nbTotalOK+nbTagOK
	
	for fr in lstFraction:
		correctoedFractionOk=correctoedFractionOk+fr
	correctoedFractionOk=correctoedFractionOk/float(len(lstFraction))
	fractionOk=0.0
	if nbTotal!=0:
		fractionOk=float(nbTotalOK)/float(nbTotal)
		
	string=chr+"\t"+str(pos)+"\t"+str(nbTotalOK)+"\t"+str(fractionOk)+"\t"+str(correctoedFractionOk)
	i=0
	for sample in dico.keys():
		string=string+"\t"+str(lstTotal[i])+"\t"+str(lstFraction[i])
		i=i+1	
	i=0
	for sample in dico.keys():
		string=string+"\t"+getString(dico,sample,index)
		i=i+1

	return string



def calculateFrequency(objreadcount, chr,start,end,outFile):
	objFile=open(outFile,"a")
	length=end-start+1
	obgGenomeSegment=GenomeSegment(length,chr,start)
	i=0
	while i < length:
		#print i, length
		pos=start+i
		string=getLineToPrint(objreadcount,i, pos, chr)
		obgGenomeSegment.addPosition(string, i)
		objFile.write(string+"\n")
		#print string
		i=i+1
	objFile.close()
	return obgGenomeSegment
	




##################################################################
#
#
#
#
#
#################################################################	
def countReadsMate(lstFile,dicoStats,chr,start,end,threshold_pval,MQ):
	dicoPos={}
	for file in lstFile:
		print file
		samfile = pysam.AlignmentFile(file, "rb")
		lstPos=[[]]*(end-start+1)
		for pileupcolumn in samfile.pileup(chr,start,end):
				position=pileupcolumn.reference_pos
				lst=[]
				if position < start:
					continue
				if position > end:
					break
				posTab=position-start
				for pReads in pileupcolumn.pileups:
					if  pReads.alignment.mapping_quality < MQ:
						lst.append(4)
					if pReads.alignment.mate_is_unmapped:
						lst.append(0)
						#lstPos[posTab].append(0)
					elif  samfile.getrname(pReads.alignment.next_reference_id) != chr:
						lst.append(3)
					else:
						rend=pReads.alignment.reference_end
						startMate=pReads.alignment.next_reference_start
						delta=abs(startMate-rend)
						mean=dicoStats[file][0]
						std=dicoStats[file][1]
						z=Z_score(delta,mean,std)
						p_value = scipy.stats.norm.sf([abs(z)])[0]
						#print pReads.alignment.next_reference_id
						#print mean, std, delta, p_value
						if p_value <  threshold_pval:
							lst.append(2)
						else:
							lst.append(1)		
				lstPos[posTab]=lst
		dicoPos[file]=lstPos	
	return dicoPos


def saveLstRegion(lstRegion, fileOut):
	objFile=open(fileOut,"a")
	for region in lstRegion:
		string=region[0]+"\t"+str(region[1])+"\t"+str(region[2])+"\n"
		objFile.write(string)
	objFile.close()
			

def main(param):
	
	dicoStats=loadStatistics(param.strConfigFile)
	
	##InitFileTo analyse
	outfile=param.outFile
	outReadCount=outfile+".rdc"
	outGoodRegion=outfile+".bed"
	objFile=open(outReadCount,"w")
	objFile.close()
	objFile=open(outGoodRegion,"w")
	objFile.close()

	lstBams=param.lstBamFiles.split(",")
	CurrStart=param.start
	CurrEnd=param.start+param.bin-1
	#print end-start
	if param.end-param.start < param.bin:
		CurrEnd=param.end
	
	while CurrEnd <=param.end:
			
		##count reads pair
		print "counting paired reads"
		hashReadCount=countReadsMate(lstBams,dicoStats,param.chr,CurrStart,CurrEnd,param.pvalMate,param.MQthreshold)
	
		## calculate some stat and create an object that represnt genome segment (save the data in file
		print " calculate frequencies"
		objGenomSegment=calculateFrequency(hashReadCount,param.chr,CurrStart,CurrEnd,outReadCount)
	
		## get the regioni
		print "smoothing count"
		objGenomSegment.smooth(param.smoothingWindows)
		print "identify regions"
		lstRegion=objGenomSegment.IdentifyGoodRegion(param.minReads, param.minFreq)
		
		## save the regions
		saveLstRegion(lstRegion,outGoodRegion)
		
		CurrStart=CurrEnd+1
		CurrEnd=CurrStart+param.bin-1
		if CurrEnd > param.end:
			CurrEnd=param.end
		if CurrEnd<=CurrStart:
			break
	

####################################################################################


parser = argparse.ArgumentParser()
parser.add_argument('--bam_files', action='store', dest='lstBamFiles', default ="", help='liste of bam file to analyse format : bam1,bam2,...,bamN',required=True)
parser.add_argument('--config', action='store', dest='strConfigFile', help='configuration file describing the mean and std of the insert per library', required=True)
parser.add_argument('--out', action='store', dest='outFile', help='output file prefix where the data will be stored ',  required=True)
parser.add_argument('--chr', action='store', dest='chr', help='chromosome to analyse',required=True)
parser.add_argument('--start', action='store', dest='start', help='start of the region to analyse',required=True, type=int)
parser.add_argument('--end', action='store', dest='end', help='end of the region to analyse\n',required=True,type=int)
parser.add_argument('--pval_mate', action='store', dest='pvalMate', help='pval threshold that two mates are in a good distance [0.0001]', default=0.0001, type=float)
parser.add_argument('--min_reads', action='store', dest='minReads', help='minimum number of reads that satisfy the pair-ends constraints required to have a "good" region [5]', default=5, type=int)
parser.add_argument('--min_freq', action='store', dest='minFreq', help='frequency threshold of reads satisfying the pair-end constraints to have a good regions [0.1]', default=0.1, type=float)
parser.add_argument('--MQ', action='store', dest='MQthreshold', help='reads with a mapping quality < MQ won\'t be considered [25]', default=25, type =int)
parser.add_argument('--smoothing_size', action='store', dest='smoothingWindows', help='size of the windows used to smooth the dataseti [5]', default=5, type=int)
parser.add_argument('--bin', action='store', dest='bin', help='number of position evaluated before storing in file (this is for performances issues) [30000]', default=30000, type=int)
param = parser.parse_args()

main(param)

