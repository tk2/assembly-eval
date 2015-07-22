#!/homes/dthybert/software/Python-2.7.5/python
import pysam
import scipy.stats
import sys
import argparse







def Z_score(val, mean,std):
        return (float(val)-float(mean))/float(std)

def isGoodRead(read,threshold_pval,dicoStats,bamFile):
	rend=read.alignment.reference_end
	startMate=read.alignment.next_reference_start
	delta=abs(startMate-rend)
	mean=dicoStats[bamFile][0]
	std=dicoStats[bamFile][1]
	z=Z_score(delta,mean,std)
	p_value = scipy.stats.norm.sf([abs(z)])[0]                              
#	print delta,mean, std, p_value 
	if p_value <  threshold_pval:
		return 0
	return 1

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
 

def getPositionInTabCoordinate(start,CurrPosition):
	return CurrPosition-start
	

def addSupport(tabSupport,start,end):
	i=start
	if start < 0 or end < 0:
		return tabSupport
	#print i,end
	while i <= end:
		if i >=len(tabSupport):
			break
		tabSupport[i]=tabSupport[i]+1
		i=i+1
	return tabSupport

def MergeList(lstOfList):
	length=len(lstOfList[0])
	#print lstOfList[0],lstOfList[1]
	i=0
	lstNew=[0]*length
	while i < length:
		j=0
		while j < len(lstOfList):
			lstNew[i]=lstNew[i]+lstOfList[j][i]
			j=j+1
		i=i+1
	return lstNew



def getConfidentRegions(tab, start, threshold):
	i=0
	CurStart=0
	CurEnd=0
	lstPosition=[]
	while i < len(tab):
		if tab[i] < threshold:
			if CurStart!=CurEnd:
				regStart=CurStart+start
				regEnd=i+start
				lstPosition.append([regStart,regEnd])
			CurStart=i
			CurEnd=i
		else:
			CurEnd=i
		i=i+1
	if CurStart!=CurEnd:
		 regStart=CurStart+start
                 regEnd=CurEnd+start
                 lstPosition.append([regStart,regEnd])
	return lstPosition



def defineRegionFile(bamFile,dicoStats,chr,start,end,threshold_pval,readLength, bin,buffer,f):
	
	samfile = pysam.AlignmentFile(bamFile, "rb")
	size=end-start+buffer
    	tabSupport=[0]*(size)
    	CurStart=start
    	CurEnd=start+bin
	if CurEnd > end:
		CurEnd=end
	while CurStart < end: #Parse the genomic region to analyse
		i=0
   		for pileupcolumn in samfile.pileup(chr,CurStart,CurEnd):#the analysis is divided in bin for memorry purpose
    			position=pileupcolumn.reference_pos
			lst=[]
        		if position < start:
        			continue
        		if position > end:
           	 		break
        		posTab=position-start
			if i % f==0: 
    				for pReads in pileupcolumn.pileups:#analyse each reads of a position
    					if pReads.alignment.mate_is_unmapped:
               					continue
            				elif  samfile.getrname(pReads.alignment.next_reference_id) != chr:
                				continue
            				else:
                				if isGoodRead(pReads,threshold_pval,dicoStats,bamFile):#in the case the read satisfy insert constraint we can take it into account
							rend=pReads.alignment.reference_end
                    					startMate=pReads.alignment.next_reference_start
                    					delta=startMate-rend
                    					if delta > 0:# take into account only whenm the mate pair is forward , this is not to count twice the relationship
                    						startTab=getPositionInTabCoordinate(start,pReads.alignment.reference_start)
                        					endMate=startMate+readLength
                       						endTab=getPositionInTabCoordinate(start,endMate)
                     						tabSupport=addSupport(tabSupport,startTab,endTab)
			i=i+1
    		CurStart=CurEnd+1
    		CurEnd=CurStart+bin
		if CurEnd > end:
			CurEnd=end
	#print tabSupport
	return tabSupport



def saveRegions(outfile,bedList,chr):
	objFile=open(outfile,"w")
	for list in bedList:		
		string=chr+"\t"+str(list[0])+"\t"+str(list[1])
		objFile.write(string+"\n") 
	objFile.close()


def saveScore(outfile,ListPosition,chr, start):
	objFile=open(outfile,"w")
	i=0
	while i < len(ListPosition):
		pos=start+i
		val=ListPosition[i]
		p=str(pos)
		v=str(val)
		string=chr+"\t"+p+"\t"+v
		objFile.write(string+"\n")
		i=i+1
	objFile.close()
		
	
def main(param):
        dicoStats=loadStatistics(param.strConfigFile)
        lstBams=param.lstBamFiles.split(",")
        lstLst=[]
        for bam in lstBams:
		print "start analysing "+bam+ " file"
		###Analyse a bam file
		lst=defineRegionFile(bam,dicoStats,param.chr,param.start,param.end,param.pvalMate,param.readLength,param.bin, param.buffer,param.frequency)
		lstLst.append(lst)
		print bam +" file treated"
	###merge all data from the different bamfile
	FinalList=MergeList(lstLst)
	
	## save the results
	outScore=param.outFile+".score"
	saveScore(outScore, FinalList, param.chr, param.start)
	outregions=param.outFile+".bed"
	bedList=getConfidentRegions(FinalList,param.start,param.threshold)
	saveRegions(outregions,bedList,param.chr)
		

parser = argparse.ArgumentParser()
parser.add_argument('--bam_files', action='store', dest='lstBamFiles', default ="", help='liste of bam file to analyse format : bam1,bam2,...,bamN',required=True)
parser.add_argument('--config', action='store', dest='strConfigFile', help='configuration file describing the mean and std of the insert per library', required=True)
parser.add_argument('--out', action='store', dest='outFile', help='output file prefix where the data will be stored ',  required=True)
parser.add_argument('--chr', action='store', dest='chr', help='chromosome to analyse',required=True)
parser.add_argument('--start', action='store', dest='start', help='start of the region to analyse',required=True, type=int)
parser.add_argument('--end', action='store', dest='end', help='end of the region to analyse\n',required=True,type=int)
parser.add_argument('--pval_mate', action='store', dest='pvalMate', help='pval threshold that two mates are in a good distance [0.01]', default=0.01, type=float)
parser.add_argument('--threshold', action='store', dest='threshold', help='coverage threshold to define a "good" region [1]', default=1, type=int)
#parser.add_argument('--min_freq', action='store', dest='minFreq', help='frequency threshold of reads satisfying the pair-end constraints to have a good regions [0.1]', default=0.1, type=float)
parser.add_argument('--bin', action='store', dest='bin', help='number of position evaluated before storing in file (this is for performances issues) [30000]', default=30000, type=int)
parser.add_argument('--read_length', action='store', dest='readLength', help='the length of the mapped read [100]', default=100, type=int)
parser.add_argument('--buffer', action='store', dest='buffer', help='the buffer size define the what is the distance after the last postion we can take into account the a the mate of the read treated.Because of the good regions can go beyond the end of the end position. Need to be at least the size of the insert [20000]', default=20000, type=int)
parser.add_argument('--f', action='store', dest='frequency', help='positon will be evaluated at every f nt', default=100, type=int)
param = parser.parse_args()

main(param)

