#Script can generate molecular formulas from FASTA file for proteins

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

listofaminoacids = []
A = {'C':3, 'H':7, 'N':1, 'O':2, 'S':0}
R = {'C':6, 'H':14,'N':4, 'O':2, 'S':0}
N = {'C':4, 'H':8, 'N':2, 'O':3, 'S':0}
D = {'C':4, 'H':7, 'N':1, 'O':4, 'S':0}
C = {'C':3, 'H':7, 'N':1, 'O':2, 'S':1}
Q = {'C':5, 'H':10,'N':2, 'O':3, 'S':0}
E = {'C':5, 'H':9, 'N':1, 'O':4, 'S':0}
G = {'C':2, 'H':5, 'N':1, 'O':2, 'S':0}
H = {'C':6, 'H':9, 'N':3, 'O':2, 'S':0}
I = {'C':6, 'H':13,'N':1, 'O':2, 'S':0}
L = {'C':6, 'H':13,'N':1, 'O':2, 'S':0}
K = {'C':6, 'H':14,'N':2, 'O':2, 'S':0}
M = {'C':5, 'H':11,'N':1, 'O':2, 'S':1}
F = {'C':9, 'H':11,'N':1, 'O':2, 'S':0}
P = {'C':5, 'H':9, 'N':1, 'O':2, 'S':0}
S = {'C':3, 'H':7, 'N':1, 'O':3, 'S':0}
T = {'C':4, 'H':9, 'N':1, 'O':3, 'S':0}
W = {'C':11,'H':12,'N':2, 'O':2, 'S':0}
Y = {'C':9, 'H':11,'N':1, 'O':3, 'S':0}
V = {'C':5, 'H':11,'N':1, 'O':2, 'S':0}



dictOfAmino = {'A':A,'R':R,'N':N,'D':D,'C':C,'Q':Q, 'E':E, 'G':G,'H':H,'I':I,'L':L,'K':K,'M':M,'F':F,'P':P,'S':S,'T':T,'W':W,'Y':Y,'V':V}

print "Note output file is appended if same file is selected twice molecular formulas \n for both runs will be present in output file"
fileName = raw_input("Protein FASTA file to generate molecular formulas for: ")
outFileName = raw_input("Output file name (include .txt): ")

fasta_file = open(fileName, "rU")
for record in SeqIO.parse(fasta_file, "fasta"):
	myseq = str(record.seq)
	analysis = ProteinAnalysis(myseq)
	listofaminoacids.append(analysis.count_amino_acids())


	
for i in listofaminoacids:
        carbonTotal = 0
        hydrogenTotal = 0
        oxygenTotal = 0
        nitrogenTotal = 0
        sulfurTotal = 0
        peptideBonds = 0
        
        for value in i:
                for amino in dictOfAmino:
                        
                        if value == amino:
                                peptideBonds = peptideBonds + i[value]
                                thisAmino = {}
                                thisAmino = dictOfAmino[amino]
                                carbonTotal = carbonTotal + (i[value]*thisAmino['C'])
                                hydrogenTotal = hydrogenTotal + (i[value]*thisAmino['H'])
                                oxygenTotal = oxygenTotal + (i[value]*thisAmino['O'])
                                nitrogenTotal = nitrogenTotal + (i[value]*thisAmino['N'])
                                sulfurTotal = sulfurTotal + (i[value]*thisAmino['S'])
                                                             

        #Correcting totals for peptide bond loss of water
        peptideBonds = peptideBonds - 1
        hydrogenTotal = hydrogenTotal -(peptideBonds*2)
        oxygenTotal = oxygenTotal - (peptideBonds*1)
        outString = "C" + str(carbonTotal) + "H" + str(hydrogenTotal) + "N" + str(nitrogenTotal) + "O" + str(oxygenTotal) + "S" + str(sulfurTotal)
        print outString
        outputFile = open(outFileName, "a")
        outputFile.write(outString)


		

fasta_file.close()
outputFile.close()
