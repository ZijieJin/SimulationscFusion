import sys
import random
import math
import numpy


def CountUpper(str):
    sum = 0
    for i in range(len(str)):
        if str[i].isupper():
            sum += 1
    return sum


def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


fastafile = open(sys.argv[1])
outdir = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
totalread = int(sys.argv[5])

chimericdict = {}
lastchimeric = ''
qualitystring = 'AFJJJJ<AA<FFJ7---AFF)7FJFJFJJJ<<JFJJJA-AFJJJJJJJJ7<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAFFF---A-<FJJJJJJJJJJJJJJJJF<<<'
currentcount = 0
weightsum = 0
SNPLevel = 0.01
CNVLevel = 0.005
for line in fastafile.readlines():
    if line.find('>') == 0:
        currentcount += 1
        #weight = math.exp(-2 * int(line[1:-1].split('_')[-1])) # M1 weight
        weight = float(line[1:-1].split('_')[-1]) * float(line[1:-1].split('_')[-1]) / (1 + math.exp(-(float(line[1:-1].split('_')[-1]) - 65)*2))
        chimericdict[currentcount] = [line[1:-1], '', weight]
        weightsum += weight
    else:
        chimericdict[currentcount][1] += line.rstrip()
fastafile.close()

factor = {}
expr = {}
for k in range(start, end+1):
    outfile1 = open(outdir + '/' + str(k) + '_1.fastq', 'a')
    outfile2 = open(outdir + '/' + str(k) + '_2.fastq', 'a')
    for i in chimericdict:
        totallength = len(chimericdict[i][1])
        if i not in factor:
            factor[i] = numpy.random.binomial(1, 0.8)
            expr[i] = int(numpy.random.binomial(1, 0.00025) * 30 * numpy.exp(numpy.random.normal(0, 0.2)) + 1 + numpy.exp(numpy.random.normal(0, 0.2)))
        numreads = int(numpy.random.binomial(totalread, chimericdict[i][2] / weightsum, 1)) * (1 + numpy.random.negative_binomial(1,0.6)) * (1 + factor[i] * numpy.random.binomial(2,0.625)) * expr[i]
        if numreads <= 0:
            continue
        for j in range(numreads):
            whilecount = 0
            while True:
                whilecount += 1
                insertsize = random.randint(min(600, totallength), min(1100, totallength))
                startpoint = random.randint(0, totallength - insertsize - 20)
                read1 = chimericdict[i][1][startpoint:startpoint + 150]
                read2 = chimericdict[i][1][startpoint + insertsize - 150:startpoint + insertsize]
                if 10 <= CountUpper(read1) + CountUpper(read2) <= 290:
                    break
                if whilecount > 10:
                    break
            if whilecount > 10:
                break
            SNPnum = numpy.random.binomial(300, SNPLevel * 4 / 3)
            CNVLength = min(numpy.random.binomial(300, CNVLevel) * (2 + numpy.random.negative_binomial(1,0.7)), 15)
            for k in range(SNPnum):
                changepos = numpy.random.randint(0, 150)
                a = numpy.random.randint(0, 2)
                newBase = numpy.random.choice(['A', 'C', 'G', 'T'], 1)[0]
                if a == 0:
                    read1 = read1[:changepos] + newBase + read1[changepos + 1:]
                else:
                    read2 = read2[:changepos] + newBase + read2[changepos + 1:]
            if CNVLength > 0:
                a = numpy.random.uniform(0, 1)
                b = numpy.random.uniform(0, 1)
                changepos = numpy.random.randint(15, 135 - CNVLength)
                if a < 0.5:
                    if b < 0.5:
                        newread1 = read1[:changepos]
                        for m in range(CNVLength):
                            newread1 += numpy.random.choice(['A', 'C', 'G', 'T'], 1)[0]
                        newread1 += read1[changepos:-CNVLength]
                        read1 = newread1
                    else:
                        newread2 = read2[:changepos]
                        for m in range(CNVLength):
                            newread2 += numpy.random.choice(['A', 'C', 'G', 'T'], 1)[0]
                        newread2 += read2[changepos:-CNVLength]
                        read2 = newread2
                else:
                    if b < 0.5:
                        read1 = read1[:changepos] + read1[changepos + CNVLength:] + chimericdict[i][1][startpoint + 150:startpoint + 150 + CNVLength]
                    else:
                        read2 = chimericdict[i][1][startpoint + insertsize - 150 - CNVLength:startpoint + insertsize - 150] + \
                                read2[:changepos - CNVLength] + read2[changepos:]
            outfile1.write('@' + chimericdict[i][0] + '_' + str(j) + 'Split\n')
            outfile2.write('@' + chimericdict[i][0] + '_' + str(j) + 'Split\n')
            outfile1.write(read1.upper() + '\n+\n' + qualitystring + '\n')
            outfile2.write(ReverseComplement(read2.upper()) + '\n+\n' + qualitystring + '\n')
    outfile2.close()
    outfile1.close()
