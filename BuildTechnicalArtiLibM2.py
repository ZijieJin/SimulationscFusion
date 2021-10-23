import sys
import random
import math


def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


def bondscore(str1, str2):
    score = 0
    str1 = str1.upper()
    str2 = str2.upper()
    if len(str1) != len(str2):
        return -1
    for i in range(len(str1)):
        if str1[i] == 'A' and str2[i] == 'T' or str1[i] == 'T' and str2[i] == 'A':
            score += 12
        if str1[i] == 'G' and str2[i] == 'C' or str1[i] == 'C' and str2[i] == 'G':
            score += 21
    return score


def poswisedist(str1, str2):
    maxsamecount = 0
    current = 0
    for i in range(len(str1)):
        if str1.upper()[i] == str2.upper()[i]:
            current += 1
            maxsamecount = max(maxsamecount, current)
        else:
            current = 0
    return len(str1) - maxsamecount


### Pars ###
LibrarySize = 30000
powerscale2 = 2


refinecdnas = {}
cdnafile = open(sys.argv[1])
chimericoutfile = open(sys.argv[2], 'w')
partialfile = open(sys.argv[3], 'w')
lastpos = ''
for line in cdnafile.readlines():
    if line.startswith('>'):
        if lastpos != '':
            if len(refinecdnas[lastpos]) < 800:
                refinecdnas.pop(lastpos)
        lastposlist = line[1:].split(' ')[2].split(':')
        lastpos = lastposlist[2] + ':' + lastposlist[3] + ':' + lastposlist[4]
        refinecdnas[lastpos] = ''
        continue
    refinecdnas[lastpos] += line.rstrip()
cdnafile.close()
if len(refinecdnas[lastpos]) < 800:
    refinecdnas.pop(lastpos)
librarycount = 0
randomcount = 0
while librarycount < LibrarySize:
    randomcount += 1
    if divmod(randomcount, 1000)[1] == 0:
        sys.stderr.write('RandomCount: ' + str(randomcount) + '. LibraryCount: ' + str(librarycount) + '\n')
    selectseqs = random.sample(refinecdnas.keys(), 2)
    seq1length = min(round(len(refinecdnas[selectseqs[0]]) * 0.75), 1000)
    seq2length = min(round(len(refinecdnas[selectseqs[1]]) * 0.75), 1000)
    if seq1length < 150 or seq2length < 150:
        continue
    #seq1start = random.randint(10, len(refinecdnas[selectseqs[0]]) - seq1length - 10)
    seq1start = refinecdnas[selectseqs[0]].rfind('GT') - seq1length
    #seq2start = random.randint(10, len(refinecdnas[selectseqs[1]]) - seq2length - 10)
    seq2start = refinecdnas[selectseqs[1]].find('AG', 10) + 2
    if seq1start < 0 or seq1start > len(refinecdnas[selectseqs[0]]) - seq1length - 10:
        continue
    if seq2start > len(refinecdnas[selectseqs[1]]) - seq2length or seq2start <= 6:
        continue
    seq1 = refinecdnas[selectseqs[0]][seq1start:seq1start + seq1length]
    seq2 = refinecdnas[selectseqs[1]][seq2start:seq2start + seq2length]
    seqinfo1 = selectseqs[0].split(':')
    seqinfo2 = selectseqs[1].split(':')
    seq1tail = refinecdnas[selectseqs[0]][seq1start+seq1length:seq1start+seq1length+6]
    brkpnt1pos = int(seqinfo1[1]) + seq1start + seq1length - 1
    direct1 = '+'
    thisseq = seq1
    seq2tail = refinecdnas[selectseqs[1]][seq2start-6:seq2start]
    brkpnt2pos = int(seqinfo2[1]) + seq2start
    direct2 = '+'
    thisseq += seq2
    totalbondscore = 0
    maxscore = -1
    for i in range(1):
        thisscore = bondscore(seq1tail[i:6], seq2[:6-i])
        totalbondscore += thisscore
        if thisscore > maxscore:
            maxscore = thisscore
            maxindex = i
    for i in range(1):
        thisscore = bondscore(seq2tail[:6-i], seq1[seq1length-6+i:])
        totalbondscore += thisscore
        if thisscore > maxscore:
            maxscore = thisscore
            maxindex = i + 6
    totalbondscore /= 2
    prob = 1 / (1 + math.exp(-(totalbondscore - 60)*3)) * 1000

    aaa = random.random()
    if aaa < prob:
        print(str(totalbondscore) + '\t' + str(aaa))
        seqpart1 = seq1
        seqpart2 = seq2
        thisread = seqpart1.upper() + seqpart2.lower()
        partial = thisread[seq1length-30:seq1length+30].upper()
        chimericoutfile.write('>' + seqinfo1[0] + ':' + str(brkpnt1pos) + ':' + direct1 + '|' + seqinfo2[0] + ':' + str(brkpnt2pos) + ':' + direct2 + '_' + str(totalbondscore) + '\n')
        chimericoutfile.write(thisread + '\n')
        partialfile.write(partial + '\t30\n')
        librarycount += 1

