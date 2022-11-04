import numpy
import random
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


_nuc = numpy.array(["A", "T", "C", "G"])


def randSeq(length):
    seqArray = _nuc[[random.randint(0, 3) for i in range(length)]]
    return(MutableSeq("".join(seqArray)))


# parameters
N = 100000
N_hic_reads = 1000
hic_insert_size = 1000
hic_read_size = 150
snp_rate = .01
recomb_rate = .001
kmer_size = 11

# "ancestral" sequence
anc = randSeq(N)
# write as "reference"
SeqIO.write([SeqRecord(anc, id='chr1')], "ref.fa", "fasta")

# 2 sets of parent haplotypes with some SNPs
parents = [[], []]
for bb in range(2):
    for ii in range(2):
        par = MutableSeq(str(anc))
        for pos in range(len(par)):
            if random.random() < snp_rate:
                new_base = _nuc[random.randint(0, 3)]
                while par[pos] == new_base:
                    new_base = _nuc[random.randint(0, 3)]
                par[pos] = new_base
        parents[bb].append(par)

# recombine parents into child haplotype pair
haps = ['', '']
for bb in range(2):
    cur_par = 0
    for pos in range(len(parents[bb][0])):
        haps[bb] += parents[bb][cur_par][pos]
        if random.random() < recomb_rate:
            cur_par = 1 - cur_par

# make gfa: two bubbles, one from 10-20% of the seq, the other 80-90%
gfa_f = open('test.gfa', 'wt')
# first node
gfa_f.write('S\t1\t{}\n'.format(haps[0][:int(.1*N)]))
# first bubble
gfa_f.write('S\tPR.2\t{}\n'.format(haps[0][int(.1*N):int(.2*N)]))
gfa_f.write('S\tPR.3\t{}\n'.format(haps[1][int(.1*N):int(.2*N)]))
# middle region
gfa_f.write('S\t4\t{}\n'.format(haps[1][int(.2*N):int(.8*N)]))
# second bubble
gfa_f.write('S\tPR.5\t{}\n'.format(haps[0][int(.8*N):int(.9*N)]))
gfa_f.write('S\tPR.6\t{}\n'.format(haps[1][int(.8*N):int(.9*N)]))
# end region
gfa_f.write('S\t7\t{}\n'.format(haps[0][int(.9*N):]))
# edges
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(1, 'PR.2'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(1, 'PR.3'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.2', 4))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.3', 4))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(4, 'PR.5'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(4, 'PR.6'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.6', 7))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.5', 7))
gfa_f.close()

# simulate hic reads
hic_f1 = open('test_hic_1.fastq', 'wt')
hic_f2 = open('test_hic_2.fastq', 'wt')
readid = 0
for rr in range(N_hic_reads):
    hap = random.randint(0, 1)
    pos = random.randint(0, N - hic_insert_size - 2*hic_read_size)
    read = haps[hap][pos:(pos+hic_read_size)]
    hic_f1.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*hic_read_size))
    pos = pos+hic_read_size+hic_insert_size
    read = haps[hap][pos:(pos+hic_read_size)]
    hic_f2.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*hic_read_size))
    readid += 1
hic_f1.close()
hic_f2.close()

# get unique parent kmers
kmers = [{}, {}]
for pos in range(N-kmer_size):
    for par in range(2):
        for phap in range(2):
            kmers[par][str(parents[par][phap][pos:(pos+kmer_size)])] = True
# select unique ones
uniq_kmers = [[], []]
for km in kmers[0]:
    if km not in kmers[1]:
        uniq_kmers[0].append(km)
for km in kmers[1]:
    if km not in kmers[0]:
        uniq_kmers[1].append(km)
# sort them
uniq_kmers[0].sort()
uniq_kmers[1].sort()
# write fasta
km_f1 = open('test_kmers_pat.fa', 'wt')
kmid = 0
for km in uniq_kmers[0]:
    km_f1.write('>{}\n{}\n'.format(kmid, km))
    kmid += 1
km_f1.close()
km_f2 = open('test_kmers_mat.fa', 'wt')
kmid = 0
for km in uniq_kmers[1]:
    km_f2.write('>{}\n{}\n'.format(kmid, km))
    kmid += 1
km_f2.close()
