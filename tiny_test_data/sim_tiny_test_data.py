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
N_par_reads = 5000
par_read_size = 150
N_porec_reads = 500
porec_read_size = 500
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
gfa_f.write('S\tPR.1.1.1.0\t{}\n'.format(haps[0][int(.1*N):int(.2*N)]))
gfa_f.write('S\tPR.1.1.1.1\t{}\n'.format(haps[1][int(.1*N):int(.2*N)]))
# middle region
gfa_f.write('S\t4\t{}\n'.format(haps[1][int(.2*N):int(.8*N)]))
# second bubble
gfa_f.write('S\tPR.2.1.1.0\t{}\n'.format(haps[0][int(.8*N):int(.9*N)]))
gfa_f.write('S\tPR.2.1.1.1\t{}\n'.format(haps[1][int(.8*N):int(.9*N)]))
# end region
gfa_f.write('S\t7\t{}\n'.format(haps[0][int(.9*N):]))
# "floating" unphased region
gfa_f.write('S\t8\t{}\n'.format(haps[1][int(.9*N):]))
# edges
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(1, 'PR.1.1.1.0'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(1, 'PR.1.1.1.1'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.1.1.1.0', 4))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.1.1.1.1', 4))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(4, 'PR.2.1.1.0'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format(4, 'PR.2.1.1.1'))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.2.1.1.0', 7))
gfa_f.write('L\t{}\t+\t{}\t+\t0M\n'.format('PR.2.1.1.1', 7))
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
pat_f = open('test.pat.fq', 'wt')
mat_f = open('test.mat.fq', 'wt')
readid = 0
for rr in range(N_par_reads):
    hap = random.randint(0, 1)
    pos = random.randint(0, N - par_read_size)
    read = parents[0][hap][pos:(pos+par_read_size)]
    pat_f.write('@r{}_pat\n{}\n+\n{}\n'.format(readid, read,
                                               '~'*par_read_size))
    hap = random.randint(0, 1)
    pos = random.randint(0, N - par_read_size)
    read = parents[1][hap][pos:(pos+par_read_size)]
    mat_f.write('@r{}_mat\n{}\n+\n{}\n'.format(readid, read,
                                               '~'*par_read_size))
    readid += 1
pat_f.close()
mat_f.close()

# simulate gfase outputs
SeqIO.write([SeqRecord(MutableSeq(haps[0][:int(.9*N)]), id='contig_0')],
            "out_phase0.fa", "fasta")
SeqIO.write([SeqRecord(MutableSeq(haps[1][:int(.9*N)]), id='contig_1')],
            "out_phase1.fa", "fasta")
SeqIO.write([SeqRecord(MutableSeq(haps[1][int(.9*N):]), id='unphased')],
            "out_unphased.fa", "fasta")

# simulate poreC data
porec_f = open('test_porec.fastq', 'wt')
readid = 0
for rr in range(N_porec_reads):
    hap = random.randint(0, 1)
    pos = random.randint(0, N - hic_insert_size - 2*porec_read_size)
    read = haps[hap][pos:(pos+porec_read_size)]
    pos = pos+porec_read_size+hic_insert_size
    read += haps[hap][pos:(pos+porec_read_size)]
    porec_f.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*len(read)))
    readid += 1
porec_f.close()
