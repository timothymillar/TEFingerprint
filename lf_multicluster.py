import os
import pysam
from multiprocessing import Pool
from subprocess import Popen, PIPE, STDOUT

os.chdir('/Users/cfltxm/Repo/leapfrog')
INFILE = "quick_test/data_full/testDanglersMappedSorted.bam"


s = pysam.view("-b", "-f", "16", INFILE, "chr11")
p = Popen(['python', 'lf_cluster.py', '--reference', 'chr11'],
          stdout=PIPE, stdin=PIPE, stderr=STDOUT)
result = p.communicate(input=s)[0]
print(result.decode())


def cluster_job(bamfile, reference, family, reverse=False):
    bam = pysam.view("-b", "-f", "16", bamfile, reference)
    proc = Popen(['python', 'lf_cluster.py',  '--reference', reference],
                 stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    result = proc.communicate(input=bam)[0]
    print(result.decode())


cluster_job(INFILE, "chr11", "family", False)

with Pool(5) as p:
    p.starmap(cluster_job, [(INFILE, "chr11", "family", False),
                            (INFILE, "chr12", "family", False)])



