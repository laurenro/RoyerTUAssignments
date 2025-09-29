import gzip
import itertools

def make_sub_vcf(infile, outfile, n=500):
    opener = gzip.open if infile.endswith('.gz') else open
    writer = gzip.open if outfile.endswith('.gz') else open

    with opener(infile, 'rt') as fin, writer(outfile, 'wt') as fout:
        # copy header lines
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
            else:
                # first variant line: write it and break to continue writing remaining (n-1) lines
                fout.write(line)
                break

        # write remaining n-1 variant lines
        for line in itertools.islice(fin, n - 1):
            fout.write(line)

# usage:
in_vcf = 'C:/Users/laure/Documents/BIOL5131/Week4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
out_vcf = "test2_500.vcf.gz"
make_sub_vcf(in_vcf, out_vcf, n=500)
print("Wrote", out_vcf, "with up to 500 variant rows")