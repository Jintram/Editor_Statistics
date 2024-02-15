# Editor_Statistics
Quantification of what various DNA editors can target

# To do:
Add explanation here how to download necessary ClinVar and ref seq Fasta files.

# How to download genome data

The genome information was downloaded via
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

using `curl --remote-name --remote-time https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz`

I used a few simple lines of Python to split these files.

```

mydirectory='/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/'
myfilepath='GCF_000001405.40_GRCh38.p14_genomic.fna'

file1=None
c=0
with open(mydirectory+myfilepath) as f:
    for line in f:
        
        
        if (line[0]=='>'):
            
            if (file1 != None):
                file1.close()
                
            print(line)
            
            outputfile=line.split(' ')[0][1:]
            file1 = open(mydirectory+outputfile+".fna", "w") 
            
        
        file1.write(line)
        
```

(Note: this also results in some files that are not necessary, I simply didn't use those.)