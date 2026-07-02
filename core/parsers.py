def load_gff(path):

    genes=[]

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            cols=line.strip().split("\t")

            gene=Gene(

                id=cols[8],

                contig=cols[0],

                start=int(cols[3]),

                end=int(cols[4]),

                strand=cols[6]

            )

            genes.append(gene)

    return genes