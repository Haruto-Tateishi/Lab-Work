# this script creates a table that contains all the info from vcf file and helps generating SFS.

import pysam

fn1 = "1kG-chr22-hg19-snpEff-prcsd-codon-reg-exon-CpG-repeat.vcf"
f1 = pysam.VariantFile(fn1)

fn2 = "1kG-chr22-hg19-table.tsv"
f2 =open(fn2, "w")

header = [
    "CHROM", "POS", "VCF_ID", "REF", "ALT", "REF_CODON", "ALT_CODON","ANN", "FTR", "REG", "SRC", "EXO_INT", "STR", "CPG", "REP", "TOTAL_AF", "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF", "TOTAL_AC", "GBR_AC", "FIN_AC", "CHS_AC", "PUR_AC", "CDX_AC", "CLM_AC", "IBS_AC", "PEL_AC", "PJL_AC", "KHV_AC", "ACB_AC", "GWD_AC", "ESN_AC", "BEB_AC", "MSL_AC", "STU_AC", "ITU_AC", "CEU_AC", "YRI_AC", "CHB_AC", "JPT_AC", "LWK_AC", "ASW_AC", "MXL_AC", "TSI_AC", "GIH_AC", "\n"
    ]

header_line = "\t".join(header)
f2.write(header_line)

for record in f1:
    print(record.pos)
    snpEff = record.info['ANN'][0].split(sep="|")
    ann = snpEff[1]
    ftr = snpEff[5]
    if type(record.info['Str']) == tuple:
        strand = "+/-"
    else:
        strand = record.info['Str']
    if type(record.info['Src']) == tuple:
        sources = record.info['Src']
        source_list = []
        for source in sources:
            source_list.append(source)
        sources = "/".join(source_list)
    else:
        sources = record.info['Src']
    list = [
        str(record.chrom), str(record.pos), record.id, record.ref, record.alts[0], record.info['RefSeq'], record.info['MutSeq'], ann, ftr, record.info['Reg'], sources, record.info['Exo/Int'], strand, record.info['CpG'], record.info['Rep'], str(record.info['AF'][0]), str(record.info['EAS_AF'][0]), str(record.info['AMR_AF'][0]), str(record.info['AFR_AF'][0]), str(record.info['EUR_AF'][0]), str(record.info['SAS_AF'][0]), str(record.info['AC'][0])
        ]
    counts = str(record).split(sep="\t")[9:]
    for count in counts:
        list.append(str(count))
    # for element in list:
    #     if type(element) != str:
    #         print(element, type(element))
    line = "\t".join(list)
    f2.write(line)

f1.close()
f2.close()