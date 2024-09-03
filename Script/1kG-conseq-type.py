# this script looks through the annotated 1kG vcf files and picks up the consequences that happen simultaneously with specific consequence.

import pysam
import gzip

def conseq_extraction(input_file, conseq_1, output_list_1, conseq_2, output_list_2):
  gunzipped = gzip.open(input_file, "r")
  in_file = pysam.VariantFile(gunzipped)
  for record in in_file:
    conseq_list = []
    try:
      vep = record.info["CSQ"]
      # vep = record.info["vep"]
    except KeyError:
      continue
    if len(vep) == 1:
      continue
    for annotation in vep:
      conseq = annotation.split(sep="|")[1]
      conseq_list.append(conseq)
    if conseq_1 in conseq_list:
      for conseq in conseq_list:
        if conseq != conseq_1 and conseq not in output_list_1:
          output_list_1.append(conseq)
          pass
        else:
          continue
    if conseq_2 in conseq_list:
      for conseq in conseq_list:
        if conseq != conseq_2 and conseq not in output_list_2:
          output_list_2.append(conseq)
          pass
        else:
          continue
  output_list_1.sort()
  output_list_2.sort()
  return output_list_1, output_list_2


def write_out(output_file, longer_conseq, longer_list, longer_len, shorter_conseq, shorter_list):
  out = open(output_file, "w")
  header = [longer_conseq, shorter_conseq, "\n"]
  header = "\t".join(header)
  out.write(header)
  for i in range(longer_len):
    try:
      out_list = [longer_list[i], shorter_list[i], "\n"]
      out_line = "\t".join(out_list)
      out.write(out_line)
    except IndexError:
      out_list = [longer_list[i], "\n"]
      out_line = "\t".join(out_list)
      out.write(out_line)


syn_out = []
mis_out = []
conseq_1 = "synonymous_variant"
conseq_2 = "missense_variant"
for i in range(1, 23):
  fn1 = f"vcf/1kG-chr{i}-vep-every.vcf.gz"
  print("start ",i,end="")
  syn_out, mis_out = conseq_extraction(fn1,conseq_1, syn_out,conseq_2, mis_out)
  print(" done")


fn1 = "1kG-vep-conseq-syn-mis.tsv"
syn_len = len(syn_out)
mis_len = len(mis_out)
if syn_len>=mis_len:
  write_out(fn1, conseq_1, syn_out, syn_len, conseq_2, mis_out)
else:
  write_out(fn1, conseq_2, mis_out, mis_len, conseq_1, syn_out)