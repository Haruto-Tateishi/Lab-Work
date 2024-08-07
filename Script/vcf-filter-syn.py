# This script fileters the vcf file so that you can get vcf file only with lives having at least one "synonymous_variant" annotation.

import pysam
import gzip

for i in range(1, 23):
  fn1 = "1kG-chr{}-filetered-pantro6.vcf.gz".format(i)
  gunzipped = gzip.open(fn1, 'r')
  f1 = pysam.VariantFile(gunzipped, 'r')
  fn2 = "1kG-chr{}-syn.vcf"
  f2 = pysam.VariantFile(f2, 'w', header =  f1.header)
  for record in f1:
    try:
      vep = record.info["CSQ"]
      # vep = record.info["vep"]
    except KeyError:
      continue
    for annotation in vep:
      conseq = annotation.split(sep="|")[1]
      if conseq == "synonymous_variant":
        f2.write(record)
        break
      else:
        pass
    pass
  f1.close()
  f2.close()