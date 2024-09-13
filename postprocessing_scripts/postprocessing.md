# Processing of functionally decorated gff3 file

Contact: katharina.hoff@uni-greifswald.de

```
add_mRNA_line.py -i decorated.gff3 -o mRNA.gff3
fix_product_names_ncbi.py -i mRNA.gff3 -o fixed_names.gff3
fix_Dbxref_attributes_in_genes.py -i fixed_names.gff3 -o fixed_dbxref.gff3
```

The contents of species-specific fixed_dbxref.gff3 files are provided at Zenodo with doi 10.5281/zenodo.13745090