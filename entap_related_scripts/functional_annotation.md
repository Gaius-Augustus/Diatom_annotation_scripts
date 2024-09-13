## Functional Annotation with EnTAP (v.1.3.0)

* Run Parameters: 
    - entap_v1.3.0_params
* Configuration:
    - entap_v1.3.0_config.ini
* Databases Searched:
    - RefSeq: complete.protein.faa.216.dmnd
    - UniProt: uniprot_sprot.fa.2.0.6.dmnd
* Contaminants Flagged:
    - Bacteria
* Dependency Versions:
    - DIAMOND: 2.1.8
    - EggNOG-mapper: 2.1.7

This is a simplified representation of how EnTAP was run using *Fistulifera pelliculosa*:
```
EnTAP --runP --entap-ini entap_v1.3.0_config.ini --run-ini entap_v1.3.0_params --out-dir entap_outfiles --input braker.aa --taxon Fistulifera_pelliculosa
```
The actual script (entap_v1.3.0.sh) loops through all the species.

## Convert BRAKER GTF to GFF with AGAT (v.1.4.0)
In order to append functional information using the Parent/ID relationship, the structural annotation needed to be converted into a GFF. Unfortunately, AGAT mishandled pre-leading introns, resulting in misplaced UTRs/exons. These were removed. 

```
agat_convert_sp_gxf2gxf.pl -g ${gff} -o ${id}.gff
```
To parallelize this process, Nextflow was leveraged, see:
* agat.sh
* agat.nf 

## Decorate GFF File with Functional Information
[NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) and [pdb_pfam_mapping.txt](https://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt) are required to decorate the GFF. Instead of using the RefSeq accession (XP*) for dbxref, we pull the GeneID. Likewise, the PFAM accession is needed instead of the name. 

```
python ncbi_gff.py entap_results.tsv Fistulifera_pelliculosa.gff pdb_pfam_mapping.txt Fistulifera_pelliculosa_ncbi.gff
```

This script will append the following to the 9th column:
* gene_biotype=protein_coding [gene]
* product [mRNA, CDS]
* Dbxref
    - GeneID [all features]
    - UniProtKB/Swiss-Prot [all features]
    - PFAM [mRNA, CDS]
    - GO [mRNA, CDS]
* note=EggNOG description [mRNA, CDS]

The full script (loop.sh) iterates through all the species.