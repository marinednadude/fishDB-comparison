## NCBI esearch

###Directory setup
```
mkdir NCBI_search
mkdir NCBI_search/names
mkdir NCBI_search/taxIDs
```

###Taxonomy DB Names Searches

ALL NAMES SEARCHES

Zack used this

```

for chunk_file in 0[124]_*; do taxa=$(cat "$chunk_file"); esearch -db taxonomy -query "$taxa" -retmax 1000000| efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element TaxId AkaTaxId Status Rank ScientificName Genus Species ModificationDate > "NCBI_taxonomy_${chunk_file}.txt"; sleep 10; done

cat NCBI_taxonomy_02_species_noNCBI_names_search* > NCBI_taxonomy_02_species_noNCBI_names_search_combined.txt

```
Sean's 

```
for chunk_file in 0[124]_*; do taxa=$(cat "$chunk_file"); esearch -db taxonomy -query "$taxa" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element TaxId AkaTaxId Status Rank ScientificName Genus Species ModificationDate > "NCBI_search/names/NCBI_taxonomy_${chunk_file}.txt"; sleep 10; done


for chunk_file in 0[124]_*; do taxa=$(cat "$chunk_file"); esearch -db taxonomy -query "$taxa" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element TaxId AkaTaxId Status Rank ScientificName Genus Species ModificationDate > "NCBI_taxonomy_${chunk_file}.txt"; sleep 10; done


cat NCBI_search/names/NCBI_taxonomy_01_noApiaID_names_search* > NCBI_search/names/NCBI_taxonomy_01_noApiaID_names_search_combined.txt
cat NCBI_search/names/NCBI_taxonomy_02_species_noNCBI_names_search* > NCBI_search/names/NCBI_taxonomy_02_species_noNCBI_names_search_combined.txt
cat NCBI_search/names/NCBI_taxonomy_04_subspecies_noNCBI_names_search* > NCBI_search/names/NCBI_taxonomy_04_subspecies_noNCBI_names_search_combined.txt

```

These will be imported back into R and brought into the master table. Then export another taxID search file for the species.

For reference, the old names search function:

```
for chunk_file in 0[124]_*; do taxa=$(cat "$chunk_file"); esearch -db nucleotide -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_search/names/NCBI_output_${chunk_file}.txt"; sleep 10; done
```

###Nucl TaxID Searches

ALL TAXID SEARCHES

```
for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db nucleotide -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_output_${chunk_file}.txt"; sleep 10; done

cat NCBI_output_03_species_singleNCBI_taxID_search_chunk_* > NCBI_output_03_species_singleNCBI_taxID_search_nucl_combined.txt


```


```
for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db nucleotide -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_search/taxIDs/NCBI_output_${chunk_file}.txt"; sleep 10; done



esearch -db nucleotide -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_search/taxIDs/NCBI_output_${chunk_file}.txt"; sleep 10; done


for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db nucleotide -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_output_${chunk_file}.txt"; sleep 10; done

for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db nuccore -query "$taxa AND mitochondrion[filter] AND 5000:5000000[SLEN]" -retmax 1000000 | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Title CreateDate Organism TaxId MolType Genome Completeness Topology Slen SubType SubName > "NCBI_output_${chunk_file}.txt"; sleep 10; done


cat NCBI_search/taxIDs/NCBI_output_03_species_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_03_species_singleNCBI_taxID_search_nucl_combined.txt
cat NCBI_search/taxIDs/NCBI_output_05_subspecies_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_05_subspecies_singleNCBI_taxID_search_nucl_combined.txt
cat NCBI_search/taxIDs/NCBI_output_06_genus_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_06_genus_singleNCBI_taxID_search_nucl_combined.txt

```

###WGS TaxID Searches

ALL TAXID SEARCHES

```
for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db sra -query "$taxa" | efetch -format runinfo > "NCBI_output_sra_${chunk_file}.txt"; sleep 10; done

for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db sra -query "$taxa" -retmax 100000 | efetch -format runinfo > "NCBI_output_sra_${chunk_file}.txt"; sleep 10; done


cat NCBI_output_sra_03_species_singleNCBI_taxID_search_chunk_* > NCBI_output_03_species_singleNCBI_taxID_search_sra_combined.txt
29d9fbc0a1617664225ea880e0e508116308	

```


```
for chunk_file in 0[356]_*; do taxa=$(cat "$chunk_file"); esearch -db sra -query "$taxa" -retmax 1000000 | efetch -format runinfo > "NCBI_search/taxIDs/NCBI_output_sra_${chunk_file}.txt"; sleep 10; done

cat NCBI_search/taxIDs/NCBI_output_sra_03_species_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_03_species_singleNCBI_taxID_search_sra_combined.txt

cat NCBI_search/taxIDs/NCBI_output_sra_05_subspecies_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_05_subspecies_singleNCBI_taxID_search_sra_combined.txt

cat NCBI_search/taxIDs/NCBI_output_sra_06_genus_singleNCBI_taxID_search_chunk_* > NCBI_search/taxIDs/NCBI_output_06_genus_singleNCBI_taxID_search_sra_combined.txt
```



