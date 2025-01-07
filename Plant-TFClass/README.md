
## Formatting Plant-TFClass data to assign plant TFs to families

The original data was provided by Romain Blanc-Mathieu, matching supplemental table 1 of the paper:

    Blanc-Mathieu R, Dumas R, Turchi L, Lucas J, Parcy F (2024) Plant-TFClass: a structural classification for plant transcription factors. Trends Plant Sci. 29(1):40-51. https://doi.org/10.1016/j.tplants.2023.06.023

We modified it as follows:

    # 1) simplify table, matching UniProt ids to families
    perl -F"\t" -lane 'printf("%s\t%s %s\n",$F[5],$F[2],$F[3])' Table_S1_withPlantsRepresentativesProt.tsv > Plant-TFClass.tsv

    # 2) retrieve FASTAs of those ids in uniprot.org -> uniprot.faa

    # 3) prepare FASTA file with Plant-TFClass annotation in header
    perl format.pl Plant-TFClass.tsv uniprot.faa > Plant-TFClass.faa
