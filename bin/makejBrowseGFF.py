sed -i 's/Parent=gene:ENSG[0-9]*;//g' /nfs/production/flicek/ensembl/havana/lucascortes/LEAP/synchronise/testFurther.gff
gff3sort $file > $file_sorted
bgzip $file
tabix $file