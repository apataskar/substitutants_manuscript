# genes.txt is a list of gene names for which protein quantification is available in PDC commons. (Quantification is downloaded from PDC commons)
open (FH,"genes.txt");
@file =<FH>;

foreach $line (@file)
{
chomp ($line);

print $line."\n";

`sed 's/CHANGE_THIS/$line/g' extract_median.R > tmp1.R`;

`R CMD BATCH tmp1.R`;

}
