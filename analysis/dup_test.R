## This script is confirming with data downloaded from Genbank that we do have 
## at least some set of individuals with identical mtDNA 

pedro_nu60 = ape::read.FASTA("data/pedro2020/nu60.fa")
pedro_nu83 = ape::read.FASTA("data/pedro2020/nu83.fa")
pedro_nu88 = ape::read.FASTA("data/pedro2020/nu88.fa")
li_NG1718 = ape::read.FASTA("data/pedro2020/NG1718_li2020.fa")

recreated = ape::read.FASTA("processed_data/papua_n.fasta")
recreated_set = recreated[c("NU60", "C001KLC-PNG387", "C001KM0-PNG419", "C001KNZ-PNG378", "C001KNP-PNG371")]

comp = c(pedro_nu60, recreated_set, pedro_nu83, pedro_nu88, li_NG1718)
ape::write.FASTA(comp, "processed_data/nu60_comparison.fa")
