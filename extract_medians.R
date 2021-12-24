#EXPE.txt is the file that contains on columns, W>F, W>Y counts and protein expression for all proteins, on rows: individual tumor entry). User has to format this file on his own computer


expe = read.table ("EXPE.txt", sep="\t", head=T)
expe_up = expe[ expe$CHANGE_THIS > 0,]
expe_down = expe[ expe$CHANGE_THIS < 0,]

lC_up = length (expe_up$id)
lC_down = length (expe_down$id)

wf_up_median= median (expe_up$wf, na.rm=T )
wf_down_median= median (expe_down$wf , na.rm=T)

wf_up_median= median (expe_up$wf, na.rm=T )
wf_down_median= median (expe_down$wf , na.rm=T)

wy_up_median= median (expe_up$wy, na.rm=T )
wy_down_median= median (expe_down$wy , na.rm=T)

wy_up_median= median (expe_up$wy, na.rm=T )
wy_down_median= median (expe_down$wy , na.rm=T)



numC = paste (wf_up_median, wf_down_median, sep="\t")

numC = paste (numC, wy_up_median, sep="\t")
numC = paste (numC, wy_down_median, sep="\t")

numC = paste (numC, lC_up, sep="\t")
numC = paste (numC, lC_down, sep="\t")

outC= paste ("CHANGE_THIS", numC, sep="\t")

write.table (outC, "GENES_WF_WY_expe_added.txt", append=T, col.names=F, quote=F)



