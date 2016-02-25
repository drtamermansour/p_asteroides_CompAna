#script_path=$(dirname "${BASH_SOURCE[0]}")

## merge the TPM for tissue specific files
targets=()
i=1
rm -f temp.*
for f in *.dataSummary_comp;do
  cat $f | tail -n+2 | awk '{print $2,$7}' | uniq > $f.gene
  cat $f | tail -n+2 | awk '{print $3,$6}' > $f.isoform
  target=${f%.dataSummary_comp}
  targets+=($target)
  if [ $i -eq 1 ];then cat $f.gene > temp.$i;else join -t" " --nocheck-order temp.$((i-1)) $f.gene > temp.$i;fi
  if [ $i -eq 1 ];then cat $f.isoform > isotemp.$i;else join -t" " --nocheck-order isotemp.$((i-1)) $f.isoform > isotemp.$i;fi
  ((i+=1))
done
#echo "geneName" "${targets[@]}" > allTissues_geneTPM
mv temp.$((i-1)) allTissues_geneTPM                         ## 726725 (7776 has no exp)
#echo "isoformName" "${targets[@]}" > allTissues_isoformTPM
mv isotemp.$((i-1)) allTissues_isoformTPM                   ## 880285 (11379 has no exp)
rm *.gene *.isoform *temp.*

echo "geneName" "${targets[@]}" > adultOnly_geneTPM
cat allTissues_geneTPM | awk '$2==0 && $3==0 && $4>0 && $5==0 && $6==0' >> adultOnly_geneTPM ## 265726
echo "geneName" "${targets[@]}" > adultOnly_isoformTPM
cat allTissues_isoformTPM | awk '$2==0 && $3==0 && $4>0 && $5==0 && $6==0' >> adultOnly_isoformTPM ## 298944

echo "geneName" "${targets[@]}" > larva_geneTPM
cat allTissues_geneTPM | awk '$2!=0 || $3!=0 || $5!=0 || $6!=0' >> larva_geneTPM ## 453222
echo "geneName" "${targets[@]}" > larva_isoformTPM
cat allTissues_isoformTPM | awk '$2!=0 || $3!=0 || $5!=0 || $6!=0' >> larva_isoformTPM ## 569961

echo "geneName" "${targets[@]}" > larvaOnly_geneTPM
cat larva_geneTPM | tail -n+2 | awk '$4==0' >> larvaOnly_geneTPM ## 106536
echo "geneName" "${targets[@]}" > larvaOnly_isoformTPM
cat larva_isoformTPM | tail -n+2 | awk '$4==0' >> larvaOnly_isoformTPM ## 141992

echo "geneName" "${targets[@]}" > unexp_geneTPM
cat allTissues_geneTPM | awk '$2==0 && $3==0 && $4==0 && $5==0 && $6==0' >> unexp_geneTPM ## 7777
echo "geneName" "${targets[@]}" > unexp_isoformTPM
cat allTissues_isoformTPM | awk '$2==0 && $3==0 && $4==0 && $5==0 && $6==0' >> unexp_isoformTPM ## 11380

echo "geneName" "${targets[@]}" > exp_geneTPM
cat allTissues_geneTPM | awk '$2!=0 || $3!=0 || $4!=0 || $5!=0 || $6!=0' >> exp_geneTPM ## 718948
echo "geneName" "${targets[@]}" > exp_isoformTPM
cat allTissues_isoformTPM | awk '$2!=0 || $3!=0 || $4!=0 || $5!=0 || $6!=0' >> exp_isoformTPM ## 868905

## merge the rawCounts for tissue specific files
#targets=()
#i=1
#rm temp.*
#for f in *.dataSummary_comp;do
#  cat $f | tail -n+2 | awk '{print $3,$5}' > $f.isoform
#  target=${f%.dataSummary_comp}
#  targets+=($target)
#  if [ $i -eq 1 ];then cat $f.isoform > isotemp.$i;else join -t" " --nocheck-order isotemp.$((i-1)) $f.isoform > isotemp.$i;fi
#  ((i+=1))
#done
#mv isotemp.$((i-1)) allTissues_isoformRaw                   ## 880285 (11379 has no exp)
#rm *.isoform *temp.*


