#step1 :ready to prepare the ref
seq=$(grep -v ">" target_genomic.fa)
len=${#seq}

echo -e "genomic_pos\tbase\tcontext" > C_context.tsv

for ((i=0;i<len;i++)); do
  base=${seq:i:1}
  if [ "$base" = "C" ]; then
    next1=${seq:i+1:1}
    next2=${seq:i+2:1}

    if [ "$next1" = "G" ]; then
      context="CG"
    elif [ "$next2" = "G" ]; then
      context="CHG"
    else
      context="CHH"
    fi

    echo -e "$((i+1))\tC\t$context" >> C_context.tsv
  fi
done

#step2 :ref convert to bisulfite
seq=$(grep -v ">" target_genomic.fa)
echo ">bisulfite_forward" > bisulfite_forward.fa
echo "$seq" | tr C T >> bisulfite_forward.fa
cat bisulfite_forward.fa

#step3 :bisulfited-ref convert to reverse strand
revcomp=$(echo "$seq" | tr C T | rev | tr ATGC TACG)
echo ">bisulfite_reverse" > bisulfite_reverse.fa
echo "$revcomp" >> bisulfite_reverse.fa
cat bisulfite_reverse.fa

#step4 :mark the reverse strand's position
amp_len=200  
echo -e "genomic_pos\tcontext\treverse_pos\tstrand\tcalc" > C_reverse.tsv

tail -n +2 C_context.tsv | while read pos base context; do
    rev_pos=$((amp_len - pos + 1))
    echo -e "$pos\t$context\t$rev_pos\treverse\tG/(G+A)"
done >> C_reverse.tsv
