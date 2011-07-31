for LINE in `ls | grep score`;
do
cat $LINE | tr -d \)\( | awk '{print $10 " " $5}' | sort -n | uniq | sort -n >> $LINE.ordered;
done

for LINE in `ls | grep ordered`;
do
cat $LINE | awk '{print $1}' | uniq -c >> $LINE.counts;
done
