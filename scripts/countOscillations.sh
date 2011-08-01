for LINE in `ls | grep ordered | grep -v count`
do 
echo "`echo $LINE | tr -d ".ordered" | tr -d "schill"`, `cat $LINE | awk '{if ($1 > 2) then num++} END {print num}'`" >> temp;
done 
cat temp | sort -n >> oscCounts;
rm temp;

