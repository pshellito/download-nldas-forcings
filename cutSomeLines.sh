#! /bin/bash

echo 'script running'
baseDir='./forcingFromNasa/'
new='_new'

for ii in `seq 1002 1249`
  do 
  echo $ii
  sed -e '203334q' $baseDir$ii.txt > $baseDir$ii$new.txt
  mv $baseDir$ii$new.txt $baseDir$ii.txt
##  echo $baseDir$ii.txt
  
done

