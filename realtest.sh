


#Clear old output
#rm test/generated/*


for entry in tests/in/*
do
    inp=`cat $entry`
    out=`./test.sh $inp`
    f=$(basename $entry)
    correct=`echo tests/out/$f`
    correctOut=`cat $correct`
    echo Testing $entry:
    diff <(echo "$out") <(echo "$correctOut")
done