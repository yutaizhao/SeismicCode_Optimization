for i in {1..11}
do
./top-stencil ../config.txt ./res.txt
cat ./res.txt >> all_res_100.txt
echo "\n"
done

