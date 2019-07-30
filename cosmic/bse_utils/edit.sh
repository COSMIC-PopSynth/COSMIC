for i in $(seq 1 200);
do
    echo "s/\[${i}\]/\[$((i-1))\]/g"
    sed -i "" "s/\[${i}\]/\[$((i-1))\]/g" zcnsts.py;
done
