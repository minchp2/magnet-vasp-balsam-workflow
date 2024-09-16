qstat -u $USER
#balsam job ls --state RUNNING
balsam job ls --by-state
echo -------------------------
for j in $(balsam job ls --state RUNNING | awk 'NR > 1{print $4}'); do echo "\n"; echo $j; tail janus_db/data/$j/vasp.out; done

