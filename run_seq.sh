rm -r log/out_Sequential.txt
echo "<== log for Sequential Run ===>" >> log/out_Sequential.txt
./bin/Sequential > log/run_seq.txt
cat log/run_seq.txt >> log/out_Sequential.txt
rm -r log/run_seq.txt