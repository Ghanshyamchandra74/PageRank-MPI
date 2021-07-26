#Script to run Programs
rm -r log/out_mpi.txt
echo "<== log for MPI ===>" >> log/out_mpi.txt
for ((n=1;n<=32;n=n*2)); do
    mpirun -np ${n} ./bin/MPI > log/1.txt
    cat log/1.txt >> log/out_mpi.txt
    rm -r log/1.txt 
done
