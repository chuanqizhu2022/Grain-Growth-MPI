mpicc main.c -o main-mpi \
&& rm -f data/*.csv figs/*.png \
&& mpirun -n 6 ./main-mpi \
&& python plot2d.py