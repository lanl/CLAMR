#!/bin/sh

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > gpuonlyruns_256/out01
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > gpuonlyruns_256/out02
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > gpuonlyruns_256/out03

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > gpuonlyruns_256/out11
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > gpuonlyruns_256/out12
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > gpuonlyruns_256/out13
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > gpuonlyruns_256/out14
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > gpuonlyruns_256/out15
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > gpuonlyruns_256/out16

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > gpuonlyruns_256/out21
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > gpuonlyruns_256/out22
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > gpuonlyruns_256/out23
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > gpuonlyruns_256/out24
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > gpuonlyruns_256/out25
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > gpuonlyruns_256/out26

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > gpuonlyruns_256/out31
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > gpuonlyruns_256/out32
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > gpuonlyruns_256/out33
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > gpuonlyruns_256/out34
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > gpuonlyruns_256/out35
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > gpuonlyruns_256/out36

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > gpuonlyruns_256/out41
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > gpuonlyruns_256/out42
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > gpuonlyruns_256/out43
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > gpuonlyruns_256/out44
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > gpuonlyruns_256/out45
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > gpuonlyruns_256/out46

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > gpuonlyruns_256/out51
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > gpuonlyruns_256/out52
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > gpuonlyruns_256/out53
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > gpuonlyruns_256/out54
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > gpuonlyruns_256/out55
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > gpuonlyruns_256/out56

./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > gpuonlyruns_256/out61
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > gpuonlyruns_256/out62
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > gpuonlyruns_256/out63
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > gpuonlyruns_256/out64
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > gpuonlyruns_256/out65
./clamr_gpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > gpuonlyruns_256/out66


./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > cpuonlyruns_256/out00
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > cpuonlyruns_256/out01
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > cpuonlyruns_256/out02
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > cpuonlyruns_256/out03

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > cpuonlyruns_256/out10
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > cpuonlyruns_256/out11
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > cpuonlyruns_256/out12
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > cpuonlyruns_256/out13
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > cpuonlyruns_256/out14
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > cpuonlyruns_256/out15
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > cpuonlyruns_256/out16

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > cpuonlyruns_256/out20
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > cpuonlyruns_256/out21
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > cpuonlyruns_256/out22
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > cpuonlyruns_256/out23
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > cpuonlyruns_256/out24
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > cpuonlyruns_256/out25
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > cpuonlyruns_256/out26

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > cpuonlyruns_256/out30
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > cpuonlyruns_256/out31
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > cpuonlyruns_256/out32
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > cpuonlyruns_256/out33
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > cpuonlyruns_256/out34
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > cpuonlyruns_256/out35
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > cpuonlyruns_256/out36

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > cpuonlyruns_256/out40
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > cpuonlyruns_256/out41
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > cpuonlyruns_256/out42
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > cpuonlyruns_256/out43
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > cpuonlyruns_256/out44
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > cpuonlyruns_256/out45
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > cpuonlyruns_256/out46

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > cpuonlyruns_256/out50
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > cpuonlyruns_256/out51
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > cpuonlyruns_256/out52
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > cpuonlyruns_256/out53
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > cpuonlyruns_256/out54
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > cpuonlyruns_256/out55
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > cpuonlyruns_256/out56

./clamr_cpuonly -N "kdtree"     -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > cpuonlyruns_256/out60
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > cpuonlyruns_256/out61
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > cpuonlyruns_256/out62
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > cpuonlyruns_256/out63
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > cpuonlyruns_256/out64
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > cpuonlyruns_256/out65
./clamr_cpuonly -N "hash_table" -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > cpuonlyruns_256/out66


./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > gpuonlyruns_512/out01
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > gpuonlyruns_512/out02
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > gpuonlyruns_512/out03

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > gpuonlyruns_512/out11
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > gpuonlyruns_512/out12
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > gpuonlyruns_512/out13
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > gpuonlyruns_512/out14
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > gpuonlyruns_512/out15
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > gpuonlyruns_512/out16

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > gpuonlyruns_512/out21
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > gpuonlyruns_512/out22
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > gpuonlyruns_512/out23
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > gpuonlyruns_512/out24
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > gpuonlyruns_512/out25
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > gpuonlyruns_512/out26

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > gpuonlyruns_512/out31
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > gpuonlyruns_512/out32
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > gpuonlyruns_512/out33
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > gpuonlyruns_512/out34
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > gpuonlyruns_512/out35
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > gpuonlyruns_512/out36

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > gpuonlyruns_512/out41
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > gpuonlyruns_512/out42
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > gpuonlyruns_512/out43
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > gpuonlyruns_512/out44
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > gpuonlyruns_512/out45
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > gpuonlyruns_512/out46

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > gpuonlyruns_512/out51
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > gpuonlyruns_512/out52
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > gpuonlyruns_512/out53
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > gpuonlyruns_512/out54
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > gpuonlyruns_512/out55
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > gpuonlyruns_512/out56

./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > gpuonlyruns_512/out61
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > gpuonlyruns_512/out62
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > gpuonlyruns_512/out63
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > gpuonlyruns_512/out64
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > gpuonlyruns_512/out65
./clamr_gpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > gpuonlyruns_512/out66


./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > cpuonlyruns_512/out00
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > cpuonlyruns_512/out01
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > cpuonlyruns_512/out02
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > cpuonlyruns_512/out03

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > cpuonlyruns_512/out10
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > cpuonlyruns_512/out11
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > cpuonlyruns_512/out12
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > cpuonlyruns_512/out13
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > cpuonlyruns_512/out14
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > cpuonlyruns_512/out15
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > cpuonlyruns_512/out16

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > cpuonlyruns_512/out20
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > cpuonlyruns_512/out21
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > cpuonlyruns_512/out22
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > cpuonlyruns_512/out23
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > cpuonlyruns_512/out24
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > cpuonlyruns_512/out25
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > cpuonlyruns_512/out26

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > cpuonlyruns_512/out30
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > cpuonlyruns_512/out31
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > cpuonlyruns_512/out32
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > cpuonlyruns_512/out33
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > cpuonlyruns_512/out34
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > cpuonlyruns_512/out35
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > cpuonlyruns_512/out36

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > cpuonlyruns_512/out40
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > cpuonlyruns_512/out41
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > cpuonlyruns_512/out42
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > cpuonlyruns_512/out43
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > cpuonlyruns_512/out44
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > cpuonlyruns_512/out45
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > cpuonlyruns_512/out46

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > cpuonlyruns_512/out50
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > cpuonlyruns_512/out51
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > cpuonlyruns_512/out52
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > cpuonlyruns_512/out53
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > cpuonlyruns_512/out54
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > cpuonlyruns_512/out55
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > cpuonlyruns_512/out56

./clamr_cpuonly -N "kdtree"     -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > cpuonlyruns_512/out60
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > cpuonlyruns_512/out61
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > cpuonlyruns_512/out62
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > cpuonlyruns_512/out63
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > cpuonlyruns_512/out64
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > cpuonlyruns_512/out65
./clamr_cpuonly -N "hash_table" -n 512 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > cpuonlyruns_512/out66

