
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > out1
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > out2
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > out3

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > out11
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > out12
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > out13
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > out14
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > out15
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > out16

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > out21
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > out22
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > out23
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > out24
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > out25
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > out26

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > out31
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > out32
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > out33
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > out34
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > out35
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > out36

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > out41
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > out42
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > out43
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > out44
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > out45
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > out46

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > out51
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > out52
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > out53
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > out54
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > out55
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > out56

./clamr_gpuonly -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > out61
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > out62
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > out63
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > out64
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > out65
./clamr_gpuonly -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > out66

