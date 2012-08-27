
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 0   > gpuonlyruns/out1
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 0   > gpuonlyruns/out2
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 0   > gpuonlyruns/out3

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 1   > gpuonlyruns/out11
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 1   > gpuonlyruns/out12
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 1   > gpuonlyruns/out13
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 1   > gpuonlyruns/out14
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 1   > gpuonlyruns/out15
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 1   > gpuonlyruns/out16

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 2   > gpuonlyruns/out21
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 2   > gpuonlyruns/out22
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 2   > gpuonlyruns/out23
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 2   > gpuonlyruns/out24
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 2   > gpuonlyruns/out25
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 2   > gpuonlyruns/out26

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 3   > gpuonlyruns/out31
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 3   > gpuonlyruns/out32
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 3   > gpuonlyruns/out33
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 3   > gpuonlyruns/out34
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 3   > gpuonlyruns/out35
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 3   > gpuonlyruns/out36

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 4   > gpuonlyruns/out41
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 4   > gpuonlyruns/out42
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 4   > gpuonlyruns/out43
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 4   > gpuonlyruns/out44
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 4   > gpuonlyruns/out45
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 4   > gpuonlyruns/out46

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 5   > gpuonlyruns/out51
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 5   > gpuonlyruns/out52
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 5   > gpuonlyruns/out53
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 5   > gpuonlyruns/out54
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 5   > gpuonlyruns/out55
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 5   > gpuonlyruns/out56

./clamr_gpucheck -n 256 -i 100 -t 1000 -P "original_order" -p "original_order" -l 6   > gpuonlyruns/out61
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "z_order"        -l 6   > gpuonlyruns/out62
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "hilbert_sort"   -l 6   > gpuonlyruns/out63
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "z_order"        -p "original_order" -l 6   > gpuonlyruns/out64
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "original_order" -l 6   > gpuonlyruns/out65
./clamr_gpucheck -n 256 -i 100 -t 1000 -P "hilbert_sort"   -p "local_hilbert"  -l 6   > gpuonlyruns/out66

