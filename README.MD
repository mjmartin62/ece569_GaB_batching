
# GPU Based QC LDPC Error Correction: Throughput optimization versus serial code baseline
This project explores the ability to speed up QC LDPC error correction stage as a forward error correction algorithm on a GPU.

### Introduction:
This project explores the implementation of GPU-based error correction for Low-Density Parity Check (LDPC) codes within a test bench framework. Leveraging the parallelization capabilities of GPUs, the study achieves a 20X throughput improvement over traditional CPU-based implementations.

### Project Highlights:
 - Efficient Parallelization: Implementation of GPU-accelerated LDPC error correction using CUDA, optimizing key kernels for speedup.
 - Structured Optimization: Utilized Design of Experiments (DOE) methodology to fine-tune parameters like block size, streams, and codewords per stream, maximizing throughput.
 - Performance Validation: Demonstrated equivalent reliability between GPU and CPU implementations through FER (Frame Error Rate) statistical comparisons.

### Key findings:
 - Significant Kernel Speedup: Optimized syndrome computation kernel achieved an 80X speedup, drastically improving test bench performance.
 - Batch Processing Efficiency: Combining CUDA streams with batched codeword processing minimized GPU idle time, enhancing throughput.
 - Overall Throughput Improvement: Final GPU implementation achieved a 20X increase in throughput compared to the baseline serial execution.

### Software Requirements:
 - CUDA 11.0

### Installation:
 - Ensure you have `gcc`, `make`, and other necessary tools installed.

### Usage:

The executable requires arguments that provide GPU configuration parameters to support automated testing and optimzation.  Command line arguments are as follows and can be used according to the example commands below:

```bash
./bin/main ~/proj/source/GaB/<H-matrix> <output file> <Stream Size> <Block Size> <Block Size> <Num Words>
 ```
 - H-matrix:    matrix in text format that connects the nodes
 - output file: output file with results
 - Stream Size: number of cuda streams
 - Block Size:  cuda threads per block
 - Num Words:   codewords per cuda stream
 - reserved:    n/a
 
Example execution from command line:
```bash
./bin/main ~/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res 35 32 35 0
 ```
Example terminal output:
```bash   
alpha		NbEr(BER)		  NbFer(FER)		  Nbtested    IterAver(Itermax)		NbUndec(Dmin)
0.02000	    373 (0.00001740)  100 (0.00604741)	  16536		  3.26(83)		        0(100000)
 ```

 Example build command from command line:
```bash
 gcc -o GaB GaB.c -lm
 ```

