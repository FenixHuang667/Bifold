# Bifold

Bifold is a C program that samples RNA sequences for two given RNA secondary structures having the same length. The sampler generates sequences are that are thermodynamically stable to the two input structures, following the Boltzmann distribution. 

The program is written in C. Please cite the software as specified at the bottom of the paper. 


### Prerequisites

GNU Automake

C Standard Library


### Installing

Use the command: 

```
./configure
make
```
The executable file will be "Bifold" in ./src. 


## Running 

### Running HamSampler 

There are several options for Bifold. Use the following command to see all options. 
```
./src/Bifold -h  
```

To specify an input file, use the command
```
./src/Bifold -i input.in 
```
The default input file is "input.in". 

To specify an output file, use the command
```
./src/Bifold -o output.out 
```
The default output file is "output.out". 

To set the number of sampled sequences to be c (a positive integer), use the command
```
./src/Bifold -m c
```
The default number of sampled sequences is c = 1000. 


### Input file style

An example of an input can be found in "input.in". The input file consists of three lines. The first line the name of the input structure pair in FASTA style. The second and thrid lines are the input secondary structure pairs. A secondary structure is presented in dot/bracket form. We recommend the length of the input structure to be below 500. Here is an example of a valid input secondary structure.  
```
.(((((((.((((((......))))))((((((.....)))))).....((((((.....)))))))))))))....
```

Currently we do not support input with multiple structures. 

### Output file

User can find c (user specificed, otherwise default c=1000) entries of sampled sequences in the output file. The included output file "output.out" is an example by running the commands
```
./src/Bifold -d 10 
```

The problem is NP-hard, and its time complexity grows exponentially with the number of exposed vertices (defined in the article). The program will run a pre-assessment on this number. If the assessment is >20, the program will stop because the runtime will be too long. 


## Contact

If you have any question or bug reports about HamSampler, please feel free to conatct the author Fenix Huang at fenixprotoss@gmail.com.  




