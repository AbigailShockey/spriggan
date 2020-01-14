### Spriggan

Note: Spriggan is still being tested

Spriggan trims and assembles short read bacterial sequence data, identifies antimicrobial resistance genes and assess assembly quality.

Spriggan uses Trimmomatic to trim reads, Shovill to assemble reads, Abricate to ID AR genes and Quast to assess assembly quality.

### Using Spriggan

The input to Spriggan is a reference genome and a file listing paths to short reads for assembly

```
usage: ./spriggan [-h] [-o output] [-t] reads reference

required arguments:
  reads             file of paths to reads
  reference         reference genome for quast
  
optional arguments:
  -h, --help    show this help message and exit
  -o            output directory, defaults to working directory
  -t            threads to use
```
