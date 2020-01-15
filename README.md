### Spriggan

Note: Spriggan is still being tested

Spriggan trims and assembles short read bacterial sequence data, identifies antimicrobial resistance genes and assess assembly quality.

Spriggan uses Trimmomatic to trim reads, Shovill to assemble reads, Abricate to ID AR genes and Quast to assess assembly quality.

Spriggan uses Docker to maintain stability, reproducibility, portability and ease of use by keeping its dependencies in a controlled environment.

### Using Spriggan

The input to Spriggan is a reference genome and a file listing paths to short reads for assembly

```
usage: ./spriggan [-h] [-o output] [-t] reads reference

required arguments:
  reads             file listing read paths; read pairs must be in same directory
  reference         reference genome for Quast
  
optional arguments:
  -h, --help    show this help message and exit
  -o            output directory, defaults to working directory
  -t            threads to use
```

### Why name it Spriggan?

Spriggan is heavily based on the [Dryad pipeline](https://github.com/k-florek/dryad). Dryads and Spriggans are both nature-based fantastical creatures from Greek and Cornish mythology, respectively. 

### Author

Abigail C. Shockey
