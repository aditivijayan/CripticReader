This is a rudimentary reader of QED data so that it can be processed by the code CRIPTIC.

This is a work in progress and currently only reads in the density information from a data file.

To get this code working, follow these steps - 

1. Clone this repo into your local directory of choice- 

    git clone --recursive https://github.com/aditivijayan/CripticReader.git
2. Do - "make". This should create an executable named "reader" in your folder.
3. Run the executable - ./reader
4. The output should be -
   Loading data from path: /Users/aditivijayan/Projects/SKA-obs/data/plt00001
   Loading data from: /Users/aditivijayan/Projects/CripticReader/data/plt00001/Level_0/Cell_D_00000
   Data from all_blocks: 2.48635e-30
   Block data: 2.48635e-30
