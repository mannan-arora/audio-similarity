#!/bin/bash

# Run make
make

# Check if make was successful
if [ $? -eq 0 ]; then
    # Run the executable
    ./similarity
    
    # Delete object & dependency files & executable
    rm *.o
    rm *.d
    rm similarity
fi