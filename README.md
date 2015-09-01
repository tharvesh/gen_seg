# gen_seg
To create segment of genome for given windows and step size
To convert input bam to bed file
To count number of reads falls in the segments

Usage: gen_seg.py [options]

Options:
  -h, --help            show this help message and exit
  -d FOLDER, --directory=FOLDER
                        Path to directory that contains uniq.sorted bam files
  -g FILE, --genome-index=FILE
                        Path to genome index
  -w INT, --window-size=INT
                        Size of the window to use
  -s INT, --step-size=INT
                        Base pairs to step before creating a new window
  -b FILE, --bedmap=FILE
                        Path to bedmap bin
  -o STRING, --outfile=STRING
                        Basename of the output file
  -p STRING, --isPair=STRING
                        Reads are paired or not(yes = paired, no = single)
