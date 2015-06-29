#!/usr/bin/python

# Importing modules

from optparse import OptionParser
import subprocess
import os
import os.path


# Method to parse options
def parse_option():
    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="folder", help="Path to directory that contains uniq.sorted bam files",
                      metavar="FOLDER")
    parser.add_option("-g", "--genome-index", dest="g_idx", help="Path to genome index", metavar="FILE")
    parser.add_option("-w", "--window-size", dest="win_size", help="Size of the window to use", metavar="INT")
    parser.add_option("-s", "--step-size", dest="step_size", help="Base pairs to step before creating a new window",
                      default=None, metavar="INT")
    (options, args) = parser.parse_args()

    bam_folder = options.folder
    g_idx = options.g_idx
    win_size = options.win_size

    if options.step_size != None:
        step_size = options.step_size
    else:
        step_size = win_size

    return (bam_folder, g_idx, win_size, step_size)


# Method to make bed file with given window and step sizes
def run_make_windows(g_idx, win_size, step_size):
    if os.path.isfile(g_idx):
        outfile = 'hg19.' + win_size + '_window_' + step_size + '_steps.bed'
        out = open(outfile, "w")
        cmd = ["bedtools", "makewindows", "-g", g_idx, "-w", win_size, "-s", step_size]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in p.stdout:
            out.write("%s" % line)
        out.close()

        return outfile
    else:
        print "Unable to find fasta index(fai) file"
        exit()


def run_bamtobed(bam_folder):
    if os.path.isdir(bam_folder):
        for fl in os.listdir(bam_folder):
            if fl.endswith('.bam'):
                outflname = os.path.splitext(fl)[0] + '.bed'
                out = open(outflname, "w")
                fl = os.path.join(bam_folder,fl)
                cmd = ["bedtools", "bamtobed", "-i", fl]
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

                for line in p.stdout:
                    out.write("%s" % line)
        return None
    else:
        print 'Unable find the given bam directory'
        exit()

def run_bedmap():
    pass


def main():
    (bam_folder, g_idx, win_size, step_size) = parse_option()

    # print bam_folder, g_idx,win_size,step_size
    windows_bed = run_make_windows(g_idx, win_size, step_size)
    run_bamtobed(bam_folder)


if __name__ == "__main__":
    main()
