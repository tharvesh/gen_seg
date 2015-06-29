#!/usr/bin/python

# Importing modules

from optparse import OptionParser
import subprocess
import os
import os.path
import itertools


# Method to parse options
def parse_option():
    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="folder", help="Path to directory that contains uniq.sorted bam files",
                      metavar="FOLDER")
    parser.add_option("-g", "--genome-index", dest="g_idx", help="Path to genome index", metavar="FILE")
    parser.add_option("-w", "--window-size", dest="win_size", help="Size of the window to use", metavar="INT")
    parser.add_option("-s", "--step-size", dest="step_size", help="Base pairs to step before creating a new window",
                      default=None, metavar="INT")
    parser.add_option("-b", "--bedmap", dest="bedmap_path", help="Path to bedmap bin", default=None, metavar="FILE")
    parser.add_option("-o", "--outfile", dest="out_base", help="Basename of the output file", default=None,
                      metavar="STRING")
    (options, args) = parser.parse_args()

    bam_folder = options.folder
    g_idx = options.g_idx
    win_size = options.win_size
    bedmap_path = options.bedmap_path
    out_base = options.out_base

    if bedmap_path == None:
        print "Path to bedmap is missing"
        exit()

    if options.step_size != None:
        step_size = options.step_size
    else:
        step_size = win_size

    return (bam_folder, g_idx, win_size, step_size, bedmap_path, out_base)


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
    # Making bed folder
    """
    Returns path to directory of the bed
    :rtype : str
    """

    if os.path.isdir(bam_folder):
        dir_to_make = os.getcwd() + "/bed"
        if not os.path.exists(dir_to_make):
            os.makedirs(dir_to_make)
        else:
            pass
        bed_dir = dir_to_make
        for fl in os.listdir(bam_folder):
            if fl.endswith('.bam'):
                # Creating output filename
                outflname = os.path.splitext(fl)[0] + '.bed'
                outflname = os.path.join(bed_dir, outflname)
                # Creata an instance of file with write mode
                out = open(outflname, "w")
                fl = os.path.join(bam_folder, fl)
                cmd = ["bedtools", "bamtobed", "-i", fl]
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

                for line in p.stdout:
                    out.write("%s" % line)
                out.close()
        return bed_dir
    else:
        print 'Unable find the given bam directory'
        exit()


def run_bedmap(bedmap_path, bed_dir, created_windows_bed_fname):
    counts = []
    created_windows_bed_fname = os.getcwd() + "/" + created_windows_bed_fname
    for bed in os.listdir(bed_dir):
        bed = bed_dir + "/" + bed
        cmd = [bedmap_path, "--fraction-map", "0.5", "--count", created_windows_bed_fname, bed]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        temp = []
        for item in p.stdout:
            temp.append(item.strip())
        counts.append(temp)
    return counts


def format_output(counts, out_base, created_w_bed_fname):
    if out_base == None:
        for count in itertools.izip_longest(*counts):
            print "\t".join(str(i) for i in count)
    else:
        created_w_bed_fname = os.getcwd() + "/" + created_w_bed_fname

        # Chunk to write first column of the text file
        first_col = []
        for line in open(created_w_bed_fname).readlines():
            line = line.strip()
            line = line.replace("\t","_")
            first_col.append(line)

        out_mat = out_base + ".matrix"
        out_txt = out_base + ".txt"
        out1 = open(out_mat,"w")
        out2 = open(out_txt,"w")

        #Writing matrix
        for count in itertools.izip_longest(*counts):
            out1.write("\t".join(str(i) for i in count) + "\n")
        out1.close()

        #Writing text
        counts.insert(0,first_col)
        for count in itertools.izip_longest(*counts):
            out2.write("\t".join(str(i) for i in count) + "\n")
        out2.close()

    return None


def main():
    (bam_folder, g_idx, win_size, step_size, bedmap_path, out_base) = parse_option()

    # print bam_folder, g_idx,win_size,step_size
    created_windows_bed_fname = run_make_windows(g_idx, win_size, step_size)
    bed_dir = run_bamtobed(bam_folder)
    counts = run_bedmap(bedmap_path, bed_dir, created_windows_bed_fname)
    format_output(counts, out_base, created_windows_bed_fname)


if __name__ == "__main__":
    main()
