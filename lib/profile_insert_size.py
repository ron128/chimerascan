'''
Created on Jan 24, 2011

@author: mkiyer
'''
import array
import logging

# local imports
import pysam

def parse_pe_sam_file(bamfh):
    pe_reads = ([], [])
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        prev_qname = qname        
        pe_reads[mate].append(read)
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def profile_insert_sizes(bamfh, min_isize, max_isize, max_samples=None):
    isize_array = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))
    count = 0
    outside_range = 0
    unmapped = 0
    isoforms = 0
    # setup debugging logging messages
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every      
    for pe_reads in parse_pe_sam_file(bamfh):
        # progress log
        debug_count += 1
        if debug_count == debug_next:
            debug_next += debug_every
            logging.debug("Processed reads: %d" % (debug_count))
            logging.debug("Unique paired reads: %d" % (count))
            logging.debug("Unmapped: %d" % (unmapped))
            logging.debug("Ambiguous (isoforms): %d" % (isoforms))
            logging.debug("Outside range: %d" % (outside_range))
        if (max_samples is not None) and count > max_samples:
            break
        # only allow mappings where there is a single
        # insert size (multiple isoforms are ambiguous)
        isizes = set()        
        for r in pe_reads[0]:
            if r.is_unmapped:
                continue
            # get insert size
            isize = r.isize
            if isize < 0: isize = -isize
            isizes.add(isize)
        # insert size must be within range
        if len(isizes) == 0:
            unmapped += 1
        elif len(isizes) > 1:
            isoforms += 1
        else:
            isize = isizes.pop()
            if (min_isize <= isize <= max_isize):
                # store in array
                isize_array[isize - min_isize] += 1
                count += 1
            else:
                outside_range += 1
    return isize_array

def get_hist_stats(a):
    # find the mode
    mode = a.index(max(a))
    # find the median
    n = sum(a)    
    half_n = n / 2.0
    median = None
    mean = 0
    count = 0
    for i,x in enumerate(a): 
        mean = mean + i*x
        count += x
        if median is None and (count >= half_n):
            median = i
    # find the mean
    mean = mean / float(n)
    # find the standard deviation
    std = 0
    for i,x in enumerate(a):
        std = std + x*((i - mean)**2)
    std = (std / float(n-1))**0.5
    return mean, median, mode, std

def profile_isize_stats(bamfh, min_isize, max_isize, max_samples=None, min_samples=1000):
    isize_array = profile_insert_sizes(bamfh, min_isize, max_isize, 
                                       max_samples=max_samples)
    # if number of samples is small, use approximation instead
    n = sum(isize_array)
    if n < min_samples:
        mean = int((max_isize - min_isize)/2.0)
        std = int((max_isize - min_isize) / 4.0)
        logging.warning("Insert size profiling yielded less than %d "
                        " samples, using approximation instead" % (min_samples))
        return mean, mean, mean, std
    mean, median, mode, std = get_hist_stats(isize_array)
    return mean, median, mode, std
    
def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option('--min-fragment-length', dest="min_fragment_length", 
                      type="int", default=0)
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--max-samples', dest="max_samples", 
                      type="int", default=None)
    parser.add_option('-o', dest="output_file", default=None) 
    options, args = parser.parse_args()
    input_bam_file = args[0]
    bamfh = pysam.Samfile(input_bam_file, "rb")
    mean, median, mode, std = profile_isize_stats(bamfh,
                                                  options.min_fragment_length,
                                                  options.max_fragment_length, 
                                                  max_samples=options.max_samples)
    bamfh.close()
    if options.output_file is not None:
        f = open(options.output_file, "w")
        print >>f, '\t'.join(map(str, [mean, median, mode, std]))

if __name__ == '__main__':
    main()