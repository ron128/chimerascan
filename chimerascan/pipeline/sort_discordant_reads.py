'''
Created on Jan 29, 2011

@author: mkiyer
'''
'''
Created on Jan 17, 2011

@author: mkiyer
'''
import logging
import subprocess

def sort_discordant_reads(input_file, output_file):
    outfh = open(output_file, "w")
    subprocess.call(["sort", "-k7,7", "-k16,16", "-k8,8g", "-k17,17g", input_file], 
                    stdout=outfh)
    outfh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <chimeras.sorted.bedpe>")
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    logging.info("Sorting discordant BEDPE file")
    logging.info("Input file: %s" % (input_file))
    logging.info("Output file: %s" % (output_file))
    sort_discordant_reads(input_file, output_file)

if __name__ == '__main__':
    main()
