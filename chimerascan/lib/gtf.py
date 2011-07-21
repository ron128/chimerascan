'''
Created on Nov 2, 2010

@author: mkiyer
'''
GTF_EMPTY_FIELD = '.'
GTF_ATTR_SEP = ';'
GTF_ATTR_TAGVALUE_SEP = ' '

class GTFFeature(object):
    '''
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.

    chr1    Cufflinks       transcript      136546  137059  1000    .       .       gene_id "VCAP_SHEZH2.657699"; transcript_id "VCAP_SHEZH2.657699.1"; FPKM "100.7219943204"; frac "1.000000"; conf_lo "80.649925"; conf_hi "120.794064"; cov "2.198209";
    '''
    __slots__ = ('seqid', 'source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 'attrs', 
                 'gene_id', 'transcript_id') 

    def __str__(self):
        line = [self.seqid,
                self.source,
                self.feature_type,
                str(self.start),
                str(self.end),
                str(self.score),
                str(self.strand),
                self.phase]
        attr_str = ' '.join('%s "%s";' % (k, v) for (k, v) in self.attrs.iteritems())
        line.append(attr_str)
        return '\t'.join(line)
        
    @staticmethod
    def from_string(line, attr_defs=None):
        f = GTFFeature()        
        # read the GTF line
        fields = line.strip().split('\t')
        f.seqid = fields[0]
        f.source = fields[1]
        f.feature_type = fields[2]
        f.start = int(fields[3])
        f.end = int(fields[4])
        f.score = 0 if (fields[5] == '.') else float(fields[5])
        strand = fields[6]
        if not (strand == '+' or strand == '-'):
            strand = GTF_EMPTY_FIELD
        f.strand = strand        
        f.phase = fields[7]        
        attrs = {}
        if fields[8] != GTF_EMPTY_FIELD:
            attr_strings = fields[8].split(GTF_ATTR_SEP)
            for a in attr_strings:
                a = a.strip()
                if len(a) == 0:
                    continue
                tag, value = a.split(GTF_ATTR_TAGVALUE_SEP, 1)
                # remove quotes
                value = value.split('"')[1]
                # apply parsing function
                if (attr_defs != None) and (tag in attr_defs) and (attr_defs[tag] != None):
                    value = attr_defs[tag](value)
                attrs[tag] = value
        f.attrs = attrs
        return f

    @staticmethod
    def parse(line_iter, attr_defs=None):
        for line in line_iter:
            # read the GTF line
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            yield GTFFeature.from_string(line, attr_defs)
