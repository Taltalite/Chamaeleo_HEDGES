from Chamaeleo.methods.default import AbstractCodingAlgorithm
from Chamaeleo.utils.indexer import connect
from Chamaeleo.methods import inherent

'''
Typical values that we use are p = 8, q = 10, s = 46, so that p + q + s = 64 bits
s bits of "salt"
low q bits of index i
p bits of previous bits  
'''


class HEDGES(AbstractCodingAlgorithm):
    def __init__(self, code_rate=0.5, p=8, q=10, s=46, need_logs=False):
        super().__init__(need_logs)
        self.s = s
        self.NPREV = p  # number of hashed previous bits
        self.NSEQBITS = q  # number of hashed sequence number bits
        self.HSALT = s  # number of hashed bits of salt

        self.VSALT = 0

        self.prevmask = (1 << self.NPREV) - 1
        self.saltmask = (1 << self.HSALT) - 1

        # just a randomly specified value
        self.MAXSEQ = 2500

        self.GC_content = [0.4, 0.6]  # set GC content: 40% ~ 60%
        self.MAXRUN = 3  # set max number of a nucleotide run
        self.NUMSTARTCHECK = 12  # set the start number of checking the dna constraints

        self.pattarr = [1 for _ in range(self.MAXSEQ + 2)]
        self.npattrn = 0
        self.pattnumber = code_rate

        self.dna_constraints_ok = []

        self.seqnomask = (1 << self.NSEQBITS) - 1

        self.__init_check__()

    def __init_check__(self):
        allow = [0.75, 0.6, 0.5, 0.333, 0.25, 0.166]
        if self.NPREV <= 0:
            raise ValueError("The parameter \"p\" is wrong, it must be greater than 0!")
        if self.NSEQBITS <= 1:
            raise ValueError("The parameter \"q\" is wrong, it must be greater than 1!")
        if self.HSALT <= 3:
            raise ValueError("The parameter \"s\" is wrong, it must be greater than 3!")
        if self.pattnumber not in allow:
            raise ValueError("The parameter \"code_rate\" is wrong, it must be 0.75, 0.6, 0.5, 0.333, 0.25 or 0.166 !")

        # set code rate
        self.setcoderate()

    def encode(self, bit_segments):
        dna_sequences = []

        # fit the bit_segments length
        # index_binary_length = int(len(str(bin(len(bit_segments)))) - 2)
        # not fit the bit_segments length
        index_binary_length = self.HSALT

        # add HSALT bits index to each bit segment
        connected_bit_segments = []
        for row in range(len(bit_segments)):
            connected_bit_segments.append(connect(row, bit_segments[row], index_binary_length))

        for segment_index, bit_segment in enumerate(connected_bit_segments):
            dna_sequence = []
            prevbits = 0
            prevcode = -1
            salt = 0
            newsalt = 0
            nucleotide_count = 0
            gc_count = 0
            run_check_count = 0

            for index, msgbit in enumerate(bit_segment):
                nbits = self.pattarr[index]
                mod = self.dnaconstraints(gc_count, nucleotide_count, run_check_count, prevcode)

                # set salt
                if index < self.HSALT:
                    salt = 0
                    newsalt = ((newsalt << 1) & self.saltmask) ^ msgbit
                elif index == self.HSALT:
                    salt = newsalt

                regout = self.digest(prevbits, index, salt, mod)
                regout = (regout + int(msgbit)) % mod

                indeed_nuc = self.dna_constraints_ok[regout]

                dna_sequence.append(inherent.index_base.get(indeed_nuc))

                prevbits = ((prevbits << nbits) & self.prevmask) | msgbit

                if indeed_nuc in (1, 2):  # code as C or G
                    gc_count += 1
                nucleotide_count += 1
                if prevcode == indeed_nuc:
                    run_check_count += 1
                else:
                    run_check_count = 1

                prevcode = indeed_nuc

            dna_sequences.append(dna_sequence)

            if self.need_logs:
                self.monitor.output(segment_index + 1, len(bit_segments))

        return dna_sequences

    def digest(self, prevbits, seq, salt, mod):
        return (self.ranhash(
            ((((seq & self.seqnomask) << self.NPREV) | prevbits) << self.HSALT) | salt
        ) % mod)

    # check if satisfied the dna constraints
    # return the value of mod
    def dnaconstraints(self, gc_count, nucleotide_count, run_check_count, prevcode):
        # initialize dna_constraints_ok
        self.dna_constraints_ok = []
        ans = 0

        def is_run_max(MAXrun, run):
            if run > MAXrun:
                return True
            else:
                return False

        if nucleotide_count < self.NUMSTARTCHECK:
            gc_count = 0.5  # pass the gc_count check
        else:
            gc_count = gc_count / nucleotide_count

        # the horrible logic tree same as the original author
        if gc_count < self.GC_content[0]:
            ans = 2
            # C and G are OK
            self.dna_constraints_ok.append(1)
            self.dna_constraints_ok.append(2)
            if is_run_max(self.MAXRUN, run_check_count):
                if prevcode == 1:  # prevcode is C
                    ans = 1
                    self.dna_constraints_ok[0] = 2  # only G is OK
                elif prevcode == 2:  # prevcode is G
                    ans = 1
                    self.dna_constraints_ok[0] = 1  # only C is OK

        elif gc_count > self.GC_content[1]:
            ans = 2
            # A and T are OK
            self.dna_constraints_ok.append(0)
            self.dna_constraints_ok.append(3)
            if is_run_max(self.MAXRUN, run_check_count):
                if prevcode == 0:  # prevcode is A
                    ans = 1
                    self.dna_constraints_ok[0] = 3  # only T is OK
                elif prevcode == 3:  # prevcode is T
                    ans = 1
                    self.dna_constraints_ok[0] = 0  # only A is OK

        else:  # no GC constraints
            ans = 4
            self.dna_constraints_ok.append(0)
            self.dna_constraints_ok.append(1)
            self.dna_constraints_ok.append(2)
            self.dna_constraints_ok.append(3)
            if is_run_max(self.MAXRUN, run_check_count):
                ans = 3
                for i in range(prevcode, 3):
                    self.dna_constraints_ok[i] = self.dna_constraints_ok[i + 1]

        return ans

    @staticmethod
    def ranhash(u):
        v = u * 3935559000370003845 + 2691343689449507681
        v ^= v >> 21
        v ^= v << 37
        v ^= v >> 4
        v *= 4768777513237032717
        v ^= v << 20
        v ^= v >> 41
        v ^= v << 5
        return v

    def setcoderate(self):
        pattrn = []
        if self.pattnumber == 0.75:  # rate 0.75
            pattrn.append(2)
            pattrn.append(1)
        elif self.pattnumber == 0.6:  # rate 0.6
            pattrn.append(2)
            for _ in range(4):
                pattrn.append(1)
        elif self.pattnumber == 0.5:  # rate 0.5
            pattrn.append(1)
        elif self.pattnumber == 0.333:  # rate 0.333
            pattrn.append(1)
            pattrn.append(1)
            pattrn.append(0)
        elif self.pattnumber == 0.25:  # rate 0.25
            pattrn.append(1)
            pattrn.append(0)
        elif self.pattnumber == 0.166:  # rate 0.166
            pattrn.append(1)
            pattrn.append(0)
            pattrn.append(0)

        self.npattrn = len(pattrn)
        for i in range(self.MAXSEQ):
            self.pattarr[i] = pattrn[i % self.npattrn]

        self.VSALT = self.vbitlen(self.HSALT)

    def vbitlen(self, nmb):
        nn = 0
        ksize = 0
        while True:
            if nn >= nmb:
                break
            if ksize >= self.MAXSEQ:
                raise ValueError("vbitlen: MAXSEQ is too small")
            nn += self.pattarr[ksize]
            ksize += 1
        return ksize
