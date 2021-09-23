
NPREV = 8 # number of hashed previous bits
NSEQBITS = 10 # number of hashed sequence number bits
HSALT = 24 # number of hashed bits of salt

seqnomask = (1 << NSEQBITS) - 1

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


if __name__ == "__main__":
    u = 0b01010110
    salt = 0b0
    input = ((((u & seqnomask) << NPREV) | 0b10) << HSALT) | salt
    v = ranhash(input)
    print(bin(v))
    print(v)
    print((v + 1) % 4)
