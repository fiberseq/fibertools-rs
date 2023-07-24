import pyft
import tqdm

bam_f = "../tmp.bam"
fiberbam = pyft.Fiberbam(bam_f)
for fiber in tqdm.tqdm(fiberbam.fetch("chr20", 0, 10_000_000)): 
    # the number of ccs passes
    fiber.ec
    # the mps start positions
    fiber.msp.starts
    # print the nuc reference starts
    fiber.nuc.reference_starts
    # lift query (fiber) positions to reference positions
    fiber.lift_query_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    # lift reference positions to query (fiber) positions
    fiber.lift_reference_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
 

for fiber in fiberbam.center("chr20", start=20_000_000, end=20_000_001, strand="-"):
    # returns the same fiber object as above; however, all the positions have been modified to be relative to the region fetched
    print(fiber.msp.reference_starts)