# import pyft
# from pyft import pyft
import pyft
import tqdm


bam_f = "../tests/data/center.bam"
fiberbam = pyft.Fiberbam(bam_f)
out_fiberbam = pyft.Fiberwriter("test.bam", bam_f)
rgn = ["chr22", 26_354_169, 26_354_170]
for fiber in tqdm.tqdm(fiberbam.fetch(*rgn)):
    # the number of ccs passes
    fiber.ec
    # the mps start positions
    fiber.msp.starts
    # the fire quality scores of the MSPs
    fiber.msp.qual
    # print the nuc reference starts
    fiber.nuc.reference_starts
    # lift query (fiber) positions to reference positions
    fiber.lift_query_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    # lift reference positions to query (fiber) positions
    fiber.lift_reference_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    out_fiberbam.write(fiber)


for fiber in fiberbam.center(rgn[0], start=rgn[1], end=rgn[2], strand="-"):
    # returns the same fiber object as above; however, all the positions have been modified to be relative to the region fetched
    # print(fiber.msp.reference_starts)
    continue


# example of reading in a footprinting table
df = pyft.read_footprint_table("../tests/data/ctcf-footprints.bed.gz", long=True)
print(df)

# read in a footprinting table and center the positions
df = pyft.read_and_center_footprint_table("../tests/data/ctcf-footprints.bed.gz")
print(df)

# read a region of a fiberbam into a pandas dataframe
df = pyft.region_to_df(fiberbam, rgn)
print(df)

# read a region of a fiberbam into a pandas dataframe and center the positions
df = pyft.region_to_centered_df(fiberbam, rgn, strand="-")
print(df)
