import pyft
import tqdm
import sys

bam_f = "../tmp.bam"
fiberdata = pyft.FiberdataFetch(bam_f, "chr20", 0, 10_000_000)
for idx, fiber in enumerate(tqdm.tqdm(fiberdata)):
    if idx < 10:
        print(fiber, file=sys.stderr)
    
    fiber.get_aligned_blocks_as_ranges()
    fiber.lift_query_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    fiber.lift_reference_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])


