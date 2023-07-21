import pyft
import tqdm

print(pyft.__all__)
bam_f = "../tests/data/center.bam"
bam_f = "chr20.bam"
fiberdata = pyft.FiberdataFetch(bam_f, "chr20", 0, 1_000_000)
for idx, fiber in enumerate(tqdm.tqdm(fiberdata)):
    if idx < 10:
        print(fiber)
    fiber.get_aligned_ranges()
