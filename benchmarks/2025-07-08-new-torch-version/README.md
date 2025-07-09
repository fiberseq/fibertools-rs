I want to upgrate to torch 0.17. On this new branch make the code changes needed to allow for this without messing up the multi threading with rayon.

Then I want you to load the conda enviroment pytorch2.6 and build with --all-features 

save the binary in the the bin subdirectory here with a name indicating the burn version and that we used pytorch

then compile only using burn features and save that binary as well by copying.

If scucessful commit these changes but dont push them.

Then go back to the branch 0.6.5 and repeat this whole process using the conda env pytorch2.2 and saving the binaries


Finally return to this benchmakring branch and use hyperfine along with this file to benchmakr all the different bianries:
/Users/mrvollger/repos/fibertools-rs/tests/data/kinetics.small.test.bam
You may need to make wrapper scripts like ./benchmarks/2025-07-08-new-torch-version/ft_0_17_1_pytorch.sh to have ft find the right pytorch libs.

