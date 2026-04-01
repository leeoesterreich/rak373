Before running the OM2BFB make sure the following:

# Create & initializing the conda virtual environment
# Note that you must use java/11.0.2-openjdk module. Other java modules do not work for BFBFinder.jar
# Make ensure BFB_vis.py, analyzing_BFBFinder.py, parsers.py, segmentation.R, is in the same directory as your script runOM2BFB.py.
# Please find all these files at https://github.com/siavashre/OM2BFB

$ module load python/ondemand-jupyter-python3.10
$ conda create --name BFBCycle python=3.7.4
$ source activate /ihome/alee/rak373/.conda/envs/BFBCycle
$ conda install numpy
$ conda install matplotlib
$ conda install pandas

# Just ensure the output directory exists before you list or write to it.
# In runOM2BFB.py, near line 793: add 
os.makedirs(args.output, exist_ok=True)
for f in os.listdir(args.output):
    if os.path.isfile(os.path.join(args.output, f)):
        os.remove(os.path.join(args.output, f))
