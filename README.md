# To download all the necessary dependencies, just run the following commands:
# These will create a python environment, source it, install all needed dependencies, and finally run the program.

python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
python3 build_phylotree.py P68871 
