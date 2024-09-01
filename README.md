# Phidl
GitHub for Phidl code, explanations, designs, and more.

Used and developed by Oliver Liebe while at the Center for Quantum Devices, at the Niels Bohr Institute of the University of Copenhagen. 

To get started, simply download the various dependencies, which are: phidl, pyyaml, and json.

After having the prerequisite packages, go into the build.ipynb Jupyter Notebooks file and run the first cell. This will create three different folders, containing the design, the marker SCON files (For Elionix marker alignment), and JSON files also containing the various marker coordinates, but more human readable. 

The first cell of build.ipynb have two main functions: build_4_chip() and build_single_test(). 

build_4_chip():
  To build the full AXL chip that was used for Oliver Liebe's Master's thesis work, with all the various routings and device specific features, such as the 4 different devices in a 2x2 grid, including test bond pads and the likes.

build_single_test():
  Builds a single device, primarily used to testing new designs. Also the function you want to use, if you are simply interested in the routing of a single device. 
