These different files of code serve a purpose detailed in the comments.
In general, I used the ROOT package for Python to extract and analyse data, as CERN stores its data from detectors and simulations in .root files to read only the necessary data very efficiently. This is because each .root file can carry information on billions of particles, and have multiple pieces of data on each one, such as transverse momentum. They store data in a structure called a TTree. You may only want to use one certain piece of data, which is why we use .root files so we can choose the branches we read, instead of wasting computational power reading the whole file.

Testing
