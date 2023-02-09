# ALENA
v1.0.0.0
Aircraft Loads Environment for Nonlinear Aeroelastics

A MATLAB-based framework for analysing aircraft using a variety of aeroelastic codes, including Nastran. 

## Getting Started

 Clone the repository using `git clone https://github.com/farg-bristol/ALENA.git`.

 Open MATLAB and navigate to the folder where you have cloned the repository. Run the function `addsandbox` to add the package folders to the MATLAB path. 

 To run the tests type `runtests('TestOMEGA')` in the MATLAB Command Prompt and press `<Enter>`. This will execute all the tests in the `TestOMEGA` class. TODO - Many more tests to be added!

 To run ALENA via the GUI, type `ALENA` the MATLAB Command Prompt and press `<Enter>`. You should see the ALENA splash page appear and the GUI will load up shortly afterwards.

 At the end of your development session run `rmsandbox` to remove the package folders from the MATLAB path.

## Dependencies

The ALENA GUI uses the [GUI Layout Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox) (GLT). 
You must install the GLT if you want to run ALENA via the GUI. 

## Change Log
All changes between versions are noted in the file "changelog.txt"
