# Calculate Irradiances from Hydrolight Lroot file

This is a Python module to obtain Irradiances (E) by integration of the radiances over the appropiate solid angles from Hydrolight Lroot file.
This module calculate Ed as the downwelling irradiance, Eu as the upwelling irradiance, El as the lateral-welling irradiance and other different type of irradiances

## Instructions

- Copy your Lroot file from Hydrolight inside the files/raw folder
- Results appear in files/csv folder
- Images appear in images folder

You will find a main.py to see how to import this module and use their basic functions:
-   calc_irradiances()
-   plot_irradiances()

## Install Dependencies

- npm install -g electron@1.8.4 orca

In case you find a "ConnectionRefusedError" when you try the fig.write_image() Plotly function, you have to allow configure orca to send requests to remote server with the following command line:

orca serve -p 32909 --plotly



