# Masters-Research-Project-Part2
## Scaled Model Aeroelastic Similarity of a Blended Wing Body through Multi-Disciplinary Optimization

> **This folder contains the following folders:**

1. **BWB Airfoil Coordinates:** This folder contains airfoil coordinates of BWB based on OpenVSP geometry. It also contains an excel to calculate and extrapolate the airfoil data.
2. **BWB Geometry:** This folder contains OpenVSP geometry and IGES files for Ribs, Spars, and Panels.
3. **BWB Nastran Mesh:** This folder contains the Nastran mesh on the BWB
4. **Extra:** This folder contains the python script for extrapolation of airfoil points.
4. **Presentation & Report:** This folder contains the final presentation and report in both PDF and editable formats.
5. **Python, XML, Airfoils:** This folder contains the Python script for generating wing aerodynamic meshes, XML for wing definition, and Airfoil files. All of them are updated and can generate a wing with any number of sections & any type of airfoils.

> **HOW TO RUN THE CODE:**

1. To run the python script successfully, all the files in "Python, XML, Airfoils" folder should always stay in the same folder. 
2. Open *"Multiple airfoil Wing Mesh.py"* in Spyder (or equivalent) and run the code. Type the name of XML file containing Wing definition data, in this case, named as *"CRM Wing definition - Multiple user defined airfoil.xml"*
3. Press Enter. 

> **Notes**
1. The XML file defines the wing geometry in numbers and the airfoil data. Hence, the name of the file which contains airoil data should be exactly the same in the XML file. 
2. All the airfoil files should be in the same folder as the XML file.
