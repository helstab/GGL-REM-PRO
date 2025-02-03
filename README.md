# GGL-REM-PRO
Geomorphic Grade Line Relative Elevation Model Python Toolbox for ArcGIS Pro

Improvements to code were done using ArcPro 3.4.

The [read me](https://github.com/helstab/GGLREM) for the ArcMap tool has a walkthrough which is almost identical to the experience in ArcPro.

The ArcPro version of the GGL-REM-PRO toolbox is compatible with model builder so you can automate these tools using model builder, an example model is shown below, note tool 1, create the centreline Feature Class, has already been run and the user will have created a valley centreline and added a route ID.

![image](https://github.com/user-attachments/assets/d90a68a2-fd28-4bff-bc18-a3c2af75f8ec)

### WARNINGS:

[1] Tools are designed to access FeatureClasses in the top level of a file geodatabase. Do not try to build your datasets within _FeatureDatasets_ as the tool will error.
