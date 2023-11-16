#-------------------------------------------------------------------------------
# WORKING UPDATE TO ARCPRO
# ArcPro version: 2.9.2
# Python version: 3.7.11
# GGLREM version: 3.1
# First release: 08/21/2023
# Latest release: 2022

# Name:        GGL REM Toolbox
# Purpose: Series of tools to build a Relative Elevation Model (REM)
# based on the Geomorphic Grade Line (GGL).
# Author:      Matt Helstab
#
# Copyright:   (c) jmhelstab 2018
# Licence:     GNU General Public License v3.0
#-------------------------------------------------------------------------------

#Import Modules
import arcpy
import os
import sys
from arcpy.sa import *
import numpy as np
import pandas as pd

#Define Toolboxs
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "GGL-REM-PRO Toolbox"
        self.alias = "GGL-REM-PRO"

        # List of tool classes associated with this toolbox
        self.tools = [Centerline, CrossSections, CenterlineStations, REM]


#Create Centerline Feature Class Tool Parameters
class Centerline(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "1. Create a Centerline Feature Class"
        self.description = "Create a polyline Feature Class in the current workspace with specific Fields and Data Types for the Create Cross Section Tool"
        self.canRunInBackground = False

    def getParameterInfo(self):
        workspaceLOC = arcpy.Parameter(
            displayName = "Input Workspace Location",
            name = "Workspace Location",
            datatype = "DEWorkspace",
            parameterType = "Required",
            direction = "Input")
        workspaceLOC.filter.list = ["File System"]

        geodatabaseLOC = arcpy.Parameter(
            displayName = "Input Project Geodatabase",
            name = "GDB Location",
            datatype = "DEWorkspace",
            parameterType = "Required",
            direction = "Input")
        geodatabaseLOC.filter.list = ["Local Database", "Remote Database"]

        centerCOORD = arcpy.Parameter(
            displayName = "Match Coordinate System to LiDAR DEM",
            name = "Coordinate System",
            datatype = "GPSpatialReference",
            parameterType = "Required",
            direction = "Input")

        nameFC = arcpy.Parameter(
            displayName = "Input Centerline Feature Class Name",
            name = "Centerline Name",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Output")

        params = [workspaceLOC, geodatabaseLOC, centerCOORD, nameFC]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        #parameters[0].value = os.path.dirname(arcpy.env.workspace)
        #parameters[1].value = arcpy.env.workspace
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        #Set Local Variables
        ws = parameters[0].valueAsText
        gdb = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        name_fc = parameters[3].valueAsText
        cl_name = "Centerline_" + name_fc
        
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        aprxMap = aprx.listMaps("Map")[0] 
        
        #Create Feature Class
        arcpy.AddMessage("Creating Feature Class")
        arcpy.CreateFeatureclass_management(gdb, cl_name, "POLYLINE", "", "", "", dem)

        #Add Route ID Field to Feature Class
        arcpy.AddMessage("Adding Route ID Field")
        arcpy.AddField_management(cl_name, "ROUTEID", "TEXT")
        
        #Add Layers
        arcpy.AddMessage("Adding Centerline Layer")
        fc_cl_path = os.path.join(gdb, cl_name)
        aprxMap.addDataFromPath(fc_cl_path)
        return


#Create Cross Section Tool Parameters
class CrossSections(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "2. Create Cross Sections and Routed Centerline"
        self.description = "Creates cross section polylines and routed centerline."
        self.canRunInBackground = False

    def getParameterInfo(self):
        #First parameter [0]
        inFC = arcpy.Parameter(
            displayName = "Input Centerline Feature Class",
            name = "Input Centerline",
            datatype = ["GPFeatureLayer", "DEFeatureClass", "DEShapefile"],
            parameterType = "Required",
            direction = "Input")

        #Second parameter [1]
        routeID = arcpy.Parameter(
            displayName = "Select Centerline Route ID",
            name = "Route ID",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        routeID.filter.type = "ValueList"
        routeID.filter.list = []

        #Third parameter [2]
        stationDirection = arcpy.Parameter(
            displayName = "Select Direction to Start Stationing From",
            name = "Station Direction",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        stationDirection.filter.type = "ValueList"
        stationDirection.filter.list = ["UPPER_LEFT", "UPPER_RIGHT", "LOWER_LEFT", "LOWER_RIGHT"]
        
        #Fourth parameter [3]
        widthVB = arcpy.Parameter(
            displayName = "Input Maximum Valley Bottom Width (m)",
            name = "Valley Bottom Width",
            datatype = "GPLong",
            parameterType = "Required",
            direction = "Input")

        

        params = [inFC, routeID, stationDirection, widthVB]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[0].value:
            with arcpy.da.SearchCursor(parameters[0].valueAsText, 'ROUTEID') as rows:
                parameters[1].filter.list = sorted(list(set([row[0] for row in rows])))
        else:
            parameters[1].filter.list = []
        return

    def updateMessages(self, parameters):
        if parameters[3].altered:
            if parameters[3].value <= 0:
                parameters[3].setErrorMessage('''Offset Value must be greater than zero.''')

        return

    def execute(self, parameters, messages):

        #Get Parameter Inputs
        fc_in = parameters[0].valueAsText
        route_id = parameters[1].valueAsText
        route_field = "ROUTEID"
        draw_dir = parameters[2].valueAsText
        o_left = parameters[3].value/2
        o_right = parameters[3].value/2
        length_id = "LOCATION"
        fc_routed = "Routed_" + route_id
        off_table = "Offset_Table_" + route_id
        merged = "Merged_" + route_id
        x_sec = "CrossSections_" + route_id
        desc = arcpy.Describe(fc_in)
        gdb = desc.path

        #Set Workspace Environment and Map Properties
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        aprxMap = aprx.listMaps("Map")[0] 

        # Process: Create Routes
        arcpy.AddMessage("Creating Routes")
        arcpy.CreateRoutes_lr(fc_in, route_field, fc_routed, "LENGTH", "", "", draw_dir, "", "", "IGNORE", "INDEX")

        #Create Table
        arcpy.AddMessage("Building Offset Table")
        arcpy.CreateTable_management(gdb, off_table)

        #Add Fields
        arcpy.AddField_management(off_table, "LOCATION", "LONG")
        arcpy.AddField_management(off_table, "OFFSET_LEFT", "LONG")
        arcpy.AddField_management(off_table, "OFFSET_RIGHT", "LONG")
        arcpy.AddField_management(off_table, route_field, "TEXT")


        #Extract values from Centerline Polyline and create variable with desired row length
        arcpy.AddMessage("Extracting Values")
        fields_centerline = ['shape_length', 'shape_Length', 'shape_LENGTH', 'Shape_length', 'Shape_Length', 'Shape_LENGTH', 'SHAPE_length', 'SHAPE_Length', 'SHAPE_LENGTH',]
        LOCATION1 = arcpy.da.SearchCursor(fc_routed, fields_centerline,).next()[0]
        LENGTH = int(LOCATION1)
        LOCATION2 = range(1,LENGTH)

        NAME =  arcpy.da.SearchCursor(fc_in, route_field,)
        NAME = [NAME] * LENGTH

        #Append Extracted Values to Offset_Table
        arcpy.AddMessage("Populating Offset Table")
        fields = ["LOCATION", "OFFSET_LEFT", "OFFSET_RIGHT", route_field]
        cursor = arcpy.da.InsertCursor(off_table, fields)
        for x in range(1, LENGTH):
            cursor.insertRow((x, o_left, o_right, route_id, ))
        del(cursor)

        #Process: Make Route Event Layers Left and Right
        arcpy.AddMessage("Creating Offset Stations")
        arcpy.MakeRouteEventLayer_lr(fc_routed,"ROUTEID",off_table,"ROUTEID POINT LOCATION", "leftoff", "OFFSET_LEFT","NO_ERROR_FIELD","NO_ANGLE_FIELD","NORMAL","ANGLE","LEFT","POINT")
        arcpy.MakeRouteEventLayer_lr(fc_routed,"ROUTEID",off_table,"ROUTEID POINT LOCATION", "rightoff", "OFFSET_RIGHT","NO_ERROR_FIELD","NO_ANGLE_FIELD","NORMAL","ANGLE","RIGHT","POINT")

        #Merge Offset Routes
        arcpy.AddMessage("Merging Offsets")
        arcpy.management.Merge(["leftoff","rightoff"], merged, "")
        #arcpy.management.CopyFeatures("leftoff", "rightoff")
        #Convert Points to Lines
        arcpy.AddMessage("Converting Offset Points to Cross Sections")
        arcpy.PointsToLine_management(merged, x_sec, "LOCATION", "LOCATION")

        #Add Layers to Map
        fc_routed_path = os.path.join(gdb, fc_routed)
        aprxMap.addDataFromPath(fc_routed_path)

        fc_x_sec_path = os.path.join(gdb, x_sec)
        aprxMap.addDataFromPath(fc_x_sec_path)

        #Delete Temporary Features
        arcpy.Delete_management(["merged", "leftoff", "rightoff"])
        return

# Create Centerline Stations Tool Parameters
class CenterlineStations(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "3. Create GGL Table and Centerline Stations"
        self.description = "Creates a Point Feature Class at each intersection of the Centerline and Cross Section polylines and appends elevation data to each point."
        self.canRunInBackground = False

    def getParameterInfo(self):
        #First parameter [0]
        inFC1 = arcpy.Parameter(
            displayName = "Input Routed Centerline Feature Class",
            name = "Input Routed Centerline",
            datatype = ["GPFeatureLayer", "DEFeatureClass"],
            parameterType = "Required",
            direction = "Input")

        #Second parameter [1]
        routeID = arcpy.Parameter(
            displayName = "Select Centerline Route ID",
            name = "Route ID",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        routeID.filter.type = "ValueList"
        routeID.filter.list = []

        #Third parameter [2]
        inFC2 = arcpy.Parameter(
            displayName = "Input Cross Section Feature Class",
            name = "Input Crosssection",
            datatype = ["GPFeatureLayer", "DEFeatureClass"],
            parameterType = "Required",
            direction = "Input")

        #Fourth parameter [3]
        inDEM = arcpy.Parameter(
            displayName = "Input LiDAR Digital Elevation Model",
            name = "Input DEM",
            datatype = ["GPRasterLayer", "GPLayer", "DEMosaicDataset", "GPMosaicLayer", "DERasterDataset", "GPRasterDataLayer"],
            parameterType = "Required",
            direction = "Input")

        params = [inFC1, routeID, inFC2, inDEM]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        if parameters[0].value:
            with arcpy.da.SearchCursor(parameters[0].valueAsText, "ROUTEID") as rows:
                parameters[1].filter.list = [row[0] for row in rows]
        else:
            parameters[1].filter.list = []
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):

        #Set Local Variables
        centerroute = parameters[0].valueAsText
        routeid = parameters[1].valueAsText
        crosssection = parameters[2].valueAsText
        raster = parameters[3].valueAsText
        table_name = "GGL_Table_" + routeid
        stations = "Stations_" + routeid

        inFeatures = [centerroute, crosssection]
        desc = arcpy.Describe(centerroute)
        gdb = desc.path
        out_table = gdb + "/" + "GGL_Table_" + routeid
        dir_loc = os.path.dirname(gdb)
        csv = table_name + ".csv"

        arcpy.env.overwriteOutput = True

        #Set Workspace Environment and Map Properties
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        aprxMap = aprx.listMaps("Map")[0] 

        arcpy.AddMessage("Intersecting Centerline and Cross Section Polylines...")
        arcpy.analysis.PairwiseIntersect(inFeatures, "xsec", "", "", "POINT")
        
        arcpy.MultipartToSinglepart_management("xsec", "xsec2")

        #Extract elevation data from DEM to Centerline Station Points
        arcpy.AddMessage("Extracting Elevation Values...")
        arcpy.sa.ExtractValuesToPoints("xsec2", raster, stations, "INTERPOLATE")
        arcpy.Delete_management("xsec")
        arcpy.Delete_management("xsec2")

        #Add Layers
        arcpy.AddMessage("Adding Stations to TOC...")
        fc_stations_path = os.path.join(gdb, stations)
        aprxMap.addDataFromPath(fc_stations_path)

        arcpy.AddMessage("Building GGL...")
        px = [row[0] for row in arcpy.da.SearchCursor(stations, "LOCATION")]
        py = [row[0] for row in arcpy.da.SearchCursor(stations, "RASTERVALU")]

        #Linear Model
        arcpy.AddMessage("LINEAR")
        polyfit_1 = np.polyfit(px, py, 1)
        p1 = np.polyval(polyfit_1, px)

        #Second Order
        arcpy.AddMessage("QUADRATIC")
        polyfit_2 = np.polyfit(px, py, 2)
        p2 = np.polyval(polyfit_2, px)

        #Third Order
        arcpy.AddMessage("THIRD ORDER POLY")
        polyfit_3 = np.polyfit(px, py, 3)
        p3 = np.polyval(polyfit_3, px)

        #Fourth Order
        arcpy.AddMessage("FOURTH ORDER POLY")
        polyfit_4= np.polyfit(px, py, 4)
        p4 = np.polyval(polyfit_4, px)

        #Fifth Order
        arcpy.AddMessage("FIFTH ORDER POLY")
        polyfit_5= np.polyfit(px, py, 5)
        p5 = np.polyval(polyfit_5, px)

        #Build Structured Array
        ##Set Data Types
        arcpy.AddMessage("Almost Done...")
        dt = {'names':['LOCATION','LIDAR', 'LINEAR', 'POLY2', 'POLY3', 'POLY4', 'POLY5'], 'formats':[np.int, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32]}

        ##Build Blank Structured Array
        poly = np.zeros(len(px), dtype=dt)

        ##Add values to Structured Array
        poly['LOCATION'] = px
        poly['LIDAR'] = py
        poly['LINEAR'] = p1
        poly['POLY2'] = p2
        poly['POLY3'] = p3
        poly['POLY4'] = p4
        poly['POLY5'] = p5

        #Convert Structured Array to Table
        if arcpy.Exists(out_table):
            arcpy.Delete_management(out_table)

        arcpy.da.NumPyArrayToTable(poly, out_table)

        arcpy.TableToTable_conversion(out_table, dir_loc, csv)

        #Join Model Output to Cross Sections and Centerline Stations Feature Classes
        arcpy.AddMessage("Joining Modeled Values to Features")

        arcpy.management.JoinField(stations, "LOCATION", out_table, "LOCATION", "LIDAR;LINEAR;POLY2;POLY3;POLY4;POLY5", "NOT_USE_FM", None)
        arcpy.management.JoinField(crosssection, "LOCATION", out_table, "LOCATION", "LIDAR;LINEAR;POLY2;POLY3;POLY4;POLY5", "NOT_USE_FM", None)

        return

class REM(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "4. Create Relative Elevation Model(s)"
        self.description = "Joins modeled GGL elevations to cross sections, converts cross sections to a raster, and lastly subtracts cross section raster from LiDAR raster to produce the REM."
        self.canRunInBackground = False

    def getParameterInfo(self):
        #First parameter
        inNAME = arcpy.Parameter(
            displayName = "Input Unique GGLREM Name",
            name = "GGLREMname",
            datatype = ["GPString"],
            parameterType = "Required",
            direction = "Input")

        #Second parameter
        inFC = arcpy.Parameter(
            displayName = "Input Cross Section Feature Class",
            name = "InputCrossSections",
            datatype = ["GPFeatureLayer", "DEFeatureClass"],
            parameterType = "Required",
            direction = "Input")

        #Third parameter
        gglLIST = arcpy.Parameter(
            displayName = "Select Values/Model to Construct Relative Eleavtion Model",
            name = "GglList",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input",
            multiValue = "True")
        gglLIST.filter.type = "ValueList"
        gglLIST.filter.list = ["Custom", "Linear Model", "Polynomial 2nd", "Polynomial 3rd", "Polynomial 4th", "Polynomial 5th"]

        #Fourth parameter
        gglCUST = arcpy.Parameter(
            displayName = "Input Custom GGL Table [ONLY IF RUNNING CUSTOM MODEL]",
            name = "CustomGglTable",
            datatype = "DETable",
            parameterType = "Optional",
            direction = "Input")

        #Fifth parameter
        gglFIELD = arcpy.Parameter(
            displayName = "Select Field with GGL Values for Detrending",
            name = "CustomGglField",
            datatype = "Field",
            parameterType = "Optional",
            direction = "Input")
        gglFIELD.filter.list =[]
        gglFIELD.parameterDependencies = [gglCUST.name]

        #Sixth parameter
        inDEM = arcpy.Parameter(
            displayName = "Input LiDAR DEM",
            name = "InputLidar",
            datatype = ["GPRasterLayer","DERasterDataset", "DERasterCatalog"],
            parameterType = "Required",
            direction = "Input")

        params = [inNAME, inFC, gglLIST, gglCUST, gglFIELD, inDEM]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        if parameters[3].value:
            parameters[4].enabled = True
        else:
            parameters[4].enabled = False
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        #Set local variables
        gglrem_name = parameters[0].valueAsText
        crosssections = parameters[1].valueAsText
        detrend = parameters[2].valueAsText
        ggl_table = parameters[3].valueAsText
        ggl_field = parameters[4].valueAsText
        lidar = parameters[5].valueAsText
        #rems = parameters[6].valueAsText
        desc = arcpy.Describe(crosssections)
        gdb = desc.path

        #Set Workspace Environment and Map Properties
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        aprxMap = aprx.listMaps("Map")[0] 

        #REM in Float Meters
        if "Custom" in detrend:
                arcpy.AddMessage("Building Cutom GGLREM")
                arcpy.CopyRows_management(ggl_table, "ggl_table_custom")
                arcpy.JoinField_management(crosssections, "LOCATION", "ggl_table_custom", "LOCATION", ggl_field)
                arcpy.PolylineToRaster_conversion(crosssections, ggl_field, gglrem_name + "_GGL_RASTER_Custom", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Custom", gglrem_name + "_GGL_REM_Custom")
                fc_Custom_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Custom_m")
                aprxMap.addDataFromPath(fc_Custom_Float_m_path)

        if "Linear Model" in detrend:
                arcpy.AddMessage("Building Linear GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "LINEAR", gglrem_name + "_GGL_RASTER_Linear", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Linear", gglrem_name + "_GGL_REM_Linear")
                fc_LINEAR_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Linear")
                aprxMap.addDataFromPath(fc_LINEAR_path)

        if "Polynomial 2nd" in detrend:
                arcpy.AddMessage("Building Quadratic GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY2", gglrem_name + "_GGL_RASTER_Poly2", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly2", gglrem_name + "_GGL_REM_Poly2")
                fc_Poly2_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly2")
                aprxMap.addDataFromPath(fc_Poly2_Float_m_path)

        if "Polynomial 3rd" in detrend:
                arcpy.AddMessage("Building 3rd Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY3", gglrem_name + "_GGL_RASTER_Poly3", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly3", gglrem_name + "_GGL_REM_Poly3")
                fc_Poly3_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly3")
                aprxMap.addDataFromPath(fc_Poly3_Float_m_path)

        if "Polynomial 4th" in detrend:
                arcpy.AddMessage("Building 4th Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY4", gglrem_name + "_GGL_RASTER_Poly4", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly4", gglrem_name + "_GGL_REM_Poly4")
                fc_Poly4_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly4")
                aprxMap.addDataFromPath(fc_Poly4_Float_m_path)

        if "Polynomial 5th" in detrend:
                arcpy.AddMessage("Building 5th Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY5", gglrem_name + "_GGL_RASTER_Poly5", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly5", gglrem_name + "_GGL_REM_Poly5")
                fc_Poly5_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly5")
                aprxMap.addDataFromPath(fc_Poly5_Float_m_path)

        arcpy.AddMessage("KEEP ASKING QUESTIONS!")

        return   
