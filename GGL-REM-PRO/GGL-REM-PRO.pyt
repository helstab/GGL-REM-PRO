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
#
# Updates by DDH (ddh@geodata.soton.ac.uk)
#   24/1/25 - Removed unused modules, variables and corrected toolbox alias, improvements to Tool 1
#-------------------------------------------------------------------------------

#Import Modules
import arcpy
import os
from arcpy.sa import *
import numpy as np

#Define Toolboxs
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "GGL-REM-PRO Toolbox"
        self.alias = "GGLREMPRO"

        # List of tool classes associated with this toolbox
        self.tools = [Centerline, CrossSections, CenterlineStations, REM]


# Create Centerline Feature Class Tool Parameters
# 24/1/25 - Now validates table name against workspace
# 24/1/25 - Removed parameter for project location as it is never used in any of the code!
# 24/1/25 - Created an output parameter and output is passed into that, allows modelbuilder connectivity.
class Centerline(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "1. Create a Centerline Feature Class"
        self.description = "Create a polyline Feature Class in the current workspace with specific Fields and Data Types for the Create Cross Section Tool"
        self.canRunInBackground = False

    def getParameterInfo(self):
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
        
        paramOut = arcpy.Parameter(
            displayName="Output Feature Class",
            name="out_features",
            datatype="GPFeatureLayer",
            parameterType="Derived",
            direction="Output")

        params = [geodatabaseLOC, centerCOORD, nameFC, paramOut]
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
        # Set Local Variables
        gdb = parameters[0].valueAsText
        dem = parameters[1].valueAsText
        name_fc = parameters[2].valueAsText
        cl_name = "Centerline_" + name_fc
        
        # Set Environment Properties
        arcpy.env.overwriteOutput = True
        arcpy.env.addOutputsToMap = False
        arcpy.env.workspace = gdb
        
        # Validate table name
        cl_name2 = arcpy.ValidateTableName(cl_name, gdb)
        if cl_name != cl_name2:
             arcpy.AddWarning("Supplied table name was invalid and has been corrected to: " + cl_name2)

        # Create Feature Class
        arcpy.AddMessage("... Creating Feature Class")
        arcpy.CreateFeatureclass_management(gdb, cl_name2, "POLYLINE", "", "", "", dem)

        # Add Route ID Field to Feature Class
        arcpy.AddMessage("... Adding Route ID Field")
        arcpy.AddField_management(cl_name2, "ROUTEID", "TEXT", None, None, 50)
        
        # Add Layers
        arcpy.AddMessage("... Adding Centerline Layer to map")
        arcpy.env.addOutputsToMap=True
        fc_cl_path = os.path.join(gdb, cl_name2)
        parameters[3].value = fc_cl_path

        arcpy.AddMessage("Now edit the centreline layer and create your valley centreline, don't forget to add a route ID!")
        return


# Create Cross Section Tool Parameters
# 24/1/25 - Created two output parameters and dropped requirement for map object, allows modelbuilder connectivity.
# 24/1/25 - Now uses AddFields to bulk create fields
# 24/1/25 - Better use of cursors
# 24/1/25 - Better filtering of parameters
# 27/1/25 - Parameter zero checks that input layer has only 1 feature in it, otherwise it rejects the layer.
# 27/1/25 - Loads route ID as a single value as layer will always have a single feature
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
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        inFC.filter.list = ["Polyline"]

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

        # Output routes [4]
        paramOut1 = arcpy.Parameter(
            displayName="Output Route Feature Class",
            name="out_Route",
            datatype="GPFeatureLayer",
            parameterType="Derived",
            direction="Output")
        
        # Output cross sections [5]
        paramOut2 = arcpy.Parameter(
            displayName="Output Cross section Feature Class",
            name="out_XS",
            datatype="GPFeatureLayer",
            parameterType="Derived",
            direction="Output")
        params = [inFC, routeID, stationDirection, widthVB, paramOut1, paramOut2]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[0].altered:
            resObj = arcpy.GetCount_management(parameters[0].valueAsText)
            n = int(resObj.getOutput(0))
            if n == 1:
                # As layer only has one feature it will have only one route ID.
                with arcpy.da.SearchCursor(parameters[0].valueAsText, 'ROUTEID') as cursor:
                    for row in cursor:
                        if row[0] not in [None, "", "<Null>"]:
                            parameters[1].value = row[0]
        return

    def updateMessages(self, parameters):
        if parameters[0].altered:
            resObj = arcpy.GetCount_management(parameters[0].valueAsText)
            n = int(resObj.getOutput(0))
            if n == 1:
                parameters[0].clearMessage()
                if parameters[1].altered == False:
                     parameters[1].setWarningMessage("Route ID is invalid, must be a non-null value.")
            else:
                 # Report error
                 parameters[0].setErrorMessage("Centreline layer must contain only 1 polyline with a route ID!")

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
        fc_routed = "Routed_" + route_id
        off_table = "Offset_Table_" + route_id
        merged = "Merged_" + route_id
        x_sec = "CrossSections_" + route_id
        desc = arcpy.Describe(fc_in)
        gdb = desc.path

        # Set Environment
        arcpy.env.overwriteOutput = True
        arcpy.env.addOutputsToMap = False
        arcpy.env.workspace = gdb

        # Process: Create Routes
        arcpy.AddMessage("... Creating Routes")
        arcpy.CreateRoutes_lr(fc_in, route_field, fc_routed, "LENGTH", "", "", draw_dir, "", "", "IGNORE", "INDEX")

        # Create Table
        arcpy.AddMessage("... Building Offset Table")
        arcpy.CreateTable_management(gdb, off_table)

        # Add Fields
        arcpy.AddMessage("... Adding fields")
        arcpy.AddFields_management(off_table, [["LOCATION", "LONG"],["OFFSET_LEFT", "LONG"],["OFFSET_RIGHT", "LONG"],[route_field, "TEXT", "", 50]])

        #Extract values from Centerline Polyline and create variable with desired row length
        arcpy.AddMessage("... Extracting Values")
        # round length to integer
        with arcpy.da.SearchCursor(fc_routed, "SHAPE@LENGTH") as cursor:
             for row in cursor:
                  LENGTH = int(row[0])
                  
        # Append Extracted Values to Offset_Table
        arcpy.AddMessage("... Populating Offset Table")
        fields = ["LOCATION", "OFFSET_LEFT", "OFFSET_RIGHT", route_field]
        with arcpy.da.InsertCursor(off_table, fields) as cursor:
            for x in range(1, LENGTH):
                cursor.insertRow((x, o_left, o_right, route_id))

        #Process: Make Route Event Layers Left and Right
        arcpy.AddMessage("... Creating Offset Stations")
        arcpy.MakeRouteEventLayer_lr(fc_routed,"ROUTEID",off_table,"ROUTEID POINT LOCATION", "leftoff", "OFFSET_LEFT","NO_ERROR_FIELD","NO_ANGLE_FIELD","NORMAL","ANGLE","LEFT","POINT")
        arcpy.MakeRouteEventLayer_lr(fc_routed,"ROUTEID",off_table,"ROUTEID POINT LOCATION", "rightoff", "OFFSET_RIGHT","NO_ERROR_FIELD","NO_ANGLE_FIELD","NORMAL","ANGLE","RIGHT","POINT")

        #Merge Offset Routes
        arcpy.AddMessage("... Merging Offsets")
        arcpy.management.Merge(["leftoff","rightoff"], merged, "")

        #Convert Points to Lines
        arcpy.AddMessage("... Converting Offset Points to Cross Sections")
        arcpy.PointsToLine_management(merged, x_sec, "LOCATION", "LOCATION")

        #Delete Temporary Features
        arcpy.AddMessage("... Deleting temporary datasets")
        arcpy.Delete_management(["merged", "leftoff", "rightoff"])

        #Add Layers to Map
        arcpy.env.addOutputsToMap = True
        fc_routed_path = os.path.join(gdb, fc_routed)
        parameters[4].value = fc_routed_path
        fc_x_sec_path = os.path.join(gdb, x_sec)
        parameters[5].value = fc_x_sec_path
        return

# Create Centerline Stations Tool Parameters
# 24/1/25 - Created output parameter and dropped requirement for map object, allows modelbuilder connectivity.
# 24/1/25 - Better filtering of parameters
# 27/1/25 - Build GGL now steps through cursor once
# 27/1/25 - Join Field tool now builds indices to improve join performance
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
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        inFC1.filter.list = ["Polyline"]

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
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        inFC2.filter.list = ["Polyline"]

        #Fourth parameter [3]
        inDEM = arcpy.Parameter(
            displayName = "Input LiDAR Digital Elevation Model",
            name = "Input DEM",
            datatype = ["GPRasterLayer", "DEMosaicDataset", "GPMosaicLayer", "DERasterDataset", "GPRasterDataLayer"],
            parameterType = "Required",
            direction = "Input")

        # Output routes [4]
        paramOut1 = arcpy.Parameter(
            displayName="Output Station Feature Class",
            name="out_Stations",
            datatype="GPFeatureLayer",
            parameterType="Derived",
            direction="Output")
        
        params = [inFC1, routeID, inFC2, inDEM, paramOut1]
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

        # Set Environment
        arcpy.env.overwriteOutput = True
        arcpy.env.addOutputsToMap = False
        arcpy.env.workspace = gdb

        # Create cross station points from interscetion of route with cross sections
        arcpy.AddMessage("... Intersecting Centerline and Cross Section Polylines")
        arcpy.analysis.PairwiseIntersect(inFeatures, "xsec", "", "", "POINT")
        arcpy.AddMessage("... Exploding into single part geometry")
        arcpy.MultipartToSinglepart_management("xsec", "xsec2")

        #Extract elevation data from DEM to Centerline Station Points
        arcpy.AddMessage("... Extracting Elevation Values")
        arcpy.sa.ExtractValuesToPoints("xsec2", raster, stations, "INTERPOLATE")
        arcpy.Delete_management("xsec")
        arcpy.Delete_management("xsec2")

        arcpy.AddMessage("... Building GGL")
        px = list()
        py = list()
        with arcpy.da.SearchCursor(stations, ["LOCATION", "RASTERVALU"]) as cursor:
             for row in cursor:
                  px.append(row[0])
                  py.append(row[1])

        #Linear Model
        arcpy.AddMessage("... LINEAR")
        polyfit_1 = np.polyfit(px, py, 1)
        p1 = np.polyval(polyfit_1, px)

        #Second Order
        arcpy.AddMessage("... QUADRATIC")
        polyfit_2 = np.polyfit(px, py, 2)
        p2 = np.polyval(polyfit_2, px)

        #Third Order
        arcpy.AddMessage("... THIRD ORDER POLY")
        polyfit_3 = np.polyfit(px, py, 3)
        p3 = np.polyval(polyfit_3, px)

        #Fourth Order
        arcpy.AddMessage("... FOURTH ORDER POLY")
        polyfit_4= np.polyfit(px, py, 4)
        p4 = np.polyval(polyfit_4, px)

        #Fifth Order
        arcpy.AddMessage("... FIFTH ORDER POLY")
        polyfit_5= np.polyfit(px, py, 5)
        p5 = np.polyval(polyfit_5, px)

        #Build Structured Array
        ##Set Data Types
        arcpy.AddMessage("... Building numpy array")
        dt = {'names':['LOCATION','LIDAR', 'LINEAR', 'POLY2', 'POLY3', 'POLY4', 'POLY5'], 'formats':[np.int64, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32]}

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
        arcpy.AddMessage("... Joining Modeled Values to Features")
        arcpy.management.JoinField(stations, "LOCATION", out_table, "LOCATION", "LIDAR;LINEAR;POLY2;POLY3;POLY4;POLY5", "NOT_USE_FM", None, "NEW_INDEXES")
        arcpy.management.JoinField(crosssection, "LOCATION", out_table, "LOCATION", "LIDAR;LINEAR;POLY2;POLY3;POLY4;POLY5", "NOT_USE_FM", None, "NEW_INDEXES")

        #Add Layers
        arcpy.env.addOutputsToMap = True
        arcpy.AddMessage("... Adding Stations to TOC")
        fc_stations_path = os.path.join(gdb, stations)
        parameters[4].value = fc_stations_path
        return

# 24/1/25 - Input parameters for cross section now filters for polyline, DEM filters for a wider range of raster formats
# 24/1/25 - Better filtering of parameters
# 24/1/25 - Created output parameter and dropped requirement for map object, allows modelbuilder connectivity.
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
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        #Second parameter
        inFC = arcpy.Parameter(
            displayName = "Input Cross Section Feature Class",
            name = "InputCrossSections",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        inFC.filter.list = ["Polyline"]

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
            datatype = ["GPRasterLayer", "DEMosaicDataset", "GPMosaicLayer", "DERasterDataset", "GPRasterDataLayer"],
            parameterType = "Required",
            direction = "Input")

        # Output routes [4]
        paramOut1 = arcpy.Parameter(
            displayName="Output REM",
            name="out_REM",
            datatype="GPRasterLayer",
            parameterType="Derived",
            direction="Output")
        
        params = [inNAME, inFC, gglLIST, gglCUST, gglFIELD, inDEM, paramOut1]
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
        desc = arcpy.Describe(crosssections)
        gdb = desc.path

        # set environment
        arcpy.env.addOutputsToMap = False
        arcpy.env.workspace = gdb
        arcpy.env.overwriteOutput = True

        #REM in Float Meters
        if "Custom" in detrend:
                arcpy.AddMessage("Building Cutom GGLREM")
                arcpy.CopyRows_management(ggl_table, "ggl_table_custom")
                arcpy.JoinField_management(crosssections, "LOCATION", "ggl_table_custom", "LOCATION", ggl_field)
                arcpy.PolylineToRaster_conversion(crosssections, ggl_field, gglrem_name + "_GGL_RASTER_Custom", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Custom", gglrem_name + "_GGL_REM_Custom")
                fc_Custom_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Custom_m")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_Custom_Float_m_path

        if "Linear Model" in detrend:
                arcpy.AddMessage("Building Linear GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "LINEAR", gglrem_name + "_GGL_RASTER_Linear", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Linear", gglrem_name + "_GGL_REM_Linear")
                fc_LINEAR_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Linear")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_LINEAR_path

        if "Polynomial 2nd" in detrend:
                arcpy.AddMessage("Building Quadratic GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY2", gglrem_name + "_GGL_RASTER_Poly2", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly2", gglrem_name + "_GGL_REM_Poly2")
                fc_Poly2_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly2")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_Poly2_Float_m_path

        if "Polynomial 3rd" in detrend:
                arcpy.AddMessage("Building 3rd Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY3", gglrem_name + "_GGL_RASTER_Poly3", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly3", gglrem_name + "_GGL_REM_Poly3")
                fc_Poly3_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly3")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_Poly3_Float_m_path

        if "Polynomial 4th" in detrend:
                arcpy.AddMessage("Building 4th Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY4", gglrem_name + "_GGL_RASTER_Poly4", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly4", gglrem_name + "_GGL_REM_Poly4")
                fc_Poly4_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly4")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_Poly4_Float_m_path

        if "Polynomial 5th" in detrend:
                arcpy.AddMessage("Building 5th Order Poly GGLREM")
                arcpy.PolylineToRaster_conversion(crosssections, "POLY5", gglrem_name + "_GGL_RASTER_Poly5", "", "", "1")
                arcpy.Minus_3d(lidar, gglrem_name + "_GGL_RASTER_Poly5", gglrem_name + "_GGL_REM_Poly5")
                fc_Poly5_Float_m_path = os.path.join(gdb, gglrem_name + "_GGL_REM_Poly5")
                arcpy.env.addOutputsToMap = True
                parameters[6].value = fc_Poly5_Float_m_path

        arcpy.AddMessage("KEEP ASKING QUESTIONS!")
        return
