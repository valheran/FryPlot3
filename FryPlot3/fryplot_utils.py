# -*- coding: utf-8 -*-
"""
/***************************************************************************
 FryPlot_Utilities
                                 A QGIS plugin
 a utility for performing autocorrelation analysis (fry plots) of point data
                              -------------------
        begin                : 2015-02-15
        git sha              : $Format:%H$
        copyright            : (C) 2015 by A Brown
        email                : email
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from __future__ import absolute_import
from builtins import object
from qgis.PyQt.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtGui import QIcon

from qgis.core import *
from qgis.gui import *

import locale
import math
import os.path

def getVectorLayerNames():
    layerMap = QgsProject.instance().mapLayers()
    layerNames = []
    for name, layer in layerMap.items():
        if layer.type() == QgsMapLayer.VectorLayer:
            layerNames.append(str(layer.name()))
    return sorted(layerNames)


def getVectorLayerByName(layerName):
    layerMap = QgsProject.instance().mapLayers()
    for name, layer in layerMap.items():
        if layer.type() == QgsMapLayer.VectorLayer and layer.name() == layerName:
            if layer.isValid():
                return layer
            else:
                return None
                
                
class FryAnalysis(object):

       
    def __init__(self):
        self.inputDict={};
        self.Xo = 0
        self.Yo = 0
    
        self.vecList=[]
        self.lenList=[]
        self.fryListX=[]
        self.fryListY=[]
        
    
    
    def listFromLayer(self, tarlayer):
        pointlist=[]
        self.tarlayer = tarlayer
        iter = self.tarlayer.getFeatures()
        spatialIndex = QgsSpatialIndex()
        maxCox = 0
        maxCoy = 0
        minCox = 99999999999999999999999.999
        minCoy = 99999999999999999999999.999
        
        for elem in iter:
            point = elem.geometry().asPoint()
            pointlist.append(point)
            #keep track of the data extents
            if maxCox < point.x():
                maxCox = point.x()
            if maxCoy < point.y():
                maxCoy = point.y()
            if minCox > point.x():
                minCox = point.x()
            if minCoy > point.y():
                minCoy = point.y()
            spatialIndex.insertFeature(elem)
        
         #find the approximate centre of the data extent
        XCentre = maxCox - ((maxCox-minCox)/2)
        YCentre = maxCoy - ((maxCoy-minCoy)/2)
        #find nearest actual point to approximate centre
        nearest = spatialIndex.nearestNeighbor(QgsPointXY(XCentre, YCentre),1)
        #set up request for Fid of nearest point
        orifeat = QgsFeatureRequest().setFilterFid(nearest[0])
        #go get the actual feature
        ori = self.tarlayer.getFeatures(orifeat)
        #extract the point as an object
        for ele in ori:
            origin = ele.geometry().asPoint()
            #print origin

        #set variables for origin using extracted nearest point
        self.Xo = origin.x()
        self.Yo = origin.y()
        return pointlist
         
    def makeFryplotlists(self, pointlist, pseudo=True):
        
        if pseudo:
            crsString = self.tarlayer.crs().authid()
            uri = f"Point?field=Length:double(10,2)&field=Azi:double(4,1)&crs={crsString}"
            pseudolayer = QgsVectorLayer(uri, "PseudoPlot", "memory")
            pseudopr = pseudolayer.dataProvider()
            fetlist=[]
        #choose a point, and iterate over all other points to obtain fry coords and vectors
        for i in pointlist:
            Xp1 = i.x()
            Yp1 = i.y()
                
            for j in pointlist:
                Xp2 = j.x()
                Yp2 = j.y()
                if (Xp1 != Xp2 and Yp1 != Yp2):
                    Xf = (Xp2 - Xp1)
                    Yf = (Yp2 - Yp1)
                    FryPoint = [Xf, Yf]
                    
                    #output using seperate lists
                    self.fryListX.append(Xf)
                    self.fryListY.append(Yf)
                    
                    
                    #calculate fry vector length and azimuth            
                    sqLength = i.sqrDist(j)
                    vecLength = math.sqrt(sqLength)
                    vecAziSign = i.azimuth(j)
                    if vecAziSign <0:
                        vecAzi= 360 + vecAziSign
                    else:
                        vecAzi=vecAziSign
                        
                    #output vector lengths to output list
                    self.lenList.append(vecLength)
                    #output vector azis to output list
                    self.vecList.append(vecAzi)
                    #add to pseudoplot layer if selected
                    if pseudo:
                        fryfet= QgsFeature()
                        #set attributes
                        fryfet.setAttributes([vecLength, vecAzi])
                        fryfet.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(self.Xo+Xf,self.Yo+Yf)))
                        fetlist.append(fryfet)
                    
                    
        if pseudo:
            pseudopr.addFeatures(fetlist)
            pseudolayer.updateExtents()
            QgsProject.instance().addMapLayer(pseudolayer)
        return self.fryListX, self.fryListY, self.lenList, self.vecList
"""
This was the original code, from early on learning python. Was refined later when used in Network Props plugin, which is used here instead. Kept for checking functionality
class FryAnalysis:

    layer = 0
    outLayer = 0
    inputDict={};
    Xo = 0
    Yo = 0
    
    vecList=[]
    lenList=[]
    fryListX=[]
    fryListY=[]
    #fryDict={};

    #some variables to help debugging
    getDataFlag = "no"
    makePseudoFryplotFlag = "no"
    makeFryplotListsFlag = "no"
    
    def __init__(self, layerName, outLayer):
        self.layerName = layerName
        self.outLayer = outLayer
        
        FryAnalysis.layer = getVectorLayerByName(layerName)
        FryAnalysis.outLayer = outLayer
    
    
    
    
    def getData(self):
        #collects up all the coordinate information in the target layer for use in calculations
        
        
        #FryAnalysis.layer = QgsVectorLayer(self.layerName)  from original 
        iter = FryAnalysis.layer.getFeatures()
        spatialIndex = QgsSpatialIndex()
        ID = 0
        maxCox = 0
        maxCoy = 0
        minCox = 99999999999999999999999.999
        minCoy = 99999999999999999999999.999
    
        for elem in iter:
            point = elem.geometry().asPoint()
            xcoord = point.x()
            ycoord = point.y()
            coords = (xcoord, ycoord)
            
            #print coords
            #print elem.id()
            #keep track of the data extents
            if maxCox < xcoord:
                maxCox = xcoord
            if maxCoy < ycoord:
                maxCoy = ycoord
            if minCox > xcoord:
                minCox = xcoord
            if minCoy > ycoord:
                minCoy = ycoord
                
            spatialIndex.insertFeature(elem)
            FryAnalysis.inputDict[ID] = coords;
            ID = ID + 1
            
        #find the approximate centre of the data extent
        XCentre = maxCox - ((maxCox-minCox)/2)
        YCentre = maxCoy - ((maxCoy-minCoy)/2)
        #find nearest actual point to approximate centre
        nearest = spatialIndex.nearestNeighbor(QgsPointXY(XCentre, YCentre),1)
        #print for debugging
        #print nearest
        #print XCentre
        #print YCentre
        #set up request for Fid of nearest point
        orifeat = QgsFeatureRequest().setFilterFid(nearest[0])
        #go get the actual feature
        ori = FryAnalysis.layer.getFeatures(orifeat)
        #extract the point as an object
        for ele in ori:
            origin = ele.geometry().asPoint()
            #print origin

        #set variables for origin using extracted nearest point
        FryAnalysis.Xo = origin.x()
        FryAnalysis.Yo = origin.y()
        #QMessageBox.about(None, 'message', 'Get Data was run')
        FryAnalysis.getDatFlag = "yes"
        
    def makePseudoFryplot(self):
        ### need to link with an output layer input
        #self.outLayer = outputLayerName
        tempLayer = QgsVectorLayer("Fry_Points?field=FryVecLen:double&field=FryVecAzi:double", "Fry Points", "memory")
        tempLayer.setCrs(FryAnalysis.layer.crs())
        prov = tempLayer.dataProvider()
        tempLayer.updateFields()
        #get layer name from dialog box input
        QgsVectorFileWriter(FryAnalysis.outLayer, "CP1250", prov.fields(),QGis.WKBPoint , tempLayer.crs(), "ESRI Shapefile", )
        #pathstring = FryAnalysis.outLayer
        #name = pathstring[(pathstring.rfind("\\")+1):pathstring.rfind('.')]
        name = os.path.basename(FryAnalysis.outLayer)
        outputLayer = QgsVectorLayer(FryAnalysis.outLayer, name,"ogr")
        outprov = outputLayer.dataProvider()



        ID=0
        fkey = 0
        resdict= {}
        for i in FryAnalysis.inputDict:
            coord= FryAnalysis.inputDict[ID];
            Xp1 = coord[0]
            Yp1 = coord[1]
            ID2 = 0
            
            for j in FryAnalysis.inputDict:
                coord = FryAnalysis.inputDict[ID2];
                Xp2 = coord[0]
                Yp2 = coord[1]
                if (Xp1 != Xp2 and Yp1 != Yp2):  #skip when central point and target point are the same
                   #co ordinate calculation
                    Xf = FryAnalysis.Xo + (Xp2 - Xp1)
                    Yf = FryAnalysis.Yo + (Yp2 - Yp1)
                    FryPoint = [Xf, Yf]
                    resdict[fkey] = FryPoint;
                    #create Qpoints for length and azimuth calc
                    centfet = QgsPointXY(Xp1, Yp1)
                    tarfet = QgsPointXY(Xp2, Yp2)
                    #calculate fry vector length and azimuth            
                    sqLength = centfet.sqrDist(tarfet)
                    vecLength = math.sqrt(sqLength)
                    vecAziSign = centfet.azimuth(tarfet)
                    if vecAziSign <0:
                        vecAzi= 360 + vecAziSign
                    else:
                        vecAzi=vecAziSign
                        
                    #output vector lengths to output list
                    #lenList[fkey] = vecLength
                    #output vector azis to output list
                    #vecList[fkey]= vecAzi
                    #create feature and write to output file
                    fryfet= QgsFeature()
                    #set attributes
                    fryfet.setAttributes([vecLength, vecAzi])
                    fryfet.setGeometry(QgsGeometry.fromPoint(QgsPointXY(Xf,Yf)))
                   
                    outprov.addFeatures([fryfet])
                    ID2 = ID2 +1
                    fkey= fkey + 1
                else:
                   ID2 = ID2 + 1
                
            ID = ID +1
        #debug output contents of results dictionary
        #for keys in resdict:
        #    print(resdict[keys]) 
        #load into canvas
        QgsProject.instance().addMapLayer(outputLayer)
 
 
 
    def makeFryplotlists(self):
       #generate fryspace co-ordinates for plotting on graph, as well as vector length and vector azimuth data
        #Iterator for fryspace Co-ordinates
        ID=0
        fkey = 0
        
        Xo = 0
        Yo = 0
        
        for i in FryAnalysis.inputDict:
            coord= FryAnalysis.inputDict[ID];
            Xp1 = coord[0]
            Yp1 = coord[1]
            ID2 = 0
            
            for j in FryAnalysis.inputDict:
                coord = FryAnalysis.inputDict[ID2];
                Xp2 = coord[0]
                Yp2 = coord[1]
                if (Xp1 != Xp2 and Yp1 != Yp2):
                    Xf = Xo + (Xp2 - Xp1)
                    Yf = Yo + (Yp2 - Yp1)
                    FryPoint = [Xf, Yf]
                    
                    #output using seperate lists
                    FryAnalysis.fryListX.append(Xf)
                    FryAnalysis.fryListY.append(Yf)
                    
                    #output frypoint to output dict this was the original idea
                    #FryAnalysis.fryDict[fkey] = FryPoint;
                    #create Qpoints for length and azimuth calc
                    centfet = QgsPointXY(Xp1, Yp1)
                    tarfet = QgsPointXY(Xp2, Yp2)
                    #calculate fry vector length and azimuth            
                    sqLength = centfet.sqrDist(tarfet)
                    vecLength = math.sqrt(sqLength)
                    vecAziSign = centfet.azimuth(tarfet)
                    if vecAziSign <0:
                        vecAzi= 360 + vecAziSign
                    else:
                        vecAzi=vecAziSign
                        
                    #output vector lengths to output list
                    FryAnalysis.lenList.append(vecLength)
                    #output vector azis to output list
                    FryAnalysis.vecList.append(vecAzi)
                    ID2 = ID2 +1
                    fkey= fkey + 1
                else:
                   ID2 = ID2 + 1
                
            ID = ID +1
            
            #print FryAnalysis.fryListX
"""