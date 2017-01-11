for name in dir():
    if not name.startswith('_'):
        del globals()[name]
import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import csv
import numpy as np
from collections import Counter
global renderer
renderer = slicer.app.layoutManager().threeDWidget(0).threeDView().renderWindow().GetRenderers().GetFirstRenderer()
renderWindow = slicer.app.layoutManager().threeDWidget(0).threeDView().renderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = slicer.app.layoutManager().threeDWidget(0).threeDView().interactor()
renderWindowInteractor.SetRenderWindow(renderWindow)
#
# MultiXplore
#

class MultiXplore(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "MultiXplore" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Saeed M. Bakhshmand (Western Uni.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# MultiXploreWidget
#

class MultiXploreWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    global parametersFormLayout,FunctionalFormLayout,MatrixCollapsibleButton,MatrixGridLayout
    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Inputs"	
    self.layout.addWidget(parametersCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    	
    #
    # Functional Area
    #
    FunctionalCollapsibleButton = ctk.ctkCollapsibleButton()
    FunctionalCollapsibleButton.text = "Functional Buttons"    	
    self.layout.addWidget(FunctionalCollapsibleButton)
    FunctionalFormLayout = qt.QGridLayout(FunctionalCollapsibleButton)
	#
			
    # Matrices Area
    MatrixCollapsibleButton = ctk.ctkCollapsibleButton()
    MatrixCollapsibleButton.text = "Matrices"
    self.layout.addWidget(MatrixCollapsibleButton)
    MatrixGridLayout = qt.QGridLayout(MatrixCollapsibleButton)
	
	
	#
	#Input tractography file selector
	#
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.inputSelector.selectNodeUponCreation = True	
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = True
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the Full-brain Tractography file." )
    parametersFormLayout.addRow("Input Fiber Bundle: ", self.inputSelector)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_fb) #fb = self.inputSelector.currentNode() add this to applybutton
	
	
    #
    # Input the Models Hierarchy Node
    #
    self.Models = slicer.qMRMLNodeComboBox()
    self.Models.nodeTypes = ["vtkMRMLModelHierarchyNode"]
    self.Models.selectNodeUponCreation = True
    self.Models.addEnabled = False
    self.Models.removeEnabled = False
    self.Models.noneEnabled = True
    self.Models.showHidden = False
    self.Models.showChildNodeTypes = False
    self.Models.setMRMLScene( slicer.mrmlScene )
    self.Models.setToolTip( "Pick the Hierarchy of Brain Surface Models." )
    parametersFormLayout.addRow("Input the Generated Models: ", self.Models)
    self.Models.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Mod)
	
	#
    # Input the Labeled Volume
    #
    self.Label = slicer.qMRMLNodeComboBox()
    self.Label.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.Label.selectNodeUponCreation = True
    self.Label.addEnabled = False
    self.Label.removeEnabled = False
    self.Label.noneEnabled = True
    self.Label.showHidden = False
    self.Label.showChildNodeTypes = False
    self.Label.setMRMLScene( slicer.mrmlScene )
    self.Label.setToolTip( "Pick the Atlas-based Segmented Brain Volume." )
    parametersFormLayout.addRow("Input the Labeled Volume: ", self.Label)
    self.Label.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_label) # missing self.OnSelect_label function
		
	
    #
    # output fiberbundle generator
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Create or pick the output FiberBundle file." )
    parametersFormLayout.addRow("Output Fiberbundle: ", self.outputSelector)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_selected_fibers)
 	
	#
    # ROI Selector
    #
    self.ROISelector = slicer.qMRMLNodeComboBox()
    self.ROISelector.nodeTypes = ["vtkMRMLAnnotationROINode"]
    self.ROISelector.selectNodeUponCreation = True
    self.ROISelector.addEnabled = True
    self.ROISelector.removeEnabled = True
    self.ROISelector.noneEnabled = True
    self.ROISelector.showHidden = False
    self.ROISelector.showChildNodeTypes = False
    self.ROISelector.setMRMLScene( slicer.mrmlScene )
    self.ROISelector.setToolTip( "Pick the ROI Interactive Box to define the region for labeling of Fibers." )
    parametersFormLayout.addRow("ROI for Labeling: ", self.ROISelector)
    self.ROISelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_ROI)
    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Load FC Matrices")
    FunctionalFormLayout.addWidget(self.applyButton,0,0,1,1)
    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)		
    self.layout.addStretch(1)
    #
    # Apply Button
    #
    self.applyButton_6 = qt.QPushButton("Load SC Matrix")
    FunctionalFormLayout.addWidget(self.applyButton_6,0,1,1,1)
    # connections
    self.applyButton_6.connect('clicked(bool)', self.onApplyButton_6)
    self.layout.addStretch(1)
	
    # Apply Button
    #
    self.applyButton_2 = qt.QPushButton("Load Headers")
    self.applyButton_2.toolTip = "Run the algorithm!!."
    #self.applyButton.enabled = True
    FunctionalFormLayout.addWidget(self.applyButton_2,0,2,1,1)
    # connections
    self.applyButton_2.connect('clicked(bool)', self.onApplyButton_2)
    self.layout.addStretch(1)

	
    #
    # Apply Button
    #
    self.applyButton_4 = qt.QPushButton("Generate Arrays")
    self.applyButton_4.toolTip = "Generate and Save Arrays!!."
    #self.applyButton_4.enabled = True
    FunctionalFormLayout.addWidget(self.applyButton_4,0,3,1,1)
    # connections
    self.applyButton_4.connect('clicked(bool)', self.onApplyButton_4)
		
	#
    # Apply Button
    #
    self.applyButton_5 = qt.QPushButton("Load Saved Array")
    self.applyButton_5.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    FunctionalFormLayout.addWidget(self.applyButton_5,1,0,1,1)
    # connections
    self.applyButton_5.connect('clicked(bool)', self.onApplyButton_5)

    #
    # Apply Button
    #
    self.applyButton_7 = qt.QCheckBox("Refresh Fiber Annotation")
    self.applyButton_7.toolTip = "Update Fiber Annotation!!."
    #self.applyButton_5.enabled = True
    FunctionalFormLayout.addWidget(self.applyButton_7,1,1,1,1)
    # connections
    self.applyButton_7.stateChanged.connect(self.onApplyButton_7)

    #
    # Apply Button
    #
    self.applyButton_8 = qt.QCheckBox("Refresh Batch Selection")
    self.applyButton_8.toolTip = "Update Batch Selection!!."
    #self.applyButton_5.enabled = True
    FunctionalFormLayout.addWidget(self.applyButton_8,1,2,1,1)

    # connections
    self.applyButton_8.stateChanged.connect(self.onApplyButton_8)

 
  def OnSelect_fb(self):
    global fb_Node
    global input
    global oldTensors
    global inPts
    global inLines
    fb_Node = self.inputSelector.currentNode()
    input = fb_Node.GetPolyData()
    oldTensors = input.GetPointData().GetTensors()
    inPts =input.GetPoints()
    inLines = input.GetLines()
    ptids = vtk.vtkIdList()
    inLines.InitTraversal()
    outarray = vtk.vtkStringArray()
    outarray.SetName('Tube Labels')
    #for lidx in range(0, input.GetNumberOfLines()):
    #    inLines.GetNextCell(ptids)
    #    for pidx in range(0, ptids.GetNumberOfIds()):
    #        outarray.InsertNextTuple1(lidx)
    #input.GetPointData().AddArray(outarray)
  

  def OnSelect_Mod(self):
    global Node_Models
    global node_names
    node_names =[]
    Node_Models = self.Models.currentNode()
    num_models = Node_Models.GetNumberOfChildrenNodes()
    for t in range(0,num_models):
       child = Node_Models.GetNthChildNode(t)
       node = child.GetModelNode()
       name = node.GetName()        
       node_names.append(name[(name.rfind('-')-2):name.__len__()])
       #print(node_names)

    	
  def OnSelect_label(self):
    global InputLabel_A
    InputLabel_A = self.Label.currentNode()
    #print('InputLabel_A')
    
  def OnSelect_selected_fibers(self):
    global Selected_fibers_Node
    Selected_fibers_Node = self.outputSelector.currentNode()
    Selected_fibers_Node.SetAndObserveTransformNodeID(fb_Node.GetTransformNodeID())
    Selected_fibers_Node.CreateDefaultDisplayNodes()
    Selected_fibers_Node.CreateDefaultStorageNode()
    Selected_fibers_Node.SetName('Selected_Fibers')
    
  def OnSelect_ROI(self):
    global bounds,markupNode
    bounds = [0,0,0,0,0,0]
    markupNode = self.ROISelector.currentNode()


  def onApplyButton(self):
    global file_name,logic,Tables,Models,Tabs,RadioButtons,RadioLayout,Button_grp,Items
    logic = MultiXploreLogic()
    file_names = qt.QFileDialog.getOpenFileNames()
    Tables = []
    Models = []
    Items = []
    RadioLayout = qt.QHBoxLayout()
    RadioButtons = []
    Button_grp = qt.QButtonGroup()
    Tabs = qt.QTabWidget()
    palette = qt.QPalette()
    palette.setColor(qt.QPalette.Background,qt.Qt.red)
    for t in range(0,file_names.__len__()):
        b_tmp,model_tmp, item_tmp = logic.run_8(file_names[t])
        Tables.append(b_tmp)
        Models.append(model_tmp)
        Items.append(item_tmp)
        RadioButtons.append(qt.QRadioButton("T"+str(1+t)))
        Button_grp.addButton(RadioButtons[t],t)
        RadioButtons[t].connect('clicked()',logic.radio_button_clicked)
        #RadioButtons[t].setPalette(palette)
        RadioLayout.addWidget(RadioButtons[t],2)
        Tabs.addTab(Tables[t],"Tab"+str(1+t))
    #logic.clearLayout(MatrixGridLayout)
    logic.run_1()
	
  def onApplyButton_6(self):
    global b2,model2
    file_name_3 = qt.QFileDialog.getOpenFileName()
    b2,model2,item2 = logic.run_8(file_name_3)
    logic.run_6()
  
  def onApplyButton_2(self):#Headers
    global file_name_2
    file_name_2 = qt.QFileDialog.getOpenFileName()
    logic.run_2()
	 
  def onApplyButton_4(self):
    global Dir_name_save
    Dir_name_save = qt.QFileDialog.getExistingDirectory() 
    logic.run_4() 
	
  def onApplyButton_5(self):
    global Directory_name
    Directory_name = qt.QFileDialog.getExistingDirectory()
    logic.run_5()
    
  def onApplyButton_7(self):#Fiber Labeling
    if self.applyButton_7.isChecked():
      logic.run_7()	  

  def onApplyButton_8(self):#Batch selection update
    if self.applyButton_8.isChecked():
      logic.Slider_Selection()
	 

class MultiXploreLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def radio_button_clicked(self):
    global selection,selectionModel,ind_2,b,model,selectionModel_2,Items,item
    selected_matrix = Button_grp.checkedId()
    if 'b' in globals():
       b.selectionModel().selectionChanged.disconnect(logic.cellSelected)
       selectionModel_2.select(selection,qt.QItemSelectionModel.Deselect)
    b = Tables[selected_matrix]
    model = Models[selected_matrix]
    item = Items[selected_matrix]
    b.selectionModel().selectionChanged.connect(logic.cellSelected)
    topleft = qt.QModelIndex()
    bottomright = model.index(r-1,c-1,qt.QModelIndex())
    topleft = model.index(0,0,qt.QModelIndex())
    selection = qt.QItemSelection(topleft,bottomright)
    selectionModel = b.selectionModel()
    ind_2 = selection.indexes()
    selectionModel_2 = b2.selectionModel()


  def run_8(self,file_name): #Sub-Function to load each FC Matrix
    b= qt.QTableView()
    model = qt.QStandardItemModel()
    b.setModel(model)
    f = open(file_name,'rb')
    csv_f = csv.reader(f)
    # To obtain number of columns and rows
    ar = np.genfromtxt (file_name, delimiter=",")
    sz = ar.shape
    global r,c
    r = sz[0]
    c = sz[1]
    model.setRowCount(r)
    model.setColumnCount(c)
    #Adjusting size of table
    if c<20:
        #b.setFixedSize((c+8)*20,(r+3)*20)
        b.setFixedHeight((r+3)*20)
    else:
        #b.setFixedSize((28)*20,)
        b.setFixedHeight((27)*20)
    item =[]
    j=0				
    for row in csv_f:
        i=0
        for field in row:
          item.append(qt.QStandardItem()) # Set text of the cell
          model.setItem(i,j,item[j*c+i])
          x = float(field)		  
          if x < -0.5:
             item[j*c+i].setBackground(qt.QColor(0,(1+x)*510,255)) # Set color of the cell
             item[j*c+i].setData(x)
          elif x < 0 and x >= -0.5:
             item[j*c+i].setBackground(qt.QColor(0,255,(-x)*510,255)) # Set color of the cell
             item[j*c+i].setData(x)
          if x >= 0 and x <= 0.5:
             item[j*c+i].setBackground(qt.QColor(x*510,255,0,255)) # Set color of the cell
             item[j*c+i].setData(x)
          elif x > 0.5:
             item[j*c+i].setBackground(qt.QColor(255,(1-x)*510,0,255)) # Set color of the cell
             item[j*c+i].setData(x)
          b.setColumnWidth(i,23)
          i=i+1
        b.setRowHeight(j,23)
        j=j+1
    return b, model, item

  def run_1(self):  # Function for loading FC Tabs
    global widget,can_btn,layout,prog_bar,slider_low,S_box_low,slider_high,S_box_high
    widget = qt.QWidget()
    layout = qt.QGridLayout()
    slider_low = qt.QSlider()
    slider_low.setRange(-100,100)
    slider_low.setValue(100)
    slider_low.setSingleStep(1) 
    S_box_low = qt.QDoubleSpinBox()
    S_box_low.setRange(-1,1)
    S_box_low.setValue(1)
    S_box_low.setSingleStep(0.01)
    slider_low.valueChanged.connect(self.slider_update)
    S_box_low.valueChanged.connect(self.spin_low_change)
    slider_high = qt.QSlider()
    slider_high.setRange(-100,100)
    slider_high.setValue(100)
    slider_high.setSingleStep(1) 
    S_box_high = qt.QDoubleSpinBox()
    S_box_high.setRange(-1,1)
    S_box_high.setValue(1)
    S_box_high.setSingleStep(0.01)
    slider_high.valueChanged.connect(self.slider_update)
    S_box_high.valueChanged.connect(self.spin_high_change)
    prog_bar = qt.QProgressBar()
    prog_bar.setRange( 0, 100 )
    can_btn = qt.QPushButton("Cancel")
    can_btn.connect('clicked(bool)', self.killLoop)
    widget.setLayout(layout)
    layout_sliders = qt.QGridLayout()
    layout_sliders.addWidget(slider_low,0,0)
    layout_sliders.addWidget(slider_high,0,1)
    layout.addLayout(layout_sliders,1,2,1,1)
    layout.addLayout(RadioLayout,0,0,1,2)
    layout.addWidget(Tabs,1,0,1,2)
    layout.addWidget(prog_bar,2,0)
    layout.addWidget(can_btn,2,1)
    layout.addWidget(S_box_low,2,2)
    layout.addWidget(S_box_high,0,2)
    MatrixGridLayout.addWidget(widget,0,0,1,1)
    #print('Count of MG:')+str(MatrixGridLayout.count())

  def run_6(self): # Adding SC Matrix
    global widget2,layout2
    widget2 = qt.QWidget()
    layout2 = qt.QGridLayout()
    widget2.setLayout(layout2)
    layout2.addWidget(b2)
    MatrixGridLayout.addWidget(widget2,1,0)

  def run_2(self): # Loading headers for FC and SC Matrices
    f = open(file_name_2,'rb')
    lst = list(csv.reader(f))
    ls = lst[0]
    sz_l = ls.__len__()
    global Lab_lst
    Lab_lst = ls
    global Lab_num
    Lab_num = lst[1]
    for k in range(sz_l):
        for u in range(0,len(Models)):		
            Models[u].setHeaderData(k,1,str(k+1)) 
            Models[u].setHeaderData(k,2,str(k+1)+'_'+ls[k])
        model2.setHeaderData(k,1,str(k+1)) 
        model2.setHeaderData(k,2,str(k+1)+'_'+ls[k])
        #print "Loaded header: "+str(ls[k])

    font = qt.QFont()
    font.setFamily('MS Shell Dlg 2')
    font.setFixedPitch(True)
    font.setPointSize(7)
    #Hheader = b.horizontalHeader()
    #Hheader.setFont(font)
    #Hheader2 = b2.horizontalHeader()
    #Hheader2.setFont(font)
	
  def run_4(self): # Calculating and saving Mem_array, Mem_points and Label_LUT 
    imageCastLabel_A = vtk.vtkImageCast()
    imageCastLabel_A.SetInputData(InputLabel_A.GetImageData())
    imageCastLabel_A.Update()
    input = fb_Node.GetPolyData()
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    InputLabel_A.GetRASToIJKMatrix(Label_A_RASToIJK)
    sp = [0]*3
    sp = imageCastLabel_A.GetOutput().GetSpacing()
    trans = vtk.vtkTransform()
    trans.Identity();
    trans.PreMultiply();
    trans.SetMatrix(Label_A_RASToIJK)
    numPts = inPts.GetNumberOfPoints()
    numLines = inLines.GetNumberOfCells()
    pts = vtk.vtkIdList()
    addLines = []
    numNewPts = int(0)
    numNewCells = int(0)
    j = int()
    label = int()
    labelDims = imageCastLabel_A.GetOutput().GetDimensions()
    LabelScala_range = imageCastLabel_A.GetOutput().GetScalarRange()
    Mem_array = np.full((numLines,2000), False, dtype=bool)
    Mem_points = np.zeros(numLines)
    Label_LUT = []
    inLines.InitTraversal()
    print(numLines)
    for inCellId in range(0,numLines):
        inLines.GetNextCell(pts);
        npts = int(pts.GetNumberOfIds());
        pIJK = 3*[0]
        pt = 3*[int()]
        p = 3*[0]
        inPtr =[]
        if npts <19:
           print('short line!!!')
           print(inCellId)
        else:
           for j in range(0,npts):
             p = inPts.GetPoint(pts.GetId(j))
             trans.TransformPoint(p,pIJK)
             pt[0]= int(pIJK[0])
             pt[1]= int(pIJK[1])
             pt[2]= int(pIJK[2])
             if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDims[0] or pt[1] >= labelDims[1] or pt[2] >= labelDims[2]:
                 continue
             inPtr = imageCastLabel_A.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
             if inPtr<1:
                continue
             if inPtr in Label_LUT:
                  Mem_array[inCellId,Label_LUT.index(inPtr)] = True
             else: 
                  Label_LUT.append(inPtr)
                  Mem_array[inCellId,Label_LUT.index(inPtr)] = True
            
        Mem_points[inCellId] = npts
	
    os.chdir(Dir_name_save)
    np.save('Mem_array.npy', Mem_array[:,0:Label_LUT.__len__()])
    np.save('Mem_points.npy',Mem_points)
    f = open("Label_LUT", "w")
    for item in Label_LUT:
         f.write("%s " % int(item))
    f.close()
    SC_Matrix = np.zeros((Label_LUT.__len__(),Label_LUT.__len__())) # Generating SC Matrix
    for i in range(0,Label_LUT.__len__()):
        for j in range(0,Label_LUT.__len__()):
             exist = Mem_array[:,i]*Mem_array[:,j]
             SC_Matrix[i,j] = exist.sum()
    np.savetxt("SC_Matrix.csv",SC_Matrix,delimiter=",") 
	
  def run_5(self):
    global Mem_array
    Mem_array = np.load(Directory_name+'/Mem_array.npy')
    global Mem_points
    Mem_points = np.load(Directory_name+'/Mem_points.npy')
    global Label_LUT
    f=open(Directory_name+'/Label_LUT', 'rb') 
    Label_LUT = map(int,f.read().split())
    

  def cellSelected(self):
    global ind,killLoopFlag
    dd = selectionModel.selectedIndexes()
    selectionModel_2.select(selection,qt.QItemSelectionModel.Deselect)
    length = dd.__len__()
    Label_l = []
    Label_list_a = []
    Label_list_b = []
    Label_num_a = []
    Label_num_b = []	
    for s in range(length): 
       ee = dd[s]
       selectionModel_2.select(ee,qt.QItemSelectionModel.Select)
       ind = [ee.row(),ee.column()]
       Label_a = Lab_lst[ind[0]]
       Label_b = Lab_lst[ind[1]]
       Label_list_a.append(Lab_lst[ind[0]])
       Label_list_b.append(Lab_lst[ind[1]])
       Label_num_a.append(int(Lab_num[ind[0]]))
       Label_num_b.append(int(Lab_num[ind[1]]))
       Label_l.append(Label_a)  
       Label_l.append(Label_b)
    Label_l = list(set(Label_l)) # Remove duplicate list items
    print "Selected Labels: "+str(Label_l)
    if Label_l.__len__()>0:
       self.Model_vis(Label_l)	      
    #Labels=[1028,2028]
    addLine = 0*(Mem_array[:,0])
    direct_links = np.full(length, False, dtype=bool)
    point_a=[]
    point_b=[]
    for s in range(length):
        exist = Mem_array[:,Label_LUT.index(Label_num_a[s])]*Mem_array[:,Label_LUT.index(Label_num_b[s])]
        addLine = addLine+(exist)
        if int(sum(exist))<1: # indirect links
           direct_links[s] = 1
           child = Node_Models.GetNthChildNode(node_names.index(Label_list_a[s]))
           node = child.GetModelNode()
           Num_pts = node.GetPolyData().GetNumberOfPoints()
           point_a.append(node.GetPolyData().GetPoint(Num_pts/2))
           child = Node_Models.GetNthChildNode(node_names.index(Label_list_b[s]))
           node = child.GetModelNode()
           Num_pts = node.GetPolyData().GetNumberOfPoints()
           point_b.append(node.GetPolyData().GetPoint(Num_pts/2))
    
    if 'Line_Actors' in globals():
      if Line_Actors.__len__()>0:
         for g in range(0,num_lns):
             renderer.RemoveActor(Line_Actors[g])
    global num_lns
    num_lns = sum(direct_links)
    Lines = []
    Line_Mappers = []
    global Line_Actors
    Line_Actors = []
    print('Number of dotted lines: '+str(num_lns))
    for g in range(0,num_lns):
       Lines.append(vtk.vtkLineSource())
       Lines[g].SetPoint1(point_a[g])
       Lines[g].SetPoint2(point_b[g])
       Lines[g].SetResolution(100)	
       Line_Mappers.append(vtk.vtkPolyDataMapper())
       Line_Mappers[g].SetInputConnection(Lines[g].GetOutputPort()) 
       Line_Actors.append(vtk.vtkActor())
       Line_Actors[g].SetMapper(Line_Mappers[g])
       renderer.AddActor(Line_Actors[g])
       Line_Actors[g].GetProperty().SetColor(1,1,0)
       Line_Actors[g].GetProperty().SetLineStipplePattern(0xf0f0)
       Line_Actors[g].GetProperty().SetLineStippleRepeatFactor(1)
       Line_Actors[g].GetProperty().SetPointSize(1)
       Line_Actors[g].GetProperty().SetLineWidth(4)
        
    
    
    numNewCells = int(sum(addLine))
    print('Number of lines: '+str(numNewCells))
    address = np.where(addLine==1)[:][0]
    numNewPts = int(sum(Mem_points[address]))
    print('Number of points: '+str(numNewPts))
    if numNewCells>0:
       Selected_fibers_Node.SetDisplayVisibility(1)
       outFibers = vtk.vtkPolyData()
       points = vtk.vtkPoints()
       points.Allocate(numNewPts)
       outFibers.SetPoints(points)
       outFibersCellArray = vtk.vtkCellArray()
       outFibersCellArray.Allocate(numNewPts+numNewCells)
       outFibers.SetLines(outFibersCellArray)
       newTensors = vtk.vtkFloatArray()
       newTensors.SetNumberOfComponents(9)
       newTensors.Allocate(9*numNewPts)
       outFibers.GetPointData().SetTensors(newTensors)
       ptId = int()
       pcount = int()
       tensor = 9*[0]
       inLines.InitTraversal()
       address1 = np.ndarray.tolist(address)
       killLoopFlag = False
       prog_bar.reset()
       for inCellId in address1:
         npts = int(Mem_points[inCellId])
         outFibersCellArray.InsertNextCell(npts,range(ptId,ptId+npts))
         st_pt = int(sum(Mem_points[0:(inCellId)]))
         points.InsertPoints(ptId,npts,st_pt,inPts)
         if (oldTensors):
             newTensors.InsertTuples(ptId,npts,st_pt,oldTensors)
         ptId = ptId+npts
         qt.QApplication.processEvents()
         prog_bar.setValue((float(ptId)/numNewPts)*100)
         qt.QApplication.processEvents()
         if (killLoopFlag):
            print('Canceled')
            break
       print('ptId: '+str(ptId))
       Selected_fibers_Node.SetAndObservePolyData(outFibers)
    else:
        Selected_fibers_Node.SetDisplayVisibility(0)
      	
  def Model_vis(self,Labels):
        len = Labels.__len__()
        print(len)
        for s in range(0,node_names.__len__()):
            child = Node_Models.GetNthChildNode(s)
            node = child.GetModelNode()
            node.SetDisplayVisibility(0)            
        for r in range(0,len):
            child = Node_Models.GetNthChildNode(node_names.index(Labels[r]))
            node = child.GetModelNode()
            d_node = node.GetDisplayNode()
            #print(node.GetName())
            node.SetDisplayVisibility(1)
            d_node.SetOpacity(0.75)		    

      
  def killLoop(self): 
      global killLoopFlag
      killLoopFlag = True
  
  def slider_update(self):
      global slid_val_low,slid_val_high
      print "slider low: " + str(slider_low.value)+"slider high:" + str(slider_high.value)
      slid_val_low = slider_low.value/float(100)
      S_box_low.setValue(slid_val_low)
      slid_val_high = slider_high.value/float(100)
      S_box_high.setValue(slid_val_high)
  def spin_low_change(self):
      print "Spin low change: " + str(S_box_low.value)
      slider_low.setValue(S_box_low.value*100) 

  def spin_high_change(self):
      print "Spin high change: " + str(S_box_high.value)
      slider_high.setValue(S_box_high.value*100)

  def Slider_Selection(self):
      selectionModel.select(selection,qt.QItemSelectionModel.Deselect)
      selectionModel_2.select(selection,qt.QItemSelectionModel.Deselect)
      for h in range(0,r*c):
        if (item[h].data()>slid_val_low) & (item[h].data()<slid_val_high):

           selectionModel.select(ind_2[h],qt.QItemSelectionModel.Select)
           selectionModel_2.select(ind_2[h],qt.QItemSelectionModel.Select)
      
  def run_7(self):#Label the fibers
    markupNode.GetRASBounds(bounds)
    line = fb_Node.GetLineDisplayNode()
    ids = vtk.vtkIdFilter()
    ids.SetInputConnection(line.GetOutputPolyDataConnection())
    ids.CellIdsOn()
    ids.FieldDataOn()
    clipper = vtk.vtkBoxClipDataSet()
    clipper.SetInputConnection(ids.GetOutputPort())
    clipper.SetBoxClip(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5])
    clipper.GenerateClippedOutputOff()
    clipper.Update()
    cc = vtk.vtkCellCenters()
    cc.SetInputConnection(clipper.GetOutputPort())

    line2 = Selected_fibers_Node.GetLineDisplayNode()
    ids2 = vtk.vtkIdFilter()
    ids2.SetInputConnection(line2.GetOutputPolyDataConnection())
    ids2.CellIdsOn()
    ids2.FieldDataOn()
    clipper2 = vtk.vtkBoxClipDataSet()
    clipper2.SetInputConnection(ids2.GetOutputPort())
    clipper2.SetBoxClip(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5])
    clipper2.GenerateClippedOutputOff()
    clipper2.Update()
    cc2 = vtk.vtkCellCenters()
    cc2.SetInputConnection(clipper2.GetOutputPort())

    cc3 = vtk.vtkCellCenters()
    visCells3 = vtk.vtkSelectVisiblePoints()
    visCells3.SetRenderer(renderer)
    cellMapper3 = vtk.vtkLabeledDataMapper()
    cellMapper3.SetInputConnection(visCells3.GetOutputPort())
    cellMapper3.SetLabelModeToLabelFieldData()
    cellMapper3.GetLabelTextProperty().SetColor(0.3, 0, 0.1)
    if 'cellLabels3' in globals():
       renderer.RemoveActor2D(cellLabels3)
    global cellLabels3
    cellLabels3 = vtk.vtkActor2D()
    cellLabels3.GetProperty().SetOpacity(0.7)
    cellLabels3.SetMapper(cellMapper3)
    
    print(bounds)
    pd_large = cc.GetOutput()
    pd_small = cc2.GetOutput()
    pd_new = vtk.vtkPolyData()
    points_new = vtk.vtkPoints()
    point_data_new = pd_new.GetPointData()
    array_new = vtk.vtkIdTypeArray()
    point_data_new.AddArray(array_new)
    point_data = pd_large.GetPointData()
    id = point_data.GetArray(0)
    points_large = pd_large.GetPoints()
    points_small = pd_small.GetPoints()
    cc.Update()
    cc2.Update()
    #print(cc.GetOutput())
    print('Small:'+str(pd_small.GetNumberOfPoints()))
    print('Large:'+str(pd_large.GetNumberOfPoints()))    
    
    for i in range(0,pd_small.GetNumberOfPoints()):
      point1 = pd_small.GetPoint(i)
      for j in range(0,pd_large.GetNumberOfPoints()):
         point2 = pd_large.GetPoint(j)
         if point1==point2:
             points_new.InsertNextPoint(point1)
             array_new.InsertNextTuple(point_data.GetArray(0).GetTuple(j))
             break
    pd_new.SetPoints(points_new)
    print('Number of new points:'+str(pd_new.GetNumberOfPoints()))
    visCells3.SetInputData(pd_new)
    renderer.AddActor2D(cellLabels3)
    #renderer.RemoveActor2D(cellLabels)
    #renderer.RemoveActor2D(cellLabels2)
    

    
	
 	

class MultiXploreTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_MultiXplore1()

  def test_MultiXplore1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = MultiXploreLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
