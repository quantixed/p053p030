#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Ternary Diagram>

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////

Menu "Macros"
	"Spindle Quantification...",  SpindleQuant()
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////

// We need to load in data that is organised in the following way:
// Condition A x Condition B where A != B Multiple files for each pair.
// So 1st level directory on disk is Condition A
// 2nd level directory on disk is Condition B
// 3rd level directory is each cell - contains one file we need, results.csv.

Function SpindleQuant()
	SetDataFolder root:
	// kill all windows and waves before we start
	CleanSlate()
	
	if(LoadSQData() == 0)
		CollectAllMeasurements()
		DetermineNamesAndLabels()
		MakeTheBoxPlots()
		MakeTheLayouts("p1",5,4)
		MakeTheScatterPlots()
		MakeTheLayouts("p2",4,2)
		DoWindow/F allp2Layout
		Legend/C/N=text0/J/F=0/S=3/A=LB "\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_CHTOG_Spindle) CHTOG_Spindle\r\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_CHTOG_Mito) CHTOG_Mito";DelayUpdate
		AppendText "\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_CLATHRIN_Spindle) CLATHRIN_Spindle\r\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_CLATHRIN_Mito) CLATHRIN_Mito";DelayUpdate
		AppendText "\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_GTSE1_Spindle) GTSE1_Spindle\r\\s(p2_GFPFKBPTACC3.GFPFKBPTACC3_GTSE1_Mito) GTSE1_Mito";DelayUpdate
		AppendText "\\s(p2_CLTAFKBPGFP.CLTAFKBPGFP_TACC3_Spindle) TACC3_Spindle\r\\s(p2_CLTAFKBPGFP.CLTAFKBPGFP_TACC3_Mito) TACC3_Mito\r"
		MakeTernaryPlots()
		MakeTheLayouts("p3",4,2)
		MakeTheBiPlots()
		MakeTheLayouts("p4",6,4)
	else
		return -1
	endif
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
// loads the data and stores in a data folder hierarchy
Function LoadSQData()
	String condADiskFolderName, condBDiskFolderName, cellDiskFolderName
	String condAList, condBList, cellList
	String condA, condB, cell
	String safeCondA, safeCondB, safeCell
	String thisFile
	Variable nCondA, nCondB, nCell
	Variable i,j,k
	
	NewPath/O/Q/M="Please find disk folder" ExpDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	PathInfo/S ExpDiskFolder
	String ExpDiskFolderName = S_Path
	
	String dataFolderName = "root:data"
	NewDataFolder/O/S $dataFolderName // make root:data: but don't put anything in it yet
	
	// get a list of all condA dirs (not full path)
	condAList = IndexedDir(expDiskFolder,-1,0)
	nCondA = ItemsInList(condAList)
	// loop condA
	for(i = 0; i < nCondA; i += 1)
		condA = StringFromList(i, condAList)
		NewPath/O/Q condAPath, ExpDiskFolderName + condA
		// get a list of all condB dirs in condA (not full path)
		condBList = IndexedDir(condAPath,-1,0)
		nCondB = ItemsInList(condBList)
		// make data folder for condA - no hyphens
		safeCondA = CleanupName(ReplaceString("-",condA,""),0)
		dataFolderName = "root:data:" + safeCondA
		NewDataFolder/O/S $dataFolderName
		
		for(j = 0; j < nCondB; j += 1)
			condB = StringFromList(j, condBList)
			NewPath/O/Q condBPath, ExpDiskFolderName + condA + ":" + condB
			// get a list of all cell dirs in condB
			cellList = IndexedDir(condBPath,-1,0)
			nCell = ItemsInList(cellList)
			// make data folder for condB
			safeCondB = CleanupName(ReplaceString("-",condB,""),0)
			dataFolderName = "root:data:" + safeCondA + ":" + safeCondB
			NewDataFolder/O/S $dataFolderName
			
			for(k = 0; k < nCell; k += 1)
				cell = StringFromList(k, cellList)
				NewPath/O/Q cellPath, ExpDiskFolderName + condA + ":" + condB + ":" + cell
				// make data folder for cell
				safeCell = CleanupName(ReplaceString("-",cell,""),0)
				dataFolderName = "root:data:" + safeCondA + ":" + safeCondB + ":" + safeCell
				NewDataFolder/O/S $dataFolderName
				// we know that we need only one file, otherwise we could get list here
				LoadWave/A/Q/J/K=0/W/O/L={0,1,0,0,0}/P=cellPath "results.csv"
				// this gives us four waves called file, type, time and measure in the exp folder
				ProcessTheRawData()
				// now store the names - disk names, then Igor names
				Make/O/N=(1,6)/T theNameWave = {{condA},{condB},{cell},{safeCondA},{safeCondB},{safeCell}}
			endfor
		endfor
	endfor
	SetDataFolder root:
	return 0
End

STATIC Function ProcessTheRawData()
	// we should be inside a data sub-sub-folder which contains the raw data
	WAVE/Z w = measure
	MatrixOp/O allData = w ^t
	// log2 of post/pre for spindle, cyto mito without or with cell correction
	Make/O/N=(1,3) dataGreenNormToPre, dataGreenNormToPreCell
	Make/O/N=(1,3) dataRedNormToPre, dataRedNormToPreCell
	// assign values - green
	dataGreenNormToPre[0][] = log(w[q + 8] / w[q]) / log(2)
	Variable bleach = w[3] / w[11] // cell value pre / post
	dataGreenNormToPreCell[0][] = log((w[q + 8] * bleach) / w[q]) / log(2)
	// assign values - red
	dataRedNormToPre[0][] = log(w[q + 12] / w[q + 4]) / log(2)
	bleach = w[7] / w[15] // cell value pre / post for red
	dataRedNormToPreCell[0][] = log((w[q + 12] * bleach) / w[q + 4]) / log(2)
	
	// fraction of signal that is spindle, cyto, mito pre or post (no cell correction because it should not matter)
	Make/O/N=(1,3) fracGreenPre, fracGreenPost
	Make/O/N=(1,3) fracRedPre, fracRedPost
	// assign values - green pre first
	fracGreenPre[0][] = w[q + 0] / sum(w,0,2)
	fracRedPre[0][] = w[q + 4] / sum(w,4,6)
	fracGreenPost[0][] = w[q + 8] / sum(w,8,10)
	fracRedPost[0][] = w[q + 12] / sum(w,12,14)
	
	return 0
End

// collects the measurements from the nested datafolders
Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	Variable nCondAFolders = CountObjectsDFR(dfr, 4)
	
	String folderNameA, folderNameB, folderNameCell
	Variable nCondBFolders, nCellFolders
	String wName, genConcStr, concStr
//	String targetList = "dataGreenNormToPre;dataGreenNormToPreCell;dataRedNormToPre;dataRedNormToPreCell;allData;"
	String targetList = "dataGreenNormToPreCell;dataRedNormToPreCell;allData;theNameWave;"
	targetList += "fracGreenPre;fracGreenPost;fracRedPre;fracRedPost;"
	
	Variable i,j,k

	for(i = 0; i < nCondAFolders; i += 1)
		folderNameA = GetIndexedObjNameDFR(dfr, 4, i)
		DFREF dfrA = root:data:$folderNameA
		nCondBFolders = CountObjectsDFR(dfrA, 4)
		
		for(j = 0; j < nCondBFolders; j += 1)
			folderNameB = GetIndexedObjNameDFR(dfrA, 4, j)
			DFREF dfrB = dfrA:$folderNameB
			nCellFolders = CountObjectsDFR(dfrB, 4)
			// store the names of the cell folders in a text wave named condA_condB_cellName
			wName = folderNameA + "_" + folderNameB + "_cellName"
			Make/O/N=(nCellFolders)/T dfr:$wName = GetIndexedObjNameDFR(dfrB, 4, p)
			genConcStr = ""
			
			for(k = 0; k < nCellFolders; k += 1)
				genConcStr += "root:data:" + folderNameA + ":" + folderNameB + ":" + GetIndexedObjNameDFR(dfrB, 4, k) + ":replaceThis;"
			endfor
			
			for(k = 0; k < ItemsInList(targetList); k += 1)
				concStr = ReplaceString("replaceThis",genConcStr,StringFromList(k,targetList))
				wName = folderNameA + "_" + folderNameB + "_" + StringFromList(k,targetList)
				Concatenate/O/NP=0 concStr, dfr:$wName
			endfor
		endfor
	endfor
	
	SetDataFolder root:
	
	return 0
End


Function DetermineNamesAndLabels()
	SetDataFolder root:data:	// relies on earlier load
	String wList = WaveList("*_theNameWave",";","")
	Concatenate/O/NP=0/KILL/T wList, allTheNames
	Duplicate/O/RMD=[][3] allTheNames,condAIgorNames
	FindDuplicates/RT=condAUniqueIgorNames condAIgornames
	Duplicate/O/RMD=[][4] allTheNames,condBIgorNames
	FindDuplicates/RT=condBUniqueIgorNames condBIgornames
	Make/O/N=(numpnts(condAUniqueIgorNames))/T condAUniqueLabels
	Make/O/N=(numpnts(condBUniqueIgorNames))/T condBUniqueLabels
	String str
	
	Variable i
	
	for(i = 0; i < numpnts(condAUniqueLabels); i += 1)
		str = condAUniqueIgorNames[i]
		FindValue/TEXT=str condAIgornames
		condAUniqueLabels[i] = allTheNames[V_Value][0]
	endfor
	for(i = 0; i < numpnts(condBUniqueLabels); i += 1)
		str = condBUniqueIgorNames[i]
		FindValue/TEXT=str condBIgornames
		condBUniqueLabels[i] = allTheNames[V_Value][1]
	endfor
	
	// build a colourwave for all B conditions
	MakeColorWave(numpnts(condBUniqueLabels),"colorWave",alpha = 127 * 257)
	
	return 0
End


Function MakeTheBoxPlots()
	SetDataFolder root:data:
	Make/O/N=(3)/T labelWave = {"Spindle","Cyto","Mito"}
	String wList = WaveList("*_data*",";","")
	wList = SortList(wList) // should sort condA then condB then dataG, dataR
	Variable nPlots = ItemsInList(wList)
	String splitStr = "", labelStr = ""
	
	String wName,plotName
	
	Variable i
	
	for(i = 0; i < nPlots; i+= 1)
		wName = StringfromList(i,wList)
		plotName = "p1_" + wName
		KillWindow/Z $plotName
		Display/N=$plotName/HIDE=1
		AppendBoxPlot/W=$plotName $wName vs labelWave
		SetAxis/W=$plotName left -2,2
		ModifyGraph/W=$plotName zero(left)=4
		if(stringmatch(wName,"*dataRed*") == 1)
			ModifyGraph/W=$plotName rgb=(235 * 257,28 * 257,36 * 257,127 * 257)
		elseif(stringmatch(wName,"*dataGreen*") == 1)	
			ModifyGraph/W=$plotName rgb=(57 * 257,181 * 257,74 * 257,127 * 257)
		endif
		ModifyBoxPlot/W=$plotName trace=$wName,whiskerMethod=4
		ModifyBoxPlot/W=$plotName trace=$wName,markers={19,-1,19},markerSizes={2,2,2}
		ModifyGraph/W=$plotName margin(left)=39,margin(right)=9
		Label/W=$plotName left "Post / Pre (Log\\B2\\M)"
		splitStr = ReplaceString("_dataRedNormToPreCell",wName,"")
		splitStr = ReplaceString("_dataGreenNormToPreCell",splitStr,"")
//		splitStr = ReplaceString("_dataRedNormToCell",splitStr,"")
//		splitStr = ReplaceString("_dataGreenNormToCell",splitStr,"")
		labelStr = "\\JR" + ReplaceString("_", splitStr, "\r")
		TextBox/W=$plotName/C/N=text0/F=0/A=RT/X=0.00/Y=0.00/E=2 labelStr
	endfor
	
	// exit back to root
	SetDataFolder root:
End

Function MakeTheScatterPlots()
	SetDataFolder root:data:
	WAVE/Z colorWave = root:colorWave
	WAVE/Z/T condBUniqueIgorNames
	DFREF dfr = GetDataFolderDFR()
	Variable nCondAFolders = CountObjectsDFR(dfr, 4)
	String folderNameA, wList, wName, plotName, correspondingWName, tName0, tName1, searchStr
	String aveName
	Variable nWaves, theRow
	
	Variable i,j
	
	for(i = 0; i < nCondAFolders; i += 1)
		folderNameA = GetIndexedObjNameDFR(dfr, 4, i)
		wList = WaveList(folderNameA + "*_dataGreenNormToPreCell",";","")
		wList = SortList(wList)
		// set up scatter plot window
		plotName = "p2_" + folderNameA
		KillWindow/Z $plotName
		Display/N=$plotName/HIDE=1
		// waves to add
		nWaves = ItemsInList(wList)
		// make wave to hold the average and sem
		aveName = folderNameA + "_aveGreen"
		Make/O/N=(numpnts(condBUniqueIgorNames),3) $aveName, $(ReplaceString("_ave",aveName,"_err"))
		Wave wAveG = $aveName
		Wave wErrG = $(ReplaceString("_ave",aveName,"_err"))
		aveName = folderNameA + "_aveRed"
		Make/O/N=(numpnts(condBUniqueIgorNames),3) $aveName, $(ReplaceString("_ave",aveName,"_err"))
		Wave wAveR = $aveName
		Wave wErrR = $(ReplaceString("_ave",aveName,"_err"))
		// set all to NaN - this way the blanks are predone
		wAveG[][] = NaN
		wErrG[][] = NaN
		wAveR[][] = NaN
		wErrR[][] = NaN
		
		for(j = 0; j < nWaves; j += 1)
			wName = StringfromList(j,wList)
			// find the condB row
			searchStr = ReplaceString(folderNameA + "_",wName,"")
			searchStr = ReplaceString("_dataGreenNormToPreCell",searchStr,"")
			FindValue/TEXT=searchStr condBUniqueIgorNames
			// row is held in V_Value
			theRow = V_Value
			// exclude TUBULIN
			if(stringmatch(wName,"*TUBULIN*") == 1)
				continue
			endif
			correspondingWName = ReplaceString("Green",wName,"Red")
			tName0 = ReplaceString("dataGreenNormToPreCell",wName,"Spindle")
			tName1 = ReplaceString("dataGreenNormToPreCell",wName,"Mito")
			AppendToGraph/W=$plotName $correspondingWName[][0]/TN=$tName0 vs $wName[][0]
			AppendToGraph/W=$plotName $correspondingWName[][2]/TN=$tName1 vs $wName[][2]
			ModifyGraph/W=$plotName mode=3
			ModifyGraph/W=$plotName marker($tName0)=19
			ModifyGraph/W=$plotName marker($tName1)=16
			ModifyGraph/W=$plotName msize=2,mrkThick=0
			ModifyGraph/W=$plotName rgb($tName0)=(colorWave[V_Value][0],colorWave[V_Value][1],colorWave[V_Value][2],colorWave[V_Value][3])
			ModifyGraph/W=$plotName rgb($tName1)=(colorWave[V_Value][0],colorWave[V_Value][1],colorWave[V_Value][2],colorWave[V_Value][3])
			// insert average
			Wave w0 = $wName
			if(numtype(sum(w0)) == 2)
				MatrixOp/O/FREE meanMat = averageCols(w0)
				MatrixOp/O/FREE semMat = sqrt(varCols(w0) / numrows(w0))
			else
				Duplicate/O/FREE w0,sumMat,countMat
				sumMat[][] = (numtype(w0[p][q]) == 2) ? 0 : w0[p][q]
				countMat[][] = (numtype(w0[p][q]) == 0) ? 1 : 0
				MatrixOp/O/FREE meanMat = sumCols(sumMat) / sumCols(countMat)
				MatrixOp/O/FREE semMat = sqrt(varCols(w0) / sumCols(countMat))
			endif
			wAveG[theRow][] = meanMat[0][q]
			wErrG[theRow][] = semMat[0][q]

			Wave w1 = $correspondingWName
			if(numtype(sum(w1)) == 0)
				MatrixOp/O/FREE meanMat = averageCols(w1)
				MatrixOp/O/FREE semMat = sqrt(varCols(w1) / numrows(w1))
			else
				Duplicate/O/FREE w1,sumMat,countMat
				sumMat[][] = (numtype(w1[p][q]) == 2) ? 0 : w1[p][q]
				countMat[][] = (numtype(w1[p][q]) == 0) ? 1 : 0
				MatrixOp/O/FREE meanMat = sumCols(sumMat) / sumCols(countMat)
				MatrixOp/O/FREE semMat = sqrt(varCols(w1) / sumCols(countMat))
			endif
			wAveR[theRow][] = meanMat[0][q]
			wErrR[theRow][] = semMat[0][q]
		endfor
		SetAxis/W=$plotName left -2,2
		SetAxis/W=$plotName bottom -2,2
		ModifyGraph/W=$plotName zero=2
		ModifyGraph width={Aspect,1}
		Label/W=$plotName left "Red Post / Pre (Log\\B2\\M)"
		Label/W=$plotName bottom "Green Post / Pre (Log\\B2\\M)"
		TextBox/W=$plotName/C/N=text0/F=0/A=RT/X=0.00/Y=0.00/E=2 folderNameA
		// add the average waves
		AppendToGraph/W=$plotName wAveR[][0]/TN=aveSpindle vs wAveG[][0]
		AppendToGraph/W=$plotName wAveR[][2]/TN=aveMito vs wAveG[][2]
		ModifyGraph/W=$plotName mode=3,marker(aveSpindle)=19,marker(aveMito)=16
		ModifyGraph/W=$plotName zColor(aveMito)={colorWave,*,*,directRGB,0}
		ModifyGraph/W=$plotName msize(aveMito)=3,mrkThick(aveMito)=0
		ModifyGraph/W=$plotName zColor(aveSpindle)={colorWave,*,*,directRGB,0}
		ModifyGraph/W=$plotName msize(aveSpindle)=3,mrkThick(aveSpindle)=0.5,useMrkStrokeRGB(aveSpindle)=1
		ModifyGraph/W=$plotName msize(aveMito)=3,mrkThick(aveMito)=0.5,useMrkStrokeRGB(aveMito)=1
		ErrorBars/W=$plotName/RGB=(0,0,0) aveSpindle XY,wave=(wErrG[*][0],wErrG[*][0]),wave=(wErrR[*][0],wErrR[*][0])
		ErrorBars/W=$plotName/RGB=(0,0,0) aveMito XY,wave=(wErrG[*][2],wErrG[*][2]),wave=(wErrR[*][2],wErrR[*][2])
	endfor
	
	// exit back to root
	SetDataFolder root:
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol)
	String prefix
	Variable nRow, nCol
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"_*",";","WIN:1")
	String modList = allList // in case we want to add filtering
	Variable nWindows = ItemsInList(modList)
	String plotName
	
	Variable i
	
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum = 1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			// LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(57, 57, 57, 57)
			ModifyLayout/W=$layoutName units=0
			ModifyLayout/W=$layoutName frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=$layoutName appendpage
				pgNum += 1
				LayoutPageAction/W=$layoutName page=(pgNum)
			endif
		endif
	endfor
End

Function MakeTernaryPlots()
	SetDataFolder root:data:
	String wList = WaveList("*_frac*pre*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave m0 = $wName
		Wave m1 = $(ReplaceString("pre",wName,"post"))
		PrepareDataForTernary(m0,m1)
	endfor
	SetDataFolder root:
End

// start with 2D waves (2 columns for each ternary axis)
// m0 is the first point m1 is the second
STATIC Function PrepareDataForTernary(m0,m1)
	Wave m0,m1
	
	Variable nRow = DimSize(m0,0)
	Variable nCol = DimSize(m0,1)
	
	if(nRow != DimSize(m1,0) || nCol != DimSize(m1,1))
		DoAlert 0, "Need two waves of the same size"
	elseif(nCol != 3)
		DoAlert 0, "Three input columns required"
	endif
	
	String mName0 = NameOfWave(m0)
	String newName = ReplaceString("_frac",mName0,"_tern")
	newName = ReplaceString("Pre",newName,"")
	String mName1 = NameOfWave(m1)
	SplitWave/O m0
	SplitWave/O m1
	
	Variable i
	
	for(i = 0; i < 3; i += 1)
		Wave loopW0 = $(mName0 + num2str(i))
		Wave loopW1 = $(mName1 + num2str(i))
		ConnectPointsInAWave(loopW0,loopW1, newName + num2str(i))
		KillWaves/Z loopW0,loopW1
	endfor
	Wave w0 = $(newName + "0")
	Wave w1 = $(newName + "1")
	Wave w2 = $(newName + "2")
	String plotName = "p3_" + newName
	TernaryDiagramModule#NewTernaryGraph(w0, w1, w2, plotName)
	ModifyGraph/W=$plotName mode=4
	ModifyGraph/W=$plotName arrowMarker(TernaryTraceData0)={_inline_,1,5,2,6,barbSharp= 1}
	if(stringmatch(newName,"*green*") == 1)
		ModifyGraph/W=$plotName rgb=(hexcolor_red(0x00A651),hexcolor_green(0x00A651),hexcolor_blue(0x00A651),32768)
	else
		ModifyGraph/W=$plotName rgb=(hexcolor_red(0xED1C24),hexcolor_green(0xED1C24),hexcolor_blue(0xED1C24),32768)
	endif
End

STATIC Function ConnectPointsInAWave(w0,w1,wName)
	Wave w0,w1
	String wName
	
	Variable nRows = numpnts(w0)
	// Make 1D wave with gaps to define paired measurements
	Make/O/FREE/N=(nRows) w2 = NaN
	Concatenate/O {w0,w1,w2}, $wName
	Wave resultW = $wName
	MatrixTranspose resultW
	Redimension/N=(3 * nRows) resultW
End

// This function assumes that the ternary plots have been made
// It will use those outputs to calculate a "biplot" version
Function MakeTheBiPlots()
	SetDataFolder root:data:
	String wList = WaveList("*_tern*0",";","")
	Variable nWaves = ItemsInList(wList)
	String wName,newName2,newName3
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave w0 = $wName
		Wave w1 = $(RemoveEnding(wName) + "2")
		newName2 = ReplaceString("_tern",wName,"_bi")
//		newName3 = RemoveEnding(newName2) + "2"
//		Make/O/N=(numpnts(w0)) $newName2, $newName3
		Make/O/N=(numpnts(w0)) $newName2
		Wave w2 = $newName2
//		Wave w3 = $newName3
		w2[] = w0[p] / (w0[p] + w1[p])
	endfor
	
	wList = WaveList("*_biGreen0",";","")
	nWaves = ItemsInList(wList)
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave w0 = $wName
		Wave w1 = $(ReplaceString("_biGreen0",wName,"_biRed0"))
		GreenRedBiPlots(w0,w1)
	endfor
	
	SetDataFolder root:
End

STATIC Function GreenRedBiPlots(w0,w1)
	Wave w0,w1
	
	String plotName = "p4_" + NameOfWave(w1)
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1 w1 vs w0 // green vs red
	ModifyGraph/W=$plotName mode=4
	ModifyGraph/W=$plotName arrowMarker={_inline_,1,3,1.5,6,barbSharp= 1}
	ModifyGraph/W=$plotName rgb=(0,0,0,32768)
	ModifyGraph/W=$plotName mrkThick=0
	SetAxis/W=$plotName left 0,1
	SetAxis/W=$plotName bottom 0,1
	ModifyGraph/W=$plotName width={Aspect,1}
	ModifyGraph/W=$plotName grid=1,mirror=1
	ModifyGraph/W=$plotName gridRGB=(48059,48059,48059)
	ModifyGraph/W=$plotName manTick(left)={1,1,0,0},manMinor(left)={3,0},manTick(bottom)={1,1,0,0},manMinor(bottom)={3,0}
//	Label/W=$plotName left "Red"
	ModifyGraph/W=$plotName axRGB(left)=(hexcolor_red(0xED1C24),hexcolor_green(0xED1C24),hexcolor_blue(0xED1C24))
	ModifyGraph/W=$plotName tlblRGB(left)=(hexcolor_red(0xED1C24),hexcolor_green(0xED1C24),hexcolor_blue(0xED1C24))
	ModifyGraph/W=$plotName alblRGB(left)=(hexcolor_red(0xED1C24),hexcolor_green(0xED1C24),hexcolor_blue(0xED1C24))
//	Label/W=$plotName bottom "Green"
	ModifyGraph/W=$plotName axRGB(bottom)=(hexcolor_red(0x00A651),hexcolor_green(0x00A651),hexcolor_blue(0x00A651))
	ModifyGraph/W=$plotName tlblRGB(bottom)=(hexcolor_red(0x00A651),hexcolor_green(0x00A651),hexcolor_blue(0x00A651))
	ModifyGraph/W=$plotName alblRGB(bottom)=(hexcolor_red(0x00A651),hexcolor_green(0x00A651),hexcolor_blue(0x00A651))
	String splitStr = ReplaceString("_biGreen0",NameOfWave(w0),"")
	splitStr = ReplaceString("_biRed0",splitStr,"")
	String labelStr = "\\JR" + ReplaceString("_", splitStr, "\r")
	TextBox/W=$plotName/C/N=text0/F=0/A=RB/X=0.00/Y=0.00/E=0 labelStr
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
// Colours are taken from Paul Tol SRON stylesheet
// Colours updated. Brighter palette for up to 6 colours, then palette of 12 for > 6
// Define colours
StrConstant SRON_1 = "0x4477aa;"
StrConstant SRON_2 = "0x4477aa;0xee6677;"
StrConstant SRON_3 = "0x4477aa;0xccbb44;0xee6677;"
StrConstant SRON_4 = "0x4477aa;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_5 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_6 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;0xaa3377;"
StrConstant SRON_7 = "0x332288;0x88ccee;0x44aa99;0x117733;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_8 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_9 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_10 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_11 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_12 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0xaa4466;0x882255;0xaa4499;"

/// @param hex		variable in hexadecimal
Function hexcolor_red(hex)
	Variable hex
	return byte_value(hex, 2) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_green(hex)
	Variable hex
	return byte_value(hex, 1) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_blue(hex)
	Variable hex
	return byte_value(hex, 0) * 2^8
End

/// @param data	variable in hexadecimal
/// @param byte	variable to determine R, G or B value
STATIC Function byte_value(data, byte)
	Variable data
	Variable byte
	return (data & (0xFF * (2^(8*byte)))) / (2^(8*byte))
End

/// @param	cond	variable for number of conditions
Function MakeColorWave(nRow, wName, [alpha])
	Variable nRow
	String wName
	Variable alpha
	
	// Pick colours from SRON palettes
	String pal
	if(nRow == 1)
		pal = SRON_1
	elseif(nRow == 2)
		pal = SRON_2
	elseif(nRow == 3)
		pal = SRON_3
	elseif(nRow == 4)
		pal = SRON_4
	elseif(nRow == 5)
		pal = SRON_5
	elseif(nRow == 6)
		pal = SRON_6
	elseif(nRow == 7)
		pal = SRON_7
	elseif(nRow == 8)
		pal = SRON_8
	elseif(nRow == 9)
		pal = SRON_9
	elseif(nRow == 10)
		pal = SRON_10
	elseif(nRow == 11)
		pal = SRON_11
	else
		pal = SRON_12
	endif
	
	Variable color
	String colorWaveFullName = "root:" + wName
	if(ParamisDefault(alpha) == 1)
		Make/O/N=(nRow,3) $colorWaveFullName
		WAVE w = $colorWaveFullName
	else
		Make/O/N=(nRow,4) $colorWaveFullName
		WAVE w = $colorWaveFullName
		w[][3] = alpha
	endif
	
	Variable i
	
	for(i = 0; i < nRow; i += 1)
		// specify colours
		color = str2num(StringFromList(mod(i, 12),pal))
		w[i][0] = hexcolor_red(color)
		w[i][1] = hexcolor_green(color)
		w[i][2] = hexcolor_blue(color)
	endfor
End

STATIC Function CleanSlate()
	String fullList = WinList("*", ";","WIN:7")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		KillDataFolder $name		
	endfor
End

// Notes:
// From illustrator text used for figures is
// Red 237,28,36 #ED1C24
// Green 0,166,81 #00A651
