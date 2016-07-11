#Module InsertLabel
#By Aron Coraor
#Last edited 7/5/16

import sys
import os

eleDict = {"C":90,"Rn":939}
eleShift = {"Rn":"BF"}
#eleDict = {}
#Example eleShift = {"Be":"H", "Pb":"Cs"}

def InsertLabel(cmlPath):
    """
    Inserts appropriate labels into a cml file. Searches through the file and
    inserts labels according to the elements of eleDict, hardcoded into InsertLabel.py,
    which is a dict of the form: {"C":90,"BF":939}. Throws a AssertionError if
    an elementType is found which is not contained in eleDict, and this element
    does not yet have a label.
    Finally, shifts any elementTypes found in the Keys of eleShift into the
    elementType of the corresponding Value. This should be used if Avogadro is
    being used to construct different OPLS types of the same element: e.g.,
    aromatic and aliphatic carbons or hydrogens.
    
    The workflow should be:
    1. Decide on element shifts (e.g. in Avogadro decide to draw all alkane Hydrogen-C's
    as Beryllium atoms, and all aromatic Hydrogen-C's as Hydrogen atoms)
    2. Hardcode/Add eleShift and eleDict information to InsertLabel.py
    (for this example: eleShift = {"Be":"H"}, eleDict = {"Be":85,"H":91})
    3. Draw the desired molecule in Avogadro and save the cml file
    4. Run InsertLabel.py from the terminal with the path to the cml file as an
    argument. The new labelled file will be "<originalFileName>_labelled.cml".
    
    Alternatively, if InsertLabel.py is called with no argument, it will insert
    labels into all .cml files in the current directory.
    
    Precondition: cmlPath is a valid path to a valid cml file. Only atoms of elements
    contained in eleDict are present. Atom labels which are present are directly
    before x3 listings.
    """
    #Open the file and store its contents
    cmlFile=open(cmlPath,'r')
    cmlTextList=list(cmlFile)
    cmlFile.close()
    
    #Update the contents of the file's list form
    for a in range(len(cmlTextList)):
        if '</atomArray>' in cmlTextList[a+2]:
            break
        text = cmlTextList[a+2]
        newText = InTextInsert(text)
        cmlTextList[a+2]=newText
    
    #Write the new contents to a new file, _labelled.cml
    newCmlFile = open(cmlPath[:-4]+"_labelled.cml",'w')
    for b in cmlTextList:
        newCmlFile.write(b)
    newCmlFile.close()


def InTextInsert(text):
    """
    Locates the elementType in a line of given text formatted as a line of a cml
    file, and inserts the correct label into the line of text before the x3 listing
    if a label is not already in that position. Uses eleDict as described in
    InsertLabel. Throws an AssertionError if an element is found which is not contained
    in eleDict and the element does not yet have a label. Returns the modified text.
    Precondition: text is a string which contains an elementType found in eleDict,
    and is properly formatted as an atom line of a cml file.
    """
    

    
    #Figure out the element and check if it needs to be switched.
    eleIndex=text.find("elementType=")
    eleSpaceIndex=text.find(" ",eleIndex)
    element = text[eleIndex+13:eleSpaceIndex-1]
    
    #Switch the element
    if element in eleShift:
        text = text[:eleIndex+13] + eleShift[element]+text[eleSpaceIndex-1:]
    
    #Check if the row is already labelled
    labelIndex= text.find("label=")
    if labelIndex != -1:
        return text
    print "text = "+`text`
    #Check to ensure the element is in eleDict
    assert element in eleDict, "Element "+`element` + " is not in dict "+`eleDict`+"."
    label = eleDict[element]
    
    #Insert the correct label in before x3
    x3Index=text.find("x3")
    text = text[:x3Index]+"label=\""+`label`+"\" "+text[x3Index:]
    return text


if __name__ == "__main__":
    
    #If you pass 1 filepath, InsertLabel on that file.
    #If you pass no filepath, InsertLabel on all files in current directory
    if len(sys.argv)==2:
        InsertLabel(sys.argv[1])
    elif len(sys.argv)==1:
        current = os.getcwd()
        pathList=os.listdir(current)
        fileList=[]
        for a in pathList:
            if os.path.isfile(a):
                fileList.append(a)
        for b in fileList:
            if b[-4:]==".cml" and b[-13:] != "_labelled.cml":
                InsertLabel(b)
    else:
        print "Improper number of arguments passed."