#Module InsertLabel
#By Aron Coraor
#Last edited 7/5/16
 
import sys
 
eleDict = {"C":90,"BF":939}
 
def InsertLabel(cmlPath):
    """
    Inserts appropriate labels into a cml file. Searches through the file and
    inserts labels according to the elements of eleDict, hardcoded into InsertLabel.py,
    which is a dict of the form: {"C":90,"BF":939}. Throws a AssertionError if
    an elementType is found which is not contained in eleDict, and this element
    does not yet have a label.
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
        if cmlTextList[a+2]==' </atomArray>\r\n':
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
     
    #Check if the row is already labelled
    labelIndex= text.find("label=")
    if labelIndex != -1:
        return text
     
    #Figure out the element and the correct label
    eleIndex=text.find("elementType=")
    eleSpaceIndex=text.find(" ",eleIndex)
    element = text[eleIndex+13:eleSpaceIndex-1]
    assert element in eleDict, "Element "+`element` + " is not in dict "+`eleDict`+"."
    label = eleDict[element]
     
    #Insert the correct label in before x3
    x3Index=text.find("x3")
    text = text[:x3Index]+"label=\""+`label`+"\" "+text[x3Index:]
    return text
 
 
if __name__ == "__main__":
    if len(sys.argv)==2:
        InsertLabel(sys.argv[1])
    else:
        print "Improper number of arguments passed."
