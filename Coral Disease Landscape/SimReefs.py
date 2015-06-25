"""
Created on Wed Jun 24 09:40:37 2015

@author: Megan
"""
execfile('C:\Users\Megan\disease-survey-methods\Coral Disease Landscape\Landscapes.py')

# Generate a test landscape and plot Band, Estimated, and Line
Ltest = Landscape([50,10],[-2.25,1.1],0.2,0)
print "Predicted Cover =" + "{:.2f}".format(Ltest.PredictCover(1000)*100)
Ltest.GenerateColonies(4000, counts=True)
print "Number of Colonies = " + "{0:02d}".format(Ltest.nCol)
print "Calculated %Cover = " + "{:.2f}".format(Ltest.PercentCover()*100.)
print Ltest.BandTransect([5,5],25,2,plot = True)
print Ltest.EstTransect([5,5],10,1,25,4,plot = True)
print Ltest.LineIntercept([5,5],25,plot = True)

#Script to generate and save landscapes that vary in
#  nCol, prev, and (eventually) clump
# and save out
#  [SFDmu,SFDsig],[xreef,yreef],prev, nCol, Tx-type,
# check here for saving generated objects: http://stackoverflow.com/questions/4529815/saving-an-object-data-persistence-in-python
